import os
import sys
import glob
import enlighten
import subprocess

from bamboo.root import gbl, loadHeader, loadBambooExtensions

# Define LorentzVector #
LV = getattr(gbl, "ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>>")

# Load extensions and HME #
loadBambooExtensions()
loadHeader("HME.h")


def runHME(path,suffix,N=None):
    import time
    path = path.replace('"',"").replace("'","")
    # Find all files that math the filename #
    files = glob.glob(path)
    chain = gbl.TChain("Events")
    for f in files:
        print ('Adding file {}'.format(f))
        chain.AddFile(f)

    # Loop over all events
    if N is None:
        N = chain.GetEntries()
    else:
        N = min(N,chain.GetEntries())
    print ("Will process {} events".format(N))
    hmes  = []
    times = []
    effs  = []
    pbar = enlighten.Counter(total=N, desc='Progress', unit='events')
    for i,event in enumerate(chain):
        pbar.update()
        if event.boosted_tag == 1:
            # Needs more inputs for boosted
            continue
        l1 = LV(event.l1_Px,event.l1_Py,event.l1_Pz,event.l1_E)
        l2 = LV(event.l2_Px,event.l2_Py,event.l2_Pz,event.l2_E)
        b1 = LV(event.j1_Px,event.j1_Py,event.j1_Pz,event.j1_E)
        b2 = LV(event.j2_Px,event.j2_Py,event.j2_Pz,event.j2_E)
        met = LV(event.met_Px,event.met_Py,event.met_Pz,event.met_E)

        start = time.time()

        hme,eff = gbl.runHME(l1,l2,b1,b2,met,event.event)

        times.append(time.time()-start)
        hmes.append(hme)
        effs.append(eff)
        
        if i >= N:
            break



    h_hme = gbl.TH1F("HME","HME",5000,0.,5000.)
    h_time= gbl.TH1F("time","time",100,0.,0.5)
    h_eff = gbl.TH1F("time","time",100,0.,100.)
    for hme,eff,time in zip(hmes,effs,times):
        h_hme.Fill(hme)
        h_time.Fill(time)
        h_eff.Fill(eff)

    C = gbl.TCanvas()
    C.SetLogy(True)

    pdfName = f"Plots/HME_{suffix}.pdf"
    C.Print(pdfName+'[')
    h_hme.GetXaxis().SetTitle("HME [GeV]")
    h_hme.GetYaxis().SetTitle("Events")
    h_hme.SetTitle(suffix)
    h_hme.Draw("hist")

    C.Print(pdfName)
    C.Clear()

    h_time.GetXaxis().SetTitle("Time [s]")
    h_time.GetYaxis().SetTitle("Events")
    h_time.SetTitle("Time")
    h_time.Draw("hist")

    C.Print(pdfName)
    C.Clear()

    h_eff.GetXaxis().SetTitle("Efficiency [%]")
    h_eff.GetYaxis().SetTitle("Events")
    h_eff.SetTitle("Efficiency")
    h_eff.Draw("hist")

    C.Print(pdfName)
    C.Clear()
    C.Print(pdfName+']')


    F = gbl.TFile(f"Plots/HME_{suffix}.root","RECREATE")
    h_hme.Write("hme")
    h_time.Write("time")
    h_eff.Write("eff")
    F.Close()


runHME(path     = sys.argv[1],
       N        = 50000,
       suffix   = sys.argv[2])

    

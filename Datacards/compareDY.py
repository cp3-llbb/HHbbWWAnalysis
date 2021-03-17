import os
import sys
import glob
import json 
from array import array
import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

def saveToJson(weight,name):
    json_dict = {'dimension': 2, 'variables': ['Eta','Pt'], 'error_type': 'absolute'}
    # Fill data #
    binning = []
    data = []
    xAxis = weight.GetXaxis()
    for i in range(1,weight.GetNbinsX()+1):
        # Get content #
        cont = weight.GetBinContent(i)
        if cont <= 0:
            continue
        # Get error #
        err = weight.GetBinError(i)

        # Get binning #
        xLow  = xAxis.GetBinLowEdge(i)
        xHigh = xAxis.GetBinUpEdge(i)
        
        if len(binning) == 0:
            binning.append(xLow)
        binning.append(xHigh)

        # Save data #
        #data.append({'bin':[xLow,xHigh],'value':cont,'error_low':err,'error_high':err}) 
        data.append({'bin':[xLow,xHigh],'value':1.,'error_low':abs(1.-cont),'error_high':abs(cont-1.)})
   
    json_dict['binning'] = {'x':[-2.5,2.5],'y':binning}
    json_dict['data'] = [{'bin':[-2.5,2.5],'values':data}]
    
    with open(f'{name}.json','w') as f:
        json.dump(json_dict,f,indent=4)
    print (f"Saved json to {name}.json")

def saveToRoot(weight,name):
    F = ROOT.TFile(f'{name}.root','RECREATE')
    weight.Write()
    F.Write()
    F.Close()

factors = {}

for era in ['2016','2017','2018']:
    #path_DY_MC          = f"/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_kinematics_DY_MC_{era}/"
    #path_DY_datadriven  = f"/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_kinematics_DY_datadrivenHT_{era}/"
    #path_DY_closure     = f"/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_kinematics_DY_closureHT_{era}/"
    path_DY_MC          = f"/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_DNN11_DY_MC_{era}/"
    path_DY_datadriven  = f"/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_DNN11_DY_datadrivenHT_{era}/"
    path_DY_closure     = f"/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_DNN11_DY_closureHT_{era}/"
    print (path_DY_MC)
    print (path_DY_datadriven)
    print (path_DY_closure)

    plotNames = [os.path.basename(name) for name in glob.glob(os.path.join(path_DY_MC,'*root'))]

    C = ROOT.TCanvas("C","C",700,800)
    C.Print(f"CompareDY/ComparisonDY_{era}.pdf[")



    for plotName in sorted(plotNames):
        print (plotName)
        C.Clear()

        f_MC = ROOT.TFile(os.path.join(path_DY_MC,plotName))
        f_datadriven = ROOT.TFile(os.path.join(path_DY_datadriven,plotName))
        f_closure = ROOT.TFile(os.path.join(path_DY_closure,plotName))

        plotName = plotName.replace('.root','')

        h_MC = f_MC.Get("DY")
        h_datadriven = f_datadriven.Get("DY")
        h_closure = f_closure.Get("DY")

#        if "resolved" in plotName:
#            xn = [0,25,50,100,150,200,300,400,500,600,1000]
#        elif "boosted" in plotName:
#            xn = [0,30,50,70,100,150,210]
        #xn = [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.]
        xn = [0.,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.9,0.95,1.]

        h_MC = h_MC.Rebin(len(xn)-1,'rebin',array('d',xn))
        h_datadriven = h_datadriven.Rebin(len(xn)-1,'rebin',array('d',xn))
        h_closure = h_closure.Rebin(len(xn)-1,'rebin',array('d',xn))

        N_MC = h_MC.Integral()
        N_datadriven = h_datadriven.Integral()
        N_closure = h_closure.Integral()

        #print (h_MC.GetNbinsX(),h_datadriven.GetNbinsX(), h_closure.GetNbinsX())

        h_MC.SetTitle("{};;Events / {:0.2f}".format(plotName,h_MC.GetXaxis().GetBinWidth(1)))
        hmax = max([h.GetMaximum() for h in [h_MC,h_datadriven,h_closure]])
        h_MC.SetMaximum(hmax*10)

        #h_MC.SetLineColor(418)
        #h_datadriven.SetLineColor(601)
        #h_closure.SetLineColor(633)
        h_MC.SetLineColor(632)
        h_datadriven.SetLineColor(600)
        h_closure.SetLineColor(1)

        h_MC.SetMarkerColor(632)
        h_datadriven.SetMarkerColor(600)
        h_closure.SetMarkerColor(1)

        h_MC.SetMarkerStyle(25)
        h_datadriven.SetMarkerStyle(8)
        h_closure.SetMarkerStyle(4)

        ratio_DY_datadriven = h_MC.Clone("ratio_DY_datadriven")
        ratio_DY_closure = h_MC.Clone("ratio_DY_closure")

        ratio_DY_datadriven.Add(h_datadriven,-1)
        ratio_DY_closure.Add(h_closure,-1)

        ratio_DY_datadriven.Divide(h_datadriven)
        ratio_DY_closure.Divide(h_closure)

#        ratio_DY_datadriven.Fit("pol1","Q")
#        ratio_DY_closure.Fit("pol1","Q")

        ratio_DY_datadriven.SetLineColor(600)
        ratio_DY_closure.SetLineColor(1)

        ratio_DY_datadriven.SetMarkerColor(600)
        ratio_DY_closure.SetMarkerColor(1)

        ratio_DY_datadriven.SetMarkerStyle(8)
        ratio_DY_closure.SetMarkerStyle(4)


        ratio_DY_datadriven.SetTitle(";;Ratio")

        C.Divide(1,2)

        pad1 = C.cd(1)
        pad1.SetLogy()
        pad1.SetTopMargin(0.1)
        pad1.SetBottomMargin(0.)
        pad1.SetLeftMargin(0.1)
        pad1.SetRightMargin(0.05)
        pad1.Draw()


        h_MC.Draw("e1p")
        h_datadriven.Draw("e1p same")
        h_closure.Draw("e1p same")

        leg_up = ROOT.TLegend(0.6,0.4,0.9,0.9)
        leg_up.AddEntry(h_MC,"#splitline{{DY (MC)}}{{Integral = {integral:0.3f}}}".format(integral=N_MC))
        leg_up.AddEntry(h_datadriven,"#splitline{{DY (data)}}{{Integral = {integral:0.3f}}}".format(integral=N_datadriven))
        leg_up.AddEntry(h_closure,"#splitline{{DY (closure)}}{{Integral = {integral:0.3f}}}".format(integral=N_closure))
        leg_up.SetBorderSize(0)
        leg_up.SetFillStyle(0)
        leg_up.Draw()

        pad2 = C.cd(2)
        pad2.SetTopMargin(0.0)
        pad2.SetBottomMargin(0.10)
        pad2.SetLeftMargin(0.1)
        pad2.SetRightMargin(0.05)
        pad2.Draw()

        ratio_DY_datadriven.Draw("e1p")
        ratio_DY_closure.Draw("e1p same")
        ratio_DY_datadriven.GetYaxis().SetRangeUser(-0.5,0.5)
        leg_down = ROOT.TLegend(0.65,0.65,0.9,1.0)
        leg_down.AddEntry(ratio_DY_datadriven,"\splitline{{#frac{{DY (MC) - DY (data)}}{{DY (data)}}}}{{Factor = {factor:0.3f}}}".format(factor=N_MC/N_datadriven))
        leg_down.AddEntry(ratio_DY_closure,"\splitline{{#frac{{DY (MC) - DY (closure)}}{{DY (closure)}}}}{{Factor = {factor:0.3f}}}".format(factor=N_MC/N_closure))
        leg_down.SetTextSize(0.04)
        leg_down.SetBorderSize(0)
        leg_down.SetFillStyle(0)
        leg_down.Draw()

        factors[plotName.replace('.root','')] = N_MC/N_closure

        C.Print(f'CompareDY/ComparisonDY_{era}.pdf','Title:{}'.format(plotName))

#        saveToJson(ratio_DY_closure,'CompareDY/'+plotName)
#        saveToRoot(ratio_DY_closure,'CompareDY/'+plotName)

        f_MC.Close()
        f_datadriven.Close()
        f_closure.Close()


    C.Print(f"CompareDY/ComparisonDY_{era}.pdf]")

with open("CompareDY/factorsDY.json",'w') as handle:
    json.dump(factors,handle,indent=4)

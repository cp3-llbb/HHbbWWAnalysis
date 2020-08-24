import os 
import sys
import copy
import urllib
import json
import traceback

from ROOT import TFile, TTree, TH1F, gROOT, TCanvas, TLegend, gStyle


class testRuns:
    def __init__(self,dict_path_dir,run_start,run_stop):
        self.dict_path_dir = dict_path_dir
        self.run_start = run_start
        self.run_stop = run_stop

        self.getJson()
        self.runAllSamples()
        input("Any key to end program")

    def getJson(self):
        #url = "https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt"
        #response = urllib.urlopen(url)
        #data = json.loads(response.read())
        with open("Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16.json") as f:
            self.golden_json = json.load(f)

    def checkRuns(self,h):
        hist_in = TH1F("hist_in","Runs inside Golden JSON",self.run_stop-self.run_start,self.run_start,self.run_stop)
        hist_out = TH1F("hist_out","Runs outside Golden JSON",self.run_stop-self.run_start,self.run_start,self.run_stop)
        for b in range(0,h.GetNbinsX()):
            run = h.GetBinLowEdge(b)
            if str(int(run)) in self.golden_json.keys():
                hist_in.SetBinContent(b,h.GetBinContent(b))
            else:
                hist_out.SetBinContent(b,h.GetBinContent(b))

        return (hist_in,hist_out)

    def listAllFiles(self,path_dir):
        list_files = []
        for subdir, dirs, files in os.walk(path_dir):
            for f in files:
                list_files.append(os.path.join(subdir, f))
        return list_files

    def getRunHisto(self,f):
        print ("Extracting run hist from file %s"%f)
        try:
            root_file = TFile.Open(f)
            tree = root_file.Get("Runs")
            tree.Draw("run>>h("+str(self.run_stop-self.run_start)+","+str(self.run_start)+","+str(self.run_stop)+")")
            h = copy.deepcopy(gROOT.FindObject("h"))
            root_file.Close()
            return h
        except Exception as e:
            print ("Exception %s for file %s"%(e,f))

    def addHisto(self,list_hist):
        hist = copy.deepcopy(list_hist[0])
        for i in range(1,len(list_hist)):
            hist.Add(list_hist[i])
        return hist
            
    def runSample(self,name,path_dir):
        list_files = self.listAllFiles(path_dir)
        list_hist = [self.getRunHisto(f) for f in list_files]
        run_hist = copy.deepcopy(self.addHisto(list_hist))
        run_hist_in,run_hist_out = self.checkRuns(run_hist)
        return (run_hist,run_hist_in,run_hist_out) 
        

    def runAllSamples(self):
        gStyle.SetOptStat(0)
        gROOT.SetBatch(True)
        dict_run_hist = dict((k,self.runSample(k,h)) for k,h in self.dict_path_dir.items())
        max_list = [tuple(h.GetMaximum() for h in hist) for hist in dict_run_hist.values()]
        max_h = max(max_list)
        gROOT.SetBatch(False)

        legend = TLegend(0.15,0.7,0.45,0.83)
        legend.SetHeader("Legend","C")
        
        c = TCanvas("C","C",2000,500)
        c.Divide(3,1)
        c.SetLogy()

        i = 0
        for k,(h,h_in,h_out) in dict_run_hist.items():
            h.SetLineColor(i+1)
            h_in.SetLineColor(i+1)
            h_out.SetLineColor(i+1)
            legend.AddEntry(h,k,"l")

            opt = "hist"
            if i!=0:
                opt += " same"

            h.SetTitle("Runs")
            h.SetMaximum(max_h[0])
            h_in.SetMaximum(max_h[1])
            h_out.SetMaximum(max_h[2])

            c.cd(1)
            h.Draw(opt)
            c.cd(2)
            h_in.Draw(opt)
            c.cd(3)
            h_out.Draw(opt)


            i+=1

        legend.Draw()
        c.Print("test.png")
        



if __name__=="__main__":
#    testRuns({"SingleMuon 2016B ver1": "/storage/data/cms/store/data/Run2016B_ver1/SingleMuon/NANOAOD/Nano1June2019_ver1-v1/",
#              "SingleMuon 2016B ver2": "/storage/data/cms/store/data/Run2016B_ver2/SingleMuon/NANOAOD/Nano1June2019_ver2-v1/"},
#             272007,
#             275376)
    testRuns({"DoubleMuon 2016B ver1": "/storage/data/cms/store/data/Run2016B_ver1/DoubleMuon/NANOAOD/Nano1June2019_ver1-v1/",
              "DoubleMuon 2016B ver2": "/storage/data/cms/store/data/Run2016B_ver2/DoubleMuon/NANOAOD/Nano1June2019_ver2-v1/"},
             272007,
             275376)




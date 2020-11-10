import os
import sys
import yaml
import copy
import glob
import json
import ROOT
import numpy as np

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

class BtagReweightingRatio:
    def __init__(self, path_on, path_off, list_hists, outputname, pdfname,era):
        self.path_on    = path_on
        self.path_off   = path_off
        self.path_off   = path_off
        self.list_hists = list_hists
        self.era        = era
        path_output = os.path.join(os.path.abspath(os.path.dirname(__file__)),'Reweighting_plots')
        if not os.path.exists(path_output):
            os.makedirs(path_output)
        self.outputname = os.path.join(path_output,outputname)
        self.pdfname = pdfname

        self.histograms_off = self.getAllHistograms(self.path_off,self.list_hists)
        self.histograms_on = self.getAllHistograms(self.path_on,self.list_hists)
        self.ratio = self.makeRatio(self.histograms_off,self.histograms_on)

        for idx,sample in enumerate(sorted(self.ratio.keys())):
            canvas,legend = self.plot(sample,self.histograms_on[sample],self.histograms_off[sample],self.ratio[sample])
            if idx == 0:
                canvas.Print(self.pdfname+'[')
            canvas.Print(self.pdfname,'Title:%s'%sample)
            self.saveToJson(sample,self.ratio[sample])
            self.saveToRoot(sample,self.ratio[sample])
        canvas.Print(self.pdfname+']')


    @staticmethod
    def getHistogram(rootfile,histname):
        f = ROOT.TFile(rootfile)
        if f.GetListOfKeys().Contains(histname):
            h = copy.deepcopy(f.Get(histname))
        else:
            print ("Could not find hist %s in %s"%(histname,rootfile))
            h = None
        f.Close()
        return h

    def getAllHistograms(self,path,hists):
        hist_dict = {}
        for f in glob.glob(os.path.join(path,'results','*.root')):
            if '__skeleton__' in f:
                continue
            hist = None
            for hname in hists:
                h = self.getHistogram(f,hname)
                if h is None:
                    raise RuntimeError('Could not find histogram from sample %s'%f)
                if hist is None:
                    hist = h.Clone()
                else:
                    hist.Add(h)
            for i in range(1,hist.GetNbinsX()+1):
                hist.SetBinError(i,1e-8)

            hist_dict[os.path.basename(f).replace('.root','')] = hist
        return hist_dict

    @staticmethod
    def makeRatio(dict_num,dict_den):
        ratio_dict = {}
        for sample in dict_num.keys():
            r = dict_num[sample].Clone("ratio")
            r.Divide(dict_den[sample])
            ratio_dict[sample] = r
        return ratio_dict

    def plot(self,sample,hist_on,hist_off,ratio):
        C = ROOT.TCanvas("c1", "c1", 600, 600)
        C.Divide(1,2)
        colors = [418,602]
        leg = ROOT.TLegend(0.65,0.60,0.90,0.80)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.05)

        C.cd(1)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetTopMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.01)
        ROOT.gPad.SetLogy()

        ROOT.gStyle.SetTitleFontSize(0.1)

        hist_on.SetTitle(sample)
        hist_on.GetYaxis().SetTitle("Event weight sum")
        hist_on.GetYaxis().SetTitleOffset(1.)
        hist_on.GetYaxis().SetLabelSize(0.05)
        hist_on.GetYaxis().SetTitleSize(0.06)
        hist_on.GetXaxis().SetTitle("")
        hist_on.GetXaxis().SetLabelSize(0.)

        hist_on.SetLineWidth(2)
        hist_on.SetLineColor(418)
        hist_off.SetLineWidth(2)
        hist_off.SetLineColor(602)

        hmax = max([hist_on.GetMaximum(),hist_off.GetMaximum()])
        hist_on.SetMaximum(hmax*1.5)

        leg.AddEntry(hist_on,'After btag weight')
        leg.AddEntry(hist_off,'Before btag weight')

        hist_on.Draw("H")
        hist_off.Draw("H same")
    
        leg.Draw()
    
        C.cd(2)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetTopMargin(0.01)
        ROOT.gPad.SetBottomMargin(0.15)
        ROOT.gPad.SetGridy()

        ratio.SetTitle("")
        ratio.SetLineWidth(2)
        ratio.SetLineColor(1)

        ratio.SetMaximum(ratio.GetMaximum()*1.1)
        ratio.SetMinimum(ratio.GetMinimum()*0.9)

        ratio.GetYaxis().SetRangeUser(0.0,1.5)


        ratio.GetYaxis().SetTitle("#frac{#sum w(before)}{#sum w(after)}")
        ratio.GetYaxis().SetLabelSize(0.05)
        ratio.GetYaxis().SetTitleSize(0.06)
        ratio.GetYaxis().SetTitleOffset(1.)

        ratio.GetXaxis().SetTitle("N(Ak4 jets)")
        ratio.GetXaxis().SetLabelSize(0.05)
        ratio.GetXaxis().SetTitleSize(0.06)

        ratio.Draw("PE")

        return C,leg


    def saveToJson(self,sample,ratio):
        json_dict = {'dimension': 1, 'variables': ['NumJets'], 'error_type': 'absolute'}
        # Fill data #
        binning = []
        data = []
        xAxis = ratio.GetXaxis()
        for i in range(1,ratio.GetNbinsX()+1):
            # Get content #
            cont = ratio.GetBinContent(i)
            if cont <= 0:
                continue
            # Get error #
            err = ratio.GetBinError(i)

            # Get binning #
            xLow  = xAxis.GetBinLowEdge(i)
            xHigh = xAxis.GetBinUpEdge(i)
            
            if len(binning) == 0:
                binning.append(xLow)
            binning.append(xHigh)

            # Save data #
            #data.append({'bin':[xLow,xHigh],'value':cont,'error_low':err,'error_high':err}) 
            data.append({'bin':[xLow,xHigh],'value':cont,'error_low':0.,'error_high':0.}) 
       
        json_dict['binning'] = {'x':binning}
        json_dict['data'] = data
        with open(self.outputname+'_%s_%s.json'%(sample,self.era),'w') as f:
            json.dump(json_dict,f,indent=4)
        print ("Saved json to",self.outputname+'_%s_%s.json'%(sample,self.era))

            
    def saveToRoot(self,sample,ratio):
        root_file = ROOT.TFile(self.outputname+'_%s_%s.root'%(sample,self.era),"recreate")
        ratio.Write("ratio")
        root_file.Write()
        root_file.Close()
        
        
        
if __name__ == "__main__":
    #---- 2016 ----#
    instance = BtagReweightingRatio(path_on     = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016NanoV7_BtagReweightingOn/',
                                    path_off    = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016NanoV7_BtagReweightingOff/',
                                    list_hists  = ['NoChannel_NoSelection_Ak4Jets_N'],
                                    outputname  = 'BtagReweightingRatio_jetN',
                                    pdfname     = 'BtagReweightingRatio_Background_2016.pdf',
                                    era         = '2016')
#    instance = BtagReweightingRatio(path_on     = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016NanoV6_BtagReweighting_On_signal/',
#                                    path_off    = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016NanoV6_BtagReweighting_Off_signal/',
#                                    list_hists  = ['NoChannel_NoSelection_Ak4Jets_N'],
#                                    outputname  = 'BtagReweightingRatio_jetN',
#                                    pdfname     = 'BtagReweightingRatio_Signal_2016.pdf',
#                                    era         = '2016')
    instance = BtagReweightingRatio(path_on     = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016NanoV7_signal_BtagReweightingOn',
                                    path_off    = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2016NanoV7_signal_BtagReweightingOff',
                                    list_hists  = ['NoChannel_NoSelection_Ak4Jets_N'],
                                    outputname  = 'BtagReweightingRatio_jetN',
                                    pdfname     = 'BtagReweightingRatio_Signal_2016.pdf',
                                    era         = '2016')

#    #---- 2017 ----#
    instance = BtagReweightingRatio(path_on     = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2017NanoV7_BtagReweightingOn/',
                                    path_off    = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2017NanoV7_BtagReweightingOff/',
                                    list_hists  = ['NoChannel_NoSelection_Ak4Jets_N'],
                                    outputname  = 'BtagReweightingRatio_jetN',
                                    pdfname     = 'BtagReweightingRatio_Background_2017.pdf',
                                    era         = '2017')
#    instance = BtagReweightingRatio(path_on     = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2017NanoV6_BtagReweighting_On_signal/',
#                                    path_off    = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2017NanoV6_BtagReweighting_Off_signal/',
#                                    list_hists  = ['NoChannel_NoSelection_Ak4Jets_N'],
#                                    outputname  = 'BtagReweightingRatio_jetN',
#                                    pdfname     = 'BtagReweightingRatio_Signal_2017.pdf',
#                                    era         = '2017')
#    #---- 2018 ----#
    instance = BtagReweightingRatio(path_on     = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2018NanoV7_BtagReweightingOn/',
                                    path_off    = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2018NanoV7_BtagReweightingOff/',
                                    list_hists  = ['NoChannel_NoSelection_Ak4Jets_N'],
                                    outputname  = 'BtagReweightingRatio_jetN',
                                    pdfname     = 'BtagReweightingRatio_Background_2018.pdf',
                                    era         = '2018')
#    instance = BtagReweightingRatio(path_on     = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2018NanoV6_BtagReweighting_On_signal/',
#                                    path_off    = '/nfs/scratch/fynu/fbury/BambooOutputHHtobbWW/full2018NanoV6_BtagReweighting_Off_signal/',
#                                    list_hists  = ['NoChannel_NoSelection_Ak4Jets_N'],
#                                    outputname  = 'BtagReweightingRatio_jetN',
#                                    pdfname     = 'BtagReweightingRatio_Signal_2018.pdf',
#                                    era         = '2018')

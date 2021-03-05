#!/usr/bin/env python
import re
import math
import glob
import csv
import yaml
import os
import sys
import logging
import copy
import argparse
import ROOT
from collections import defaultdict

def get_options():
    p = argparse.ArgumentParser(description='Skimmer Ntuples Check')
    p.add_argument("-v", "--verbose", action="store_true", default=False,
                   help="increase output verbosity")
    p.add_argument('--config', action='store', required=True, type=str,
                   help = 'Name of Config File')
    p.add_argument('--era', action='store', required=True, type=str,
                   help = 'Era')
    p.add_argument('--category', action='store', required=True, type=str,
                   help = 'resolved2b2j || resolved2b1j || resolved1b2j || resolved2b0j')
    p.add_argument('--channel', action='store', required=True, type=str,
                   help = 'El || Mu')
    p.add_argument('--dostack', action='store_true', required=False,
                   help = 'to make stack plot')
    p.add_argument('--donorm', action='store_true', required=False,
                   help = 'to make normalised plot')


    args = p.parse_args()
    return args

def load_input (file):
    with open(file, 'r') as inf:
        return yaml.safe_load(inf)

def add_hists(histList, node):
    hNode = histList[0].Clone()
    hNode.Reset()
    for h in histList:
        hNode.Add(h, 1.0)
    hNode.SetName(node)
    #print('Node: {}, nEntries: {}'.format(node,hNode.GetEntries()))
    return hNode

def make_stack(histDict, nodeList, partName):
    if not os.path.exists(partName):
        os.mkdir(partName)
    for var, nodehist in histDict.items():
        print('Making Stack historam of {}'.format(var))
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetTextFont(42)
        d = ROOT.TCanvas("d", "", 1200, 900)
        pad = ROOT.TPad('pad-%s' % var, '', 0, 0, 1, 1)
        #pad.SetLogy()
        pad.SetGrid()
        pad.Draw()
        pad.cd()  
        legend = ROOT.TLegend(0.62, 0.70, 0.82, 0.88)
        legend.SetFillColor(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.03)
        hs = ROOT.THStack()
        signalHists = []
        for i,node in enumerate(nodeList):
            #if node == 'HH':
            hsig = ROOT.TH1F()
            if 'HH' in node :
                hsig = nodehist[i].Clone()
                hsig.Scale(20000)
                hsig.SetLineColor(2*i+5)
                hsig.SetLineWidth(2)
                legend.AddEntry(hsig, node, "l")
                signalHists.append(hsig)
                hsig.SetDirectory(0)
            else :
                legend.AddEntry(nodehist[i], node, "f")
                nodehist[i].SetFillColor(2*i + 6)
                hs.Add(nodehist[i])
        hs.Draw("HIST")
        for i,hist in enumerate(signalHists):
            hist.Draw("SAME HIST")
            d.Update()
        legend.Draw()
        d.SaveAs(partName+'/'+var+'.root')


def make_norm(histDict, nodeList, partName, colors):
    if not os.path.exists(partName):
        os.mkdir(partName)
    for var, nodehist in histDict.items():
        print('Making Normalised historam of {}'.format(var))
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetTextFont(42)
        d = ROOT.TCanvas("d", "", 1200, 900)
        pad = ROOT.TPad('pad-%s' % var, '', 0, 0, 1, 1)
        #pad.SetLogy()
        pad.SetGrid()
        pad.Draw()
        pad.cd()  
        legend = ROOT.TLegend(0.62, 0.70, 0.82, 0.88)
        legend.SetFillColor(0)
        legend.SetBorderSize(0)
        legend.SetTextSize(0.03)
        for i,node in enumerate(nodeList):
            nodehist[i].Scale(1/nodehist[i].Integral())
            #nodehist[i].SetLineColor(i + 1)
            nodehist[i].SetLineColor(colors[i])
            nodehist[i].SetLineWidth(2)
            if i == 0 : 
                nodehist[i].Draw('HIST')
            else : 
                nodehist[i].Draw('HIST SAME')
            legend.AddEntry(nodehist[i], node, "l")
            d.Update()
        #hs.Draw("HIST")
        #hsig.Draw("SAME HIST")
        legend.Draw()
        #d.SaveAs(partName+var+'.pdf')
        d.SaveAs(partName+'/'+var+'.root')
        

def main():
    # Enable multi-threading
    ROOT.ROOT.EnableImplicitMT(30)
    args      = get_options()
    logging.basicConfig(level=logging.INFO,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)

    pwd = os.getcwd()    
    configPath = os.path.join(pwd, args.config)
    if not os.path.exists(configPath):
        logging.info('Config file does not exist')
    #infoDict  = load_input('sampleListSL.yml')
    infoDict  = load_input(configPath)
    basepath  = str(infoDict.get('sampleDir'))
    channels  = infoDict.get('channelList')
    variablesInfo = infoDict.get('variableList')
    eras      = list(infoDict.get('sampleDict').keys())
    xsecs     = infoDict.get('xSecDict')
    logging.info('Eras: %s'%(eras))
    lumiDict = {'2016':35922, '2017':41529.152060112, '2018':59740.565201546}
    colors = infoDict.get('colors')
    print(colors)
    catSampleDict = infoDict.get('sampleDict').get(args.era)

    nodes = []
    histsToStack = defaultdict(list)
    for key, val in catSampleDict.items():
        if not key.startswith(args.category+'_'+args.channel):
            continue
        logging.info('Category,Channel,Node: %s'%(key))
        node = key.split('_')[-1]
        nodes.append(node)
        histsPerNode = defaultdict(list)
        for sample in val:
            if str(sample).startswith('XX'):
                continue
            ntuple = os.path.join(basepath,sample)
            tfile = ROOT.TFile.Open(ntuple, "READ")
            #evTree = tfile.Get("Events")
            #print(tfile.GetListOfKeys)
            rdf = ROOT.RDataFrame("Events", tfile)
            rdfRuns = ROOT.RDataFrame("Runs", tfile)
            eventWtSum = rdfRuns.Sum("genEventSumw").GetValue()
            lumi = lumiDict.get(args.era)
            xsec = xsecs.get(ntuple.split('/')[-1])
            lumiweight = lumi * xsec / eventWtSum
            rdf = rdf.Define("lweight","{}".format(lumiweight)).Define("weight","lweight * total_weight")
            for varinfo in variablesInfo:
                var   = varinfo.split('/')[0]
                nbins = int(varinfo.split('/')[1])
                xmin  = float(varinfo.split('/')[2])
                xmax  = float(varinfo.split('/')[3])
                hvarPoint  = rdf.Histo1D((var, "", nbins, xmin, xmax), var,"weight")
                hvar = hvarPoint.GetValue()
                histsPerNode[str(var)].append(hvar)
            logging.info('Samples: {} | \tnEntries: {}'.format(ntuple,rdf.Count().GetValue()))
        for varinfo in variablesInfo:
            var = varinfo.split('/')[0]
            htemp = add_hists(histsPerNode.get(str(var)), node)
            histsToStack[str(var)].append(htemp)

    
    fname_stack = args.category+'_'+args.channel+'_'+args.era+'_stack'
    fname_norm  = args.category+'_'+args.channel+'_'+args.era+'_norm'
    #print(histsToStack)
    #print(nodes)
    
    if args.dostack : 
        make_stack(histsToStack,nodes,fname_stack)
    if args.donorm : 
        make_norm(histsToStack,nodes,fname_norm,colors)
        

if __name__ == "__main__":
    main()

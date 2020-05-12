import os.path
from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule
from bamboo import treefunctions as op    
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.plots import Plot, EquidistantBinning, SummedPlot


def plotNumber(contName,cont,sel,Nmax,xTitle):
    return [Plot.make1D(contName,
                        op.rng_len(cont),
                        sel,
                        EquidistantBinning(Nmax,0.,Nmax),
                        xTitle=xTitle)]

def plotJetReco(contName,cont,sel,partName):
    plots = []

    plots.append(Plot.make1D(contName+"_deltaR",
                             op.map(cont,lambda m : op.deltaR(m.p4,m.genJet.p4)),
                             sel,
                             EquidistantBinning(100,0.,0.2),
                             xTitle="#Delta R(%s_{gen},%s_{reco})"%(partName,partName)))
    plots.append(Plot.make1D(contName+"_deltaR2",
                             op.map(cont,lambda m : op.pow(op.deltaR(m.p4,m.genJet.p4),2)),
                             sel,
                             EquidistantBinning(100,0.,0.04),
                             xTitle="#Delta R^{2}(%s_{gen},%s_{reco})"%(partName,partName)))

    plots.append(Plot.make1D(contName+"_deltaPtRel",
                             op.map(cont,lambda m : (m.pt-m.genJet.pt)/m.pt),
                             sel,
                             EquidistantBinning(100,-2.,2.),
                             xTitle="P_{T}(%s_{reco})-P_{T}(%s_{gen})/P_{T}(%s_{reco})"%(partName,partName,partName)))
    return plots

def plotMatching(contName,list_cont,list_sel,partName):
    plots = []

    dr_plots = []
    dr2_plots = []
    ptrel_plots = []
    TF_plots = []
    for idx,(cont,sel) in enumerate(zip(list_cont,list_sel)):
        dr_plots.append(Plot.make1D(contName+str(idx)+"_deltaR",
                                 op.map(cont,lambda m : op.deltaR(m[0].p4,m[1].p4)), 
                                 sel,
                                 EquidistantBinning(100,0.,0.2),
                                 xTitle="#Delta R(%s_{gen},%s_{reco})"%(partName,partName)))
        dr2_plots.append(Plot.make1D(contName+str(idx)+"_deltaR2",
                                 op.map(cont,lambda m : op.pow(op.deltaR(m[0].p4,m[1].p4),2)), 
                                 sel,
                                 EquidistantBinning(100,0.,0.04),
                                 xTitle="#Delta R^{2}(%s_{gen},%s_{reco})"%(partName,partName)))

        ptrel_plots.append(Plot.make1D(contName+str(idx)+"_deltaPtRel",
                                 op.map(cont,lambda m : (m[1].pt-m[0].pt)/m[1].pt),
                                 sel,
                                 EquidistantBinning(100,-2.,2.),
                                 xTitle="P_{T}(%s_{reco})-P_{T}(%s_{gen})/P_{T}(%s_{reco})"%(partName,partName,partName)))
        TF_plots.append(Plot.make2D(contName+str(idx)+"_TF",
                                 (op.map(cont,lambda m : m[0].p4.E()),op.map(cont,lambda m : m[1].p4.E()-m[0].p4.E())),
                                 sel,
                                 [EquidistantBinning(100,0.,500.),EquidistantBinning(400,-200.,200.)],
                                 xTitle="E^{parton}(e^{-})",
                                 yTitle="#Delta E = E^{reco}(%s)-E^{parton}(%s)"%(partName,partName)))

    plots.append(SummedPlot(contName+"_deltaR",dr_plots,xTitle="#Delta R(%s_{gen},%s_{reco})"%(partName,partName)))
    plots.append(SummedPlot(contName+"_deltaR2",dr2_plots,xTitle="#Delta R^{2}(%s_{gen},%s_{reco})"%(partName,partName)))
    plots.append(SummedPlot(contName+"_deltaPtRel",ptrel_plots,xTitle="P_{T}(%s_{reco})-P_{T}(%s_{gen}/P_{T}(%s_{reco})"%(partName,partName,partName)))
    plots.append(SummedPlot(contName+"_TF",TF_plots,xTitle="#Delta E = E^{reco}(%s)-E^{patron}(%s)"%(partName,partName)))
    
    return plots



class TransferFunction(NanoAODHistoModule):
    def addArgs(self, parser):
        super(TransferFunction, self).addArgs(parser)
        parser.add_argument("--backend", type=str, default="dataframe", help="Backend to use, 'dataframe' (default) or 'lazy'")

    def prepareTree(self, tree, sample=None, sampleCfg=None):
        from bamboo.treedecorators import NanoAODDescription, nanoRochesterCalc, nanoJetMETCalc
        # JEC's Recommendation for Full RunII: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # JER : -----------------------------: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
        # Get base aguments #
        era = sampleCfg['era']
        self.is_MC = self.isMC(sample) # no confusion between boolean is_MC and method isMC()
        metName = "METFixEE2017" if era == "2017" else "MET"
        tree,noSel,be,lumiArgs = super(TransferFunction,self).prepareTree(tree          = tree, 
                                                                          sample        = sample, 
                                                                          sampleCfg     = sampleCfg, 
                                                                          description   = NanoAODDescription.get("v5", year=(era if era else "2016"),
                                                                          isMC          = self.is_MC,
                                                                          systVariations= [ nanoRochesterCalc, (nanoJetMETCalc_METFixEE2017 if era == "2017" else nanoJetMETCalc) ]), 
                                                                                # will do Jet and MET variations, and the Rochester correction
                                                                          lazyBackend   = self.args.backend == "lazy" or self.args.onlypost)
    
        from bamboo.analysisutils import configureJets ,configureRochesterCorrection, configureType1MET 

        ## Check distributed option #
        isNotWorker = (self.args.distributed != "worker") 

        # Check era #
        if era != "2016" and era != "2017" and era != "2018":
            raise RuntimeError("Unknown era {0}".format(era))

        # Rochester and JEC corrections (depends on era) #     
        cachJEC_dir = '/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/cacheJEC'

        if era == "2016":
            # Rochester corrections #
            configureRochesterCorrection(variProxy  = tree._Muon,
                                         paramsFile = os.path.join(os.path.dirname(__file__), "data", "RoccoR2016.txt"),
                                         isMC       = self.is_MC,
                                         backend    = be, 
                                         uName      = sample)
            # JetMET treatment #
            if self.is_MC:   # if MC -> needs smearing
                configureJets(variProxy             = tree._Jet, 
                              jetType               = "AK4PFchs",
                              jec                   = "Summer16_07Aug2017_V20_MC",
                              smear                 = "Summer16_25nsV1_MC",
                              jesUncertaintySources = ["Total"],
                              mayWriteCache         = isNotWorker,
                              isMC                  = self.is_MC,
                              backend               = be, 
                              uName                 = sample,
                              cachedir              = cachJEC_dir)
                configureType1MET(variProxy             = tree._MET, 
                                  jec                   = "Summer16_07Aug2017_V20_MC",
                                  smear                 = "Summer16_25nsV1_MC",
                                  jesUncertaintySources = ["Total"],
                                  mayWriteCache         = isNotWorker,
                                  isMC                  = self.is_MC,
                                  backend               = be,
                                  uName                 = sample,
                                  cachedir              = cachJEC_dir)

            else:                   # If data -> extract info from config 
                jecTag = None
                if "2016B" in sample or "2016C" in sample or "2016D" in sample:
                    jecTag = "Summer16_07Aug2017BCD_V11_DATA"
                elif "2016E" in sample or "2016F" in sample:
                    jecTag = "Summer16_07Aug2017EF_V11_DATA"
                elif "2016G" in sample or "2016H" in sample:
                    jecTag = "Summer16_07Aug2017GH_V11_DATA"
                configureJets(variProxy             = tree._Jet, 
                              jetType               = "AK4PFchs",
                              jec                   = jecTag,
                              mayWriteCache         = isNotWorker,
                              isMC                  = self.is_MC,
                              backend               = be, 
                              uName                 = sample,
                              cachedir              = cachJEC_dir)
                configureType1MET(variProxy         = tree._MET, 
                                  jec               = jecTag,
                                  mayWriteCache     = isNotWorker,
                                  isMC              = self.is_MC,
                                  backend           = be, 
                                  uName             = sample,
                                  cachedir          = cachJEC_dir)

        #----- Pile up -----#
        puWeightsFile = None
        # Select sfTaf #
        if era == "2016":
            self.sfTag="94X"
        elif era == "2017":
            self.sfTag="94X"
        elif era == "2018":
            self.sfTag="102X"
        
        # Get MC PU weight file #
        if self.is_MC:
            if "pufile" not in sampleCfg:
                raise KeyError("Could not find 'pufile' entry for sample %s in the YAML file"%sampleCfg["sample"])
            puWeightsFile = os.path.join(os.path.dirname(__file__), "data", "pileup",sampleCfg["pufile"])
            if not os.path.exists(puWeightsFile):
                raise RuntimeError("Could not find pileup file %s"%puWeightsFile)
            from bamboo.analysisutils import makePileupWeight
            PUWeight = makePileupWeight(puWeightsFile, tree.Pileup_nTrueInt, systName="pileup",nameHint=f"puweightFromFile{sample}".replace('-','_'))                                                      
            noSel = noSel.refine("puWeight", weight=PUWeight)


        return tree,noSel,be,lumiArgs


    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        era = sampleCfg['era']   
        plots = []
        
        ##########################################################
        #                       Leptons                          #
        ##########################################################
        #----- Gen particles -----# 
        gen_e_minus = op.select(t.GenPart,lambda g : op.AND( g.pdgId==11 , g.statusFlags & ( 0x1 << 13), g.pt>10, op.abs(g.eta)<2.5 ))
        gen_e_plus = op.select(t.GenPart,lambda g : op.AND( g.pdgId==-11 , g.statusFlags & ( 0x1 << 13), g.pt>10, op.abs(g.eta)<2.5 ))
        gen_mu_minus = op.select(t.GenPart,lambda g : op.AND( g.pdgId==13 , g.statusFlags & ( 0x1 << 13), g.pt>10, op.abs(g.eta)<2.4  ))
        gen_mu_plus = op.select(t.GenPart,lambda g : op.AND( g.pdgId==-13 , g.statusFlags & ( 0x1 << 13), g.pt>10, op.abs(g.eta)<2.4  ))

        plots.extend(plotNumber("N_e_minus_gen",gen_e_minus,noSel,10,"N(e^{-}) gen level"))
        plots.extend(plotNumber("N_e_plus_gen",gen_e_plus,noSel,10,"N(e^{+}) gen level"))
        plots.extend(plotNumber("N_mu_minus_gen",gen_mu_minus,noSel,10,"N(\mu^{-}) gen level"))
        plots.extend(plotNumber("N_mu_plus_gen",gen_mu_plus,noSel,10,"N(\mu^{+}) gen level"))

        #----- Reco particles -----# 
        reco_e_minus = op.select(t.Electron, lambda ele : op.AND(ele.charge == -1 , ele.pdgId==11 , ele.pt>10 , op.abs(ele.eta)<2.5,
                                                                 op.OR(ele.genPartFlav==1, ele.genPartFlav==15, ele.genPartFlav==22)))
        reco_e_plus = op.select(t.Electron, lambda ele : op.AND(ele.charge == +1 , ele.pdgId==-11 , ele.pt>10 , op.abs(ele.eta)<2.5,
                                                                 op.OR(ele.genPartFlav==1, ele.genPartFlav==15, ele.genPartFlav==22)))
        reco_mu_minus = op.select(t.Muon, lambda mu : op.AND(mu.charge == -1 , mu.pdgId==13 , mu.pt>10 , op.abs(mu.eta)<2.4, 
                                                                 op.OR(mu.genPartFlav==1, mu.genPartFlav==15)))
        reco_mu_plus = op.select(t.Muon, lambda mu : op.AND(mu.charge == +1 , mu.pdgId==-13 , mu.pt>10 , op.abs(mu.eta)<2.4 ,
                                                                 op.OR(mu.genPartFlav==1, mu.genPartFlav==15)))

        plots.extend(plotNumber("N_e_minus_reco",reco_e_minus,noSel,10,"N(e^{-}) reco level"))
        plots.extend(plotNumber("N_e_plus_reco",reco_e_plus,noSel,10,"N(e^{+}) reco level"))
        plots.extend(plotNumber("N_mu_minus_reco",reco_mu_minus,noSel,10,"N(\mu^{-}) reco level"))
        plots.extend(plotNumber("N_mu_plus_reco",reco_mu_plus,noSel,10,"N(\mu^{+}) reco level"))

        #---- Matching -----#
        lambda_lepton_match = lambda l_gen,l_reco : op.AND(op.deltaR(l_gen.p4,l_reco.p4)<0.1, op.abs(l_gen.pt-l_reco.pt)/l_gen.pt<0.2)

        match_e_minus = op.combine((gen_e_minus,reco_e_minus), pred=lambda_lepton_match)
        match_e_plus = op.combine((gen_e_plus,reco_e_plus), pred=lambda_lepton_match)
        match_mu_minus = op.combine((gen_mu_minus,reco_mu_minus), pred=lambda_lepton_match)
        match_mu_plus = op.combine((gen_mu_plus,reco_mu_plus), pred=lambda_lepton_match)
        
        SelEMinusMatch = noSel.refine("SelEMinusMatch",cut=[op.rng_len(match_e_minus)>=1])
        SelEPlusMatch = noSel.refine("SelEPlusMatch",cut=[op.rng_len(match_e_plus)>=1])
        SelMuMinusMatch = noSel.refine("SelMuMinusMatch",cut=[op.rng_len(match_mu_minus)>=1])
        SelMuPlusMatch = noSel.refine("SelMuPlusMatch",cut=[op.rng_len(match_mu_plus)>=1])

        plots.extend(plotNumber("N_e_minus_match",match_e_minus,noSel,10,"N(e^{-}) match level"))
        plots.extend(plotNumber("N_e_plus_match",match_e_plus,noSel,10,"N(e^{+}) match level"))
        plots.extend(plotNumber("N_mu_minus_match",match_mu_minus,noSel,10,"N(\mu^{-}) match level"))
        plots.extend(plotNumber("N_mu_plus_match",match_mu_plus,noSel,10,"N(\mu^{+}) match level"))

        plots.extend(plotMatching("e",[match_e_minus,match_e_plus],[SelEMinusMatch,SelEPlusMatch],"e^{#pm}"))
        plots.extend(plotMatching("mu",[match_mu_minus,match_mu_plus],[SelMuMinusMatch,SelMuPlusMatch],"#mu^{#pm}"))


        ##########################################################
        #                       B jets                           #
        ##########################################################
        #----- RecoJets -----#
        #recoJet_b = op.select(t.Jet, lambda j : op.AND(j.partonFlavour == 5, j.pt>10, op.abs(j.eta)<2.4 ,(op.abs(j.pt-j.genJet.pt)/j.genJet.pt)<1 ))
        #recoJet_bbar = op.select(t.Jet, lambda j : op.AND(j.partonFlavour == -5, j.pt>10, op.abs(j.eta)<2.4 ,(op.abs(j.pt-j.genJet.pt)/j.genJet.pt)<1  ))
        recoJet_b = op.select(t.Jet, lambda j : op.AND(j.partonFlavour == 5, j.pt>10, op.abs(j.eta)<2.4, j.genJet.isValid))
        recoJet_bbar = op.select(t.Jet, lambda j : op.AND(j.partonFlavour == -5, j.pt>10, op.abs(j.eta)<2.4, j.genJet.isValid))
        plots.extend(plotNumber("N_jetb_reco",recoJet_b,noSel,5,"N(b) reco jet"))
        plots.extend(plotNumber("N_jetbbar_reco",recoJet_bbar,noSel,5,"N(#bar{b}) reco jet"))

        SelBReco = noSel.refine("SelBReco",cut=[op.rng_len(recoJet_b)>=1])
        SelBBarReco = noSel.refine("SelBBarReco",cut=[op.rng_len(recoJet_bbar)>=1])

#        plots.extend(plotJetReco("bjet",recoJet_b,SelBReco,"b jet"))
#        plots.extend(plotJetReco("bbarjet",recoJet_bbar,SelBBarReco,"#bar{b} jet"))
#
#
#       #----- Parton -----#
#        gen_b = op.select(t.GenPart,lambda g : op.AND( g.pdgId==5 , g.statusFlags & ( 0x1 << 13) , g.pt>10, op.abs(g.eta)<2.4))
#        gen_bbar = op.select(t.GenPart,lambda g : op.AND( g.pdgId==-5 , g.statusFlags & ( 0x1 << 13) , g.pt>10, op.abs(g.eta)<2.4))
#        plots.extend(plotNumber("N_b_gen",gen_b,noSel,5,"N(b) gen quark"))
#        plots.extend(plotNumber("N_bbar_gen",gen_bbar,noSel,5,"N(#bar{b}) gen quark"))
#
#        #----- Matching -----#
#        lambda_bjet_match = lambda b_gen,b_jet : op.deltaR(b_gen.p4,b_jet.genJet.p4)<0.2
#
#        match_b = op.combine((gen_b,recoJet_b),pred=lambda_bjet_match)
#        match_bbar = op.combine((gen_bbar,recoJet_bbar),pred=lambda_bjet_match)
#
#        plots.extend(plotNumber("N_b_match",match_b,noSel,5,"N(b) match level"))
#        plots.extend(plotNumber("N_bbar_match",match_bbar,noSel,5,"N(#bar{b}) match level"))
#
#        SelBMatch = noSel.refine("SelBMatch",cut=[op.rng_len(match_b)>=1])
#        SelBBarMatch = noSel.refine("SelBBarMatch",cut=[op.rng_len(match_bbar)>=1])
#
#        plots.extend(plotMatching("b",[match_b,match_bbar],[SelBMatch,SelBBarMatch],"b"))

        return plots

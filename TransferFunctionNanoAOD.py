import os.path
import sys
from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, DataDrivenBackgroundHistogramsModule
from bamboo import treefunctions as op    
from bamboo.analysisutils import makeMultiPrimaryDatasetTriggerSelection
from bamboo.plots import Plot, EquidistantBinning, SummedPlot

sys.path.append(os.path.dirname(os.path.abspath(__file__))) # Add scripts in this directory

from BaseHHtobbWW import BaseNanoHHtobbWW

def plotNumber(contName,cont,sel,Nmax,xTitle):
    return [Plot.make1D(contName,
                        op.rng_len(cont),
                        sel,
                        EquidistantBinning(Nmax,0.,Nmax),
                        xTitle=xTitle)]

def plotMatching(contName,cont,sel,partName,bregCorr=False):
    plots = []

    plots.append(Plot.make1D(contName+"_deltaR",
                             op.map(cont,lambda m : op.deltaR(m[0].p4,m[1].p4)), 
                             sel,
                             EquidistantBinning(100,0.,0.2),
                             xTitle="#Delta R(%s_{gen},%s_{reco})"%(partName,partName)))
    plots.append(Plot.make1D(contName+"_deltaPtRel",
                             op.map(cont,lambda m : (m[1].pt-m[0].pt)/m[1].pt),
                             sel,
                             EquidistantBinning(100,-2.,2.),
                             xTitle="P_{T}(%s_{reco})-P_{T}(%s_{gen})/P_{T}(%s_{reco})"%(partName,partName,partName)))
    plots.append(Plot.make2D(contName+"_Pt2D",
                             [op.map(cont,lambda m : m[0].pt),op.map(cont,lambda m : m[1].pt)],
                             sel,
                             [EquidistantBinning(1000,0,1000),EquidistantBinning(1000,0,1000)],
                             xTitle="P_{T}(gen)",
                             yTitle="P_{T}(reco)"))
    plots.append(Plot.make2D(contName+"_TF",
                             (op.map(cont,lambda m : m[0].p4.E()),op.map(cont,lambda m : m[1].p4.E()-m[0].p4.E())),
                             sel,
                             [EquidistantBinning(3000,0.,3000.),EquidistantBinning(5000,-500.,500.)],
                             xTitle="E^{parton}(e^{-})",
                             yTitle="#Delta E = E^{reco}(%s)-E^{parton}(%s)"%(partName,partName)))
    if bregCorr:
        def getCorrBp4(j):
            return  op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (j.pt*j.bRegCorr, j.eta, j.phi, j.mass*j.bRegCorr)) 
        plots.append(Plot.make2D(contName+"_TF_bregCorr",
                                 (op.map(cont,lambda m : m[0].p4.E()),op.map(cont,lambda m : getCorrBp4(m[1]).E()-m[0].p4.E())),
                                 sel,
                                 [EquidistantBinning(3000,0.,3000.),EquidistantBinning(5000,-500.,500.)],
                                 xTitle="E^{parton}(e^{-})",
                                 yTitle="#Delta E = E^{reco}(%s)-E^{parton}(%s)"%(partName,partName)))


    return plots

def plotRecoForGen(sel,gen,reco,lambda_match,name):
    n_reco = [op.select(reco, lambda r: op.switch(op.rng_len(gen)>i, lambda_match(gen[i],r), op.c_bool(False))) for i in range(10)]

    plots = SummedPlot("n_reco_{}_per_gen".format(name),
                       [Plot.make1D("n_reco_{}_per_gen_{}".format(name,i),
                        op.rng_len(nre),
                        sel,
                        EquidistantBinning(10,0,10),
                        xTitle="N reco {} for gen {}".format(name,i))  
                                 for i,nre in enumerate(n_reco)],
                       xTitle = "N reco {} for each gen".format(name))

    return plots 

def fatjetPlots(sel,fatjets,name):
    plots = []
    plots.append(Plot.make2D(name+'_Pt2D',
                             [op.map(fatjets, lambda fatjet : fatjet.subJet1.pt),op.map(fatjets, lambda fatjet : fatjet.subJet2.pt)],
                             sel,
                             [EquidistantBinning(1000,0.,1000.),EquidistantBinning(1000,0.,1000.)],
                             xTitle = "subjet 1 P_{T}",
                             yTitle = "subjet 2 P_{T}"))
    plots.append(Plot.make2D(name+'_E2D',
                             [op.map(fatjets, lambda fatjet : fatjet.subJet1.p4.E()),op.map(fatjets, lambda fatjet : fatjet.subJet2.p4.E())],
                             sel,
                             [EquidistantBinning(1000,0.,1000.),EquidistantBinning(1000,0.,1000.)],
                             xTitle = "subjet 1 E",
                             yTitle = "subjet 2 E"))
    plots.append(Plot.make2D(name+'_subjet1VSFatjet_E2D',
                             [op.map(fatjets, lambda fatjet : fatjet.subJet1.p4.E()),op.map(fatjets, lambda fatjet : fatjet.p4.E())],
                             sel,
                             [EquidistantBinning(1000,0.,1000.),EquidistantBinning(1000,0.,1000.)],
                             xTitle = "subjet 1 E",
                             yTitle = "fatjet E"))
    plots.append(Plot.make2D(name+'_subjet2VSFatjet_E2D',
                             [op.map(fatjets, lambda fatjet : fatjet.subJet2.p4.E()),op.map(fatjets, lambda fatjet : fatjet.p4.E())],
                             sel,
                             [EquidistantBinning(1000,0.,1000.),EquidistantBinning(1000,0.,1000.)],
                             xTitle = "subjet 2 E",
                             yTitle = "fatjet E"))
    plots.append(Plot.make2D(name+'_Eta2D',
                             [op.map(fatjets, lambda fatjet : fatjet.subJet1.p4.eta()),op.map(fatjets, lambda fatjet : fatjet.subJet2.p4.eta())],
                             sel,
                             [EquidistantBinning(500,-2.5,2.5),EquidistantBinning(500,-2.5,2.5)],
                             xTitle = "subjet 1 #eta",
                             yTitle = "subjet 2 #eta"))
    plots.append(Plot.make2D(name+'_Phi2D',
                             [op.map(fatjets, lambda fatjet : fatjet.subJet1.p4.phi()),op.map(fatjets, lambda fatjet : fatjet.subJet2.p4.phi())],
                             sel,
                             [EquidistantBinning(320,-3.2,3.2),EquidistantBinning(320,-3.2,3.2)],
                             xTitle = "subjet 1 #phi",
                             yTitle = "subjet 2 #phi"))

    return plots



class TransferFunction(BaseNanoHHtobbWW,DataDrivenBackgroundHistogramsModule):
    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(TransferFunction,self).prepareObjects(t, noSel, sample, sampleCfg, 'DL')
        noSel = self.beforeJetselection(noSel)
        era = sampleCfg['era']   
        plots = []

        if not self.is_MC:
            return plots

        lambda_is_matched = lambda lep : op.OR(lep.genPartFlav==1,  # Prompt muon or electron
                                               lep.genPartFlav==15) # From tau decay
                                               #lep.genPartFlav==22) # From photon conversion (only available for electrons)

        plots.append(Plot.make1D("totalWeight",
                                 noSel.weight,
                                 noSel,
                                 EquidistantBinning(100,-int(10**5),int(10**5)),
                                 xTitle="total weight"))
        plots.append(Plot.make1D("genWeight",
                                 t.genWeight,
                                 noSel,
                                 EquidistantBinning(100,-int(10**5),int(10**5)),
                                 xTitle="gen weight"))
        
        ##########################################################
        #                       Leptons                          #
        ##########################################################
        #----- Gen particles -----# 
        gen_e  = op.select(t.GenPart,lambda g : op.AND(op.abs(g.pdgId) == 11 , 
                                                       g.statusFlags & ( 0x1 << 13), 
                                                       g.pt >= 10, 
                                                       op.abs(g.eta)<2.5 ))
        gen_m = op.select(t.GenPart,lambda g : op.AND(op.abs(g.pdgId) == 13 , 
                                                      g.statusFlags & ( 0x1 << 13), 
                                                      g.pt >= 10, 
                                                      op.abs(g.eta)<2.4 ))

        #----- Reco particles -----# 
        lambda_reco_e = lambda el : op.AND(op.abs(el.pdgId) == 11 , 
                                           el.pt >= 10 , 
                                           op.abs(el.eta) <2.5,
                                           lambda_is_matched(el))
        lambda_reco_m = lambda mu : op.AND(op.abs(mu.pdgId) == 13 , 
                                           mu.pt >= 10 , 
                                           op.abs(mu.eta) < 2.4, 
                                           lambda_is_matched(mu))

        reco_e = op.select(t.Electron, lambda_reco_e)
        reco_m = op.select(t.Muon, lambda_reco_m)

        fakeable_e = op.select(self.electronsFakeSel, lambda_reco_e)
        fakeable_m = op.select(self.muonsFakeSel, lambda_reco_m)

        tight_e = op.select(self.electronsTightSel, lambda_reco_e)
        tight_m = op.select(self.muonsTightSel, lambda_reco_m)
                                                            
        #---- Matching -----#
        lambda_lepton_match = lambda gen,reco : op.AND(op.deltaR(gen.p4,reco.p4)<0.1)

        match_e = op.combine((gen_e,reco_e), pred=lambda_lepton_match)
        match_m = op.combine((gen_m,reco_m), pred=lambda_lepton_match)

        plots.extend(plotMatching("e",match_e,noSel,"e^{#pm}"))
        plots.extend(plotMatching("m",match_m,noSel,"#mu^{#pm}"))

        plots.append(plotRecoForGen(noSel,gen_e,reco_e,lambda_lepton_match,"e"))
        plots.append(plotRecoForGen(noSel,gen_e,reco_m,lambda_lepton_match,"m"))

        match_fakeable_e = op.combine((gen_e,fakeable_e), pred=lambda_lepton_match)
        match_fakeable_m = op.combine((gen_m,fakeable_m), pred=lambda_lepton_match)

        plots.extend(plotMatching("e_fakeable",match_fakeable_e,noSel,"e^{#pm}"))
        plots.extend(plotMatching("m_fakeable",match_fakeable_m,noSel,"#mu^{#pm}"))

        plots.append(plotRecoForGen(noSel,gen_e,fakeable_e,lambda_lepton_match,"e_fakeable"))
        plots.append(plotRecoForGen(noSel,gen_e,fakeable_m,lambda_lepton_match,"m_fakeable"))

        match_tight_e = op.combine((gen_e,tight_e), pred=lambda_lepton_match)
        match_tight_m = op.combine((gen_m,tight_m), pred=lambda_lepton_match)

        plots.extend(plotMatching("e_tight",match_tight_e,noSel,"e^{#pm}"))
        plots.extend(plotMatching("m_tight",match_tight_m,noSel,"#mu^{#pm}"))

        plots.append(plotRecoForGen(noSel,gen_e,tight_e,lambda_lepton_match,"e_tight"))
        plots.append(plotRecoForGen(noSel,gen_e,tight_m,lambda_lepton_match,"m_tight"))

        ##########################################################
        #                       Ak4 B jets                       #
        ##########################################################
        #----- RecoJets -----#
        recoJet_b = op.select(t.Jet, lambda j : op.AND(op.abs(j.partonFlavour) == 5,
                                                       j.pt>=20,
                                                       op.abs(j.eta)<2.4,
                                                       j.genJet.isValid))
        #----- Parton -----#
        gen_b = op.select(t.GenPart,lambda g : op.AND(op.abs(g.pdgId) == 5 ,
                                                      g.statusFlags & ( 0x1 << 13) ,
                                                      g.pt>=20,
                                                      op.abs(g.eta)<2.4))
        #----- Matching -----#
        lambda_jet_match = lambda gen,reco : op.AND(op.deltaR(gen.p4,reco.p4)<0.2)
        match_b = op.combine((gen_b, recoJet_b), pred=lambda_jet_match)
        plots.append(plotRecoForGen(noSel,gen_b,recoJet_b,lambda_jet_match,"ak4b"))
        plots.extend(plotMatching("ak4b",match_b,noSel,"b",bregCorr=True))

        ##########################################################
        #                       Ak4 C jets                       #
        ##########################################################
        #----- RecoJets -----#
        recoJet_c = op.select(t.Jet, lambda j : op.AND(op.abs(j.partonFlavour) == 4,
                                                       j.pt>=20,
                                                       op.abs(j.eta)<2.4,
                                                       j.genJet.isValid))
        #----- Parton -----#
        gen_c = op.select(t.GenPart,lambda g : op.AND(op.abs(g.pdgId) == 4 ,
                                                      g.statusFlags & ( 0x1 << 13) ,
                                                      g.pt>=20,
                                                      op.abs(g.eta)<2.4))
        #----- Matching -----#
        match_c = op.combine((gen_c, recoJet_c), pred=lambda_jet_match)
        plots.append(plotRecoForGen(noSel,gen_c,recoJet_c,lambda_jet_match,"ak4c"))
        plots.extend(plotMatching("ak4c",match_c,noSel,"c"))

        ##########################################################
        #                     Ak4 lightjets                      #
        ##########################################################
        #----- RecoJets -----#
        recoJet_l = op.select(t.Jet, lambda j : op.AND(op.OR(op.abs(j.partonFlavour) == 1,
                                                             op.abs(j.partonFlavour) == 2, 
                                                             op.abs(j.partonFlavour) == 3),
                                                       j.pt>=20,
                                                       op.abs(j.eta)<2.4,
                                                       j.genJet.isValid))
        #----- Parton -----#
        gen_l = op.select(t.GenPart,lambda g : op.AND(op.OR(op.abs(g.pdgId) == 1,
                                                            op.abs(g.pdgId) == 2,
                                                            op.abs(g.pdgId) == 3),
                                                      g.statusFlags & ( 0x1 << 13) ,
                                                      g.pt>=20,
                                                      op.abs(g.eta)<2.4))
        #----- Matching -----#
        match_l = op.combine((gen_l, recoJet_l), pred=lambda_jet_match)
        plots.append(plotRecoForGen(noSel,gen_l,recoJet_l,lambda_jet_match,"ak4l"))
        plots.extend(plotMatching("ak4l",match_l,noSel,"l"))

        ##########################################################
        #                        Ak8 jets                        #
        ##########################################################
        #----- RecoJets -----#
        recoFatJet = op.select(t.FatJet, lambda j : op.AND(j.pt>=100,
                                                           op.abs(j.eta)<2.4,
                                                           op.AND(j.subJet1.isValid, j.subJet1.pt >= 20. , op.abs(j.subJet1.eta)<=2.4,
                                                                  j.subJet2.isValid, j.subJet2.pt >= 20. , op.abs(j.subJet2.eta)<=2.4),
                                                           op.AND(j.msoftdrop >= 30, j.msoftdrop <= 210),
                                                           j.tau2/j.tau1 <= 0.75))
        #recoFatJet_b = op.select(recoFatJet, lambda j : op.abs(j.partonFlavour) == 5)
        #recoFatJet_c = op.select(recoFatJet, lambda j : op.abs(j.partonFlavour) == 4)
        #recoFatJet_l = op.select(recoFatJet, lambda j : op.OR(op.abs(j.partonFlavour) == 1,
        #                                                   op.abs(j.partonFlavour) == 2, 
        #                                                   op.abs(j.partonFlavour) == 3))

        plots.extend(fatjetPlots(noSel,recoFatJet,"fatjet"))

        return plots

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

def plotMatching(contName,cont,sel,partName,bregCorr=False):
    plots = []

    plots.append(Plot.make1D(contName+"_deltaR",
                          op.map(cont,lambda m : op.deltaR(m[0].p4,m[1].p4)), 
                          sel,
                          EquidistantBinning(100,0.,0.2),
                          xTitle="#Delta R(%s_{gen},%s_{reco})"%(partName,partName)))
    plots.append(Plot.make1D(contName+"_deltaR2",
                         op.map(cont,lambda m : op.pow(op.deltaR(m[0].p4,m[1].p4),2)), 
                         sel,
                         EquidistantBinning(100,0.,0.04),
                         xTitle="#Delta R^{2}(%s_{gen},%s_{reco})"%(partName,partName)))
    plots.append(Plot.make1D(contName+"_deltaPtRel",
                       op.map(cont,lambda m : (m[1].pt-m[0].pt)/m[1].pt),
                       sel,
                       EquidistantBinning(100,-2.,2.),
                       xTitle="P_{T}(%s_{reco})-P_{T}(%s_{gen})/P_{T}(%s_{reco})"%(partName,partName,partName)))
    plots.append(Plot.make2D(contName+"_TF",
                          (op.map(cont,lambda m : m[0].p4.E()),op.map(cont,lambda m : m[1].p4.E()-m[0].p4.E())),
                          sel,
                          [EquidistantBinning(1000,0.,1000.),EquidistantBinning(1000,-500.,500.)],
                          xTitle="E^{parton}(e^{-})",
                          yTitle="#Delta E = E^{reco}(%s)-E^{parton}(%s)"%(partName,partName)))
    if bregCorr:
        def getCorrBp4(j):
            return  op.construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (j.pt*j.bRegCorr, j.eta, j.phi, j.mass*j.bRegCorr)) 
        plots.append(Plot.make2D(contName+"_TF_bregCorr",
                              (op.map(cont,lambda m : m[0].p4.E()),op.map(cont,lambda m : getCorrBp4(m[1]).E()-m[0].p4.E())),
                              sel,
                              [EquidistantBinning(1000,0.,1000.),EquidistantBinning(1000,-500.,500.)],
                              xTitle="E^{parton}(e^{-})",
                              yTitle="#Delta E = E^{reco}(%s)-E^{parton}(%s)"%(partName,partName)))


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
        reco_e_minus = op.select(t.Electron, lambda ele : op.AND(ele.charge == -1 ,
                                                                 ele.pdgId == 11 , 
                                                                 ele.pt>10 , 
                                                                 op.abs(ele.eta)<2.5,
                                                                 lambda_is_matched(ele)))
        reco_e_plus = op.select(t.Electron, lambda ele : op.AND(ele.charge == +1 , 
                                                                ele.pdgId == -11 , 
                                                                ele.pt>10 , 
                                                                op.abs(ele.eta)<2.5,
                                                                lambda_is_matched(ele)))
        reco_mu_minus = op.select(t.Muon, lambda mu : op.AND(mu.charge == -1 , 
                                                             mu.pdgId == 13 , 
                                                             mu.pt>10 , 
                                                             op.abs(mu.eta)<2.4, 
                                                             lambda_is_matched(mu)))
        reco_mu_plus = op.select(t.Muon, lambda mu : op.AND(mu.charge == +1, 
                                                            mu.pdgId == -13,
                                                            mu.pt>10, 
                                                            op.abs(mu.eta)<2.4,
                                                            lambda_is_matched(mu)))
                                                            

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

        has_match_e_minus   = noSel.refine("has_match_e_minus",cut=[op.rng_len(match_e_minus)>0])
        has_match_e_plus    = noSel.refine("has_match_e_plus",cut=[op.rng_len(match_e_plus)>0])
        has_match_mu_minus  = noSel.refine("has_match_mu_minus",cut=[op.rng_len(match_mu_minus)>0])
        has_match_mu_plus   = noSel.refine("has_match_mu_plus",cut=[op.rng_len(match_mu_plus)>0])
       
        plots.extend(plotNumber("N_e_minus_match",match_e_minus,noSel,10,"N(e^{-}) match level"))
        plots.extend(plotNumber("N_e_plus_match",match_e_plus,noSel,10,"N(e^{+}) match level"))
        plots.extend(plotNumber("N_mu_minus_match",match_mu_minus,noSel,10,"N(\mu^{-}) match level"))
        plots.extend(plotNumber("N_mu_plus_match",match_mu_plus,noSel,10,"N(\mu^{+}) match level"))

        plots.extend(plotMatching("e_minus",match_e_minus,has_match_e_minus,"e^{-}"))
        plots.extend(plotMatching("e_plus",match_e_plus,has_match_e_plus,"e^{+}"))
        plots.extend(plotMatching("mu_minus",match_mu_minus,has_match_mu_minus,"#mu^{-}"))
        plots.extend(plotMatching("mu_plus",match_mu_plus,has_match_mu_plus,"#mu^{+}"))


        ##########################################################
        #                       B jets                           #
        ##########################################################
        #----- RecoJets -----#
        recoJet_b = op.select(t.Jet, lambda j : op.AND(op.abs(j.partonFlavour) == 5,
                                                       j.pt>10,
                                                       op.abs(j.eta)<2.4,
                                                       j.genJet.isValid))
        #----- Parton -----#
        gen_b = op.select(t.GenPart,lambda g : op.AND(op.abs(g.pdgId) == 5 ,
                                                      g.statusFlags & ( 0x1 << 13) ,
                                                      g.pt>10,
                                                      op.abs(g.eta)<2.4))
        #----- Matching -----#
        match_b = op.combine((gen_b, recoJet_b), pred=(lambda bg,br : op.deltaR(bg.p4, br.p4) < 0.2))

        plots.extend(plotNumber("N_jetb_reco",recoJet_b,noSel,5,"N(b) reco jet"))
        plots.extend(plotNumber("N_b_gen",gen_b,noSel,5,"N(b) gen quark"))
        plots.extend(plotNumber("N_b_match",match_b,noSel,5,"N(b) match level"))
        plots.extend(plotMatching("b",match_b,noSel,"b",bregCorr=True))

        return plots

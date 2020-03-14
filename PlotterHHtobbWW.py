import os
import sys

from bamboo.analysismodules import HistogramsModule
from bamboo import treefunctions as op
from bamboo.plots import Plot, EquidistantBinning, SummedPlot

sys.path.append('/home/ucl/cp3/fbury/bamboodev/HHbbWWAnalysis/') # Add scripts in this directory -- TODO : make cleaner
from BaseHHtobbWW import BaseNanoHHtobbWW
from plotDef import *
from scalefactorsbbWW import ScaleFactorsbbWW

#===============================================================================================#
#                                 PlotterHHtobbWW                                               #
#===============================================================================================#
class PlotterNanoHHtobbWW(BaseNanoHHtobbWW,HistogramsModule):
    """ Plotter module: HH->bbW(->e/µ nu)W(->e/µ nu) histograms from NanoAOD """
    def __init__(self, args):
        super(PlotterNanoHHtobbWW, self).__init__(args)

    def definePlots(self, t, noSel, sample=None, sampleCfg=None): 
        noSel = super(PlotterNanoHHtobbWW,self).prepareObjects(t, noSel, sample, sampleCfg)

        plots = []

        era = sampleCfg['era']

        isMC = self.isMC(sample)

        yieldPlots = makeYieldPlots()

        #############################################################################
        #                             Scalefactors                                  #
        #############################################################################
        if isMC:
            # Initialize scalefactors class #
            SF = ScaleFactorsbbWW()    
            ####  Muons ####
           # muTightIDSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, self.sfTag), "id_tight"), combine="weight", systName="muid")
           # muTightISOSF = SF.get_scalefactor("lepton", ("muon_{0}_{1}".format(era, self.sfTag), "iso_tight_id_medium"), combine="weight", systName="muiso")
            muLooseId = SF.get_scalefactor("lepton", 'muon_loose_{}'.format(era), combine="weight", systName="mu_loose")
            muTightMVA = SF.get_scalefactor("lepton", 'muon_tightMVA_{}'.format(era), combine="weight", systName="mu_tightmva")
            ####  Electrons ####
           # elMVA80noisoSF = SF.get_scalefactor("lepton", ("electron_{0}_{1}".format(era,self.sfTag), "id_mva80noiso"), systName="elid") # Only applied on fakeable but not tight lepton
            elLooseRecoPtLt20 = SF.get_scalefactor("lepton", ('electron_loosereco_{}'.format(era) , 'electron_loosereco_ptgt20'), combine="weight", systName="el_looserecoptlt20")
            elLooseRecoPtGt20 = SF.get_scalefactor("lepton", ('electron_loosereco_{}'.format(era) , 'electron_loosereco_ptlt20'), combine="weight", systName="el_looserecoptgt20")
            elLooseId = SF.get_scalefactor("lepton", 'electron_looseid_{}'.format(era) , combine="weight", systName="el_looseid")
            elLooseEff = SF.get_scalefactor("lepton", 'electron_looseeff_{}'.format(era) , combine="weight", systName="el_looseeff")
            elTightMVA = SF.get_scalefactor("lepton", 'electron_tightMVA_{}'.format(era) , combine="weight", systName="el_tightmva")
            #### Triggers (from TTH) ####
            ttH_doubleMuon_trigSF = op.systematic(op.c_float(1.010), name="ttH_doubleMuon_trigSF", up=op.c_float(1.020), down=op.c_float(1.000))
            ttH_doubleElectron_trigSF = op.systematic(op.c_float(1.020), name="ttH_doubleElectron_trigSF", up=op.c_float(1.040), down=op.c_float(1.000))
            ttH_electronMuon_trigSF = op.systematic(op.c_float(1.020), name="ttH_electronMuon_trigSF", up=op.c_float(1.030), down=op.c_float(1.010))
            #### Ak4 Btagging ####
            DeepJetTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepFlavB}
            DeepJetMediumSF = SF.get_scalefactor("jet", ("btag_"+era+"_"+self.sfTag, "DeepJet_medium"), additionalVariables=DeepJetTag_discriVar, systName="deepjet") # For RESOLVED
            #### Ak8 Btagging ####
            DeepCSVTag_discriVar = {"BTagDiscri": lambda j : j.btagDeepB}
            DeepCSVMediumSF = SF.get_scalefactor("jet", ("subjet_btag_"+era+"_"+self.sfTag, "DeepCSV_medium"), additionalVariables=DeepCSVTag_discriVar, systName="deepcsv") # For BOOSTED (btag on subjet)

        ##############################################################################
        #                               Dilepton                                    #
        #############################################################################
        # Make dilepton fakeable pairs wrt to channel #
        lambdaOS  = lambda l1,l2 : l1.charge != l2.charge           # Opposite sign charges
        lambdaPTCut = lambda l1,l2: op.OR( l1.pt > 25 , l2.pt > 25) # PT cut 
        lambdaFakeableDilepton = lambda l1, l2 : op.AND( lambdaOS(l1,l2) , lambdaPTCut(l1,l2) )
        OsElElDilepton = op.combine(self.electronsFakeSel, N=2, pred=lambdaFakeableDilepton)
        OsMuMuDilepton = op.combine(self.muonsFakeSel, N=2, pred=lambdaFakeableDilepton)
        OsElMuDilepton = op.combine((self.electronsFakeSel, self.muonsFakeSel), pred=lambdaFakeableDilepton)

        # Applied SF #
        MuonSF = lambda mu : [muLooseId(mu)] 
        ElectronSF = lambda el : [elLooseId(el) , elLooseEff(el), op.switch(el.pt>20 , elLooseRecoPtGt20(el) , elLooseRecoPtLt20(el))] 

        ElElLooseSF = [ttH_doubleElectron_trigSF]+ElectronSF(OsElElDilepton[0][0])+ElectronSF(OsElElDilepton[0][1]) if isMC else None
        MuMuLooseSF = [ttH_doubleMuon_trigSF]+MuonSF(OsMuMuDilepton[0][0])+MuonSF(OsMuMuDilepton[0][1]) if isMC else None
        ElMuLooseSF = [ttH_electronMuon_trigSF]+ElectronSF(OsElMuDilepton[0][0])+MuonSF(OsElMuDilepton[0][1]) if isMC else None

        ElElTightSF = [elTightMVA(OsElElDilepton[0][0]),elTightMVA(OsElElDilepton[0][1])] if isMC else None
        MuMuTightSF = [muTightMVA(OsMuMuDilepton[0][0]),muTightMVA(OsMuMuDilepton[0][1])] if isMC else None
        ElMuTightSF = [elTightMVA(OsElMuDilepton[0][0]),muTightMVA(OsElMuDilepton[0][1])] if isMC else None

        # At least 2 fakeables leptons #
        ElElHasFakeableDilepton = noSel.refine("ElElHasFakeableDilepton",
                                        cut = [op.rng_len(OsElElDilepton)>0],
                                        weight = ElElLooseSF)
        MuMuHasFakeableDilepton = noSel.refine("MuMuHasFakeableDilepton",
                                        cut = [op.rng_len(OsMuMuDilepton)>0],
                                        weight = MuMuLooseSF)
        ElMuHasFakeableDilepton = noSel.refine("ElMuHasFakeableDilepton",
                                        cut = [op.rng_len(OsElMuDilepton)>0],
                                        weight = ElMuLooseSF)

        # Yield #
        plots.append(yieldPlots.addYield(ElElHasFakeableDilepton,"ElElHasFakeableDilepton","OS fakeable leptons (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDilepton,"MuMuHasFakeableDilepton","OS fakeable leptons (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDilepton,"ElMuHasFakeableDilepton","OS fakeable leptons (channel $e^{\pm}\mu^{\mp}$)"))
    
        # Mll cut + Z peak exclusion (charge already done in previous selection) in **preselected** leptons # 
        OsElElPresel = op.combine(self.electronsPreSel, N=2, pred=lambdaOS)
        OsMuMuPresel = op.combine(self.muonsPreSel, N=2, pred=lambdaOS)
        OsElMuPresel = op.combine((self.electronsPreSel, self.muonsPreSel), pred=lambdaOS)
        lambda_lowMllCut    = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.invariant_mass(dilep[0].p4, dilep[1].p4)<12.))
            # If any dilepton preselected below 12GeV, returns False
        lambda_outZ         = lambda dileptons: op.NOT(op.rng_any(dileptons, lambda dilep : op.in_range(80.,op.invariant_mass(dilep[0].p4, dilep[1].p4),100.)))
            # If any dilepton preselected within Z peak, returns False

        ElElHasFakeableDileptonPreMllCutOutZ = ElElHasFakeableDilepton.refine("ElElHasFakeableDileptonPreMllCutOutZ", 
                                cut=[lambda_lowMllCut(OsElElPresel), lambda_outZ(OsElElPresel)])
        MuMuHasFakeableDileptonPreMllCutOutZ = MuMuHasFakeableDilepton.refine("MuMuHasFakeableDileptonPreMllCutOutZ", 
                                cut=[lambda_lowMllCut(OsMuMuPresel),lambda_outZ(OsMuMuPresel)])
        ElMuHasFakeableDileptonPreMllCut     = ElMuHasFakeableDilepton.refine("ElMuHasFakeableDileptonPreMllCut", 
                                cut=[lambda_lowMllCut(OsElMuPresel)])

        # Yield #
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonPreMllCutOutZ,"ElElHasFakeableDileptonPreMllCutOutZ","OS fakeable leptons + $M_{ll}$ cut (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonPreMllCutOutZ,"MuMuHasFakeableDileptonPreMllCutOutZ","OS fakeable leptons + $M_{ll}$ cut (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonPreMllCut,"ElMuHasFakeableDileptonPreMllCut","OS fakeable leptons + $M_{ll}$ cut (channel $e^{\pm}\mu^{\mp}$)"))

        # Tight cut veto : no more than two tight leptons in the event #
        ElElHasFakeableDileptonPreMllCutOutZTightVeto = ElElHasFakeableDileptonPreMllCutOutZ.refine("ElElHasFakeableDileptonPreMllCutOutZTightVeto", 
                                cut=[op.rng_len(self.muonsTightSel)+op.rng_len(self.electronsTightSel) <= 2])
        MuMuHasFakeableDileptonPreMllCutOutZTightVeto = MuMuHasFakeableDileptonPreMllCutOutZ.refine("MuMuHasFakeableDileptonPreMllCutOutZTightVeto", 
                                cut=[op.rng_len(self.muonsTightSel)+op.rng_len(self.electronsTightSel) <= 2])
        ElMuHasFakeableDileptonPreMllCutTightVeto     = ElMuHasFakeableDileptonPreMllCut.refine("ElMuHasFakeableDileptonPreMllCutTightVeto", 
                                cut=[op.rng_len(self.muonsTightSel)+op.rng_len(self.electronsTightSel) <= 2])

        # Yield #
        plots.append(yieldPlots.addYield(ElElHasFakeableDileptonPreMllCutOutZTightVeto,"ElElHasFakeableDileptonPreMllCutOutZTightVeto","OS fakeable leptons + $M_{ll}$ cut + Two tights veto (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasFakeableDileptonPreMllCutOutZTightVeto,"MuMuHasFakeableDileptonPreMllCutOutZTightVeto","OS fakeable leptons + $M_{ll}$ cut + Two tights veto (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasFakeableDileptonPreMllCutTightVeto,"ElMuHasFakeableDileptonPreMllCutTightVeto","OS fakeable leptons + $M_{ll}$ cut + Two tights veto (channel $e^{\pm}\mu^{\mp}$)"))

        # Tight selected OS pairs #
        if isMC:
            lambda_is_matched = lambda dilep : op.AND(dilep[0].genPartFlav==1,dilep[1].genPartFlav==1)
        else:
            lambda_is_matched = lambda dilep : op.c_bool(True)



        ElElHasTightDileptonPreMllCutOutZTightVeto = ElElHasFakeableDileptonPreMllCutOutZTightVeto.refine("ElElHasTightDileptonPreMllCutOutZTightVeto",
                                            cut=[self.lambda_electronTightSel(OsElElDilepton[0][0]),
                                                 self.lambda_electronTightSel(OsElElDilepton[0][1]),
                                                 lambda_is_matched(OsElElDilepton[0])],
                                            weight=ElElTightSF)
        MuMuHasTightDileptonPreMllCutOutZTightVeto = MuMuHasFakeableDileptonPreMllCutOutZTightVeto.refine("MuMuHasTightDileptonPreMllCutOutZTightVeto",
                                            cut=[self.lambda_muonTightSel(OsMuMuDilepton[0][0]),
                                                 self.lambda_muonTightSel(OsMuMuDilepton[0][1]),
                                                 lambda_is_matched(OsMuMuDilepton[0])],
                                            weight=MuMuTightSF)
        ElMuHasTightDileptonPreMllCutTightVeto = ElMuHasFakeableDileptonPreMllCutTightVeto.refine("ElMuHasTightDileptonPreMllCutTightVeto",
                                            cut=[self.lambda_electronTightSel(OsElMuDilepton[0][0]),
                                                 self.lambda_muonTightSel(OsElMuDilepton[0][1]),
                                                 lambda_is_matched(OsElMuDilepton[0])],
                                            weight=ElMuTightSF)

        # Yield #
        plots.append(yieldPlots.addYield(ElElHasTightDileptonPreMllCutOutZTightVeto,"ElElHasTightDileptonPreMllCutOutZTightVeto","OS tight leptons + $M_{ll}$ cut + Two tights veto (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasTightDileptonPreMllCutOutZTightVeto,"MuMuHasTightDileptonPreMllCutOutZTightVeto","OS tight leptons + $M_{ll}$ cut + Two tights veto (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasTightDileptonPreMllCutTightVeto,"ElMuHasTightDileptonPreMllCutTightVeto","OS tight leptons + $M_{ll}$ cut + Two tights veto (channel $e^{\pm}\mu^{\mp}$)"))


        # Trigger Plots #
        plots.extend(triggerPlots(sel         = ElElHasTightDileptonPreMllCutOutZTightVeto,
                                  triggerDict = self.triggersPerPrimaryDataset,
                                  suffix      = "ElElHasTightDileptonPreMllCutOutZTightVeto",
                                  channel     = "ElEl"))
        plots.extend(triggerPlots(sel         = MuMuHasTightDileptonPreMllCutOutZTightVeto,
                                  triggerDict = self.triggersPerPrimaryDataset,
                                  suffix      = "MuMuHasTightDileptonPreMllCutOutZTightVeto",
                                  channel     = "MuMu"))
        plots.extend(triggerPlots(sel         = ElMuHasTightDileptonPreMllCutTightVeto,
                                  triggerDict = self.triggersPerPrimaryDataset,
                                  suffix      = "ElMuHasTightDileptonPreMllCutTightVeto",
                                  channel     = "ElMu"))

        # Dilepton channel plots #
        plots.extend(channelPlot(sel        = ElElHasTightDileptonPreMllCutOutZTightVeto,
                                 DilepElEl  = OsElElDilepton,
                                 DilepMuMu  = OsMuMuDilepton,
                                 DilepElMu  = OsElMuDilepton,
                                 suffix     = "ElElHasTightDileptonPreMllCutOutZTightVeto",
                                 channel    = "ElEl"))
        plots.extend(channelPlot(sel        = MuMuHasTightDileptonPreMllCutOutZTightVeto,
                                 DilepElEl  = OsElElDilepton,
                                 DilepMuMu  = OsMuMuDilepton,
                                 DilepElMu  = OsElMuDilepton,
                                 suffix     = "MuMuHasTightDileptonPreMllCutOutZTightVeto",
                                 channel    = "MuMu"))
        plots.extend(channelPlot(sel        = ElMuHasTightDileptonPreMllCutTightVeto,
                                 DilepElEl  = OsElElDilepton,
                                 DilepMuMu  = OsMuMuDilepton,
                                 DilepElMu  = OsElMuDilepton,
                                 suffix     = "ElMuHasTightDileptonPreMllCutTightVeto",
                                 channel    = "ElMu"))


        #############################################################################
        #                              AK4 Jets                                     #
        #############################################################################
        # Jet Selection #
        ElElHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets = ElElHasTightDileptonPreMllCutOutZTightVeto.refine("ElElHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets", 
                                cut=[op.rng_len(self.ak4Jets)>=2])
        MuMuHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets = MuMuHasTightDileptonPreMllCutOutZTightVeto.refine("MuMuHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets", 
                                cut=[op.rng_len(self.ak4Jets)>=2])
        ElMuHasTightDileptonPreMllCutTightVetoTwoAk4Jets     = ElMuHasTightDileptonPreMllCutTightVeto.refine("ElMuHasTightDileptonPreMllCutTightVetoTwoAk4Jets", 
                                cut=[op.rng_len(self.ak4Jets)>=2])

        # Yield 
        plots.append(yieldPlots.addYield(ElElHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets,
                                         "ElElHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets",
                                         "OS leptons + Ak4 Jets $\geq 2$ (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets,
                                         "MuMuHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets",
                                         "OS leptons + Ak4 Jets $\geq 2$ (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasTightDileptonPreMllCutTightVetoTwoAk4Jets,
                                         "ElMuHasTightDileptonPreMllCutTightVetoTwoAk4Jets",
                                         "OS leptons + Ak4 Jets $\geq 2$ (channel $e^{\pm}\mu^{\mp}$)"))

        # Plots #
        HasTightDileptonTwoAk4JetsChannelList = [
                             {'channel':'ElEl','sel':ElElHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets,'dilepton':OsElElDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets'},
                             {'channel':'MuMu','sel':MuMuHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets,'dilepton':OsMuMuDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets'},
                             {'channel':'ElMu','sel':ElMuHasTightDileptonPreMllCutTightVetoTwoAk4Jets,    'dilepton':OsElMuDilepton[0],'suffix':'HasTightDileptonPreMllCutTightVetoTwoAk4Jets'},
                                                  ]
        commonItems = ['channel','sel','suffix']
        HasTightDileptonTwoAk4JetsBaseDictJetsPlots = {'leadjet':self.ak4Jets[0],'subleadjet':self.ak4Jets[1],'lead_is_b':False,'sublead_is_b':False}
        HasTightDileptonTwoAk4JetsBaseDictJetsN = {'objName':'Ak4Jets','objCont':self.ak4Jets,'Nmax':10,'xTitle':'N(Ak4 jets)'}
                                                    
        for channelDict in HasTightDileptonTwoAk4JetsChannelList:
            # Dilepton #
            plots.extend(makeDileptonPlots(**channelDict,isMC=isMC))
            # Ak4 Jets #
            plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**HasTightDileptonTwoAk4JetsBaseDictJetsN))
            plots.extend(makeAk4JetsPlots(**{k:channelDict[k] for k in commonItems},**HasTightDileptonTwoAk4JetsBaseDictJetsPlots))
            # MET #
            plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))

        # Check plots #
        ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets = ElElHasFakeableDileptonPreMllCutOutZTightVeto.refine("ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets",
                                cut=[op.rng_len(self.ak4Jets)>=2])
        MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets = MuMuHasFakeableDileptonPreMllCutOutZTightVeto.refine("MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets",
                                cut=[op.rng_len(self.ak4Jets)>=2])
        ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets = ElMuHasFakeableDileptonPreMllCutTightVeto.refine("ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets",
                                cut=[op.rng_len(self.ak4Jets)>=2])

        # Dilepton channel plots #
        plots.extend(makeDileptonCheckPlots(self        = self,
                                            sel         = ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets,
                                            suffix      = "ElElHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets",
                                            dilepton    = OsElElDilepton[0],
                                            channel     = 'ElEl'))
        plots.extend(makeDileptonCheckPlots(self        = self,
                                            sel         = MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets,
                                            suffix      = "MuMuHasFakeableDileptonPreMllCutOutZTightVetoTwoAk4Jets",
                                            dilepton    = OsMuMuDilepton[0],
                                            channel     = 'MuMu'))
        plots.extend(makeDileptonCheckPlots(self        = self,
                                            sel         = ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets,
                                            suffix      = "ElMuHasFakeableDileptonPreMllCutTightVetoTwoAk4Jets",
                                            dilepton    = OsElMuDilepton[0],
                                            channel     = 'ElMu'))

        #############################################################################
        #                              Ak4 Btagging                                 #
        #############################################################################
        # Resolved selection #
        # Inclusive resolved #
#        ElElHasTightDileptonPreMllCutOutZTightVetoInclusiveResolved = ElElHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
#                                "ElElHasTightDileptonPreMllCutOutZTightVetoInclusiveResolved", 
#                                cut=[op.rng_len(self.ak4BJets) >= 1])
#        MuMuHasTightDileptonPreMllCutOutZTightVetoInclusiveResolved = MuMuHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
#                                "MuMuHasTightDileptonPreMllCutOutZTightVetoInclusiveResolved", 
#                                cut=[op.rng_len(self.ak4BJets) >= 1])
#        ElMuHasTightDileptonPreMllCutTightVetoInclusiveResolved = ElMuHasTightDileptonPreMllCutTightVetoTwoAk4Jets.refine("
#                                ElMuHasTightDileptonPreMllCutTightVetoInclusiveResolved", 
#                                cut=[op.rng_len(self.ak4BJets) >= 1])
        # Applied SF #
        Ak4BtaggingOneBtag  = [DeepJetMediumSF(self.ak4BJets[0])] if isMC else None
        Ak4BtaggingTwoBtags = [DeepJetMediumSF(self.ak4BJets[0]),DeepJetMediumSF(self.ak4BJets[1])] if isMC else None

        # Exclusive Resolved No Btag #
        ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag = ElElHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
                                "ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag", 
                                cut=[op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0])
        MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag = MuMuHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
                                "MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag", 
                                cut=[op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0])
        ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedNoBtag = ElMuHasTightDileptonPreMllCutTightVetoTwoAk4Jets.refine(
                                "ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedNoBtag", 
                                cut=[op.rng_len(self.ak4BJets)==0,op.rng_len(self.ak8BJets)==0])
        # Exclusive resolved 1 Btag #
        ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag = ElElHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
                                "ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag", 
                                cut=[op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0],
                                weight=Ak4BtaggingOneBtag)
        MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag = MuMuHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
                                "MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag", 
                                cut=[op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0],
                                weight=Ak4BtaggingOneBtag)
        ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedOneBtag = ElMuHasTightDileptonPreMllCutTightVetoTwoAk4Jets.refine(
                                "ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedOneBtag", 
                                cut=[op.rng_len(self.ak4BJets)==1,op.rng_len(self.ak8BJets)==0],
                                weight=Ak4BtaggingOneBtag)
        # Exclusive resolved 2 Btags #
        ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags = ElElHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
                                "ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags", 
                                cut=[op.rng_len(self.ak4BJets)==2,op.rng_len(self.ak8BJets)==0],
                                weight=Ak4BtaggingTwoBtags)
        MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags = MuMuHasTightDileptonPreMllCutOutZTightVetoTwoAk4Jets.refine(
                                "MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags", 
                                cut=[op.rng_len(self.ak4BJets)==2,op.rng_len(self.ak8BJets)==0],
                                weight=Ak4BtaggingTwoBtags)
        ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags = ElMuHasTightDileptonPreMllCutTightVetoTwoAk4Jets.refine(
                                "ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags", 
                                cut=[op.rng_len(self.ak4BJets)==2,op.rng_len(self.ak8BJets)==0],
                                weight=Ak4BtaggingTwoBtags)

                                
        # Yield 
#        plots.append(yieldPlots.addYield(ElElHasTightDileptonPreMllCutOutZTightVetoInclusiveResolved,
#                                         "ElElHasTightDileptonPreMllCutOutZTightVetoInclusiveResolved",
#                                         "OS leptons + Inclusive Resolved Jets (channel $e^+e^-$)"))
#        plots.append(yieldPlots.addYield(MuMuHasTightDileptonPreMllCutOutZTightVetoInclusiveResolved,
#                                         "MuMuHasTightDileptonPreMllCutOutZTightVetoInclusiveResolved",
#                                         "OS leptons + Inclusive Resolved Jets (channel $\mu^+\mu^-$)"))
#        plots.append(yieldPlots.addYield(ElMuHasTightDileptonPreMllCutTightVetoInclusiveResolved,
#                                         "ElMuHasTightDileptonPreMllCutTightVetoInclusiveResolved",
#                                         "OS leptons + Inclusive Resolved Jets (channel $e^{\pm}\mu^{\mp}$)"))

        plots.append(yieldPlots.addYield(ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,
                                         "ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag",
                                         "OS leptons + Exclusive Resolved Jets (0 bjet) (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,
                                         "MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag",
                                         "OS leptons + Exclusive Resolved Jets (0 bjet) (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedNoBtag,
                                         "ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedNoBtag",
                                         "OS leptons + Exclusive Resolved Jets (0 bjet) (channel $e^{\pm}\mu^{\mp}$)"))

        plots.append(yieldPlots.addYield(ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,
                                         "ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag",
                                         "OS leptons + Exclusive Resolved Jets (1 bjet) (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,
                                         "MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag",
                                         "OS leptons + Exclusive Resolved Jets (1 bjet) (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedOneBtag,
                                         "ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedOneBtag",
                                         "OS leptons + Exclusive Resolved Jets (1 bjet) (channel $e^{\pm}\mu^{\mp}$)"))

        plots.append(yieldPlots.addYield(ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags,
                                         "ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags",
                                         "OS leptons + Exclusive Resolved Jets (2 bjets) (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags,
                                         "MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags",
                                         "OS leptons + Exclusive Resolved Jets (2 bjets) (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags,
                                         "ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags",
                                         "OS leptons + Exclusive Resolved Jets (2 bjets) (channel $e^{\pm}\mu^{\mp}$)"))

        # Dilepton plots #
        HasTightDileptonResolvedJetsChannelListLeptons = [
#             {'channel':'ElEl','sel':ElElHasTightDileptonPreMllCutOutZTightVetoInclusiveResolved,'dilepton':OsElElDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoInclusiveResolved'},
#             {'channel':'MuMu','sel':MuMuHasTightDileptonPreMllCutOutZTightVetoInclusiveResolved,'dilepton':OsMuMuDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoInclusiveResolved'},
#             {'channel':'ElMu','sel':ElMuHasTightDileptonPreMllCutTightVetoInclusiveResolved,    'dilepton':OsElMuDilepton[0],'suffix':'HasTightDileptonPreMllCutTightVetoInclusiveResolved'},
             {'channel':'ElEl','sel':ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,'dilepton':OsElElDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag'},
             {'channel':'MuMu','sel':MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,'dilepton':OsMuMuDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag'},
             {'channel':'ElMu','sel':ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedNoBtag,    'dilepton':OsElMuDilepton[0],'suffix':'HasTightDileptonPreMllCutTightVetoExclusiveResolvedNoBtag'},
             {'channel':'ElEl','sel':ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,'dilepton':OsElElDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag'},
             {'channel':'MuMu','sel':MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,'dilepton':OsMuMuDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag'},
             {'channel':'ElMu','sel':ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedOneBtag,    'dilepton':OsElMuDilepton[0],'suffix':'HasTightDileptonPreMllCutTightVetoExclusiveResolvedOneBtag'},
             {'channel':'ElEl','sel':ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags,'dilepton':OsElElDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags'},
             {'channel':'MuMu','sel':MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags,'dilepton':OsMuMuDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags'},
             {'channel':'ElMu','sel':ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags,    'dilepton':OsElMuDilepton[0],'suffix':'HasTightDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags'},
                                                            ]

        commonItems = ['channel','sel','suffix']
        HasTightDileptonResolvedJetsBaseDictJetsN = {'objName':'Ak4BJets','objCont':self.ak4BJets,'Nmax':10,'xTitle':'N(Ak4 Bjets)'}

        for channelDict in HasTightDileptonResolvedJetsChannelListLeptons:
            # Dilepton #
            plots.extend(makeDileptonPlots(**channelDict,isMC=isMC))
            # Number of jets #
            plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**HasTightDileptonResolvedJetsBaseDictJetsN))
            # MET #
            plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))
 

        # Dijet plots #
        HasTightDileptonResolvedJetsChannelListJets = [
            {'channel':'ElEl','sel':ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,'suffix':'ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag','leadjet':self.ak4LightJets[0],'subleadjet':self.ak4LightJets[1],'lead_is_b':False,'sublead_is_b':False},
            {'channel':'MuMu','sel':MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,'suffix':'MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag','leadjet':self.ak4LightJets[0],'subleadjet':self.ak4LightJets[1],'lead_is_b':False,'sublead_is_b':False},
            {'channel':'ElMu','sel':ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedNoBtag,'suffix':'ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedNoBtag','leadjet':self.ak4LightJets[0],'subleadjet':self.ak4LightJets[1],'lead_is_b':False,'sublead_is_b':False},

            {'channel':'ElEl','sel':ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,'suffix':'ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag','leadjet':self.ak4BJets[0],'subleadjet':self.ak4LightJets[0],'lead_is_b':True,'sublead_is_b':False},
            {'channel':'MuMu','sel':MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,'suffix':'MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag','leadjet':self.ak4BJets[0],'subleadjet':self.ak4LightJets[0],'lead_is_b':True,'sublead_is_b':False},
            {'channel':'ElMu','sel':ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedOneBtag,'suffix':'ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedOneBtag','leadjet':self.ak4BJets[0],'subleadjet':self.ak4LightJets[0],'lead_is_b':True,'sublead_is_b':False},

            {'channel':'ElEl','sel':ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags,'suffix':'ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags','leadjet':self.ak4BJets[0],'subleadjet':self.ak4BJets[1],'lead_is_b':True,'sublead_is_b':True},
            {'channel':'MuMu','sel':MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags,'suffix':'MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags','leadjet':self.ak4BJets[0],'subleadjet':self.ak4BJets[1],'lead_is_b':True,'sublead_is_b':True},
            {'channel':'ElMu','sel':ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags,'suffix':'ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags','leadjet':self.ak4BJets[0],'subleadjet':self.ak4BJets[1],'lead_is_b':True,'sublead_is_b':True},
                                                          ]
        for channelDict in HasTightDileptonResolvedJetsChannelListJets:
            plots.extend(makeAk4JetsPlots(**channelDict))

        #############################################################################
        #                              AK8 Jets                                     #
        #############################################################################
        # Jet Selection #
        ElElHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet = ElElHasTightDileptonPreMllCutOutZTightVeto.refine("ElElHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet", 
                                cut=[op.rng_len(self.ak8Jets)>=1])
        MuMuHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet = MuMuHasTightDileptonPreMllCutOutZTightVeto.refine("MuMuHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet", 
                                cut=[op.rng_len(self.ak8Jets)>=1])
        ElMuHasTightDileptonPreMllCutTightVetoOneAk8Jet     = ElMuHasTightDileptonPreMllCutTightVeto.refine("ElMuHasTightDileptonPreMllCutTightVetoOneAk8Jet", 
                                cut=[op.rng_len(self.ak8Jets)>=1])

        # Yield 
        plots.append(yieldPlots.addYield(ElElHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet,
                                         "ElElHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet",
                                         "OS leptons + Ak8 Jets $\geq 1$ (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet,
                                         "MuMuHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet",
                                         "OS leptons + Ak8 Jets $\geq 1$ (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasTightDileptonPreMllCutTightVetoOneAk8Jet,
                                         "ElMuHasTightDileptonPreMllCutTightVetoOneAk8Jet",
                                         "OS leptons + Ak8 Jets $\geq 1$ (channel $e^{\pm}\mu^{\mp}$)"))

        # Plots #
        HasTightDileptonOneAk8JetJetsChannelList = [
                             {'channel':'ElEl','sel':ElElHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet,'dilepton':OsElElDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoOneAk8Jet'},
                             {'channel':'MuMu','sel':MuMuHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet,'dilepton':OsMuMuDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoOneAk8Jet'},
                             {'channel':'ElMu','sel':ElMuHasTightDileptonPreMllCutTightVetoOneAk8Jet,    'dilepton':OsElMuDilepton[0],'suffix':'HasTightDileptonPreMllCutTightVetoOneAk8Jet'},
                                                  ]
        commonItems = ['channel','sel','suffix']
        HasTightDileptonOneAk8JetsBaseDictFatjetsN = {'objName':'Ak8Jets','objCont':self.ak8Jets,'Nmax':10,'xTitle':'N(Ak8 jets)'}

                                                    
        for channelDict in HasTightDileptonOneAk8JetJetsChannelList:
            # Dilepton #
            plots.extend(makeDileptonPlots(**channelDict,isMC=isMC))
            # Ak8 Jets #
            plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**HasTightDileptonOneAk8JetsBaseDictFatjetsN))
            plots.extend(makeAk8Plots(fatjet=self.ak8Jets[0],**{k:channelDict[k] for k in commonItems}))
            # MET #
            plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))


        #############################################################################
        #                             Boosted jets                                  #
        #############################################################################
        # Boosted selection #
        ElElHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted = ElElHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet.refine("ElElHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted", 
                                cut=[op.rng_len(self.ak8BJets)>=1])
        MuMuHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted = MuMuHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet.refine("MuMuHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted", 
                                cut=[op.rng_len(self.ak8BJets)>=1])
        ElMuHasTightDileptonPreMllCutTightVetoInclusiveBoosted     = ElMuHasTightDileptonPreMllCutTightVetoOneAk8Jet.refine("ElMuHasTightDileptonPreMllCutTightVetoInclusiveBoosted", 
                                cut=[op.rng_len(self.ak8BJets)>=1])

#        ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveBoosted = ElElHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet.refine("ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveBoosted", 
#                                cut=[op.rng_len(self.ak8BJets)>=1, op.OR(op.rng_len(self.ak4Jets)<=1 , op.rng_len(self.ak4BJets)==0)])
#        MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveBoosted = MuMuHasTightDileptonPreMllCutOutZTightVetoOneAk8Jet.refine("MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveBoosted", 
#                                cut=[op.rng_len(self.ak8BJets)>=1, op.OR(op.rng_len(self.ak4Jets)<=1 , op.rng_len(self.ak4BJets)==0)])
#        ElMuHasTightDileptonPreMllCutTightVetoExclusiveBoosted     = ElMuHasTightDileptonPreMllCutTightVetoOneAk8Jet.refine("ElMuHasTightDileptonPreMllCutTightVetoExclusiveBoosted", 
#                                cut=[op.rng_len(self.ak8BJets)>=1, op.OR(op.rng_len(self.ak4Jets)<=1 , op.rng_len(self.ak4BJets)==0)])

        # Yield 
        plots.append(yieldPlots.addYield(ElElHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted,
                                         "ElElHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted",
                                         "OS leptons + Inclusive Boosted Jets (channel $e^+e^-$)"))
        plots.append(yieldPlots.addYield(MuMuHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted,
                                         "MuMuHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted",
                                         "OS leptons + Inclusive Boosted Jets (channel $\mu^+\mu^-$)"))
        plots.append(yieldPlots.addYield(ElMuHasTightDileptonPreMllCutTightVetoInclusiveBoosted,
                                         "ElMuHasTightDileptonPreMllCutTightVetoInclusiveBoosted",
                                         "OS leptons + Inclusive Boosted Jets (channel $e^{\pm}\mu^{\mp}$)"))
#        plots.append(yieldPlots.addYield(ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveBoosted,
#                                         "ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveBoosted",
#                                         "OS leptons + Exclusive Boosted Jets (channel $e^+e^-$)"))
#        plots.append(yieldPlots.addYield(MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveBoosted,
#                                         "MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveBoosted",
#                                         "OS leptons + Exclusive Boosted Jets (channel $\mu^+\mu^-$)"))
#        plots.append(yieldPlots.addYield(ElMuHasTightDileptonPreMllCutTightVetoExclusiveBoosted,
#                                         "ElMuHasTightDileptonPreMllCutTightVetoExclusiveBoosted",
#                                         "OS leptons + Exclusive Boosted Jets (channel $e^{\pm}\mu^{\mp}$)"))

        # Plots #
        HasTightDileptonBoostedJetsChannelList = [
             {'channel':'ElEl','sel':ElElHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted,'dilepton':OsElElDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted'},
             {'channel':'MuMu','sel':MuMuHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted,'dilepton':OsMuMuDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted'},
             {'channel':'ElMu','sel':ElMuHasTightDileptonPreMllCutTightVetoInclusiveBoosted,    'dilepton':OsElMuDilepton[0],'suffix':'HasTightDileptonPreMllCutTightVetoInclusiveBoosted'},
#             {'channel':'ElEl','sel':ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveBoosted,'dilepton':OsElElDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoExclusiveBoosted'},
#             {'channel':'MuMu','sel':MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveBoosted,'dilepton':OsMuMuDilepton[0],'suffix':'HasTightDileptonPreMllCutOutZTightVetoExclusiveBoosted'},
#             {'channel':'ElMu','sel':ElMuHasTightDileptonPreMllCutTightVetoExclusiveBoosted,    'dilepton':OsElMuDilepton[0],'suffix':'HasTightDileptonPreMllCutTightVetoExclusiveBoosted'},
                                                  ]
        commonItems = ['channel','sel','suffix']
        HasTightDileptonBoostedJetsBaseDictJetsN = {'objName':'Ak8BJets','objCont':self.ak8BJets,'Nmax':10,'xTitle':'N(Ak8 Bjets)'}

        for channelDict in HasTightDileptonBoostedJetsChannelList:
            # Dilepton #
            plots.extend(makeDileptonPlots(**channelDict,isMC=isMC))
            # Ak8 Jets #
            plots.append(objectsNumberPlot(**{k:channelDict[k] for k in commonItems},**HasTightDileptonBoostedJetsBaseDictJetsN))
            plots.extend(makeAk8Plots(fatjet=self.ak8BJets[0],**{k:channelDict[k] for k in commonItems}))
            # MET #
            plots.extend(makeMETPlots(**{k:channelDict[k] for k in commonItems}, met=self.corrMET))



        #############################################################################
        #                     High-level combinations                               #
        #############################################################################
        HasOsPreMllCutHighLevelVariablesChannelList = [
            # Resolved No Btag #
            {'channel'      : 'ElEl',
             'met'          : self.corrMET,
             'l1'           : OsElElDilepton[0][0],
             'l2'           : OsElElDilepton[0][1],
             'j1'           : self.ak4Jets[0],
             'j2'           : self.ak4Jets[1],
             'sel'          : ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag ,
             'suffix'       : 'ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag'},
            {'channel'      : 'MuMu',
             'met'          : self.corrMET,
             'l1'           : OsMuMuDilepton[0][0],
             'l2'           : OsMuMuDilepton[0][1],
             'j1'           : self.ak4Jets[0],
             'j2'           : self.ak4Jets[1],
             'sel'          : MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag,
             'suffix'       : 'MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedNoBtag'},
            {'channel'      : 'ElMu',
             'met'          : self.corrMET,
             'l1'           : OsElMuDilepton[0][0],
             'l2'           : OsElMuDilepton[0][1],
             'j1'           : self.ak4Jets[0],
             'j2'           : self.ak4Jets[1],
             'sel'          : ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedNoBtag,
             'suffix'       : 'ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedNoBtag'},
            # Resolved One Btag #
            {'channel'      : 'ElEl',
             'met'          : self.corrMET,
             'l1'           : OsElElDilepton[0][0],
             'l2'           : OsElElDilepton[0][1],
             'j1'           : self.ak4BJets[0],
             'j2'           : self.ak4Jets[0],
             'sel'          : ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag,
             'suffix'       : 'ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag'},
            {'channel'      : 'MuMu',
             'met'          : self.corrMET,
             'l1'           : OsMuMuDilepton[0][0],
             'l2'           : OsMuMuDilepton[0][1],
             'j1'           : self.ak4BJets[0],
             'j2'           : self.ak4Jets[0],
             'sel'          : MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag, 
             'suffix'       : 'MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedOneBtag'},
            {'channel'      : 'ElMu',
             'met'          : self.corrMET,
             'l1'           : OsElMuDilepton[0][0],
             'l2'           : OsElMuDilepton[0][1],
             'j1'           : self.ak4BJets[0],
             'j2'           : self.ak4Jets[0],
             'sel'          : ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedOneBtag, 
             'suffix'       : 'ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedOneBtag'},
            # Resolved Two Btags  #
            {'channel'      : 'ElEl',
             'met'          : self.corrMET,
             'l1'           : OsElElDilepton[0][0],
             'l2'           : OsElElDilepton[0][1],
             'j1'           : self.ak4BJets[0],
             'j2'           : self.ak4BJets[1],
             'sel'          : ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags, 
             'suffix'       : 'ElElHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags'},
            {'channel'      : 'MuMu',
             'met'          : self.corrMET,
             'l1'           : OsMuMuDilepton[0][0],
             'l2'           : OsMuMuDilepton[0][1],
             'j1'           : self.ak4BJets[0],
             'j2'           : self.ak4BJets[1],
             'sel'          : MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags, 
             'suffix'       : 'MuMuHasTightDileptonPreMllCutOutZTightVetoExclusiveResolvedTwoBtags'},
            {'channel'      : 'ElMu',
             'met'          : self.corrMET,
             'l1'           : OsElMuDilepton[0][0],
             'l2'           : OsElMuDilepton[0][1],
             'j1'           : self.ak4BJets[0],
             'j2'           : self.ak4BJets[1],
             'sel'          : ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags,
             'suffix'       : 'ElMuHasTightDileptonPreMllCutTightVetoExclusiveResolvedTwoBtags'},
            # Boosted #
            {'channel'      : 'ElEl',
             'met'          : self.corrMET,
             'l1'           : OsElElDilepton[0][0],
             'l2'           : OsElElDilepton[0][1],
             'j1'           : self.ak8BJets[0].subJet1,
             'j2'           : self.ak8BJets[0].subJet2,
             'sel'          : ElElHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted, 
             'suffix'       : 'ElElHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted'},
            {'channel'      : 'MuMu',
             'met'          : self.corrMET,
             'l1'           : OsMuMuDilepton[0][0],
             'l2'           : OsMuMuDilepton[0][1],
             'j1'           : self.ak8BJets[0].subJet1,
             'j2'           : self.ak8BJets[0].subJet2,
             'sel'          : MuMuHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted, 
             'suffix'       : 'MuMuHasTightDileptonPreMllCutOutZTightVetoInclusiveBoosted'},
            {'channel'      : 'ElMu',
             'met'          : self.corrMET,
             'l1'           : OsElMuDilepton[0][0],
             'l2'           : OsElMuDilepton[0][1],
             'j1'           : self.ak8BJets[0].subJet1,
             'j2'           : self.ak8BJets[0].subJet2,
             'sel'          : ElMuHasTightDileptonPreMllCutTightVetoInclusiveBoosted, 
             'suffix'       : 'ElMuHasTightDileptonPreMllCutTightVetoInclusiveBoosted'},
                                      ]

        for channelDict in HasOsPreMllCutHighLevelVariablesChannelList:
            plots.extend(makeHighLevelQuantities(**channelDict))

        return plots


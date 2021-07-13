import os
import sys
from copy import copy, deepcopy
from bamboo import treefunctions as op
from itertools import chain
import random

###################################################################################################################
# ModelPathDict is from Plotter. Keep the order same for function List in Stage-2
# Apply the JPA categorization BDT on each event using all possible permutations of jets
# Take the maximum BDT response
# Keep the max score for each BDT-Stage 1
# Pass them to the 2nd BDT {Multiclass: 7 classes for resolved, 2 for boosted}
# Return the index of maximum output, selectionObject and the final combination of Jets                                 
###################################################################################################################

###################################################################################################################
#               Jet-Parton Assignment BDTs (input variables for each of the JPA categories)                       #
###################################################################################################################
# Resolved 2b2Wj :: JPA_4jet

def bJetCorrPT(j):
    return j.pt*j.bRegCorr

def evaluateJPA_2b2Wj(lepton, muons, electrons, ak4jets, jets, bJetsL, bJetsM, met, model, HLL):
    #invars = [op.c_float(0.)]*14
    invars = [jets[0].btagDeepFlavB,                                                  # bjet1_btagCSV
              bJetCorrPT(jets[1]),                                                    # bjet2_ptReg
              jets[1].btagDeepFlavB,                                                  # bjet2_btagCSV
              jets[1].qgl,                                                            # bjet2_qgDiscr
              op.max(jets[2].btagDeepFlavB,jets[3].btagDeepFlavB),                    # maxwjetbtagCSV
              bJetCorrPT(jets[3]),                                                    # wjet2_ptReg
              jets[3].qgl,                                                            # wjet2_qgDiscr
              HLL.Wjj_simple(jets[2].p4, jets[3].p4).M(),                             # HadW_mass
              op.deltaR(HLL.Wjj_simple(jets[2].p4, jets[3].p4), lepton.p4),           # dR_HadW_lep
              op.max(HLL.dR_HadW_bjet(jets[0].p4, jets[2].p4, jets[3].p4),            # max_dR_HadW_bjet
                     HLL.dR_HadW_bjet(jets[1].p4, jets[2].p4, jets[3].p4)),
              op.deltaR(jets[2].p4, jets[3].p4),                                      # dR_wjet1wjet2 
              (HLL.bJetCorrP4(jets[0]) + HLL.Wjj_simple(jets[2].p4, jets[3].p4)).M(), # mTop1 = (WhadP4 + bJet1P4).M()
              op.rng_len(bJetsM),                                                     # nBJetMedium
              (HLL.bJetCorrP4(jets[0]) + HLL.bJetCorrP4(jets[1])).M()                 # Hbb_massReg
          ]
    return model(*invars, defineOnFirstUse=False)[0]

# Resolved 2b1Wj :: JPA_MissingWjet
def evaluateJPA_2b1Wj(lepton, muons, electrons, ak4jets, jets, bJetsL, bJetsM, met, model, HLL):
    #invars = [op.c_float(0.)]*13
    invars = [jets[0].btagDeepFlavB,                                                  # bjet1_btagCSV
              bJetCorrPT(jets[1]),                                                    # bjet2_ptReg
              jets[1].btagDeepFlavB,                                                  # bjet2_btagCSV
              jets[1].qgl,                                                            # bjet2_qgDiscr
              op.abs(lepton.eta - jets[1].eta),                                       # dEta_bjet2_lep
              bJetCorrPT(jets[2]),                                                    # wjet1_ptReg
              jets[2].btagDeepFlavB,                                                  # wjet1_btagCSV
              jets[2].qgl,                                                            # wjet1_qgDiscr
              op.abs(lepton.eta - jets[2].eta),                                       # dEta_wjet1_lep
              op.deltaR(jets[0].p4, jets[1].p4),                                      # dR_bjet1bjet2 
              op.rng_len(ak4jets),                                                    # nJets
              op.rng_len(bJetsL),                                                     # nBJetLoose
              op.rng_len(bJetsM)                                                      # nBJetMedium
          ]
    return model(*invars, defineOnFirstUse=False)[0]

# Resolved 2b0Wj :: JPA_MissingAllWjet
def evaluateJPA_2b0Wj(lepton, muons, electrons, ak4jets, jets, bJetsL, bJetsM, met, model, HLL):
    #invars = [op.c_float(0.)]*9
    invars = [bJetCorrPT(jets[0]),                                                    # bjet1_ptReg
              jets[0].btagDeepFlavB,                                                  # bjet1_btagCSV
              bJetCorrPT(jets[1]),                                                    # bjet2_ptReg
              jets[1].btagDeepFlavB,                                                  # bjet2_btagCSV
              jets[1].qgl,                                                            # bjet2_qgDiscr
              op.abs(lepton.eta - jets[1].eta),                                       # dEta_bjet2_lep
              (HLL.bJetCorrP4(jets[0]) + HLL.bJetCorrP4(jets[1])).M(),                # Hbb_massReg
              op.rng_len(ak4jets),                                                    # nJets
              op.rng_len(bJetsM)                                                      # nBJetMedium
          ]
    return model(*invars, defineOnFirstUse=False)[0]
    
# Resolved 1b2Wj :: JPA_MissingBjet
def evaluateJPA_1b2Wj(lepton, muons, electrons, ak4jets, jets, bJetsL, bJetsM, met, model, HLL):
    #invars = [op.c_float(0.)]*12
    invars = [bJetCorrPT(jets[0]),                                                    # bjet1_ptReg
              jets[0].btagDeepFlavB,                                                  # bjet1_btagCSV
              bJetCorrPT(jets[1]),                                                    # wjet1_ptReg
              jets[1].btagDeepFlavB,                                                  # wjet1_btagCSV
              jets[2].btagDeepFlavB,                                                  # wjet2_btagCSV
              jets[2].qgl,                                                            # wjet2_qgDiscr
              HLL.Wjj_simple(jets[1].p4, jets[2].p4).M(),                             # HadW_mass
              op.deltaR(HLL.Wjj_simple(jets[1].p4, jets[2].p4), lepton.p4),           # dR_HadW_lep
              op.deltaR(jets[1].p4, jets[2].p4),                                      # dR_wjet1wjet2 
              (HLL.bJetCorrP4(jets[0]) + HLL.Wjj_simple(jets[1].p4, jets[2].p4)).M(), # mTop1 = (WhadP4 + bJet1P4).M()
              op.rng_len(ak4jets),                                                    # nJets
              op.rng_len(bJetsM)                                                      # nBJetMedium
          ]
    return model(*invars, defineOnFirstUse=False)[0]
    
# Resolved 1b1Wj :: JPA_MissingBjet_MissingWJet
def evaluateJPA_1b1Wj(lepton, muons, electrons, ak4jets, jets, bJetsL, bJetsM, met, model, HLL):
    #invars = [op.c_float(0.)]*12
    invars = [bJetCorrPT(jets[0]),                                                    # bjet1_ptReg
              jets[0].btagDeepFlavB,                                                  # bjet1_btagCSV
              jets[0].qgl,                                                            # bjet1_qgDiscr
              op.abs(lepton.eta - jets[0].eta),                                       # dEta_bjet1_lep              
              bJetCorrPT(jets[1]),                                                    # wjet1_ptReg
              jets[1].btagDeepFlavB,                                                  # wjet1_btagCSV
              jets[1].qgl,                                                            # wjet1_qgDiscr
              op.abs(lepton.eta - jets[1].eta),                                       # dEta_wjet1_lep
              op.rng_len(ak4jets),                                                    # nJets
              op.rng_len(bJetsL),                                                     # nBJetLoose
              op.rng_len(bJetsM),                                                     # nBJetMedium
              op.rng_len(muons) + op.rng_len(electrons)                               # nLep
          ]
    return model(*invars, defineOnFirstUse=False)[0]
    
# Resolved 1b1Wj :: JPA_MissingBjet_MissingAllWJet
def evaluateJPA_1b0Wj(lepton, muons, electrons, ak4jets, jets, bJetsL, bJetsM, met, model, HLL):
    #invars = [op.c_float(0.)]*8
    invars = [lepton.pt,                                                              # leptonPt                                         
              bJetCorrPT(jets[0]),                                                    # bjet1_ptReg
              jets[0].btagDeepFlavB,                                                  # bjet1_btagCSV
              jets[0].qgl,                                                            # bjet1_qgDiscr
              op.abs(lepton.eta - jets[0].eta),                                       # dEta_bjet1_lep
              op.rng_len(ak4jets),                                                    # nJets
              op.rng_len(bJetsL),                                                     # nBJetLoose
              op.rng_len(bJetsM)                                                      # nBJetMedium
          ]
    return model(*invars, defineOnFirstUse=False)[0]
    
# Boosted Hbb2Wj :: JPA_Boosted
def evaluateJPA_Hbb2Wj(lepton, muons, electrons, fatJets, jets, bJetsL, bJetsM, met, model, HLL):
    #invars = [op.c_float(0.)]*14
    invars = [jets[0].btagDeepB,                                                      # wjet1_btagCSV
              jets[0].qgl,                                                            # wjet1_qgDiscr
              jets[1].btagDeepB,                                                      # wjet2_btagCSV
              jets[1].qgl,                                                            # wjet2_qgDiscr
              (fatJets[0].subJet1.p4 + fatJets[0].subJet2.p4).Pt(),                   # Hbb_Pt
              HLL.Wjj_simple(jets[0].p4, jets[1].p4).M(),                             # HadW_mass
              #HLL.comp_cosThetaS(jets[0].p4, jets[1].p4),                             # HadW_cosTheta
              op.extMethod("HHbbWWJPA::cosThetaS", returnType="float")(jets[0].p4, jets[1].p4),
              #op.c_float(random.choice([0,1])),
              op.deltaR(HLL.Wjj_simple(jets[0].p4, jets[1].p4), lepton.p4),           # dR_HadW_lep
              op.min(HLL.dR_HadW_bjet(fatJets[0].subJet1.p4, jets[0].p4, jets[1].p4), # min_dR_HadW_bjet
                     HLL.dR_HadW_bjet(fatJets[0].subJet2.p4, jets[0].p4, jets[1].p4)),
              op.deltaR(jets[0].p4, jets[1].p4),                                      # dR_wjet1wjet2
              #HLL.HWW_simple(jets[0].p4, jets[1].p4, lepton.p4, met).M()              # Hww_mass 
              #(jets[0].p4+jets[1].p4+lepton.p4+met.p4).M()                            # Hww_mass 
              op.extMethod("HHbbWWJPA::HWW_SimpleMassWithRecoNeutrino", returnType="float")(jets[0].p4, jets[1].p4, lepton.p4, met.p4)
    ]
    
    return model(*invars, defineOnFirstUse=False)[0]
    
# Boosted Hbb1Wj :: JPA_Boosted_MissingWjet
def evaluateJPA_Hbb1Wj(lepton, muons, electrons, fatJets, jets, bJetsL, bJetsM, met, model, HLL):
    #invars = [op.c_float(0.)]*14
    invars = [op.switch(fatJets[0].subJet1.btagDeepB > fatJets[0].subJet2.btagDeepB,  # dEta_bjet1_lep 
                        op.abs(lepton.eta-fatJets[0].subJet1.eta),
                        op.abs(lepton.eta-fatJets[0].subJet2.eta)),        
              bJetCorrPT(jets[0]),                                                    # wjet1_ptReg
              jets[0].btagDeepB,                                                      # wjet1_btagCSV
              jets[0].qgl,                                                            # wjet1_qgDiscr
              op.deltaR(lepton.p4, jets[0].p4),                                       # dR_wjet1_lep
              (fatJets[0].subJet1.p4 + fatJets[0].subJet2.p4).Pt(),                   # Hbb_Pt
              op.rng_len(bJetsM)                                                      # nBJetMedium
          ]

    return model(*invars, defineOnFirstUse=False)[0]
    
###################################################################################################################
#                            Jet-Parton Assignment BDTs (asiign categories and jets)                              #
###################################################################################################################
def toJetIndices(jpaJetList):
    #maxJPAJets = tojetIndices(JPAJetList)[maxIdx]
    # use as `t.Jet[maxJPAJets[0]]` etc.
    print('--- Finding the indices of JPA jets')
    idx_per_entry = []
    for jpaJets in jpaJetList:
        N = len(jpaJets.cont.ranges) ## a hack, should not be needed after https://gitlab.cern.ch/cp3-cms/bamboo/-/issues/76 is implemented
        idx_per_entry.append(op.initList("std::vector<std::size_t>", "std::size_t", (jpaJets[i].idx for i in range(N))))
    return op.construct("std::vector<std::vector<std::size_t>>", op.initList("std::vector<std::vector<std::size_t>>", "std::vector<std::size_t>", idx_per_entry))


def makeOddEvenEvaluator(isOdd, oddPath, evenPath, mvaType="TMVA"):
    model_odd  = op.mvaEvaluator(oddPath , mvaType=mvaType)
    model_even = op.mvaEvaluator(evenPath, mvaType=mvaType)
    methName = model_even.evaluate._name
    assert model_odd.evaluate._name == methName
    switchModelEval = op.switch(isOdd, model_odd.evaluate._objStb, model_even.evaluate._objStb)
    return op.MVAEvaluator(getattr(switchModelEval, methName))

# ----------------------------------------- Resolved --------------------------------------- #

def findJPACategoryResolved (self, selObj, lepton, muons, electrons, jets, bJetsL, bJetsM, met, modelPathDict, event, HLL, nodeList, plot_yield=False):
    JPAfuncDict = {'f1':evaluateJPA_2b2Wj, 
                   'f2':evaluateJPA_2b1Wj, 
                   'f3':evaluateJPA_1b2Wj, 
                   'f4':evaluateJPA_2b0Wj, 
                   'f5':evaluateJPA_1b1Wj, 
                   'f6':evaluateJPA_1b0Wj}

    JPAMaxScoreList = []
    bestCombo_per_cat = []
    
    combo2_1b0W_1Wj = op.combine(jets, N=2, samePred=lambda j1,j2 : j1.idx != j2.idx)
    combo2_2b0Wj    = op.combine(jets, N=2, pred=lambda j1,j2 : j1.pt > j2.pt,  samePred=lambda j1,j2 : j1.idx != j2.idx)
    combo3_1b2Wj    = op.combine(jets, N=3, pred=lambda j1,j2,j3 : j2.pt > j3.pt,  samePred=lambda j1,j2 : j1.idx != j2.idx)
    combo3_2b1Wj    = op.combine(jets, N=3, pred=lambda j1,j2,j3 : j1.pt > j2.pt,  samePred=lambda j1,j2 : j1.idx != j2.idx)
    combo4          = op.combine(jets, N=4, pred=lambda j1,j2,j3,j4 : op.AND(j1.pt > j2.pt, j3.pt > j4.pt), samePred=lambda j1,j2 : j1.idx != j2.idx)

    funckeys = [k for k in JPAfuncDict.keys()]
    
    for idx, func in enumerate(funckeys):
        node        = nodeList[idx]
        modelpaths  = modelPathDict.get(node)
        model = makeOddEvenEvaluator(event%2, modelpaths[1], modelpaths[0], mvaType="TMVA")
        lambdaFunc = lambda jetCombo : JPAfuncDict[func](lepton, muons, electrons, jets, jetCombo, bJetsL, bJetsM, met, model, HLL)
        
        if idx == 0:
            best = op.rng_max_element_by(combo4, lambdaFunc)
            maxScore = op.switch(best.idx != -1, best.idx.op.this.result.second, op.c_float(-1.)) 
            #best.idx.op.this.canDefine=False
        elif idx == 1:
            best = op.rng_max_element_by(combo3_2b1Wj, lambdaFunc)            ## hack: index of best is first in a pair, with the maximum value as second
            maxScore = best.idx.op.this.result.second
        elif idx == 2:
            best = op.rng_max_element_by(combo3_1b2Wj, lambdaFunc)            ## hack: index of best is first in a pair, with the maximum value as second
            maxScore = best.idx.op.this.result.second
        elif idx == 3:
            best = op.rng_max_element_by(combo2_2b0Wj, lambdaFunc)            ## hack: index of best is first in a pair, with the maximum value as second
            maxScore = best.idx.op.this.result.second
        elif idx == 4:
            best = op.rng_max_element_by(combo2_1b0W_1Wj, lambdaFunc)            ## hack: index of best is first in a pair, with the maximum value as second
            maxScore = best.idx.op.this.result.second
        else:
            best = op.rng_max_element_by(combo2_1b0W_1Wj, lambdaFunc)            ## hack: index of best is first in a pair, with the maximum value as second
            maxScore = best.idx.op.this.result.second
            
        JPAMaxScoreList.append(op.pow((1.0 + op.sqrt((1 - maxScore)/(1 + maxScore))), -1))
        #JPAMaxScoreList.append(maxScore)
        bestCombo_per_cat.append(best)

    evtCat = makeOddEvenEvaluator(event%2, modelPathDict.get('evCat')[1], modelPathDict.get('evCat')[0], mvaType="TMVA")
    
    evtCatOutList = evtCat(*JPAMaxScoreList)
    maxIdx = op.rng_max_element_index(evtCatOutList)

    newSelObj  = copy(selObj)
    selObjJPAjetsIdxDict = {}

    for i, node in enumerate(nodeList):
        outSelObj = copy(newSelObj)
        outSelObj.selName += '%s'%node
        outSelObj.yieldTitle += " in %s node"%node 
        outSelObj.refine(cut = (maxIdx == i)) 
        if i < 6:
            selObjJPAjetsIdxDict[node] = [outSelObj, bestCombo_per_cat[i]]
        else:
            selObjJPAjetsIdxDict[node] = [outSelObj, None]

    return JPAMaxScoreList, evtCatOutList, selObjJPAjetsIdxDict
    
# ----------------------------------- Boosted -------------------------------------- #
def findJPACategoryBoosted (self, selObj, lepton, muons, electrons, fatJets, jets, bJetsL, bJetsM, met, modelPathDict, event, HLL, nodeList, plot_yield=False):
    JPAfuncDict = {'f1':evaluateJPA_Hbb2Wj, 
                   'f2':evaluateJPA_Hbb1Wj}

    JPAMaxScoreList = []
    bestCombo_per_cat = []
    
    combo2     = op.combine(jets, N=2, pred=lambda j1,j2 : j1.pt > j2.pt, samePred=lambda j1,j2 : j1.idx != j2.idx)
    fakeCombo2 = op.combine(jets, N=2, pred=lambda j1,j2 : j1.pt >= j2.pt, samePred=None)

    funckeys = [k for k in JPAfuncDict.keys()]
    for idx, func in enumerate(funckeys):
        node        = nodeList[idx]
        modelpaths  = modelPathDict.get(node)
        model = makeOddEvenEvaluator(event%2, modelpaths[1], modelpaths[0], mvaType="TMVA")
        lambdaFunc = lambda jetCombo : JPAfuncDict[func](lepton, muons, electrons, fatJets, jetCombo, bJetsL, bJetsM, met, model, HLL)

        if idx == 0:
            best = op.rng_max_element_by(combo2, lambdaFunc)
            maxScore = op.switch(best.idx != -1, best.idx.op.this.result.second, op.c_float(-1.))
        else:
            best = op.rng_max_element_by(fakeCombo2, lambdaFunc)
            #best = op.rng_max_element_by(combo2, lambdaFunc)
            maxScore = best.idx.op.this.result.second
            #maxScore = op.switch(best.idx != -1, best.idx.op.this.result.second, op.c_float(-1.))

        JPAMaxScoreList.append(op.pow((1.0 + op.sqrt((1 - maxScore)/(1 + maxScore))), -1))
        #JPAMaxScoreList.append(maxScore)
        bestCombo_per_cat.append(best)

    evtCat = makeOddEvenEvaluator(event%2, modelPathDict.get('evCat')[1], modelPathDict.get('evCat')[0], mvaType="TMVA")
    JPAL2outList = evtCat(*JPAMaxScoreList)
    maxIdx = op.rng_max_element_index(JPAL2outList)

    newSelObj  = copy(selObj)
    selObjJPAjetsIdxDict = {}
    for i, node in enumerate(nodeList):
        outSelObj = copy(newSelObj)
        outSelObj.selName += '%s'%node
        outSelObj.yieldTitle += " in %s node"%node 
        outSelObj.refine(cut = [maxIdx == op.c_int(i)])
        #if plot_yield:
        #    outSelObj.makeYield(self.yieldPlots)
        if i < 2:
            selObjJPAjetsIdxDict[node] = [outSelObj, bestCombo_per_cat[i]]
        else:
            selObjJPAjetsIdxDict[node] = [outSelObj, None]

    return JPAMaxScoreList, JPAL2outList, selObjJPAjetsIdxDict
    

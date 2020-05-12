import math 
from bamboo import treefunctions as op

def METFilter(flags, era, isMC):
    # from https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    #
    # from  https://gitlab.cern.ch/cms-hh-bbww/cms-hh-to-bbww/blob/master/Legacy/event_selection.md
    #   primary vertex filter ("Flag_goodVertices"): data and MC 2016, 2017, 2018.
    #   beam halo filter ("Flag_globalSuperTightHalo2016Filter"): data and MC 2016, 2017, 2018.
    #   HBHE noise filter ("Flag_HBHENoiseFilter"): data and MC 2016, 2017, 2018.
    #   HBHEiso noise filter ("Flag_HBHENoiseIsoFilter"): data and MC 2016, 2017, 2018.
    #   ECAL TP filter ("Flag_EcalDeadCellTriggerPrimitiveFilter"): data and MC 2016, 2017, 2018.
    #   Bad PF Muon Filter ("Flag_BadPFMuonFilter"): data and MC 2016, 2017, 2018.
    #   ECAL bad calibration filter ("Flag_ecalBadCalibReducedMINIAODFilter"): data and MC 2017, 2018.
    #   ECAL endcap bad SC noise filter ("Flag_eeBadScFilter"): data 2016, 2017, 2018.
    cuts = [
           flags.goodVertices,
           flags.globalSuperTightHalo2016Filter,
           flags.HBHENoiseFilter,
           flags.HBHENoiseIsoFilter,
           flags.EcalDeadCellTriggerPrimitiveFilter,
           flags.BadPFMuonFilter
            ]
#    if "2017" in era or "2018" in era:
#        cuts.append(flags.ecalBadCalibReducedMINIAODFilter) # Only 2017-2018 : both MC and data
        # Could not find in NanoAOD : TODO 
    if not isMC:
        cuts.append(flags.eeBadScFilter) # Only data 

    return cuts


class METcorrection(object):
    # https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection.h
    def __init__(self,rawMET,pv,sample,era,isMC):
        if(era=='2016'):
            if isMC:
                xcorr = (0.195191,   0.170948)
                ycorr = (0.0311891, -0.787627)
            else:
                if '2016B' in sample:
                    xcorr = (0.0478335,  0.108032)
                    ycorr = (-0.125148, -0.355672)
                elif '2016C' in sample: 
                    xcorr = ( 0.0916985, -0.393247)
                    ycorr = (-0.151445,  -0.114491)
                elif '2016D' in sample:
                    xcorr = ( 0.0581169, -0.567316)
                    ycorr = (-0.147549,  -0.403088)
                elif '2016E' in sample:
                    xcorr = ( 0.065622, -0.536856)
                    ycorr = (-0.188532, -0.495346)
                elif '2016F' in sample:
                    xcorr = ( 0.0313322, -0.39866)
                    ycorr = (-0.16081,   -0.960177)
                elif '2016G' in sample:
                    xcorr = (-0.040803,   0.290384)
                    ycorr = (-0.0961935, -0.666096)
                else:
                    xcorr = (-0.0330868, 0.209534)
                    ycorr = (-0.141513, -0.816732)


        elif(era=='2017'):
            if isMC:
                xcorr = (0.217714, -0.493361)
                ycorr = (-0.177058, 0.336648)
            #these are the corrections for v2 MET recipe (currently recommended for 2017)
            else: 
                if '2017B' in sample:
                    xcorr = ( 0.19563, -1.51859)
                    ycorr = (-0.306987, 1.84713)
                elif '2017C' in sample:
                    xcorr = ( 0.161661, -0.589933)
                    ycorr = (-0.233569,  0.995546)
                elif '2017D' in sample:
                    xcorr = ( 0.180911, -1.23553)
                    ycorr = (-0.240155,  1.27449)
                elif '2017E' in sample:
                    xcorr = ( 0.149494, -0.901305)
                    ycorr = (-0.178212,  0.535537)
                else:
                    xcorr = ( 0.165154, -1.02018)
                    ycorr = (-0.253794, -0.75776)
        elif(era=='2018'):
            if isMC:
                xcorr = (-0.296713,  0.141506)
                ycorr = (-0.115685, -0.0128193)
            else:
                if '2018A' in sample:
                    xcorr= (-0.362865,  1.94505)
                    ycorr= (-0.0709085, 0.307365)
                elif'2018B' in sample:
                    xcorr = (-0.492083, 2.93552)
                    ycorr = (-0.17874,  0.786844)
                elif '2018C' in sample:
                    xcorr = (-0.521349, 1.44544)
                    ycorr = (-0.118956, 1.96434)
                else:
                    xcorr = (-0.531151,  1.37568)
                    ycorr = (-0.0884639, 1.57089)
        
        METxcorr=xcorr[0] *pv.npvs+xcorr[1]
        METycorr=ycorr[0] *pv.npvs+ycorr[1]
        
        corrMETx=rawMET.pt*op.cos(rawMET.phi) +METxcorr
        corrMETy=rawMET.pt*op.sin(rawMET.phi) +METycorr
        
        self.pt=op.sqrt(corrMETx**2 +corrMETy**2)
        atan=op.atan(corrMETy/corrMETx)
        self.phi=op.multiSwitch((corrMETx> 0,atan),(corrMETy> 0,atan+math.pi),atan-math.pi)


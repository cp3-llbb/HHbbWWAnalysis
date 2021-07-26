#!/bin/bash

# Electrons #
python extractElectronScaleFactorsFromRoot.py --suffix 2016 ../data/ScaleFactors_POG/POG/2016/2016LegacyReReco_ElectronMVA90_Fall17V2.root 
python extractElectronScaleFactorsFromRoot.py --suffix 2017 ../data/ScaleFactors_POG/POG/2017/2017_ElectronMVA90.root 
python extractElectronScaleFactorsFromRoot.py --suffix 2018 ../data/ScaleFactors_POG/POG/2018/2018_ElectronMVA90.root 

# Muons #
python extractMuonScaleFactorsFromRoot.py -s 2016_RunBCDEFGH_ID ../data/ScaleFactors_POG/POG/2016/RunBCDEFGH_SF_ID_DESY.root
python extractMuonScaleFactorsFromRoot.py -s 2016_RunBCDEFGH_ISO ../data/ScaleFactors_POG/POG/2016/RunBCDEFGH_SF_ISO_DESY.root

#python extractMuonScaleFactorsFromRoot.py -s 2017_RunBCDEF_ID ../data/ScaleFactors_POG/POG/2017/RunBCDEF_SF_ID_syst.root
#python extractMuonScaleFactorsFromRoot.py -s 2017_RunBCDEF_ISO ../data/ScaleFactors_POG/POG/2017/RunBCDEF_SF_ISO_syst.root
python extractMuonScaleFactorsFromRoot.py -s 2017_RunBCDEF_ID ../data/ScaleFactors_POG/POG/2017/RunBCDEF_SF_ID_DESY_2017.root
python extractMuonScaleFactorsFromRoot.py -s 2017_RunBCDEF_ISO ../data/ScaleFactors_POG/POG/2017/RunBCDEF_SF_ISO_DESY_2017.root

#python extractMuonScaleFactorsFromRoot.py -s 2018_RunABCD_ID ../data/ScaleFactors_POG/POG/2018/RunABCD_SF_ID.root
#python extractMuonScaleFactorsFromRoot.py -s 2018_RunABCD_ISO ../data/ScaleFactors_POG/POG/2018/RunABCD_SF_ISO.root
python extractMuonScaleFactorsFromRoot.py -s 2018_RunABCD_ID ../data/ScaleFactors_POG/POG/2018/RunABCD_SF_ID_DESY_2018.root
python extractMuonScaleFactorsFromRoot.py -s 2018_RunABCD_ISO ../data/ScaleFactors_POG/POG/2018/RunABCD_SF_ISO_DESY_2018.root

# ID -> NUM_TightID_DEN_genTracks_pt_abseta (2016/2017) and NUM_TightID_DEN_TrackerMuons_pt_abseta (2018)
# ISO -> NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta

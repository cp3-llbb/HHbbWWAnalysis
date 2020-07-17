#!/bin/bash
# Electron Loose efficiency #  
python extractElectronScaleFactorsFromRoot.py --suffix LooseEff2016 ../data/ScaleFactors_ttH/2016/ElectronLooseEfficiency/TnP_loose_ele_2016.root  &&
python extractElectronScaleFactorsFromRoot.py --suffix LooseEff2017 ../data/ScaleFactors_ttH/2017/ElectronLooseEfficiency/TnP_loose_ele_2017.root  &&
python extractElectronScaleFactorsFromRoot.py --suffix LooseEff2018 ../data/ScaleFactors_ttH/2018/ElectronLooseEfficiency/TnP_loose_ele_2018.root  &&

# Electron Reco efficiency #  
python extractElectronScaleFactorsFromRoot.py --suffix RecoEff2016PtGt20 ../data/ScaleFactors_ttH/2016/ElectronRecoEfficiency/el_scaleFactors_gsf_ptGt20.root   &&
python extractElectronScaleFactorsFromRoot.py --suffix RecoEff2016PtLt20 ../data/ScaleFactors_ttH/2016/ElectronRecoEfficiency/el_scaleFactors_gsf_ptLt20.root   &&
python extractElectronScaleFactorsFromRoot.py --suffix RecoEff2017PtGt20 ../data/ScaleFactors_ttH/2017/ElectronRecoEfficiency/el_scaleFactors_gsf_ptGt20.root   &&
python extractElectronScaleFactorsFromRoot.py --suffix RecoEff2017PtLt20 ../data/ScaleFactors_ttH/2017/ElectronRecoEfficiency/el_scaleFactors_gsf_ptLt20.root   &&
python extractElectronScaleFactorsFromRoot.py --suffix RecoEff2018 ../data/ScaleFactors_ttH/2018/ElectronRecoEfficiency/el_scaleFactors_gsf.root  &&

# Electron Loose MVA #  
python extractElectronScaleFactorsFromRoot.py --suffix LooseMVA2016 ../data/ScaleFactors_ttH/2016/ElectronLooseMVASF/TnP_loosettH_ele_2016.root   &&
python extractElectronScaleFactorsFromRoot.py --suffix LooseMVA2017 ../data/ScaleFactors_ttH/2017/ElectronLooseMVASF/TnP_loosettH_ele_2017.root   &&
python extractElectronScaleFactorsFromRoot.py --suffix LooseMVA2018 ../data/ScaleFactors_ttH/2018/ElectronLooseMVASF/TnP_loosettH_ele_2018.root  &&
  
# Electron Tight MVA #  
python extractTTHTightMVAScaleFactorsFromRoot.py --SF_file ../data/ScaleFactors_ttH/2016/ElectronTightMVASF/egammaEffi.txt_EGM2D.root --err_eta ../data/ScaleFactors_ttH/2016/ElectronTightMVASF/SFttbar_2016_ele_eta.root --err_pt ../data/ScaleFactors_ttH/2016/ElectronTightMVASF/SFttbar_2016_ele_pt.root  --suffix ElectronTight2016  &&
python extractTTHTightMVAScaleFactorsFromRoot.py --SF_file ../data/ScaleFactors_ttH/2017/ElectronTightMVASF/egammaEffi.txt_EGM2D.root --err_eta ../data/ScaleFactors_ttH/2017/ElectronTightMVASF/SFttbar_2017_ele_eta.root --err_pt ../data/ScaleFactors_ttH/2017/ElectronTightMVASF/SFttbar_2017_ele_pt.root  --suffix ElectronTight2017  &&
python extractTTHTightMVAScaleFactorsFromRoot.py --SF_file ../data/ScaleFactors_ttH/2018/ElectronTightMVASF/egammaEffi.txt_EGM2D.root --err_eta ../data/ScaleFactors_ttH/2018/ElectronTightMVASF/SFttbar_2018_ele_eta.root --err_pt ../data/ScaleFactors_ttH/2018/ElectronTightMVASF/SFttbar_2018_ele_pt.root  --suffix ElectronTight2018  &&
  
  
# Muon Loose MVA #  
python extractElectronScaleFactorsFromRoot.py --suffix Loose2016 ../data/ScaleFactors_ttH/2016/MuonLooseSF/TnP_loose_muon_2016.root   &&
python extractElectronScaleFactorsFromRoot.py --suffix Loose2017 ../data/ScaleFactors_ttH/2017/MuonLooseSF/TnP_loose_muon_2017.root   &&
python extractElectronScaleFactorsFromRoot.py --suffix Loose2018 ../data/ScaleFactors_ttH/2018/MuonLooseSF/TnP_loose_muon_2018.root   &&
  
  
# Rename (because scripts give the electron name ... anyway let's do it !)  
mv Electron_EGamma_SF2D_Loose2016.json Muon_EGamma_SF2D_Loose2016.json  &&
mv Electron_EGamma_SF2D_Loose2017.json Muon_EGamma_SF2D_Loose2017.json  &&
mv Electron_EGamma_SF2D_Loose2018.json Muon_EGamma_SF2D_Loose2018.json  &&

# Muon Tight MVA #  &&
python extractTTHTightMVAScaleFactorsFromRoot.py --SF_file ../data/ScaleFactors_ttH/2018/MuonTightMVASF/egammaEffi.txt_EGM2D.root --err_eta ../data/ScaleFactors_ttH/2018/MuonTightMVASF/SFttbar_2018_muon_eta.root --err_pt ../data/ScaleFactors_ttH/2018/MuonTightMVASF/SFttbar_2018_muon_pt.root  --suffix MuonTight2018  &&
python extractTTHTightMVAScaleFactorsFromRoot.py --SF_file ../data/ScaleFactors_ttH/2017/MuonTightMVASF/egammaEffi.txt_EGM2D.root --err_eta ../data/ScaleFactors_ttH/2017/MuonTightMVASF/SFttbar_2017_muon_eta.root --err_pt ../data/ScaleFactors_ttH/2017/MuonTightMVASF/SFttbar_2017_muon_pt.root  --suffix MuonTight2017  &&
python extractTTHTightMVAScaleFactorsFromRoot.py --SF_file ../data/ScaleFactors_ttH/2016/MuonTightMVASF/egammaEffi.txt_EGM2D.root --err_eta ../data/ScaleFactors_ttH/2016/MuonTightMVASF/SFttbar_2016_muon_eta.root --err_pt ../data/ScaleFactors_ttH/2016/MuonTightMVASF/SFttbar_2016_muon_pt.root  --suffix MuonTight2016  

# Single trigger SF #
python extractTTHSingleLeptonTriggerScaleFactorsFromRoot.py -s SingleElectron_2016 --SF_file ../data/ScaleFactors_ttH/2016/SingleElectronTriggerSF/Electron_Ele25WPTight_eff.root
python extractTTHSingleLeptonTriggerScaleFactorsFromRoot.py -s SingleElectron_2017 --SF_file ../data/ScaleFactors_ttH/2017/SingleElectronTriggerSF/Electron_Ele32orEle35_eff.root
python extractTTHSingleLeptonTriggerScaleFactorsFromRoot.py -s SingleElectron_2018 --SF_file ../data/ScaleFactors_ttH/2018/SingleElectronTriggerSF/Electron_Run2018_Ele32orEle35.root

python extractTTHSingleLeptonTriggerScaleFactorsFromRoot.py -s SingleMuon_2016 --SF_file ../data/ScaleFactors_ttH/2016/SingleMuonTriggerSF/Muon_Mu22OR_eta2p1_eff.root
python extractTTHSingleLeptonTriggerScaleFactorsFromRoot.py -s SingleMuon_2017 --SF_file ../data/ScaleFactors_ttH/2017/SingleMuonTriggerSF/Muon_IsoMu24orIsoMu27_eff.root
python extractTTHSingleLeptonTriggerScaleFactorsFromRoot.py -s SingleMuon_2018 --SF_file ../data/ScaleFactors_ttH/2018/SingleMuonTriggerSF/Muon_Run2018_IsoMu24orIsoMu27.root

# Fake Rates
python extractTTHFakeRates.py --era 2016 --file ../data/FakeRates/2016/FR_lep_ttH_mva_2016_CERN_2019Jul08.root 
python extractTTHFakeRates.py --era 2017 --file ../data/FakeRates/2017/FR_lep_ttH_mva_2017_CERN_2019Jul08.root 
python extractTTHFakeRates.py --era 2018 --file ../data/FakeRates/2018/FR_lep_ttH_mva_2018_CERN_2019Jul08.root 


from itertools import chain
import json
import yaml
import subprocess
import os.path
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

#triggers_per_era = {}
#triggers_per_era["2016B1"] = [
#    "HLT_IsoMu24", "HLT_IsoTkMu24", "HLT_Mu50", "HLT_Mu45_eta2p1",
#    "HLT_Ele27_WPTight_Gsf", "HLT_Ele105_CaloIdVT_GsfTrkIdT", "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon175",
#    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", "HLT_Mu30_TkMu11",
#    "HLT_Mu30_TkMu11",
#    "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
#    "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL",
#    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW"
#    ]
#triggers_per_era["2016B2"] = triggers_per_era["2016B1"] + [ "HLT_TkMu50" ]
#triggers_per_era["2016C"] = triggers_per_era["2016B2"]
#triggers_per_era["2016D"] = triggers_per_era["2016B2"]
#triggers_per_era["2016D"] = triggers_per_era["2016B2"]
#triggers_per_era["2016E"] = triggers_per_era["2016B2"]
#triggers_per_era["2016F1"] = triggers_per_era["2016B2"]
#triggers_per_era["2016F2"] = triggers_per_era["2016B2"] + [ "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL" ]
#triggers_per_era["2016G"] = triggers_per_era["2016F2"]
#triggers_per_era["2016H"] = [ trig for trig in triggers_per_era["2016F2"] if trig not in 
#    [ "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", "HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL", "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL" ]
#    ] + [ "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL", "HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ" ]
###
#triggers_per_era["2017B"] = [
#    "HLT_IsoMu24", "HLT_IsoMu24_eta2p1", "HLT_IsoMu27", "HLT_Mu50",
#    "HLT_Ele35_WPTight_Gsf", "HLT_Photon200",
#    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
#    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
#    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_DoubleEle33_CaloIdL_MW",
#    "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter", "Flag_ecalBadCalibFilterV2"
#    ]
#triggers_per_era["2017C1"] = triggers_per_era["2017B"]+[
#    "HLT_OldMu100", "HLT_TkMu100",
#    "HLT_Ele115_CaloIdVT_GsfTrkIdT",
#    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
#    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"
#    ]
#triggers_per_era["2017C2"] = triggers_per_era["2017C1"]+[
#    "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
#    "HLT_Mu37_TkMu27",
#    "HLT_Mu27_Ele37_CaloIdL_MW", "HLT_Mu37_Ele27_CaloIdL_MW"
#    ]
#triggers_per_era["2017D"] = triggers_per_era["2017C2"]
#triggers_per_era["2017E"] = triggers_per_era["2017C2"]
#triggers_per_era["2017F"] = triggers_per_era["2017C2"]
#triggers_per_era["2018ABC"] = [
#    "HLT_IsoMu24", "HLT_IsoMu27", "HLT_Mu50", "HLT_OldMu100", "HLT_TkMu100",
#    "HLT_Ele32_WPTight_Gsf", "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Photon200",
#    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_DoubleEle25_CaloIdL_MW",
#    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", "HLT_Mu37_TkMu27",
#    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
#    "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
#    "HLT_Mu27_Ele37_CaloIdL_MW", "HLT_Mu37_Ele27_CaloIdL_MW"
#    ]
#triggers_per_era["2018D"] = triggers_per_era["2018ABC"]
#
#dasQueries = {
#    "2017C1" : "/{pd}/piedavid-Run2017C-31Mar2018-v1_TopNanoAODv6-1-1_2017-9721c24ccc7f925c513e24ff74941177/USER instance=prod/phys03 run between [299368, 302025]",
#    "2017C2" : "/{pd}/piedavid-Run2017C-31Mar2018-v1_TopNanoAODv6-1-1_2017-9721c24ccc7f925c513e24ff74941177/USER instance=prod/phys03 run between [302026, 302029]",
#    }
#dasQueries.update(dict(
#    ("2017{0}".format(era), "/{{pd}}/piedavid-Run2017{era}-31Mar2018-v1_TopNanoAODv6-1-1_2017-9721c24ccc7f925c513e24ff74941177/USER instance=prod/phys03".format(era=era))
#    for era in ("BDEF")))
#
#dasQueries.update({
#    "2016B1": "/{pd}/schoef-TopNanoAODv6-1-2-4_{pd}_Run2016B_ver2-ba7e3129b1ff910ad1abce6981b33804/USER run between [272007, 274953] instance=prod/phys03",
#    "2016B2": "/{pd}/schoef-TopNanoAODv6-1-2-4_{pd}_Run2016B_ver2-ba7e3129b1ff910ad1abce6981b33804/USER run between [274954, 275376] instance=prod/phys03",
#    "2016F1": "/{pd}/schoef-TopNanoAODv6-1-2-4_{pd}_Run2016F-ba7e3129b1ff910ad1abce6981b33804/USER run between [277772, 278272] instance=prod/phys03",
#    "2016F2": "/{pd}/schoef-TopNanoAODv6-1-2-4_{pd}_Run2016F-ba7e3129b1ff910ad1abce6981b33804/USER run between [278273, 278808] instance=prod/phys03"
#    })
#dasQueries.update(dict(
#    ("2016{0}".format(era), "/{{pd}}/schoef-TopNanoAODv6-1-2-4_{{pd}}_Run2016{era}-ba7e3129b1ff910ad1abce6981b33804/USER instance=prod/phys03".format(era=era))
#    for era in ("C", "D", "E", "G", "H")))
#
#dasQueries.update({
#    "2018ABC": "/{pd}/palencia-TopNanoAODv6-1-2_2018-2667fe41f354e79b08df8a25806ccf17/USER instance=prod/phys03",
#    "2018D": "/{pd}/palencia-TopNanoAODv6-1-2_2018-831765d0aa9cd559fee11ff659127d4e/USER instance=prod/phys03"
#    })
#

triggers_per_era = {}
triggers_per_era["2016"] = ["IsoMu22","IsoTkMu22","IsoMu22_eta2p1","IsoTkMu22_eta2p1","IsoMu24","IsoTkMu24",
                            "Ele27_WPTight_Gsf","Ele25_eta2p1_WPTight_Gsf","Ele27_eta2p1_WPLoose_Gsf",
                            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL","Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ","Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL","Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
                            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL","Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ","Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL","Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ"]
                 


golden = {
    "2016": os.path.join(os.path.dirname(__file__), "Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt"),
    "2017": os.path.join(os.path.dirname(__file__), "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt"),
    "2018": os.path.join(os.path.dirname(__file__), "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt")
    }

storageBase = "/storage/data/cms"

if __name__ == "__main__":
    import argparse
    from argparse import ArgumentParser
    parser = ArgumentParser(description="Check for triggers in the files, and split eras as needed")
    parser.add_argument("--year", action="append", default=[], help="Data-taking year(s) to consider (multiple arguments possible)")
    #parser.add_argument("--siteinfo", required=True, help="siteinfo.yml path")
    parser.add_argument("--localprefix", default="/store/user/piedavid/topNanoAOD/v6-1-1/", help="Prefix in the local storage, for replicated LFNs")
    parser.add_argument("--output", "-o", default=".", help="Directory to write the output files to")
    args = parser.parse_args()

#    with open(args.siteinfo) as yF:
#        siteinfo = yaml.load(yF)
#    lfnprefixes = set(chain.from_iterable(pxs for uCfg in siteinfo["users"].values() for pxs in uCfg["prefix"].values()))

    from cppyy import gbl
    for year in args.year:
        with open(golden[year]) as jsf:
            goldenLumis = json.load(jsf)

        if year == "2018":
            primary_datasets = ["SingleMuon", "EGamma", "DoubleMuon", "MuonEG"]
        else:
            primary_datasets = ["SingleMuon", "SingleElectron", "DoubleMuon", "DoubleEG", "MuonEG"]
        for era, eraTriggers in triggers_per_era.items():
            if era.startswith(year):
                for pd in primary_datasets:
                    dasout = subprocess.check_output(["dasgoclient", "-query", "file dataset={0}".format(dasQueries[era].format(pd=pd))])
                    lfnList = [ ln.strip() for ln in dasout.decode().split("\n") if ln.strip() ]
                    notWithPrefix = [ lfn for lfn in lfnList if not any(lfn.startswith(pref) for pref in lfnprefixes) ]
                    if notWithPrefix:
                        raise RuntimeError("Some LFNs don't start with one of the known prefixes: {0}".format(notWithPrefix))
                    pfnList = [ "".join((storageBase, args.localprefix, lfn[len(next(pref for pref in lfnprefixes if lfn.startswith(pref))):])) for lfn in lfnList ]
                    assert len(pfnList) == len(lfnList)
                    logger.info(f"{len(pfnList)} files for {pd}_{era}")
                    nGood = 0
                    canSkip = []
                    for lfn,pfn in zip(lfnList, pfnList):
                        f = gbl.TFile.Open(pfn, "READ")
                        tEvt = f.Get("Events")
                        lvs = set(lv.GetName() for lv in tEvt.GetListOfLeaves())
                        trig_missing = [ trig for trig in eraTriggers if trig not in lvs ]
                        if not trig_missing:
                            nGood += 1
                        else:
                            runs_in_file = []
                            tRuns = f.Get("Runs")
                            for i in range(tRuns.GetEntries()):
                                tRuns.GetEntry(i)
                                runs_in_file.append(str(tRuns.run))
                            logger.info(f"File {pfn} for {pd}_{era} is missing {len(trig_missing):d} branches")
                            if not any(run in goldenLumis for run in runs_in_file):
                                logger.info(f"Contains events from runs {', '.join(runs_in_file)}, none of which on golden, so the file can be skipped")
                                canSkip.append((lfn, pfn))
                            else:
                                logger.info(f"Contains events from runs {', '.join(runs_in_file)}, missing branches are {', '.join(trig_missing)}")
                    logger.info(f"{nGood:d}/{len(pfnList):d} files passed the check")
                    outName = os.path.join(args.output, f"{pd}_Run{era}.dat")
                    skippedName = os.path.join(args.output, f"{pd}_Run{era}_skipped.dat")
                    if canSkip and len(canSkip)+nGood == len(pfnList):
                        if os.path.exists(outName):
                            logger.error(f"Output file {outName} exists, will not overwrite")
                        else:
                            with open(outName, "w") as outF:
                                outF.write("\n".join(pfn for pfn in pfnList if not any(ap == pfn for al,ap in canSkip)))
                            logger.info(f"All failing files can be skipped, new PFN list written to {outName}")
                            with open(skippedName, "w") as outF:
                                outF.write("\n".join(lfn for lfn,pfn in canSkip))
                    elif not os.path.exists(outName):
                        with open(outName, "w") as outF:
                            outF.write("\n".join(pfn for pfn in pfnList if pfn not in canSkip))
                        logger.info(f"New PFN list written to {outName}")

## ## NOTE the above recovers correctly for
## samples/TopNanoAODv6p1/2017/DoubleEG_Run2017C1.dat
## samples/TopNanoAODv6p1/2017/DoubleEG_Run2017E.dat
## samples/TopNanoAODv6p1/2017/SingleElectron_Run2017C1.dat
## samples/TopNanoAODv6p1/2017/SingleElectron_Run2017E.dat
## ## There is one file remaining then, /storage/data/cms/store/user/piedavid/topNanoAOD/v6-1-1/2017/DoubleEG/TopNanoAODv6-1-1_2017/200615_080356/0000/tree_221.root
## ## it misses some dimuon and muoneg triggers, but only contains the run 302027 in the C2 range, which is not in the golden json
## samples/TopNanoAODv6p1/2017/DoubleEG_Run2017C2.dat

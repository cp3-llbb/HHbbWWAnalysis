import os
import sys
import re
import glob
import collections
import logging
import argparse

class TemplateLatex:
    def __init__(self,dirpath,logplot=False):
        self.dirpath = dirpath
        self.texfile = os.path.join(self.dirpath,"plots.tex")
        self.filesdict = {}
        self.logplot = logplot

        # Ordered dicts so that plots keep that order #
        self.channels = collections.OrderedDict([
                        ("ElEl", "$e^+e^-$"),
                        ("MuMu", "$\mu^+\mu^-$"),
                        ("ElMu", "$e^{\pm}\mu^{\mp}"),
                        ])
            
        self.selections = collections.OrderedDict([
                        ('ExclusiveResolvedNoBtag_' , "Exclusive Resolved Jets (no Btag)"),
                        ('ExclusiveResolvedOneBtag_' , "Exclusive Resolved Jets (1 Btag)"),
                        ('ExclusiveResolvedTwoBtags_' , "Exclusive Resolved Jets (2 Btags)"),
                        ('InclusiveBoostedJets_' , "Inclusive Boosted Jets"),
                        ])
        self.variables = collections.OrderedDict([
                        # Triggers #
                        ("_trigger_singleMuon" , "Single Muon Trigger"),
                        ("_trigger_singleElectron" , "Single Electron Trigger"),
                        ("_trigger_doubleMuon" , "Double Muon Trigger"),
                        ("_trigger_doubleElectron" , "Double Electron Trigger"),
                        ("_trigger_MuonElectron" , "Muon-Electron Trigger"),
                        # MET #
                        ("_met_pt" , "MET $P_T$"),
                        ("_met_phi" , "MET $\phi$"),
                        # Leptons #
                        ("_firstlepton_pt" , "First lepton $P_T$"),
                        ("_firstlepton_eta" , "First lepton $\eta$"),
                        ("_firstlepton_phi" , "First lepton $\phi$"),
                        ("_secondlepton_pt" , "Second lepton $P_T$"),
                        ("_secondlepton_eta" , "Second lepton $\eta$"),
                        ("_secondlepton_phi" , "Second lepton $\phi$"),
                        # Dilepton #
                        ("_dilepton_pt" , "Dilepton $P_T$"),
                        ("_dilepton_deltaPhi" , "Dilepton $\Delta \phi$"),
                        ("_dilepton_deltaR" , "Dilepton $\Delta R$"),
                        ("_dilepton_invariantMass" , "Dilepton invariant mass"),
                        # Resolved Jets #
                        ("_leadjet_pt" , "Leading jet $P_T$"),
                        ("_leadjet_eta" , "Leading jet $\eta$"),
                        ("_leadjet_phi" , "Leading jet $\phi$"),
                        ("_subleadjet_pt" , "Subleading jet $P_T$"),
                        ("_subleadjet_eta" , "Subleading jet $\eta$"),
                        ("_subleadjet_phi" , "Subleading jet $\phi$"),
                        ("_leadbjet_pt" , "Leading jet $P_T$"),
                        ("_leadbjet_eta" , "Leading jet $\eta$"),
                        ("_leadbjet_phi" , "Leading jet $\phi$"),
                        ("_subleadbjet_pt" , "Subleading jet $P_T$"),
                        ("_subleadbjet_eta" , "Subleading jet $\eta$"),
                        ("_Ak4Jets_N" , "Number of Ak4 jets"),
                        # Resolved Dijet #
                        ("_dijet_pt" , "Dijet $P_T$"),
                        ("_dijet_deltaPhi" , "Dijet $\Delta \phi$"),
                        ("_dijet_deltaR" , "Dijet $\Delta R$"),
                        ("_dijet_invariantMass" , "Dijet invariant mass"),
                        ("_Ak4BJets_N" , "Number of Ak4 bjets"),
                        # Fatjet #
                        ("_fatjet_pt" , "Fatjet $P_T$"),
                        ("_fatjet_eta" , "Fatjet $\eta$"),
                        ("_fatjet_phi" , "Fatjet $\phi$"),
                        ("_fatjet_mass" , "Fatjet mass"),
                        ("_fatjet_softdropmass" , "Fatjet SoftDrop mass"),
                        ("_Ak8Jets_N" , "Number of Ak8 jets"),
                        ("_Ak8BJets_N" , "Number of Ak8 bjets"),
                        ("_boostedjets_N" , "Number of boosted bjets (fatjets)"),
                        # Highlevel variables 
                        ("_highlevelvariable_DilepMETdeltaPhi", "Dilepton-MET $\Delta \phi$"),
                        ("_highlevelvariable_DilepMETdeltapt", "Dilepton-MET $P_T$"),
                        ("_highlevelvariable_mTll", "Dilepton-MET tranverse mass"),
                        ("_highlevelvariable_mTlljj", "Dilepton-Dijet-MET tranverse mass"),
                        ("_highlevelvariable_HT2.", "$H_{T2}$ (scalar sum of the magnitudes of the HH decay momentas)"),
                        ("_highlevelvariable_HT2R.", "$H_{T2}^R$ (ratio of $H_{T2}$ and scalar sum of the transverse momenta HH decay products"),

                     #   ("_ElectronFatjet_DeltaR" , "\Delta R(Electron, Fatjet)"),
                     #   ("_ElectronJet_DeltaR" , "\Delta R(Electron, Jet)"),
                     #   ("_MuonFatjet_DeltaR" , "\Delta R(Muon, Fatjet)"),
                     #   ("_MuonJet_DeltaR" , "\Delta R(Muon, Fatjet)"),
                        ])
                # "_" at beginning is useful to distinguish leading and subleading
                # "." at end is useful to distinguish "HT2" and "HT2R"
        self.content = ""

        self.ClassifyPlots()
        self.makeHeader()
        self.makeFrames()
        self.makeFooter()
        self.dumpToFile()

    def ClassifyPlots(self):
        self.channeldict = {}
        for channel in self.channels.keys():
            logging.info("="*60)
            logging.info("Channel : %s"%channel)
            selectiondict = {} 
            for selection in self.selections:
                logging.info("... "+"-"*40)
                logging.info("... Selection : "+selection)
                variabledict = {}
                for pdf in glob.glob(os.path.join(self.dirpath,channel+"*"+selection+"*.pdf")):
                    pdf = os.path.basename(pdf)
                    log = "logy" in pdf
                    if self.logplot ^ log: # xor
                        continue
                    passed = False
                    for variable in self.variables:
                        if variable in pdf:
                            passed = True
                            var = variable
                            break
                    if passed:
                        logging.info("...... Variable : %s"%var)
                        logging.debug(pdf)
                        variabledict[var] = pdf
                selectiondict[selection] = variabledict
            self.channeldict[channel] = selectiondict

    def makeHeader(self): 
        self.content += r"""
\documentclass{beamer}
\mode<presentation> {
    \usetheme{Madrid}
    \usecolortheme{beaver}
    \setbeamertemplate{navigation symbols}{}
    \setbeamerfont{frametitle}{size=\small}
}
\usepackage{graphicx}
\usepackage[english]{babel}
\usepackage[absolute,overlay]{textpos}
\begin{document}
                        """

    def makeFooter(self):
        self.content += "\n\end{document}"

    def makeFrames(self):
        for sel,sel_name in self.selections.items():
            for var, var_name in self.variables.items():
                valid_plot = False # If empty frame, will not print it
                frame = r"\begin{frame}"
                frame += "\n\t"
                frame += r"\frametitle{Selection : %s, variable : %s}"%(sel_name,var_name) 
                frame += "\n\t"
                frame += r"\begin{figure}"
                for chan,chan_name in self.channels.items():
                    frame += "\n\t\t"
                    try:
                        frame += "\includegraphics[width=0.45\linewidth]{%s}"%(self.channeldict[chan][sel][var])
                        valid_plot = True
                    except Exception as e:
                        logging.debug('Could not find plot with channel "%s", selection "%s" and variable "%s" due to Exception "%s"'%(chan,sel,var,e))
                frame += r"""
    \end{figure}
\end{frame}
                        """
                
                if valid_plot:
                    self.content += frame
                
    def dumpToFile(self):
        with open(self.texfile,"w") as f:
            f.write(self.content)

        logging.info("Tex file written at %s"%self.texfile)
        


if __name__ == "__main__":
    # Start logging #
    logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')
    # Argparse #
    parser = argparse.ArgumentParser(description='Utility to generate Latex slides on the spot for the HH->bbWW analysis')
    parser.add_argument('-d','--dirname', action='store', required=False, type=str, default='',
                        help='Path of the directory where the plots are')
    parser.add_argument('--log', action='store_true', required=False, default=False,
                        help='Wether to use the log plots (must be present)')
    parser.add_argument('-v','--verbose', action='store_true', required=False, default=False,
                        help='Show DEGUG logging')
    opt = parser.parse_args()

    # Change verbose level if not requested #
    if not opt.verbose:
        logging.getLogger().setLevel(logging.INFO)

    instance = TemplateLatex(opt.dirname)

    

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

        self.channels = collections.OrderedDict([
                        ("ElEl", "$e^+e^-$"),
                        ("MuMu", "$\mu^+\mu^-$"),
                        ("ElMu", "$e^{\pm}\mu^{\mp}"),
                        ])
            
        self.selections = collections.OrderedDict([
                        ('ExclusiveResolvedJets_' , "Resolved Jets"),
                        ('ExclusiveResolvedJets_inclusive_' , "Resolved Jets (inclusive)"),
                        ('ExclusiveResolvedJets_bjets_' , "Resolved Jets (2 bjets)"),
                        ('ExclusiveResolvedJets_mixed_' , "Resolved Jets (bjet + lightjet)"),
                        ('BoostedJets_' , "Boosted Jets"),
                        ])
        self.variables = collections.OrderedDict([
                        ("_firstlepton_pt" , "First lepton $P_T$"),
                        ("_firstlepton_eta" , "First lepton $\eta$"),
                        ("_firstlepton_phi" , "First lepton $\phi$"),
                        ("_secondlepton_pt" , "Second lepton $P_T$"),
                        ("_secondlepton_eta" , "Second lepton $\eta$"),
                        ("_secondlepton_phi" , "Second lepton $\phi$"),
                        ("_dilepton_invariantMass" , "Dilepton invariant mass"),
                        ("_leadjet_pt" , "Leading jet $P_T$"),
                        ("_leadjet_eta" , "Leading jet $\eta$"),
                        ("_leadjet_phi" , "Leading jet $\phi$"),
                        ("_subleadjet_pt" , "Subleading jet $P_T$"),
                        ("_subleadjet_eta" , "Subleading jet $\eta$"),
                        ("_subleadjet_phi" , "Subleading jet $\phi$"),
                        ("_bjet_pt" , "Bjet $P_T$"),
                        ("_bjet_eta" , "Bjet $\eta$"),
                        ("_bjet_phi" , "Bjet $\phi$"),
                        ("_lightjet_pt" , "Lightjet $P_T$"),
                        ("_lightjet_eta" , "Lightjet $\eta$"),
                        ("_lightjet_phi" , "Lightjet $\phi$"),
                        ("_fatjet_pt" , "Fatjet $P_T$"),
                        ("_fatjet_eta" , "Fatjet $\eta$"),
                        ("_fatjet_phi" , "Fatjet $\phi$"),
                        ("_fatjet_mass" , "Fatjet mass"),
                        ("_fatjet_softdropmass" , "Fatjet SoftDrop mass"),
                        ("_dijet_invariantMass" , "Dijet invariant mass"),
                     #   ("_ElectronFatjet_DeltaR" , "\Delta R(Electron, Fatjet)"),
                     #   ("_ElectronJet_DeltaR" , "\Delta R(Electron, Jet)"),
                     #   ("_MuonFatjet_DeltaR" , "\Delta R(Muon, Fatjet)"),
                     #   ("_MuonJet_DeltaR" , "\Delta R(Muon, Fatjet)"),
                        ])
                # "_" at beginning is useful to distinguish leading and subleading
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
}
\usepackage{graphicx}
\usepackage[english]{babel}
\begin{document}
                        """

    def makeFooter(self):
        self.content += "\n\end{document}"

    def makeFrames(self):
        for sel,sel_name in self.selections.items():
            for var, var_name in self.variables.items():
                valid_plot = False # If empty frame, will not print it
                frame = ""
                frame += "\n"
                frame += "%"+"-"*80
                frame += "\n"
                frame += r"\begin{frame}"
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
                frame += "\n\t"
                frame += r"\end{figure}"
                frame += "\n"
                frame += r"\end{frame}"
                frame += "\n"
                
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

    

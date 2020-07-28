import os
import sys
import re
import glob
import collections
import logging
import argparse
import subprocess
from functools import reduce
import operator

class TemplateLatex:
    def __init__(self,dirpaths,logplot=False):
        self.dirpaths = [os.path.abspath(dirpath) for dirpath in dirpaths]
        self.texfile = os.path.join(os.path.dirname(self.dirpaths[0]),"plots.tex")
        self.filesdict = {}
        self.logplot = logplot

        # Ordered dicts so that plots keep that order #
        self.channels = collections.OrderedDict([
                        ("ElEl", "$e^+e^-$"),
                        ("MuMu", "$\mu^+\mu^-$"),
                        ("ElMu", "$e^{\pm}\mu^{\mp}"),
                        ])
            
        self.selections = collections.OrderedDict([
                        ('TwoAk4Jets_' , "At least two Ak4 Jets (before Btagging)"),
                        ('OneAk8Jet_' , "At least one Ak8 Jet (before Btagging)"),
                        ('TwoAk4JetsExclusiveResolvedNoBtag' , "Exclusive Resolved Jets (no Btag)"),
                        ('TwoAk4JetsExclusiveResolvedOneBtag' , "Exclusive Resolved Jets (1 Btag)"),
                        ('TwoAk4JetsExclusiveResolvedTwoBtags' , "Exclusive Resolved Jets (2 Btags)"),
                        ('OneAk8JetInclusiveBoosted' , "Inclusive Boosted Jets"),
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
                        ("_highlevelvariable_DilepJetInvariantMass", "Dilepton + Jets invariant mass"),
                        ("_highlevelvariable_DilepMETpt", "Dilepton-MET $P_T$"),
                        ("_highlevelvariable_DilepMETdeltaPhi", "Dilepton-MET $|\Delta \phi|$"),
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
                pdfs = [glob.glob(os.path.join(dirpath,channel+"*"+selection+"*.pdf")) for dirpath in self.dirpaths]
                pdfs = reduce(operator.concat, pdfs)
                for pdf in pdfs:
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
\graphicspath{
                        """
        for dirpath in self.dirpaths:
            self.content += "{%s/},"%os.path.basename(dirpath)
        self.content += "}"
        self.content += r"""
\begin{document}
                        """
        self.content += "\n"

    def makeFooter(self):
        self.content += "\n\end{document}\n"

    def makeFrames(self):
        for sel,sel_name in self.selections.items():
            for var, var_name in self.variables.items():
                valid_plot = False # If empty frame, will not print it
                frame = r"\begin{frame}"
                frame += "\n\t"
                frame += r"\frametitle{Selection : %s, variable : %s}"%(sel_name,var_name) 
                frame += "\n\t"
                frame += r"\begin{figure}"
                for i,(chan,chan_name) in enumerate(self.channels.items()):
                    frame += "\n\t\t"
                    try:
                        frame += "\includegraphics[width=0.45\linewidth]{%s}"%(self.channeldict[chan][sel][var])
                        if i == 1:
                            frame += r"\\" # Back to line after two first plots
                        valid_plot = True
                    except Exception as e:
                        logging.debug('Could not find plot with channel "%s", selection "%s" and variable "%s" due to Exception "%s"'%(chan,sel,var,e))
                frame += r"""
    \end{figure}
\end{frame}
                        """
                frame += "\n%"+"-"*80+"\n"
                
                if valid_plot:
                    self.content += frame
                
    def dumpToFile(self):
        with open(self.texfile,"w") as f:
            f.write(self.content)

        logging.info("Tex file written at %s"%self.texfile)

    def producePDF(self):
        cwd = os.getcwd()
        wd = os.path.dirname(self.texfile)
        try:
            os.chdir(wd)
        except:
            logging.critical("Could not go to directory %s, will abort"%wd)
            sys.exit()
        try:
            with open('plots_compilation.log', 'w') as output:
                process= subprocess.Popen(['pdflatex', os.path.basename(self.texfile)], # We moved to the directory
                                          stdout=output,
                                          stderr=output)
        except:
            logging.critical("Could not produce the pdf from the Tex file, try manually 'pdflatex %s'"%self.texfile)
        logging.info("Generated pdf file %s"%(self.texfile.replace(".tex",".pdf")))
        logging.info("Log available at %s"%(os.path.join(wd,'plots_compilation.log')))
        os.chdir(cwd)
        
    def compileYield(self):
        # Load tex file and edit content in string #
        text = r"""
\documentclass{report}
    \usepackage{booktabs}
    \usepackage{graphicx}
    \usepackage[a4paper,bindingoffset=0.5cm,left=0cm,right=1cm,top=2cm,bottom=2cm,footskip=0.25cm]{geometry}

\begin{document}
                """
        wd = os.path.dirname(self.texfile)
        for dirpath in self.dirpaths:
            with open(os.path.join(dirpath,"yields.tex"),"r") as f:
                while True: 
                    line = f.readline() 
                    if not line:
                        break
                    if line.startswith(r"\begin{document}"):
                        continue
                    if line.startswith(r"\begin{tabular}"):
                        text += r"\resizebox{\textwidth}{!}{"+"\n"
                    text += line
        text += "\n}\n"+r"\end{document}"

        # Save new tex file #
        with open(os.path.join(wd,"yieldsTable.tex"),"w") as f:
            f.write(text)
        # Move to directory and compile #
        cwd = os.getcwd()
        try:
            os.chdir(wd)
        except:
            logging.critical("Could not go to directory %s"%wd)
        try:
            with open('yields_compilation.log', 'w') as output:
                process= subprocess.Popen(['pdflatex', "yieldsTable.tex"], # We moved to the directory
                                          stdout=output,
                                          stderr=output)
        except:
            logging.critical("Could not produce the pdf from the Tex file, try manually 'pdflatex yieldsTable.tex'")
        logging.info("Generated pdf file %s"%(os.path.join(wd,"yieldsTable.pdf")))
        logging.info("Log available at %s"%(os.path.join(wd,'yields_compilation.log')))
        os.chdir(cwd)


             
    
if __name__ == "__main__":
    # Start logging #
    logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')
    # Argparse #
    parser = argparse.ArgumentParser(description='Utility to generate Latex slides on the spot for the HH->bbWW analysis')
    parser.add_argument('-d','--dirnames', action='store', required=False, type=str, default='', nargs='+',
                        help='Path of the directory where the plots are (can be several separated by spaces')
    parser.add_argument('-y','--yields', action='store_true', required=False, default=False,
                        help='Wether to compile the yields.tex into yields.pdf')
    parser.add_argument('--pdf', action='store_true', required=False, default=False,
                        help='Wether to produce the PDF ')
    parser.add_argument('--log', action='store_true', required=False, default=False,
                        help='Wether to use the log plots (must be present)')
    parser.add_argument('-v','--verbose', action='store_true', required=False, default=False,
                        help='Show DEGUG logging')
    opt = parser.parse_args()

    # Change verbose level if not requested #
    if not opt.verbose:
        logging.getLogger().setLevel(logging.INFO)

    instance = TemplateLatex(opt.dirnames)
    if opt.yields:
        instance.compileYield()
    if opt.pdf:
        instance.producePDF()

    

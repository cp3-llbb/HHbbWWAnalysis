import os
import sys
import re
import glob
import collections
import logging
import argparse
import subprocess

class TemplateLatex:
    def __init__(self,dirpath,logplot=False):
        self.dirpath = os.path.abspath(dirpath)
        self.texfile = os.path.join(self.dirpath,"plots.tex")
        self.filesdict = {}
        self.logplot = logplot

        # Ordered dicts so that plots keep that order #
        self.channels = collections.OrderedDict([
                        ("El", "$e^{\pm}$"),
                        ("Mu", "$\mu^{\pm}$"),
                        ])
            
        self.selections = collections.OrderedDict([
<<<<<<< HEAD
                        ('Ak4JetsLooseExclusiveResolved1b2j' , "At least three Ak4 Jets (One is b-jet)"),
                        ('Ak4JetsLooseExclusiveResolved2b1j' , "At least three Ak4 Jets (Two is b-jet)"),
                        ('Ak4JetsTightExclusiveResolved1b3j' , "At least four Ak4 Jets (one is b-jet)"),
                        ('Ak4JetsTightExclusiveResolved2b2j' , "At least four Ak4 Jets (Two is b-jet)"),
                        ('Ak8BJetsHbbBoostedWtoJ' , "At least one Ak8 Jet with b-tag and another Ak4 Jet"),
                        ('Ak8BJetsHbbBoostedWtoJJ' , "At least one Ak8 Jet with b-tag and another 2 Ak4 Jets"),
=======
                        ('TwoAk4Jets_' , "At least two Ak4 Jets (before Btagging)"),
                        ('OneAk8Jet_' , "At least one Ak8 Jet (before Btagging)"),
                        ('TwoAk4JetsExclusiveResolvedNoBtag' , "Exclusive Resolved Jets (no Btag)"),
                        ('TwoAk4JetsExclusiveResolvedOneBtag' , "Exclusive Resolved Jets (1 Btag)"),
                        ('TwoAk4JetsExclusiveResolvedTwoBtags' , "Exclusive Resolved Jets (2 Btags)"),
                        ('OneAk8JetInclusiveBoosted' , "Inclusive Boosted Jets"),
>>>>>>> 351efd47214d22a9f931b61cf7c288851194d8c8
                        ])
        self.variables = collections.OrderedDict([
                        # Triggers #
                        ("_trigger_singleMuon" , "Single Muon Trigger"),
                        ("_trigger_singleElectron" , "Single Electron Trigger"),
                        # MET #
                        ("_met_pt" , "MET $P_T$"),
                        ("_met_phi" , "MET $\phi$"),
                        # Leptons #
                        ("_lepton_pt" , "Lepton $P_T$"),
                        ("_lepton_eta" , "Lepton $\eta$"),
                        ("_lepton_phi" , "Lepton $\phi$"),
                        # Single Jets #
                        ("_eta_leadAk4B" , "Leading Ak4 b-jet $\eta$"),
                        ("_eta_subLeadAk4B" , "Sub-leading Ak4 b-jet $\eta$"),
                        ("_eta_leadAk4" , "Leading Ak4 jet $\eta$"),
                        ("_eta_subLeadAk4" , "Sub-leading Ak4 $\eta$"),
                        ("_phi_leadAk4B" , "Leading Ak4 b-jet $\phi$"),
                        ("_phi_subLeadAk4B" , "Sub-leading Ak4 b-jet $\phi$"),
                        ("_phi_leadAk4" , "Leading Ak4 jet $\phi$"),
                        ("_phi_subLeadAk4" , "Sub-leading Ak4 $\phi$"),
                        ("_pT_leadAk4B" , "Leading Ak4 b-jet $P_T$"),
                        ("_pT_subLeadAk4B" , "Sub-leading Ak4 b-jet $P_T$"),
                        ("_pT_leadAk4" , "Leading Ak4 jet $P_T$"),
                        ("_pT_subLeadAk4" , "Sub-leading Ak4 $P_T$"),
                        ("_Ak4BJets_N" , "No of b-tagged Ak4 jet"),
                        ("_Ak4LightJets_N" , "No of all Ak4 jets"),
                        # Dijet Correlation#
                        ("_DeltaPhi_leadAk4B_subLeadAk4B" , "$\Delta\phi$ of lead and sublead b-jet"),
                        ("_DeltaPhi_leadAk4B_leadAk4" , "$\Delta\phi$ of lead b-jet and lead Ak4 jet"),
                        ("_DeltaPhi_leadAk4B_subLeadAk4" , "$\Delta\phi$ of lead b-jet and sublead Ak4 jet"),
                        ("_DeltaPhi_subLeadAk4B_leadAk4" , "$\Delta\phi$ of sublead b-jet and lead Ak4 jet"),
                        ("_DeltaPhi_leadAk4_subLeadAk4" , "$\Delta\phi$ of lead and sublead Ak4 jet"),
                        ("_DeltaR_leadAk4B_subLeadAk4B" , "$\Delta$R of lead and sublead b-jet"),
                        ("_DeltaR_leadAk4B_leadAk4" , "$\Delta$R of lead b-jet and lead Ak4 jet"),
                        ("_DeltaR_leadAk4B_subLeadAk4" , "$\Delta$R of lead b-jet and sublead Ak4 jet"),
                        ("_DeltaR_subLeadAk4B_leadAk4" , "$\Delta$R of sublead b-jet and lead Ak4 jet"),
                        ("_DeltaR_leadAk4_subLeadAk4" , "$\Delta$R of lead and sublead Ak4 jet"),
                        ("_DiJetPT_leadAk4B_subLeadAk4B" , "Dijet $P_T$ of lead and sublead b-jet"),
                        ("_DiJetPT_leadAk4B_leadAk4" , "Dijet $P_T$ of lead b-jet and lead Ak4 jet"),
                        ("_DiJetPT_leadAk4_subLeadAk4" , "Dijet $P_T$ of lead and sublead Ak4 jet"),
                        ("_InvariantMass_leadAk4B_leadAk4" , "Invariant mass of lead b-jet and lead Ak4 jet"),
                        ("_InvariantMass_leadAk4B_subLeadAk4B" , "Invariant mass of lead and sublead b-jet"),
                        ("_InvariantMass_leadAk4_subLeadAk4" , "Invariant mass of lead and sublead Ak4 jet"),
                        # Fatjet #
                        ("_pT_Ak8bJet" , "Fatjet $P_T$"),
                        ("_eta_Ak8bJet" , "Fatjet $\eta$"),
                        ("_phi_Ak8bJet" , "Fatjet $\phi$"),
                        ("_mass_Ak8bJet" , "Fatjet mass"),
                        ("_massSD_Ak8bJet" , "Fatjet SoftDrop mass"),
                        ("_Ak8Jets_N" , "Number of Ak8 jets"),
                        ("_Ak8BJets_N" , "Number of Ak8 bjets"),
                        ("_DeltaPhi_Ak8bJet_Ak4Jet" , "$\Delta\phi$ between Ak8 and Ak4 jet"),
                        ("_DijetPT_Ak8bJet_Ak4Jet" , "Dijet $P_T$ of Ak8 and Ak4 jet"),
                        ("_DeltaR_Ak8bJet_Ak4Jet" , "$\Delta$R between Ak8 and Ak4 jet"),
                        # Highlevel variables 
                        ("_highlevelvariable_DeltaPhi_lep_leadAk4B" , "$\Delta\phi$ of lepton and lead Ak4 b-jet"),
                        ("_highlevelvariable_DeltaPhi_lep_leadAk4" , "$\Delta\phi$ of lepton and lead Ak4 jet"),
                        ("_highlevelvariable_DeltaPhi_lep_subLeadAk4" , "$\Delta\phi$ of lepton and sublead Ak4 jet"),
                        ("_highlevelvariable_DeltaPhi_lep_subLeadAk4B" , "$\Delta\phi$ of lepton and sublead Ak4 b-jet"),
                        ("_highlevelvariable_DeltaR_lep_leadAk4B" , "$\Delta$R of lepton and lead Ak4 b-jet"),
                        ("_highlevelvariable_DeltaR_lep_leadAk4" , "$\Delta$R of lepton and lead Ak4 jet"),
                        ("_highlevelvariable_DeltaR_lep_subLeadAk4" , "$\Delta$R of lepton and sublead Ak4 jet"),
                        ("_highlevelvariable_DeltaR_lep_subLeadAk4B" , "$\Delta$R of lepton and sublead Ak4 b-jet"),
                        ("_highlevelvariable_HT2_l3jmet", "$H_{T2}$ (scalar sum of the magnitudes of the HH decay momentas)"),
                        ("_highlevelvariable_HT2R_l3jmet", "$H_{T2}^R$ (ratio of $H_{T2}$ and scalar sum of the transverse momenta HH decay products"),

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
                        frame += "\includegraphics[width=0.49\linewidth]{%s}"%(self.channeldict[chan][sel][var])
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
        try:
            os.chdir(self.dirpath)
            try:
                with open('plots_compilation.log', 'w') as output:
                    process= subprocess.Popen(['pdflatex', os.path.basename(self.texfile)], # We moved to the directory
                                              stdout=output,
                                              stderr=output)
            except:
                logging.critical("Could not produce the pdf from the Tex file, try manually 'pdflatex %s'"%self.texfile)
            logging.info("Generated pdf file %s"%(self.texfile.replace(".tex",".pdf")))
            logging.info("Log available at %s"%(os.path.join(self.dirpath,'plots_compilation.log')))
        except:
            logging.critical("Could not go to directory %s"%self.dirpath)
        os.chdir(cwd)
        
    def compileYield(self):
        # Load tex file adn edit content in string #
        text = r"""
\documentclass{report}
    \usepackage{booktabs}
    \usepackage{graphicx}
    \usepackage[a4paper,bindingoffset=0.5cm,left=0cm,right=1cm,top=2cm,bottom=2cm,footskip=0.25cm]{geometry}
                """
<<<<<<< HEAD
        with open(os.path.join(self.dirpath,"yields.tex"),"r") as f:
            while True: 
                line = f.readline() 
                if not line:
                    break
                if line.startswith(r"\begin{tabular}"):
                    text += r"\begin{document}"+"\n"
                    text += r"\resizebox{\textwidth}{!}{"+"\n"
                text += line
        text += "\n}\n"+r"\end{document}"
=======
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
            text += "\n}\n"
        text += r"\end{document}"
>>>>>>> 351efd47214d22a9f931b61cf7c288851194d8c8

        # Save new tex file #
        with open(os.path.join(self.dirpath,"yieldsTable.tex"),"w") as f:
            f.write(text)
        # Move to directory and compile #
        cwd = os.getcwd()
        try:
            os.chdir(self.dirpath)
            try:
                with open('yields_compilation.log', 'w') as output:
                    process= subprocess.Popen(['pdflatex', "yieldsTable.tex"], # We moved to the directory
                                              stdout=output,
                                              stderr=output)
            except:
                logging.critical("Could not produce the pdf from the Tex file, try manually 'pdflatex yieldsTable.tex'")
            logging.info("Generated pdf file %s"%(os.path.join(self.dirpath,"yieldsTable.pdf")))
            logging.info("Log available at %s"%(os.path.join(self.dirpath,'yields_compilation.log')))
        except:
            logging.critical("Could not go to directory %s"%self.dirpath)
        os.chdir(cwd)


             
    
if __name__ == "__main__":
    # Start logging #
    logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')
    # Argparse #
    parser = argparse.ArgumentParser(description='Utility to generate Latex slides on the spot for the HH->bbWW analysis')
    parser.add_argument('-d','--dirname', action='store', required=False, type=str, default='',
                        help='Path of the directory where the plots are')
    parser.add_argument('-y','--yields', action='store_true', required=False, default=True,
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

    instance = TemplateLatex(opt.dirname)
    if opt.yields:
        instance.compileYield()
    if opt.pdf:
        instance.producePDF()

    

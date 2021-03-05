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
        if self.logplot:
            self.texfile = self.texfile.replace('.tex','_log.tex')

        # Ordered dicts so that plots keep that order #
        self.channels = collections.OrderedDict([
                        ("HH_cat_e_", "$e^{\pm}$"),
                        ("HH_cat_m_", "$\mu^{\pm}$"),
                        ("HH_cat_ee_", "$e^+e^-$"),
                        ("HH_cat_mm_", "$\mu^+\mu^-$"),
                        ("HH_cat_em_", "$e^{\pm}\mu^{\mp}"),
                        ("HH_cat_ss_", "$e^+e^- + \mu^+\mu^-$"),
                        ])
            
        self.selections = collections.OrderedDict([
                        ('alljet_' , "Resolved (1 btag + 2 btag) + Boosted (1 btag)"),
                        ('resolved_' , "Exclusive Resolved Jets (1 btag + 2 btag)"),
                        ('boosted_' , "Inclusive Boosted Jets (1 btag)"),
                        ('resolved0b_' , "Exclusive Resolved Jets (0 btag)"),
                        ('resolved1b_' , "Exclusive Resolved Jets (1 btag)"),
                        ('resolved2b_' , "Exclusive Resolved Jets (2 btag)"),
                        ('boosted0b_' , "Inclusive Boosted Jets (0 btag)"),
                        ('boosted1b_' , "Inclusive Boosted Jets (1 btag)"),
                        ('resolved2b2Wj_' , "Exclusive Resolved Jets 2b2Wj"),
                        ('resolved2b1Wj_' , "Exclusive Resolved Jets 2b1Wj"),
                        ('resolved2b0Wj_' , "Exclusive Resolved Jets 2b0Wj"),
                        ('resolved1b2Wj_' , "Exclusive Resolved Jets 1b2Wj"),
                        ('resolved1b1Wj_' , "Exclusive Resolved Jets 1b1Wj"),
                        ('resolved1b0Wj_' , "Exclusive Resolved Jets 1b0Wj"),
                        ('resolved0b0_' , "Exclusive Resolved Jets 0b"),
                        ])
        self.variables = collections.OrderedDict([
                        ("VBFnode", "DNN VBF node"),
                        ("GGFnode", "DNN GGF node"),
                        ("TTnode", "DNN TT node"),
                        ("DYnode", "DNN DY node"),
                        ("STnode", "DNN ST node"),
                        ("TTVXnode", "DNN TTVX node"),
                        ("Hnode", "DNN H node"),
                        ("VVVnode", "DNN VVV node"),
                        ("Rarenode", "DNN Rare node"),
                        ("WJetsnode", "DNN WJets node"),
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
                #if self.logplot:
                base_pdf = channel+"*"+selection+"*.pdf"
                if self.logplot:
                    base_pdf = base_pdf.replace('.pdf','_logy.pdf')
                pdfs = [glob.glob(os.path.join(dirpath,base_pdf)) for dirpath in self.dirpaths]
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
    \setbeamerfont{frametitle}{size=\footnotesize}
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
             
    
if __name__ == "__main__":
    # Start logging #
    logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')
    # Argparse #
    parser = argparse.ArgumentParser(description='Utility to generate Latex slides on the spot for the HH->bbWW analysis')
    parser.add_argument('-d','--dirnames', action='store', required=False, type=str, default='', nargs='+',
                        help='Path of the directory where the plots are (can be several separated by spaces')
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

    instance = TemplateLatex(opt.dirnames,opt.log)
    if opt.pdf:
        instance.producePDF()

    

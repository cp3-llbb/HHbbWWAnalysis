# HHbbWWAnalysis

## Code organization 
- BaseHHtobbWW.py : contains all the objects selections (leptons, dileptons, Ak4 jets, Ak8 jets and Btagging), the trigger selections and all the corrections (JEC,JER,MET,Rochester)
- SkimmerHHtobbWW.py : inherits the object selection from BaseHHtobbWW.py, produces a skim for the synchronization
- PlotterHHtobbWW.py : inherits the object selection from BaseHHtobbWW.py, produces plots at differents stages of the event selection
- plotDef.py : contains all the plots definitions, is called from PlotterHHtobbWW.py
- METScripts.py : contains the scripts to filter and correct the MET
- scalefactorsutils.py : utils to build the whole SF dict
- scalefactorsbbWW.py : class that build all the SF paths (saved in full dict with scalefactorsutils.py) and produces SF passed to PlotterHHtobbWW.py
- HHtobbWW.py : deprecated version of the code, kept for comparison

## How to run 

### Skimmer

If you want to produce the skim for the synchronization, the command is :

bambooRun -m SkimmerHHtobbWW.py:SkimmerNanoHHtobbWW analysis2016_synchro.yml -o Synchronization --outputTreeName syncTree 

-> Will produce a tree with the variables defined within the script and name "syncTree"

### Plotter 

The plotter can be run with :

bambooRun (--distributed=driver) -m PlotterHHtobbWW.py:PlotterNanoHHtobbWW analysis2016.yml -o folderName

Note : like that nothing will be produced, you need to specify the lepton selection 
Arguments :
- `--Preselected` : will use the preselected leptons to build dilepton pairs (and procede with the highest-Pt one)
- `--Fakeable` : will use the fakeable leptons to build dilepton pairs (and procede with the highest-Pt one), additionaly will perform several cuts for the fakeable selection
- `--Tight`: will use the tight leptons to build dilepton pairs (and procede with the highest-Pt one), the fakeable dilepton pair needs to pass the tight criteria for both leptons (and if MC be gen latched)
- `FakeExtrapolation` : the dilepton pairs are the one that failed the tight criterium for at least one of the lepton, but also they need both to be gen matched (apart from that condition : Tight + FakeExtrapolation = Fakeable)
- `OnlyYield`: Additional parameter to only plot the yield and no other plots






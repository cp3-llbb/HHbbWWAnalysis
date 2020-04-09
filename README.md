# HHbbWWAnalysis

## Code organization 
- BaseHHtobbWW.py : contains all the objects definition (leptons, dileptons, Ak4 jets, Ak8 jets and Btagging), the trigger selections and all the corrections (JEC,JER,MET,Rochester, SF)
- SkimmerHHtobbWW.py : inherits the object selection from BaseHHtobbWW.py, produces a skim for the synchronization
- PlotterHHtobbWW.py : inherits the object selection from BaseHHtobbWW.py, produces plots at differents stages of the event selection
- plotDef.py : contains all the plots definitions, is called from PlotterHHtobbWW.py
- selectionDef.py : contains all the selection, uses a class object encapsulating the selection
- METScripts.py : contains the scripts to filter and correct the MET
- scalefactorsutils.py : utils to build the whole SF dict
- scalefactorsbbWW.py : class that build all the SF paths (saved in full dict with scalefactorsutils.py) and produces SF passed to PlotterHHtobbWW.py
- HHtobbWW.py : deprecated version of the code, kept for comparison

## Development rules 
It is proposed that the different part  object definition must be done in very specific parts of the framework :
- Object selection (using combine and select) should only be done in BaseHHtobbWW.py so that they are all available for the dauther classes
- All the selections should be done in selectionDef.py and encapsulated in SelectionObject class objects. This was the workflow is more easily manageable and the selection steps can be called either from the Plotter or the Skimmer
- All plots must be defined in plotDef.py with functions that return a list of plots. This is to make sure that all the plots are easily accounted for
- The Plotter will call the selection and plot functions in a sensible manner (aka, define the plot right after the selection and not at the end, this is to optimize the RooDataFrame workflow)
- The Skimmer will call the selection functions and define the variables dictionnary
- All arguments must be defined in BaseHHtobbWW.py to be available throughout the framework)

## How to run 

### Arguments 
Use --help for more details
Lepton arguments:
- `--Preselected` : will use the preselected leptons to build dilepton pairs (and procede with the highest-Pt one)
- `--Fakeable` : will use the fakeable leptons to build dilepton pairs (and procede with the highest-Pt one), additionaly will perform several cuts for the fakeable selection
- `--Tight`: will use the tight leptons to build dilepton pairs (and procede with the highest-Pt one), the fakeable dilepton pair needs to pass the tight criteria for both leptons (and if MC be gen latched)
- `--FakeExtrapolation` : the dilepton pairs are the one that failed the tight criterium for at least one of the lepton, but also they need both to be gen matched (apart from that condition : Tight + FakeExtrapolation = Fakeable)
Jet arguments :
- `--Ak4` : Selection with at least 2 Ak4 jets 
- `--Ak8` : Selection with at least 1 Ak8 jet
- `--Resolved0Btag` : Exclusive resolved selection with no btagged Ak4 jets
- `--Resolved1Btag` : Exclusive resolved selection with one btagged Ak4 jets
- `--Resolved2Btag` : Exclusive resolved selection with two btagged Ak4 jets
- `--Boosted` : Inclusive boosted selection with one Ak8 jets having a least one btagged subjet
Plotter arguments :
- `--OnlyYield`: Additional parameter to only plot the yield and no other plots
Skimmer arguments :
- --Synchronization : will produce the syncTree without triggers, flags or corrections
- --Channel : for ntuple production

### Skimmer

If you want to produce the skim for the synchronization, the command is :
```
bambooRun (--distributed=driver) -m SkimmerHHtobbWW.py:SkimmerNanoHHtobbWW analysis2016_synchro.yml -o folderName --outputTreeName treeName  (--Synchronization) (lepton arguments) (jet arguments) (--Channel)
```

To produce sync tree : Only specify `--Synchronization`, other arguments will be ignored

To produce analysis ntuples, note that only one selection at a time can be used so far
-> Must specify :
  - One of the lepton argument
  - One of the jet variable 
  - A channel : --Channel [ElEl/MuMu/ElMu]
Enforced in the script to have one and only one for each

### Plotter 

The plotter can be run with :
```
bambooRun (--distributed=driver) -m PlotterHHtobbWW.py:PlotterNanoHHtobbWW analysis2016.yml -o folderName (lepton arguments) (jet arguments) 
```

Note : If no arguments is provided, will plot all jet and lepton selection (**this means a lot of histograms** : set a high enough time and memory, each systematic variation will produce many histogram)
One or more of each argument can be used






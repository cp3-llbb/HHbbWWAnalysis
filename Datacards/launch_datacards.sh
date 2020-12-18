#!/bin/bash

python produceDataCards.py --yaml dataset_fit_TTHIDLoose_DNN08_2016.yml
python produceDataCards.py --yaml dataset_fit_TTHIDLoose_DNN08_2017.yml
python produceDataCards.py --yaml dataset_fit_TTHIDLoose_DNN08_2018.yml
python produceDataCards.py --yaml dataset_fit_TTHIDLoose_ZPeak_DNN08_2016.yml
python produceDataCards.py --yaml dataset_fit_TTHIDLoose_ZPeak_DNN08_2017.yml
python produceDataCards.py --yaml dataset_fit_TTHIDLoose_ZPeak_DNN08_2018.yml

plotIt -i /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_DNN08_2016/plotit -o /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_DNN08_2016/plotit/plots_2016 -e 2016 /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_DNN08_2016/plotit/plots.yml
plotIt -i /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_DNN08_2017/plotit -o /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_DNN08_2017/plotit/plots_2017 -e 2017 /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_DNN08_2017/plotit/plots.yml
plotIt -i /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_DNN08_2018/plotit -o /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_DNN08_2018/plotit/plots_2018 -e 2018 /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_DNN08_2018/plotit/plots.yml

plotIt -i /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_ZPeak_DNN08_2016/plotit -o /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_ZPeak_DNN08_2016/plotit/plots_2016 -e 2016 /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_ZPeak_DNN08_2016/plotit/plots.yml
plotIt -i /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_ZPeak_DNN08_2017/plotit -o /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_ZPeak_DNN08_2017/plotit/plots_2017 -e 2017 /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_ZPeak_DNN08_2017/plotit/plots.yml
plotIt -i /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_ZPeak_DNN08_2018/plotit -o /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_ZPeak_DNN08_2018/plotit/plots_2018 -e 2018 /home/users/f/b/fbury/bamboodev/HHbbWWAnalysis/Datacards/dataset_fit_TTHIDLoose_ZPeak_DNN08_2018/plotit/plots.yml


python compileDatacardPlots.py --pdf -d dataset_fit_TTHIDLoose_DNN08_2016/plotit/plots_2016
python compileDatacardPlots.py --pdf -d dataset_fit_TTHIDLoose_DNN08_2017/plotit/plots_2017
python compileDatacardPlots.py --pdf -d dataset_fit_TTHIDLoose_DNN08_2018/plotit/plots_2018
python compileDatacardPlots.py --pdf -d dataset_fit_TTHIDLoose_ZPeak_DNN08_2016/plotit/plots_2016
python compileDatacardPlots.py --pdf -d dataset_fit_TTHIDLoose_ZPeak_DNN08_2017/plotit/plots_2017
python compileDatacardPlots.py --pdf -d dataset_fit_TTHIDLoose_ZPeak_DNN08_2018/plotit/plots_2018

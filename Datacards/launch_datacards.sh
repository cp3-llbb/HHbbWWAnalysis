#!/bin/bash

python produceDataCards.py --yaml dataset_fit_TTHIDLoose_DNN11_mc_noSyst_2017.yml &
python produceDataCards.py --yaml dataset_fit_TTHIDLoose_DNN11_mc_noSyst_2018.yml &
python produceDataCards.py --yaml dataset_fit_TTHIDLoose_DNN11_datadriven_noSyst_2017.yml &
python produceDataCards.py --yaml dataset_fit_TTHIDLoose_DNN11_datadriven_noSyst_2018.yml &

python produceDataCards.py --yaml dataset_fit_TTHIDLoose_DNN11_mc_2017.yml &
python produceDataCards.py --yaml dataset_fit_TTHIDLoose_DNN11_mc_2018.yml &
python produceDataCards.py --yaml dataset_fit_TTHIDLoose_DNN11_datadriven_2017.yml &
python produceDataCards.py --yaml dataset_fit_TTHIDLoose_DNN11_datadriven_2018.yml &

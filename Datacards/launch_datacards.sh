#!/bin/bash

python produceDataCards.py --yaml dataset_DYCR_boosted_2016.yml  &
python produceDataCards.py --yaml dataset_DYCR_boosted_2017.yml  &
python produceDataCards.py --yaml dataset_DYCR_boosted_2018.yml  &
python produceDataCards.py --yaml dataset_DYCR_resolved_2016.yml &
python produceDataCards.py --yaml dataset_DYCR_resolved_2017.yml &
python produceDataCards.py --yaml dataset_DYCR_resolved_2018.yml &
python produceDataCards.py --yaml dataset_SR_for_DY_boosted_2016.yml    &
python produceDataCards.py --yaml dataset_SR_for_DY_boosted_2017.yml    &
python produceDataCards.py --yaml dataset_SR_for_DY_boosted_2018.yml    &
python produceDataCards.py --yaml dataset_SR_for_DY_resolved_2016.yml   &
python produceDataCards.py --yaml dataset_SR_for_DY_resolved_2017.yml   &
python produceDataCards.py --yaml dataset_SR_for_DY_resolved_2018.yml   &
python produceDataCards.py --yaml dataset_TTCR_2016.yml          &

python produceDataCards.py --yaml dataset_fit_POGID_2016.yml &

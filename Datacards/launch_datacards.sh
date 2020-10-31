#!/bin/bash

python produceDataCards.py --yaml dataset_DYCR_boosted_2016.yml  &
python produceDataCards.py --yaml dataset_DYCR_boosted_2017.yml  &
python produceDataCards.py --yaml dataset_DYCR_boosted_2018.yml  &
python produceDataCards.py --yaml dataset_DYCR_resolved_2016.yml &
python produceDataCards.py --yaml dataset_DYCR_resolved_2017.yml &
python produceDataCards.py --yaml dataset_DYCR_resolved_2018.yml &
python produceDataCards.py --yaml dataset_SR_boosted_2016.yml    &
python produceDataCards.py --yaml dataset_SR_boosted_2017.yml    &
python produceDataCards.py --yaml dataset_SR_boosted_2018.yml    &
python produceDataCards.py --yaml dataset_SR_resolved_2016.yml   &
python produceDataCards.py --yaml dataset_SR_resolved_2017.yml   &
python produceDataCards.py --yaml dataset_SR_resolved_2018.yml   &
python produceDataCards.py --yaml dataset_TTCR_2016.yml          &

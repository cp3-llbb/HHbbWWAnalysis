#!/bin/bash
era=$1
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_1_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_2_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_3_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_4_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_5_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_6_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_7_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_8_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_8a_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_9_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_10_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_11_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_12_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_1_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_2_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_3_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_4_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_5_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_6_$era.yml
cp datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_SM_$era.yml datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_7_$era.yml

sed -i 's/SM/1/g'        datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_1_$era.yml
sed -i 's/SM/2/g'        datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_2_$era.yml
sed -i 's/SM/3/g'        datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_3_$era.yml
sed -i 's/SM/4/g'        datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_4_$era.yml
sed -i 's/SM/5/g'        datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_5_$era.yml
sed -i 's/SM/6/g'        datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_6_$era.yml
sed -i 's/SM/7/g'        datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_7_$era.yml
sed -i 's/SM/8/g'        datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_8_$era.yml
sed -i 's/SM/8a/g'       datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_8a_$era.yml
sed -i 's/SM/9/g'        datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_9_$era.yml
sed -i 's/SM/10/g'       datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_10_$era.yml
sed -i 's/SM/11/g'       datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_11_$era.yml
sed -i 's/SM/12/g'       datacard_fit_TTHIDLoose_DNN11_datadriven_BM_classic_12_$era.yml
sed -i 's/SM/cluster1/g' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_1_$era.yml
sed -i 's/SM/cluster2/g' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_2_$era.yml
sed -i 's/SM/cluster3/g' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_3_$era.yml
sed -i 's/SM/cluster4/g' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_4_$era.yml
sed -i 's/SM/cluster5/g' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_5_$era.yml
sed -i 's/SM/cluster6/g' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_6_$era.yml
sed -i 's/SM/cluster7/g' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_7_$era.yml

sed -i 's/datacard_cluster1/datacard_1/' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_1_$era.yml
sed -i 's/datacard_cluster2/datacard_2/' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_2_$era.yml
sed -i 's/datacard_cluster3/datacard_3/' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_3_$era.yml
sed -i 's/datacard_cluster4/datacard_4/' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_4_$era.yml
sed -i 's/datacard_cluster5/datacard_5/' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_5_$era.yml
sed -i 's/datacard_cluster6/datacard_6/' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_6_$era.yml
sed -i 's/datacard_cluster7/datacard_7/' datacard_fit_TTHIDLoose_DNN11_datadriven_BM_new_7_$era.yml

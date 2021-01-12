import os
import sys
from copy import copy
import logging
logger = logging.getLogger(__name__) 
from bamboo import treefunctions as op

# strongly dependent on nSlices
def makeXvalEvaluator(evt, modelFile1,modelFile2,modelFile3,modelFile4,modelFile5, input_names, output_names, mvaType="Tensorflow"):
    model_1  = op.mvaEvaluator(modelFile1 , mvaType=mvaType, otherArgs=(input_names, output_name))
    model_2  = op.mvaEvaluator(modelFile2 , mvaType=mvaType, otherArgs=(input_names, output_name))
    model_3  = op.mvaEvaluator(modelFile3 , mvaType=mvaType, otherArgs=(input_names, output_name))
    model_4  = op.mvaEvaluator(modelFile4 , mvaType=mvaType, otherArgs=(input_names, output_name))
    model_5  = op.mvaEvaluator(modelFile5 , mvaType=mvaType, otherArgs=(input_names, output_name))
    methName = model_1.evaluate._name
    assert model_2.evaluate._name == methName
    assert model_3.evaluate._name == methName
    assert model_4.evaluate._name == methName
    assert model_5.evaluate._name == methName
    switchModelEval = op.multiSwitch(
        (evt%5 == 0, model_1.evaluate._objStb),
        (evt%5 == 1, model_2.evaluate._objStb),
        (evt%5 == 2, model_3.evaluate._objStb),
        (evt%5 == 3, model_4.evaluate._objStb),
        model_5.evaluate._objStb)
    return op.MVAEvaluator(getattr(switchModelEval, methName))



# input varList
def getInputVars_DNN (lepton, jets, met, HLL):
    return [ lepton.pt,
             lepton.eta,
             lepton.phi]

# output DNN
def getOutput_DNN (evt, modelFile1,modelFile2,modelFile3,modelFile4,modelFile5, input_names, output_names,
                   lepton, jets, met, HLL):
    DNNScoreEvaluator = makeXvalEvaluator(evt, modelFile1,modelFile2,modelFile3,modelFile4,modelFile5, input_names, output_names, mvaType="Tensorflow")
    return DNNScoreEvaluator(*getInputVars_DNN(lepton, jets, met, HLL))


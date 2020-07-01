from bamboo import treefunctions as op
from bamboo.treeoperations import ScaleFactorWithSystOp

def makeBtagRatioReweighting(jsonFile, numJets , variation="Nominal", systName=None, nameHint=None):
    """ Construct a btag ratio for MC, based on the weights in a JSON file

    :param btagRatioFile: path of the JSON file with weights (binned in NumJets)
    :param numJets : expression to get the number of selected jets
    :param systName: name of the associated systematic nuisance parameter
    """
    paramVType = "Parameters::value_type::value_type"
    args = op.construct("Parameters", (op.initList("std::initializer_list<{0}>".format(paramVType), paramVType,
        (op.initList(paramVType, "float", (op.extVar("int", "BinningVariable::NumJets"), numJets)),)),))
    wFun = op.define("ILeptonScaleFactor", 'const ScaleFactor <<name>>{{"{0}"}};'.format(jsonFile), nameHint=nameHint)
    expr = wFun.get(args, op.extVar("int", variation))
    if systName and variation == "Nominal":
        expr._parent = ScaleFactorWithSystOp(expr._parent, systName)
    return expr

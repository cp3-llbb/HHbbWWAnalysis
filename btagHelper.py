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

def makeBTagCalibrationReader(taggerName, csvFileName, wp=None, sysType="central", otherSysTypes=None, measurementType=None, sel=None, uName=None):
    if otherSysTypes is None:
        otherSysTypes = []
    if measurementType is None: ## BTV recommendation, see
        measurementType = {"B": "comb", "C": "comb", "UDSG": "incl"}
    elif isinstance(measurementType, str): ## if string, use it for all
        measurementType = {fl: measurementType for fl in ("B", "C", "UDSG")}
    calibName = sel._fbe.symbol(f'const BTagCalibration <<name>>{{"{taggerName}", "{csvFileName}"}};', nameHint=f"bTagCalib_{taggerName}")
    uName ="bTagReader_{0}".format("".join(c for c in uName if c.isalnum()))
    readerName = sel._fbe.symbol('BTagCalibrationReader <<name>>{{BTagEntry::OP_{0}, "{1}", {{ {2} }} }}; // for {3}'.format(wp.upper(), sysType, ", ".join(f'"{sv}"' for sv in otherSysTypes), uName), nameHint=f"bTagReader_{uName}")
    from bamboo.root import gbl
    calibHandle = getattr(gbl, calibName)
    readerHandle = getattr(gbl, readerName)
    for flav,measType in measurementType.items():
        readerHandle.load(calibHandle, getattr(gbl.BTagEntry, f"FLAV_{flav}"), measType)
    import bamboo.treefunctions as op
    return op.extVar("BTagCalibrationReader", readerName)

class BtagSF:
    def _nano_getPt(jet):
        import bamboo.treeproxies as _tp
        if isinstance(jet._parent, _tp.AltCollectionProxy):
            bs = jet._parent._base
            return bs.brMap["pt"].wrapped.result[jet.idx] ## use nominal always
        else:
            return jet.pt
    def _nano_getEta(jet):
        import bamboo.treefunctions as op
        return op.abs(jet.eta)
    def _nano_getJetFlavour(jet):
        import bamboo.treefunctions as op
        return op.extMethod("BTagEntry::jetFlavourFromHadronFlavour")(jet.hadronFlavour)

    def __init__(self, reader, getPt=None, getEta=None, getJetFlavour=None):
        self.reader = reader
        self.getPt = getPt if getPt is not None else BtagSF._nano_getPt
        self.getEta = getEta if getEta is not None else BtagSF._nano_getEta
        self.getJetFlavour = getJetFlavour if getJetFlavour is not None else BtagSF._nano_getJetFlavour

    def _evalFor(self, var, jet):
        import bamboo.treefunctions as op
        return self.reader.eval_auto_bounds(op._tp.makeConst(var, "std::string"), self.getJetFlavour(jet), self.getEta(jet), self.getPt(jet))

    def __call__(self, jet, systVars=None, systName=None):
        import bamboo.treefunctions as op
        nom = self._evalFor("central", jet)
        if systVars is None:
            return nom
        else:
            if systName is None:
                raise RuntimeError("A name for the systematic uncertainty group is needed")
            return op.systematic(nom, name=systName, **{ var : self._evalFor(var, jet) for var in systVars })

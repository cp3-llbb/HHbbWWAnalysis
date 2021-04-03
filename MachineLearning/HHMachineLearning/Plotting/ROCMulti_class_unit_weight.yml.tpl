ROC_multi:
  tree: tree
  classes:
    - Ewk
    - GGF
    - H
    - Top
    - VBF
    - WJets
  prob_branches:
    - output_Ewk
    - output_GGF
    - output_H
    - output_Top
    - output_VBF
    - output_WJets
  labels:
    - P(Ewk)
    - P(GGF)
    - P(H)
    - P(Top)
    - P(VBF)
    - P(WJets)
  colors:
    - '#610596'
    - '#288a24'
    - '#06b894'
    - '#cc7a16'
    - '#8f0a1e'
    - '#d95564'
  weight : '1'
  title : Multiclass_unit
  cut : ''
  selector :
    'Ewk'   : 'Ewk'
    'GGF'   : 'GGF'
    'H'     : 'H'
    'Top'   : 'Top'
    'VBF'   : 'VBF'
    'WJets' : 'WJets'

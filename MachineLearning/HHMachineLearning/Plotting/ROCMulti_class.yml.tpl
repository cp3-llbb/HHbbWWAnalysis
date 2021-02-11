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
<<<<<<< HEAD
    - "#288a24"
    - "#06b894"
    - "#610596"
    - "#99053d"
    - "#8f0a1e"
    - "#d95564"

=======
    - '#610596'
    - '#288a24'
    - '#06b894'
    - '#cc7a16'
    - '#8f0a1e'
    - '#d95564'
>>>>>>> eed57b2aa0f195c36370ef8af9adb2772972dcc4
  weight : '1'
  title : Multiclass
  cut : ''
  selector :
    'Ewk'   : 'Ewk'
    'GGF'   : 'GGF'
    'H'     : 'H'
    'Top'   : 'Top'
    'VBF'   : 'VBF'
    'WJets' : 'WJets'

ROC_multi:
  tree: tree
  classes:
    - DY
    - GGF
    - H
    - Rare
    - ST
    - TT
    - TTVX
    - VVV
  prob_branches:
    - output_DY
    - output_GGF
    - output_H
    - output_Rare
    - output_ST
    - output_TT
    - output_TTVX
    - output_VVV
  labels:
    - P(DY)
    - P(GGF)
    - P(H)
    - P(Rare)
    - P(ST)
    - P(TT)
    - P(TTVX)
    - P(VVV)
  colors:
    - '#1a83a1'
    - '#228B22'
    - '#06b894'
    - '#610596'
    - '#99053d'
    - '#cc7a16'
    - '#174704'
    - '#ccbf45'
  weight : '1'
  title : Multiclass_unit
  cut : ''
  selector :
    'DY'   : 'DY'
    'GGF'  : 'GGF'
    'H'    : 'H'
    'Rare' : 'Rare'
    'ST'   : 'ST'
    'TT'   : 'TT'
    'TTVX' : 'TTVX'
    'VVV'  : 'VVV'

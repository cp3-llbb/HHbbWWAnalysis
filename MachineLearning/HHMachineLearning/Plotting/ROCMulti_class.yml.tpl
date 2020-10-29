ROC_multi:
  tree: tree
  classes:
    - DY
    - H
    - HH
    - Rare
    - ST
    - TT
    - TTVX
    - VVV
    - WJets
  prob_branches:
    - output_DY
    - output_H
    - output_HH
    - output_Rare
    - output_ST
    - output_TT
    - output_TTVX
    - output_VVV
    - output_WJets
  labels:
    - P(DY)
    - P(H)
    - P(HH)
    - P(Rare)
    - P(ST)
    - P(TT)
    - P(TTVX)
    - P(VVV)
    - P(WJets)
  colors:
    - dodgerblue
    - mediumaquamarine
    - red
    - darkviolet
    - firebrick
    - darkorange
    - darkgreen
    - y
    - salmon 

  weight : 'event_weight'
  title : Multiclass
  cut : ''
  selector :
    'DY' : 'DY'
    'H' : 'H'
    'HH' : 'HH'
    'Rare' : 'Rare'
    'ST' : 'ST'
    'TT' : 'TT'
    'TTVX' : 'TTVX'
    'VVV': 'VVV'
    'WJets': 'WJets'
 


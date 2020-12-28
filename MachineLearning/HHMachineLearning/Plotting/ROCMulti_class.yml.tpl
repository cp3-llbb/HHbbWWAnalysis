ROC_multi:
  tree: tree
  classes:
    - DY
    - GGF
    - ST
  prob_branches:
    - output_DY
    - output_GGF
    - output_ST
  labels:
    - P(DY)
    - P(GGF)
    - P(ST)
  colors:
    - "#1a83a1"
    - "#288a24"
    - "#cc7a16"

  weight : 'event_weight'
  title : Multiclass
  cut : ''
  selector :
    'DY'    : 'DY'
    'GGF'    : 'GGF'
    'ST'   : 'ST'
 


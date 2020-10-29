multiclass:
  filename: 
  tree: tree
  weight: 'event_weight'
  name: multiclass 
  bins: 50
  xmin: 0
  xmax: 1
  title: ''
  xlabel: DNN output
  ylabel: events
  list_variable:
    - output_DY
    - output_H
    - output_HH
    - output_Rare
    - output_ST
    - output_TT
    - output_TTVX
    - output_VVV
    - output_WJets
  list_color:
    - 861
    - 845
    - 633
    - 877
    - 894
    - 805
    - 419
    - 397
    - 625
  list_legend:
    - node DY
    - node H
    - node HH
    - node Rare
    - node ST
    - node TT
    - node TTVX
    - node VVV
    - node WJets
  list_cut : '1'
  legend_pos:
    - 0.4
    - 0.60
    - 0.75
    - 0.92
  norm : True

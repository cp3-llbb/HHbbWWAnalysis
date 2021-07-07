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
    - output_GGF
    - output_H
    - output_Rare
    - output_ST
    - output_TT
    - output_TTVX
    - output_VVV
  list_color:
    - '#1a83a1'
    - '#228B22'
    - '#06b894'
    - '#610596'
    - '#99053d'
    - '#cc7a16'
    - '#174704'
    - '#ccbf45'
  list_legend:
    - node DY
    - node GGF
    - node H
    - node Rare
    - node ST
    - node TT
    - node TTVX
    - node VVV
  list_cut : 
    - output_DY > output_GGF && output_DY > output_H && output_DY > output_Rare && output_DY > output_ST && output_DY > output_TT && output_DY > output_TTVX && output_DY > output_VVV 
    - output_GGF > output_DY && output_GGF > output_H && output_GGF > output_Rare && output_GGF > output_ST && output_GGF > output_TT && output_GGF > output_TTVX && output_GGF > output_VVV 
    - output_H > output_DY && output_H > output_GGF && output_H > output_Rare && output_H > output_ST && output_H > output_TT && output_H > output_TTVX && output_H > output_VVV 
    - output_Rare > output_DY && output_Rare > output_GGF && output_Rare > output_H && output_Rare > output_ST && output_Rare > output_TT && output_Rare > output_TTVX && output_Rare > output_VVV 
    - output_ST > output_DY && output_ST > output_GGF && output_ST > output_H && output_ST > output_Rare && output_ST > output_TT && output_ST > output_TTVX && output_ST > output_VVV 
    - output_TT > output_DY && output_TT > output_GGF && output_TT > output_H && output_TT > output_Rare && output_TT > output_ST && output_TT > output_TTVX && output_TT > output_VVV 
    - output_TTVX > output_DY && output_TTVX > output_GGF && output_TTVX > output_H && output_TTVX > output_Rare && output_TTVX > output_ST && output_TTVX > output_TT && output_TTVX > output_VVV 
    - output_VVV > output_DY && output_VVV > output_GGF && output_VVV > output_H && output_VVV > output_Rare && output_VVV > output_ST && output_VVV > output_TT && output_VVV > output_TTVX
  legend_pos:
    - 0.4
    - 0.60
    - 0.75
    - 0.92
  norm : True

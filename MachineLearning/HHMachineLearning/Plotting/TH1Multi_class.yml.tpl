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
    - output_GGF
    - output_H
    - output_Rare
    - output_ST
    - output_TT
    - output_VBF
    - output_WJets
  list_color:
    - "#288a24"
    - "#06b894"
    - "#610596"
    - "#99053d"
    - "#cc7a16"
    - "#8f0a1e"
    - "#d95564"
  list_legend:
    - node GGF
    - node H
    - node Rare
    - node ST
    - node TT
    - node VBF
    - node WJets
  list_cut : '1'
  legend_pos:
    - 0.4
    - 0.60
    - 0.75
    - 0.92
  norm : True

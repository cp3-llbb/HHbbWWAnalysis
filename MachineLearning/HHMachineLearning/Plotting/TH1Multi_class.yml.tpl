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
    - output_ST
  list_color:
    - "#1a83a1"
    - "#cc7a16"
    - "#288a24"
  list_legend:
    - node DY
    - node GGF
    - node ST
  list_cut : '1'
  legend_pos:
    - 0.4
    - 0.60
    - 0.75
    - 0.92
  norm : True

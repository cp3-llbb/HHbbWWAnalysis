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
    - output_Ewk
    - output_GGF
    - output_H
    - output_Top
    - output_VBF
    - output_WJets
  list_color:
    - '#610596'
    - '#288a24'
    - '#06b894'
    - '#cc7a16'
    - '#8f0a1e'
    - '#d95564'
  list_legend:
    - node Ewk
    - node GGF
    - node H
    - node Top
    - node VBF
    - node WJets
  list_cut : 
    - output_Ewk > output_GGF && output_Ewk > output_H && output_Ewk > output_Top && output_Ewk > output_VBF && output_Ewk > output_WJets
    - output_GGF > output_Ewk && output_GGF > output_H && output_GGF > output_Top && output_GGF > output_VBF && output_GGF > output_WJets
    - output_H > output_GGF && output_H > output_Ewk && output_H > output_Top && output_H > output_VBF && output_H > output_WJets
    - output_Top > output_GGF && output_Top > output_H && output_Top > output_Top && output_Ewk > output_VBF && output_Top > output_WJets
    - output_VBF > output_GGF && output_VBF > output_H && output_VBF > output_Top && output_VBF > output_VBF && output_Ewk> output_WJets
    - output_WJets > output_GGF && output_WJets > output_H && output_WJets > output_Top && output_WJets > output_VBF && output_WJets > output_Ewk
  legend_pos:
    - 0.4
    - 0.60
    - 0.75
    - 0.92
  norm : True

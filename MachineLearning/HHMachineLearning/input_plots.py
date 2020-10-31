import sys
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import parameters

def InputPlots(data,list_inputs):
    logging.info("\tPlotting %s"%inp)
    pp = PdfPages('inputs.pdf')
    for inp in list_inputs:
        data[inp] = data[inp].astype(float)
        xmin = min(0.,data[inp].quantile(0.05))
        xmax = data[inp].quantile(0.95)
        sw = {node:data[data['tag']==node]['event_weight'].sum() for node in parameters.nodes}
        order = [_[0] for _ in sorted(sw.items(), key=lambda x : x[1])]
        order.reverse()
        w = {node:data[data['tag']==node]['event_weight'] for node in order}
        x = {node:data[data['tag']==node][inp] for node in order}
        c = {node:parameters.node_colors[node] for node in order}
        l = [node for node in order]

        # Hist plot #
        fig,axs = plt.subplots(1,1,figsize=(7,6))
        axs.hist(x          = list(x.values()),
                 bins       = 50,
                 range      = (xmin,xmax),
                 color      = list(c.values()),
                 histtype   = 'step',
                 density    = True,                    
                 weights    = list(w.values()),
                 label      = l)
        axs.set_xlabel(inp,fontsize=18)
        axs.set_ylabel("Events",fontsize=18)
        axs.legend(loc="upper right",fontsize=18)
        pp.savefig(fig)


        # Stack plot #
        fig,axs = plt.subplots(1,1,figsize=(7,6))
        axs.hist(x          = [v for k,v in x.items() if k != 'HH'],
                 bins       = 50,
                 range      = (xmin,xmax),
                 color      = [v for k,v in c.items() if k != 'HH'],
                 histtype   = 'barstacked',
                 stacked    = True,
                 weights    = [v for k,v in w.items() if k != 'HH'],
                 label      = [v for v in l if v != 'HH'])
        axs.hist(x          = x['HH'],
                 bins       = 50,
                 range      = (xmin,xmax),
                 color      = c['HH'],
                 histtype   = 'step',
                 weights    = w['HH']*(0.5*sum([v for k,v in sw.items() if k!='HH'])/sw['HH']),
                 label      = 'HH')
        axs.set_xlabel(inp,fontsize=18)
        axs.set_ylabel("Events",fontsize=18)
        axs.legend(loc="upper right",fontsize=18)
        pp.savefig(fig)
    pp.close()
    logging.info('Produced inputs.pdf')

import yaml
import os
import pandas as pd
import matplotlib as mpl
import seaborn as sns

from .utils import *

def rm_color_cats(palette, order, cats):
    if cats:
        keys = palette.keys()
        pop_list = []
        for key in keys:
            if key not in cats:
                pop_list.append(key)
        for p in pop_list:
            del palette[p]
        order = [o for o in order if o in cats]
    return palette, order

#######################################
############# Color palettes

def get_novelty_colors(cats=None):
    palette = {'FSM': "#61814B",
               'ISM': "#8EDE95", 
               'NIC': "#356CA1", 
               'NNC': "#C8773C",  
               'Intergenic': "darkred", 
               'Genic': "#B5B5B5", 
               'Fusion': "#4F4F4F", 
               'Antisense': "#6E5353"}
    order = ['FSM', 'ISM', 'NIC', 'NNC',
             'Intergenic', 'Genic', 'Fusion', 'Antisense']

    palette, order = rm_color_cats(palette, order, cats)
    return palette, order
    
def get_population_colors(cats=None):
    palette = {'ITU': '#db72f2',
                 'PEL': '#ff3a33',
                 'HAC': '#4cb33e',
                 'AJI': '#46bff0',
                 'LWK': '#A09136',
                 'YRI': '#DFBD00',
                 'CEU': '#347eed',
                 'MPC': '#eb9d0c'}

    order = list(palette.keys())
    order.sort()

    palette, order = rm_color_cats(palette, order, cats)
    return palette, order

def get_sex_colors(cats=None):
    palette = {'Female': '#ff69b4',
               'Male': '#1e90ff'}

    order = list(palette.keys())
    order.sort()

    palette, order = rm_color_cats(palette, order, cats)
    return palette, order

def get_cell_line_colors(cats=None):
    palette = {'panc1': '#2b9b81',
         'hepg2': '#fe9b00',
         'k562': '#f4c40f'}

    order = list(palette.keys())
    order.sort()

    palette, order = rm_color_cats(palette, order, cats)
    return palette, order

def get_rdna_region_colors(cats=None):
    palette = {'5_ETS': '#590003',
               '18S': '#a50706',
               'ITS1': '#427c8f',
               '5.8S': '#d98326',
               'ITS2': '#255a57',
               '28S': '#60d395',
               '3_ETS': '#d9d6d5'}

    order = get_rdna_regions()

    palette, order = rm_color_cats(palette, order, cats)
    return palette, order

def get_shade_colors(color, order):
    c_dict = {}
    min_color = '#FFFFFF'
    cmap = mpl.colors.LinearSegmentedColormap.from_list('temp', [color, min_color], N=len(order)+1)
    for i, cat in enumerate(order):
        c_dict[cat] = mpl.colors.to_hex(cmap(i))

    return c_dict, order

def init_plot_settings(font_scale=2,
                       aspect='square'):
    """
    Initialize default plotting settings
    """
    sns.set_context('paper', font_scale=2)
    mpl.rcParams['font.family'] = 'Helvetica'
    mpl.rcParams['pdf.fonttype'] = 42

    if aspect=='square':
        mpl.rcParams['figure.figsize'] = (5,5)
    elif aspect=='rectangle':
        mpl.rcParams['figure.figsize'] = (7,5)

def get_dataset_colors(cats):
    """
    Get the color associated with each dataset (library)
    """
    df = load_meta()
    df = df[['dataset', 'sample_color_hex_code']]
    df.set_index('dataset', inplace=True)
    d = df.to_dict(orient='index')
    df.reset_index(inplace=True)
    palette = {}
    for key, item in d.items():
        palette[key] = item['sample_color_hex_code']
    order = df.dataset.tolist()

    palette, order = rm_color_cats(palette, order, cats)
    return palette, order

    return c_dict, order

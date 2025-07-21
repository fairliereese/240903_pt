import yaml
import os
import pandas as pd
import matplotlib as mpl
import seaborn as sns
from matplotlib.colors import to_rgb, to_hex
import colorsys


from .utils import *

def mute_color(hex_color, factor=1.4):
    """Lighten a hex color by increasing its lightness in HLS space."""
    r, g, b = to_rgb(hex_color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    l = min(1.0, l * factor)  # Lighten
    r_muted, g_muted, b_muted = colorsys.hls_to_rgb(h, l, s)
    return to_hex((r_muted, g_muted, b_muted))

def adjust_lightness(color, amount=0.5):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], max(0, min(1, amount * c[1])), c[2])

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

def get_annot_colors(cats=None):
    palette = {'GENCODE': "#6c5a52",
               'PODER': "#a0ad55",
               'Enhanced GENCODE': '#a7a9ab',
               'Enhanced\nGENCODE':'#a7a9ab'}
    order = ['GENCODE', 'PODER', 'Enhanced GENCODE', 'Enhanced\nGENCODE']

    palette, order = rm_color_cats(palette, order, cats)
    return palette, order

def get_afr_colors(cats=None):
    palette = {'AFR': "#D5AC4B",
               'OOA': "#259069"}
    order = ['AFR', 'OOA']

    palette, order = rm_color_cats(palette, order, cats)
    return palette, order

def get_eur_colors(cats=None):
    palette = {'European': "#466995",
               'Non-European': "A53860"}
    order = ['European', 'Non-European']

    palette, order = rm_color_cats(palette, order, cats)
    return palette, order

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

def my_theme(base_size=11, w=4, h=3):
    """
    Custom plotnine theme with:
    - White background
    - Clean styling
    - Axes and ticks retained

    Parameters:
    - base_size: Base font size

    Returns:
    - plotnine.theme object
    """
    return (
        theme_minimal(base_size=base_size)
        + theme(
            # White background
            panel_background=element_rect(fill='white', color=None),
            plot_background=element_rect(fill='white', color=None),

            # Remove grid lines
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            panel_border=element_blank(),

            # Keep axis lines & ticks (don't blank them)
            axis_line=element_line(color='black'),
            axis_ticks=element_line(color='black'),

            plot_title=element_text(hjust=0.5, family='Helvetica'),
            axis_title_x=element_text(hjust=0.5, family='Helvetica'),
            axis_title_y=element_text(hjust=0.5, margin={'t':0, 'r':-2, 'b':0, 'l':0}, family='Helvetica'),
            
            # Styling text
            legend_title=element_blank(),
            axis_title=element_text(size=base_size + 1, family='Helvetica'),
            legend_text=element_text(size=base_size-2, family='Helvetica'),
            axis_text=element_text(size=base_size, color='black', family='Helvetica'),
            strip_text_x=element_text(size=base_size-1),
            strip_text_y=element_text(size=base_size-1),
            figure_size=(w, h),  # Controls plot dimensions (width x height in inches)
            plot_margin=0.05      # Shrinks surrounding white space
        )
    )

def clean_figure(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis="x", rotation=45)
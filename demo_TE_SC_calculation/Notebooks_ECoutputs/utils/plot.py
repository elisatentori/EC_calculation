import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import ScalarFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, to_rgba

import numpy as np

from . import colormaps as maps

#=================================================================================================#
from matplotlib import font_manager, rcParams
font_file = "/home/tentori/.local/avenir_ff/AvenirLTStd-Roman.otf"
font_file_b = "/home/tentori/.local/avenir_ff/AvenirLTStd-Black.otf"
font_file_c = "/home/tentori/.local/avenir_ff/AvenirLTStd-Book.otf"
font_manager.fontManager.addfont(font_file)
font_manager.fontManager.addfont(font_file_b)
font_manager.fontManager.addfont(font_file_c)

# Imposta il font predefinito su Avenir
rcParams['font.family'] = "Avenir LT Std"

DIM = 22

plt.rcParams.update({
    'font.size': DIM,
    'axes.labelsize': DIM,
    'axes.titlesize': DIM,
    'xtick.labelsize': DIM,
    'ytick.labelsize': DIM
})
#=================================================================================================#

import os
def Set_Dir_Plots(path):
    if not os.path.exists(path):
        os.mkdir(path)

coldhot_cmap   = maps.create_cmaphot()
coldhot_cmap_r = maps.create_cmaphot_r()

colors = ['#2F7FC3','#E62A08','#464646','#FFD700','#32CD32','#8A2BE2']


#=================================================================================================#
# Axes formatter for plots 

def set_format(ax, axis_ticks = 'both', pwr_x_min=-1, pwr_x_max=1, pwr_y_min=-1, pwr_y_max=1,  cbar = None, pwr_cbar_min=-1, pwr_cbar_max=1,  DIM = 30):
    
    import seaborn as sns
    
    sns.despine(ax=ax, trim=False)
    ax.set_facecolor('none')
    
    # - - -  TICKS
    ax.tick_params(axis=axis_ticks, which='major', labelsize=DIM)
    
    # - - -  FORMATTER x axis
    formatter_x = ScalarFormatter(useMathText=True)   
    formatter_x.set_scientific(True)
    formatter_x.set_powerlimits((pwr_x_min, pwr_x_max))
    ax.xaxis.set_major_formatter(formatter_x)
    ax.xaxis.offsetText.set_fontsize(DIM-10)
    
    from matplotlib.transforms import ScaledTranslation
    dx, dy = 15/72, 15/72
    offset = ScaledTranslation(dx, dy, ax.figure.dpi_scale_trans)
    ax.xaxis.offsetText.set_transform(ax.xaxis.offsetText.get_transform() + offset)

    # - - -  FORMATTER y axis
    formatter_y = ScalarFormatter(useMathText=True)    
    formatter_y.set_scientific(True) 
    formatter_y.set_powerlimits((pwr_y_min, pwr_y_max))
    ax.yaxis.set_major_formatter(formatter_y);
    ax.yaxis.offsetText.set_fontsize(DIM-10)
    
    if cbar:
        # - - -  FORMATTER cbar
        formatter_cbar = ScalarFormatter(useMathText=True)   
        formatter_cbar.set_scientific(True)
        formatter_cbar.set_powerlimits((pwr_cbar_min, pwr_cbar_max))
        cbar.ax.yaxis.set_major_formatter(formatter_cbar); 
        cbar.ax.yaxis.offsetText.set_fontsize(DIM-10)
        cbar.ax.xaxis.set_major_formatter(formatter_cbar); 
        cbar.ax.xaxis.offsetText.set_fontsize(DIM-10)
        
        # Move the offset text to the top of the colorbar
        dx, dy = 0.8, 0.3  # Adjust dy for vertical and dx for horizontal shifts
        cbar_offset = ScaledTranslation(dx, dy, cbar.ax.figure.dpi_scale_trans)
        cbar.ax.yaxis.offsetText.set_transform(cbar.ax.yaxis.offsetText.get_transform() + cbar_offset)
        
#=================================================================================================#
# utils for raster plot

def neuronID_to_cluster(neuron_ids, channel, cluster):
    max_chan = np.max(channel) + 1
    lookup = np.full(max_chan, -1, dtype=int)
    lookup[channel] = cluster
    return lookup[neuron_ids]
      
def alternating_colormap(n_clusters, color1='black', color2='crimson'):
    """Crea una ListedColormap alternata con due colori contrastanti."""
    colors = [color1 if i % 2 == 0 else color2 for i in range(n_clusters)]
    return mcolors.ListedColormap(colors)


#=================================================================================================#
# Raster plot

def rasterplot(st, channel, cluster, title=None, tmin=0, tmax=120, cmap='viridis',figsize=(15,5), outf : str = None, show_plot = True):
    spike_times = st[0]
    neuron_ids  = st[1].astype(int)
    cluster_ids = neuronID_to_cluster(neuron_ids, channel, cluster)
    condition   = np.logical_and(spike_times>=tmin, spike_times<=tmax)

    # Colormap
    n_clusters = len(np.unique(cluster))
    if cmap is None:
        cmap = alternating_colormap(n_clusters)
    else:
        cmap = cmap
        
    plt.subplots(figsize=figsize)
    plt.scatter(spike_times[condition], neuron_ids[condition], s=0.1, c=cluster_ids[condition], cmap=cmap)
    plt.xlabel("time (s)")
    plt.ylabel("channel ID")
    #plt.xlim(tmin,tmax)
    if title:
        plt.title(title)
    plt.colorbar(label="cluster ID")

    if outf:
        plt.savefig(outf, bbox_inches='tight')
        if not show_plot:
            plt.close()
    else:
        plt.show()


#=================================================================================================#
# Channel map with colorbar

def plot_map(pos, car_array, cbar_label='cluster ID', title='recording channels map', cmap='viridis', 
             figsize=(15, 7), ax=None, outf : str = None, show_plot = True):

    # Colormap
    n_clusters = len(np.unique(car_array))
    if cmap is None:
        cmap = alternating_colormap(n_clusters)
    else:
        cmap = cmap
        
    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure() 
        
    # channels map
    im = ax.scatter(pos[:,0],pos[:,1],c=car_array,s=1,marker='s',cmap=cmap)
    xt=ax.set_xlabel('x (mm)'); ax.set_ylabel('y (mm)')

    if title:
        ax.set_title(title)
    fig.colorbar(im, ax=ax, label=cbar_label)
    
    if ax is None:
        if outf:
            ax.savefig(outf, bbox_inches='tight')
            if not show_plot:
                plt.close()
        else:
            plt.show()

#=================================================================================================#
# plot matrix – aspect='auto'

def plot_mat_aspect(mat, vmin=None, vmax=None, cmap='viridis', title=None, xlabel='target', ylabel='source', 
                    cbarlabel='spike count', invert_y: bool = True, xticklabels: list = None, yticklabels: list = None, 
                    figsize=(15, 10), ax = None, outf: str = None, show_plot: bool = True):
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure() 
        
    if vmin is None and vmax is None:
        vmin, vmax = np.percentile(mat, [5, 97.5])
        
    im = ax.imshow(mat, aspect='auto', vmin=vmin, vmax=vmax, cmap=cmap)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title)

    if xticklabels is not None:
        ax.set_xticks(np.arange(len(xticklabels)))
        ax.set_xticklabels(xticklabels, rotation=45)

    if yticklabels is not None:
        ax.set_yticks(np.arange(len(yticklabels)))
        ax.set_yticklabels(yticklabels)

    if invert_y:
        ax.invert_yaxis()

    fig.colorbar(im, ax=ax, label=cbarlabel)

    if ax is None:
        if outf:
            plt.savefig(outf, bbox_inches='tight')
            
            if not show_plot:
                plt.close()
        else:
            plt.show()
        
# ------------------------------------------------------------------------------------ #

def plot_mat(mat, title=None, xlabel='target', ylabel='source', cmap='viridis', cbarlabel='spike count', invert_y: bool = True, figsize=(12, 10), ax = None, outf : str = None, show_plot = True):

    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure() 
        
    if vmin is None and vmax is None:
        vmin, vmax = np.percentile(mat, [5, 97.5])
                
    im = ax.imshow(mat, vmin=vmin, vmax=vmax, cmap=cmap)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if title:
        ax.set_title(title)
        
    if invert_y:
        ax.invert_yaxis()
    
    fig.colorbar(im, ax=ax, label=cbarlabel)
    
    if outf:
        plt.savefig(outf, bbox_inches='tight')
    
    if ax is None:
        if outf:
            plt.savefig(outf, bbox_inches='tight')
            
            if not show_plot:
                plt.close()
        else:
            plt.show()




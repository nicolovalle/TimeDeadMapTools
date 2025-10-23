import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from matplotlib.cm import ScalarMappable
import matplotlib.gridspec as gridspec
import numpy as np
import os

from mylogger import *

log_file = 'QApy.log'
LHCOrbitNS = 88924.6
# === Canvas Configuration ===
view_file_path = "1dview_coordinates.txt"
qc_view_file_path = "1dview_coordinates_qc.txt"
output_dir = "canvas"
save_each = False
os.makedirs(output_dir, exist_ok=True)

#____________________________
def LOG(severity, *message):
    logger = Logger(log_file)
    logger.set_highlight_keyword(True)
    logger.log(severity,*message)

#########################################################
def make_canvas1(
        lane_dead_time = [0.5,1]*(3816//2),
        stave_dead_time = [0.5,1] * (192 // 2),
        number_of_fully_dead = [3,]*3816,
        stave_recovery_rate = [0.5,1] * (192 // 2),
        text1 = 'k#First line#r#Second line#g#Third Line', # Format: 'color#Senetence
        text2 = 'b#Text2',
        gaps = [1,]*100000,
        dead_fraction = [{'both':list(range(1000))},{'IB': [999] * 1000, 'OB': [999] * 1000}],
        dead_fraction_rolling = [{'both':list(range(1000))},{'IB': [999] * 1000, 'OB': [999] * 1000}],
        words = [[20,]*100, [30,]*100],
        worst_ob = [1,0]*(3816//2),
        WorstOBstep = 20,
        WorstIBstep = 20,
        worst_ib = [1,0]*(3816//2),
        last = [1,0]*(3816//2)
        ):

    
    LOG(INFO,f'Creating canvas 1')

    # === Load bins from file ===
    bins = []
    bins2 = []
    with open(view_file_path) as f:
        for line in f:
            parts = list(map(float, line.strip().split()))
            if len(parts) != 9:
                continue
            _, x0, x1, x2, x3, y0, y1, y2, y3 = parts
            x = [x0, x1, x2, x3]
            y = [y0, y1, y2, y3]
            bins.append((x, y))
    
    with open(qc_view_file_path) as f:
        for line in f:
            parts = list(map(float, line.strip().split()))
            if len(parts) != 7:
                continue
            _, x0, x1, x2, y0, y1, y2 = parts
            x = [x0/10, x1/10, x2/10]  # saved as cm
            y = [y0/10, y1/10, y2/10]  # saved as cm
            bins2.append((x,y))

    
    # === Create canvas ===
    fig, axes = plt.subplots(4, 4, figsize=(20, 20))
    unused = list(range(len(fig.axes))) # updated along the loop over axes
    axes = axes.flatten()
    WorstIBtime = sum([gaps[i]*LHCOrbitNS*1.e-9 for i in range(WorstIBstep)])
    WorstOBtime = sum([gaps[i]*LHCOrbitNS*1.e-9 for i in range(WorstOBstep)])
    plot_title = ['Lane dead time after trg ramp','Number of dead chips in the lane', 'INFO', 'QUALITY',
                  'Orbit gaps','Gap distribution', 'Dead fraction', 'Dead fraction rolling avg',
                  f'Number of words (tot: {sum(words[1])/1000:.1f}k)', f'Worst IB: step {WorstIBstep} = {WorstIBtime:.1f} s', f'Worst OB: step {WorstOBstep} = {WorstOBtime:.1f} s', 'Lanes with single dead chips',
                  'Stave dead time','Recovery per hour']

    x_axis_title = ['','','','',
                    'step','seconds','step','time (min)',
                    'words','','','',
                    'x (cm) (IBx3)','x (cm) (IBx3)']
    y_axis_title = ['','','','',
                    'orbits','counts','','',
                    'counts','','','',
                    'y (cm) (IBx3)','y (cm) (IBx3)']
    
    # === TH2Poly-like plots: indices [0, 1, 9, 10, 11] ===
    poly_indices = [0, 1, 9, 10, 11]
    for i, ax_idx in enumerate(poly_indices):
        unused.remove(ax_idx)
        ax = axes[ax_idx]
        poly_data = bins #poly_sets[i]
        patches_list = []
        values = []

    
        for x, y in poly_data:
            xy = np.column_stack([x, y])
            polygon = plt.Polygon(xy, closed=True)
            patches_list.append(polygon)
            
        if ax_idx == 0:
            values = lane_dead_time
        if ax_idx == 1:
            values = number_of_fully_dead
        if ax_idx == 9:
            values = worst_ib
        if ax_idx == 10:
            values = worst_ob
        if ax_idx == 11:
            values = last

        vmin = min(values)
        vmax = max(values)
        if 'Lane dead time' in plot_title[ax_idx] and vmin < vmax: # and vmin > 0:
            norm = LogNorm(vmin=1.e-5, vmax=1)
        else:    
            norm = Normalize(vmin=vmin, vmax=vmax)

        cmap = 'coolwarm'
        if vmin == vmax:
            cmap += '_r'
        if 'single dead chips' in plot_title[ax_idx]:
            cmap = 'Accent'
        pc = PatchCollection(patches_list, cmap=cmap, edgecolor='k', linewidth=0.1, norm=norm)
        pc.set_array(np.array(values))
        ax.add_collection(pc)
        ax.autoscale_view()
        if np.min(np.array(values)) == np.max(np.array(values)):
            pc.set_clim(np.min(np.array(values)),np.min(np.array(values)))
        
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(plot_title[ax_idx])
        if ax_idx < 8:
            fig.colorbar(pc, ax=ax, shrink=0.8)# .set_label("Value")

        if save_each:
            fig_single, ax_single = plt.subplots()
            pc_single = PatchCollection(patches_list, cmap="coolwarm", edgecolor='k', linewidth=0.1, norm=norm)
            pc_single.set_array(np.array(values))
            ax_single.add_collection(pc_single)
            ax_single.autoscale_view()
            ax_single.set_aspect('equal')
            ax_single.set_xticks([])
            ax_single.set_yticks([])
            if ax_idx < 8:
                fig_single.colorbar(pc_single, ax=ax_single, shrink=0.8).set_label("Value")
            fig_single.savefig(f"{output_dir}/poly_plot_{ax_idx+1}.png")
            plt.close(fig_single)

    # === TH2Poly-like triangular style: indices [12, 13 ] ===
    poly2_indices = [12, 13]
    for i, ax_idx in enumerate(poly2_indices):
        unused.remove(ax_idx)
        ax = axes[ax_idx]
        poly2_data = bins2
        patches_list = []
        values = []

        for x,y in poly2_data:
            xy = np.column_stack([x,y])
            polygon = plt.Polygon(xy, closed = True)
            patches_list.append(polygon)

        if ax_idx == 12:
            values = stave_dead_time
        if ax_idx == 13:
            values = stave_recovery_rate

        vmin = min(values)
        try:
            vminpositive = min([i for i in values if i > 0])
        except:
            vminpositive = 0
        vmax = max(values)

        if 'Stave dead time' in plot_title[ax_idx] and vmin < vmax: # and vmin > 0:
            norm = LogNorm(vmin=1.e-5, vmax=1)
        elif 'Recover' in plot_title[ax_idx] and vmin < vmax: # and vmin > 0:
            
            norm = LogNorm(vmin=0.8*vminpositive if vminpositive > 0 else 1.e-2, vmax=vmax)
        else:
            norm = Normalize(vmin=vmin, vmax=vmax)

        cmap = 'coolwarm'
        if vmin == vmax:
            cmap += '_r'
        pc = PatchCollection(patches_list, cmap=cmap, edgecolor='k', linewidth=0.1, norm=norm)
        pc.set_array(np.array(values))
        ax.add_collection(pc)
        ax.autoscale_view()
        if np.min(np.array(values)) == np.max(np.array(values)):
            pc.set_clim(np.min(np.array(values)),np.min(np.array(values)))
        
        ax.set_aspect('equal')
        #ax.set_xticks([])
        #ax.set_yticks([])
        ax.set_title(plot_title[ax_idx])
        ax.set_xlabel(x_axis_title[ax_idx])
        ax.set_ylabel(y_axis_title[ax_idx])
        fig.colorbar(pc, ax=ax, shrink=0.8)

        
    # === Text Info: indices [2, 3] ===
    texts = [text1, text2]
    for idx, text in zip([2, 3], texts):
        unused.remove(idx)
        ax = axes[idx]
        ax.axis('off')
        dec = text.split('#')
        colors = dec[::2]
        sentences = dec[1::2]
        for i in range(len(sentences)):
            ax.text(0.05, 1-0.06-0.055*i, sentences[i], ha='left', va='center', fontsize=11, color=colors[i])

        if save_each:
            fig_text, ax_text = plt.subplots()
            ax_text.axis('off')
            for i in range(len(sentences)):
                fig_ax.text(0.05, 1-0.08*i, sentences[i], ha='left', va='center', fontsize=10, color=colors[i])
            fig_text.savefig(f"{output_dir}/text_{idx+1}.png")
            plt.close(fig_text)
    
    # === Multigraphs at indices [4, 6, 7] ===
    multi_indices = [4, 6, 7]
    for idx in multi_indices:
        unused.remove(idx)
        ax = axes[idx]
        if idx == 4:
            x_ = list(range(len(gaps)))
            data1 = data2 = [x_,gaps]
            markersize = 1
        if idx == 6:
            if 'both' in dead_fraction[0]:
                data1 = [dead_fraction[0]['both'],dead_fraction[1]['IB']]
                data2 = [dead_fraction[0]['both'],dead_fraction[1]['OB']]
            else:
                data1 = [dead_fraction[0]['IB'],dead_fraction[1]['IB']]
                data2 = [dead_fraction[0]['OB'],dead_fraction[1]['OB']]
            markersize = 1
        if idx == 7:
            if 'both' in dead_fraction_rolling[0]:
                data1 = [dead_fraction_rolling[0]['both'],dead_fraction_rolling[1]['IB']]
                data2 = [dead_fraction_rolling[0]['both'],dead_fraction_rolling[1]['OB']]
            else:
                data1 = [dead_fraction_rolling[0]['IB'],dead_fraction_rolling[1]['IB']]
                data2 = [dead_fraction_rolling[0]['OB'],dead_fraction_rolling[1]['OB']]
            markersize = 3
        x1 = np.array(data1[0])
        y1 = np.array(data1[1]) 
        x2 = np.array(data2[0])
        y2 = np.array(data2[1])
        data1.clear()
        data2.clear()

        try:      
            ax.plot(x1, y1, label='IB', marker='o', markersize=markersize)
            ax.plot(x2, y2, label='OB', marker='s', markersize=markersize)
            ax.set_title(plot_title[idx])
            ax.set_xlabel(x_axis_title[idx])
            ax.set_ylabel(y_axis_title[idx])
            if 'Dead fraction' in plot_title[idx]:# and min(y1) > 0 and min(y2) > 0:
                ax.set_yscale('log')
            if not np.array_equal(y1,y2):
                ax.legend()
     
            if save_each:
                fig_mg, ax_mg = plt.subplots()
                ax_mg.plot(x1, y1, label='IB', marker='o', markersize=markersize)
                ax_mg.plot(x2, y2, label='OB', marker='s', markersize=markersize)
                ax_mg.set_title(plot_title[idx])
                ax_mg.set_xlabel(x_axis_title[idx])
                ax_mg.set_ylabel(y_axis_title[idx])
                if 'Dead fraction' in plot_title[idx] and min(y1) > 0 and min(y2) > 0:
                    ax.set_yscale('log')
                if not np.array_equal(y1,y2):
                    ax_mg.legend()
                fig_mg.savefig(f"{output_dir}/multigraph_{idx+1}.png")
                plt.close(fig_mg)
        except Exception as e:
            LOG(WARNING,f'Exception while creating plot in slot {idx}: {e}')
    
    # === 1D histograms ===
    hist_indices = [5,8]
    for i in hist_indices:
        unused.remove(i)
        ax = axes[i]
        if i == 5:
            data = np.array(gaps) * (LHCOrbitNS * 1.e-9)
        if i == 8:
            data = np.array(words[1])
        ax.hist(data, bins=20, color='lightcoral', edgecolor='black')
        ax.set_title(plot_title[i])
        ax.set_xlabel(x_axis_title[i])
        ax.set_ylabel(y_axis_title[i])
        ax.set_yscale('log')

        if save_each:
            fig_hist, ax_hist = plt.subplots()
            ax_hist.hist(data, bins=20, color='lightcoral', edgecolor='black')
            ax_hist.set_title(plot_title[i])
            ax_hist.set_xlabel(x_axis_title[i])
            ax_hist.set_ylabel(y_axis_title[i])
            ax_hist.set_yscale('log')
            fig_hist.savefig(f"{output_dir}/histogram_{i+1}.png")
            plt.close(fig_hist)

    # === Scatter plot ====
    #hist_indices = [8,]
    #for i in hist_indices:
    #    unused.remove(i)
    #    ax = axes[i]
    #    if i == 8:
    #       datax = words[0] # chips
    #       datay = words[1] # words
    #   ax.hist2d(datax, datay, bins=30, cmap='RdPu')
    #   ax.set_title(plot_title[i])
    #   ax.set_xlabel(x_axis_title[i])
    #   ax.set_ylabel(y_axis_title[i])
    #   ax.set_yscale('log')
    #   ax.set_xscale('log')


    # remove empty pads
    for idx in unused:
        fig.delaxes(axes[idx])
    
    # === Save full canvas ===
    outpath = os.path.join(output_dir, 'full_canvas1.png')
    LOG(INFO,f'Saving to {outpath}')
    plt.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
    #LOG(INFO,'Done')

###############################################################
def make_canvas2( 
        lanemap={k: np.arange(100)/100 for k in range(10)},
        idx = [[2,8], [10,20], [20,30]],
        offset_sec = 0, # this defines the first value on the x axis
        run = '?',
        spec0 = '', # a string for the canvas header
        spec1 = ['','',''], # three strings to append to the plot titles
        spec2 = '' # a string to append to the file name
        ):


    LOG(INFO,f'Creating canvas 2 {spec2}')
    
    # Extract the keys and arrays from the lanemap
    keys = list(lanemap.keys())  # x-axis values (keys of the lanemap dictionary)
    cansize = (20,10)
    candpi = 500
    LOG(INFO,f'Evaluating sampling the map for display')
    scale_legend = "No sampling"
    if max([keys[i] - keys[i-1] for i in range(1,len(keys))]) > 17000: # approx 1.5 seconds
        LOG(INFO,f'There are large gaps, original size {len(keys)} will be kept')
        arrays = [lanemap[k] for k in keys]  # Corresponding y-values (arrays)
    else:
        maxbinpp = candpi*cansize[0]/3
        scale0 = int(len(keys) / maxbinpp)-1
        if scale0 < 2:
            LOG(INFO,f'Not so many points, original size {len(keys)} will be kept')
            arrays = [lanemap[k] for k in keys]  # Corresponding y-values (arrays)
        else:
            keys = keys[::scale0]
            LOG(INFO,f'Sampling every {scale0} entry of the original map. Sampled map has {len(keys)} entries')
            arrays = [lanemap[k] for k in keys]  # Corresponding y-values (arrays)
            scale_legend = f"Sampled 1:{scale0}"
    
    

    # Convert arrays into a 2D numpy array
    fulldata = np.array(arrays)


    #fig, axes = plt.subplots(1, 3, figsize=(14, 7))
    #axes = axes.flatten()
    #mpl.rcParams['savefig.dpi'] = 'figure'
    fig = plt.figure(figsize=cansize)
    gs = gridspec.GridSpec(2, 3, height_ratios=[1,0.02], hspace=0.2)
    axes = [fig.add_subplot(gs[0,i]) for i in range(3)]

    colors = ['#ffffff', '#5773e0', '#7b98ef', '#9FBEFF', '#C8D8F1', '#F3C6B3', '#E0664F', '#A70E2A']
    new_cmap = mcolors.ListedColormap(colors)

 
    plot_titles = [f'Inner Barrel{spec1[0]}', f'Layer 3,4{spec1[1]}', f'Layer 5,6{spec1[2]}']

    def getlabel(b,l): # b is barrel, l is the lane in barrel
        NStaves = [ 12, 16, 20, 24, 30, 42, 48 ]
        if b == 0:
            if l < 12*9:
                return f'L0_{l//9}'
            elif l < 28*9:
                return f'L1_{(l-12*9)//9}'
            else:
                return f'L2_{(l-28*9)//9}'
        elif b == 1:
            if l < 24*16:
                return f'L3_{l//16}'
            else:
                return f'L4_{(l-24*16)//16}'
        else:
            if l < 42*28:
                return f'L5_{l//28}'
            else:
                return f'L6_{(l-42*28)//28}'
            

    for i in [0,1,2]:

        hline_step = [9,16,28]

        data = (fulldata.T)[idx[i][0]:idx[i][1]]
        zero_row = np.zeros((1,data.shape[1]))
        data = np.vstack([zero_row,data])
        # Prepare the bin edges for pcolormesh
        x_edges = np.array(keys)
        x_edges = (x_edges - x_edges[0]) * LHCOrbitNS * 1.e-9 / 60 + offset_sec / 60
        x_edges = np.append(x_edges, x_edges[-1]+2/60) #adding 2 seconds at the end for visualization of the last step

        xlabel = 'time (min)'
        if x_edges[-1]-x_edges[0] < 5:
            xlabel = 'time (sec)'
            x_edges = x_edges * 60

        # y_edges should have data.shape[0] + 1 elements (rows + 1)
        y_edges = np.arange(idx[i][1]-idx[i][0]+1)
        y_edges = np.insert(y_edges, 0, 0-hline_step[i])


        # Select the second axis for the plot
        ax = axes[i]
        
        # Plot the data using pcolormesh (respecting x and y bin edges)
        c = ax.pcolormesh(x_edges, y_edges, data, cmap=new_cmap, shading='auto', vmin=0, vmax=1)
        #TO BE TRIED: c = ax.pcolorfast(x_edges, y_edges, data, cmap=new_cmap, vmin=0, vmax=1)


        
        for y in range(0, data.shape[0]-1, hline_step[i]):
            ax.axhline(y, color='grey', linestyle='--', linewidth=0.2)                   
            ax.text(x_edges[0], y+hline_step[i]/2, getlabel(i,y)+' ', color='black', va='center', ha='right', fontsize=6 if i<2 else 5)
        for x in x_edges:
            ax.vlines(x, 0-hline_step[i], 0, color='grey', linewidth=0.2)
        ax.axhline(0, color='k', linestyle='-', linewidth=0.5)   
        
 

        ax.set_title(plot_titles[i])
        ax.set_xlabel(xlabel)
        ax.set_yticklabels([])
        ax.tick_params(axis='y', which='both', size=0)  # Removes y-axis ticks


    # Create a scalar mappable for the shared colorbar
    sm = ScalarMappable(cmap=new_cmap, norm=plt.Normalize(vmin=0, vmax=1))
    sm.set_array([])  # Required for older matplotlib versions
    cbar_ax = fig.add_subplot(gs[1, 1])
    cbar_ax.tick_params(size=0)
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='horizontal')
    ticks = np.linspace(0, 1, 2)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{t:.2f}" for t in ticks])  # or custom labels
    #cbar.set_label('z value')
    cbar.set_ticks(np.linspace(0, 1, len(colors)))

    
    fig.suptitle(f'Run {run}{spec0}', fontsize=13, y=0.985)

    # Save the figure
    name_spec = f'_{spec2}' if spec2 else ''
    outpath = os.path.join(output_dir, f'full_canvas2{name_spec}.png')
    LOG(INFO,f'Saving to {outpath}')
    #plt.tight_layout()
    fig.subplots_adjust(
        left=0.04,
        right=0.98,
        top=0.925,
        bottom=0.035,
        wspace=0.1
    )

    plt.figtext(0.02, 0.02, scale_legend, ha="left", va="bottom", fontsize=8, color="gray")
        
    fig.savefig(outpath, dpi=candpi)
    #fig.show()

    # Close the figure to release resources
    plt.close(fig)


def make_canvas4( 
        lane_dead_time = [0.1*i for i in range(9)]*(3816//9),
        idx = [[0, 108], [108, 252], [252, 432], [432, 816], [816, 1296], [1296, 2472], [2472, 3816]],
        run = '?'
        ):

    LOG(INFO,f'Creating canvas 4')
    
    # Convert arrays into a 2D numpy array
    fulldata = np.array(lane_dead_time)


    fig, axes = plt.subplots(7, 1, figsize=(5, 15))
    axes = axes.flatten()

    vsteps = [9,9,9,16,16,28,28]
    for i in [0,1,2,3,4,5,6]:
        X = np.arange(idx[i][1]-idx[i][0]+1)
        Y = fulldata[idx[i][0]:idx[i][1]]
        Y = np.append(Y,0)

        ax = axes[i]

        ax.step(X,Y, where='post', color='blue', linewidth=0.2)
        ax.set_xlim(X[0],X[-1])
        ax.set_ylim(1e-4,1.2)
        ax.set_yscale('log')
        ax.set_ylabel(f'Layer {i}')
        ax.set_xticklabels([])
        ax.tick_params(axis='x', which='both', size=0)

        
        for x in range(0, len(X)-1, vsteps[i]):
            ax.axvline(x, color='grey', linestyle='--', linewidth=0.2)
            ax.text(x+1, 1.3, x//vsteps[i], color='k', fontsize=4)

        if i==0:
            ax.text(108, 4, f'Lane dead time after trigger ramp - run {run}', color='k', ha='right', fontsize=7)

    outpath = os.path.join(output_dir, 'full_canvas4.png')
    plt.tight_layout()
    LOG(INFO,f'Saving to {outpath}')
    fig.savefig(outpath, dpi=200)
    plt.close(fig)
        
       
def make_canvas5(dead0 = [[i for i in range(5)], [10-i for i in range(5)]],
                 dead1 = [[i for i in range(6)], [12-i for i in range(6)]],
                 dead2 = [[],[]],
                 dead3 = [[],[]],
                 dead4 = [[],[]],
                 dead5 = [[],[]],
                 dead6 = [[],[]],
                 reco0 = [[],[]],
                 reco1 = [[],[]],
                 reco2 = [[],[]],
                 reco3 = [[],[]],
                 reco4 = [[],[]],
                 reco5 = [[],[]],
                 reco6 = [[],[]],
                 window_size_reco = 1,
                 run = '?'
                 ):


     LOG(INFO,f'Creating canvas 5')
     
     fig, axes = plt.subplots(3,2,figsize=(8,9))
     
     axes = axes.flatten()
     plot_title = [f'Dead fraction IB - run {run}',f'Dead fraction OB - run {run}',f'Approx reco/hour, IB - run {run}',f'Approx reco/hour, OB - run {run}',f'IB', f'OB']
    

     markersize = 3
     markers = ['o', 's', 'h', 'd']
      
     # Dead time IB
     idx = 0
     ax = axes[idx]
     
     X=[np.array(dead0[0])/60, np.array(dead1[0])/60, np.array(dead2[0])/60]
     Y=[np.array(dead0[1]), np.array(dead1[1]), np.array(dead2[1])]
     
    
     for i in range(3):        
         try:             
             ax.plot(X[i],Y[i], label=f'L{i}', marker = markers[i], markersize=markersize)
             ax.set_title(plot_title[idx])
             ax.set_yscale('log')
             ax.legend()
             ax.set_xlabel('min')
         except Exception as e:
             
             LOG(WARNING,f'Exception while creating {plot_title[idx]}, L{i}: {e}')
   

     # Dead time OB
     idx = 1
     ax = axes[idx]
     X=[np.array(dead3[0])/60, np.array(dead4[0])/60, np.array(dead5[0])/60, np.array(dead6[0])/60]
     Y=[np.array(dead3[1]), np.array(dead4[1]), np.array(dead5[1]), np.array(dead6[1])]
    
     for i,l in enumerate([3,4,5,6]):
         try:
             ax.plot(X[i],Y[i], label=f'L{l}', marker = markers[i], markersize=markersize)
             ax.set_title(plot_title[idx])
             ax.legend()
             ax.set_yscale('log')
             ax.set_xlabel('min')
         except Exception as e:
             LOG(WARNING,f'Exception while creating {plot_title[idx]}, L{l}: {e}')


     # Reco rate IB
     idx = 2
     ax = axes[idx]
     
     X=[np.array(reco0[0])/60, np.array(reco1[0])/60, np.array(reco2[0])/60]
     Y=[np.array(reco0[1])*3600, np.array(reco1[1])*3600, np.array(reco2[1])*3600]
     
    
     for i in range(3):         
         try:            
             ax.plot(X[i],Y[i], label=f'L{i}', marker = markers[i], markersize=markersize)
             ax.set_title(plot_title[idx])
             ax.legend()
             ax.set_yscale('log')
         except Exception as e:            
             LOG(WARNING,f'Exception while creating {plot_title[idx]}, L{i}: {e}')
   

     # Reco rate OB
     idx = 3
     ax = axes[idx]
     X=[np.array(reco3[0])/60, np.array(reco4[0])/60, np.array(reco5[0])/60, np.array(reco6[0])/60]
     Y=[np.array(reco3[1])*3600, np.array(reco4[1])*3600, np.array(reco5[1])*3600, np.array(reco6[1])*3600]
     
    
     for i,l in enumerate([3,4,5,6]):
         try:
             ax.plot(X[i],Y[i], label=f'L{l}', marker = markers[i], markersize=markersize)
             ax.set_title(plot_title[idx])
             ax.legend()
             ax.set_yscale('log')
         except Exception as e:
             LOG(WARNING,f'Exception while creating {plot_title[idx]}, L{l}: {e}')

     # Scatter IB
     idx = 4
     ax = axes[idx]
     try:
         ax.scatter(dead0[1], [3600*r for r in reco0[1]],  marker=markers[0], s=12, label='L0', alpha=0.7)
         ax.scatter(dead1[1], [3600*r for r in reco1[1]],  marker=markers[1], s=12, label='L1', alpha=0.7)
         ax.scatter(dead2[1], [3600*r for r in reco2[1]],  marker=markers[2], s=12, label='L2', alpha=0.7)
         
         ax.set_xlabel("dead fraction")
         ax.set_ylabel("recovery rate")
         ax.set_title(plot_title[idx])
         ax.set_xscale('log')
         ax.set_yscale('log')
         #ax.grid(True)

         xmin, xmax = ax.get_xlim()
         ymin, ymax = ax.get_ylim()

         xref = (sum(dead0[1])+sum(dead1[1])+sum(dead2[1]))/(len(dead0[1])+len(dead1[1])+len(dead2[1]))
         yref = 3600*(sum(reco0[1])+sum(reco1[1])+sum(reco2[1]))/(len(reco0[1])+len(reco1[1])+len(reco2[1]))
         kLin = yref/xref
         kQuad = yref/(xref*xref)
         xmaxLin = min(xmax, ymax/kLin)
         xminLin = max(xmin, ymin/kLin)
         xmaxQuad = min(xmax, np.sqrt(ymax/kQuad))
         xminQuad = max(xmin, np.sqrt(ymin/kQuad))

         if xmaxLin > xmin and xmaxQuad > xmin:
             ax.plot([xminLin, xmaxLin], [kLin*xminLin, kLin*xmaxLin], linestyle='--', linewidth=0.5, color='grey', label='y = x', zorder=0)
             ax.plot([xminQuad, xmaxQuad], [kQuad*xminQuad*xminQuad, kQuad*xmaxQuad*xmaxQuad], linestyle='--', linewidth=0.5, color='grey', label='y = x', zorder=0) # valid only in log-log scale
         else:
             LOG(WARNING,f'Could not draw the reference line for scatter IB plot')

     except Exception as e:
         LOG(WARNING,f'Exception while creating correlation plot {plot_title[idx]}: {e}')


     # Scatter OB
     idx = 5
     ax = axes[idx]
     try:
         ax.scatter(dead3[1], [3600*r for r in reco3[1]],  marker=markers[0], s=12, label='L3', alpha=0.7)
         ax.scatter(dead4[1], [3600*r for r in reco4[1]],  marker=markers[1], s=12, label='L4', alpha=0.7)
         ax.scatter(dead5[1], [3600*r for r in reco5[1]],  marker=markers[2], s=12, label='L5', alpha=0.7)
         ax.scatter(dead6[1], [3600*r for r in reco6[1]],  marker=markers[3], s=12, label='L6', alpha=0.7)
         
         ax.set_xlabel("dead fraction")
         ax.set_ylabel("recovery rate")
         ax.set_title(plot_title[idx])
         ax.set_xscale('log')
         ax.set_yscale('log')
         #ax.grid(True)

         xmin, xmax = ax.get_xlim()
         ymin, ymax = ax.get_ylim()

         xref = (sum(dead3[1])+sum(dead4[1])+sum(dead5[1])+sum(dead6[1]))/(len(dead3[1])+len(dead4[1])+len(dead5[1])+len(dead6[1]))
         yref = 3600*(sum(reco3[1])+sum(reco4[1])+sum(reco5[1])+sum(reco6[1]))/(len(reco3[1])+len(reco4[1])+len(reco5[1])+len(reco6[1]))
         kLin = yref/xref
         kQuad = yref/(xref*xref)
         xmaxLin = min(xmax, ymax/kLin)
         xminLin = max(xmin, ymin/kLin)
         xmaxQuad = min(xmax, np.sqrt(ymax/kQuad))
         xminQuad = max(xmin, np.sqrt(ymin/kQuad))

         if xmaxLin > xmin and xmaxQuad > xmin:
             ax.plot([xminLin, xmaxLin], [kLin*xminLin, kLin*xmaxLin], linestyle='--', linewidth=0.5, color='grey', label='y = x', zorder=0)
             ax.plot([xminQuad, xmaxQuad], [kQuad*xminQuad*xminQuad, kQuad*xmaxQuad*xmaxQuad], linestyle='--', linewidth=0.5, color='grey', label='y = x', zorder=0) # valid only in log-log scale
         else:
             LOG(WARNING,f'Could not draw the reference line for scatter OB plot')
     except Exception as e:
         LOG(WARNING,f'Exception while creating correlation plot {plot_title[idx]}: {e}')

         

     # ===== Save full canvas ======
     outpath = os.path.join(output_dir,'full_canvas5.png')
     LOG(INFO,f'Saving to {outpath}')
     plt.tight_layout()
     fig.savefig(outpath, dpi = 150)
     plt.close(fig)
    
             

##########################################

def make_canvas22( 
        lanemap={k: np.arange(100)/100 for k in range(10)},
        idx = [[2,8], [10,20], [20,30]],
        offset_sec = 0, # first value on the x axis
        run = '?',
        spec0 = '',
        spec1 = ['','',''],
        spec2 = '',
        output_dir='.',

        # --- performance knobs ---
        final_dpi=500,                     # same dpi for figure and save → single render
        png_compress_level=1,              # faster PNG writing (slightly larger file)
        antialiased=False,                 # faster for huge images/lines
):
    mpl.use("Agg")
    plt.ioff()

    LOG(ERROR,f'canvas22')

    # ===== Data prep (no manual rebinning) =====
    keys   = np.asarray(list(lanemap.keys()))
    arrays = [lanemap[k] for k in keys]
    fulldata = np.asarray(arrays)         # shape: (nx, ny?) depends on your inputs
    F = fulldata.T                        # we index on rows below; cache transpose once

    # Time axis in minutes, with your offset and final 2s pad
    # (we’ll still pass extent to imshow so the scale is preserved)
    x_edges = (keys - keys[0]) * LHCOrbitNS * 1.0e-9 / 60.0 + offset_sec / 60.0
    x_edges = np.append(x_edges, x_edges[-1] + 2/60.0)
    x_is_sec = (x_edges[-1] - x_edges[0]) < 5
    x_label = 'time (sec)' if x_is_sec else 'time (min)'
    if x_is_sec:
        x_edges = x_edges * 60.0

    # ===== Figure / axes =====
    fig = plt.figure(figsize=(20, 10), dpi=final_dpi)
    gs = gridspec.GridSpec(2, 3, height_ratios=[1, 0.02], hspace=0.2, figure=fig)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]

    # Custom colormap (unchanged)
    colors = ['#ffffff', '#5773e0', '#7b98ef', '#9FBEFF', '#C8D8F1', '#F3C6B3', '#E0664F', '#A70E2A']
    new_cmap = mcolors.ListedColormap(colors)

    plot_titles = [f'Inner Barrel{spec1[0]}', f'Layer 3,4{spec1[1]}', f'Layer 5,6{spec1[2]}']
    hline_step = [9, 16, 28]

    def getlabel(b, l):
        if b == 0:
            if l < 12*9:
                return f'L0_{l//9}'
            elif l < 28*9:
                return f'L1_{(l-12*9)//9}'
            else:
                return f'L2_{(l-28*9)//9}'
        elif b == 1:
            if l < 24*16:
                return f'L3_{l//16}'
            else:
                return f'L4_{(l-24*16)//16}'
        else:
            if l < 42*28:
                return f'L5_{l//28}'
            else:
                return f'L6_{(l-42*28)//28}'

    # Precompute X/Y extents for imshow so axes show your physical bin edges
    # imshow expects extents as (xmin, xmax, ymin, ymax)
    xmin, xmax = x_edges[0], x_edges[-1]

    for i in range(3):
        # Slice rows once (no manual rebinning)
        row0, row1 = idx[i]
        data = F[row0:row1, :]                          # shape: (ny, nx)
        # Your original code prepended a zero row and shifted y-edges. Preserve visually:
        data = np.vstack([np.zeros((1, data.shape[1])), data])

        # Build y_edges compatible with the above data (len = rows+1)
        y_edges = np.arange(row1 - row0 + 1)
        y_edges = np.insert(y_edges, 0, 0 - hline_step[i])
        ymin, ymax = y_edges[0], y_edges[-1]

        ax = axes[i]

        # ---- FAST path: single raster image instead of millions of quads ----
        # Use 'nearest' to avoid costly filtering; keep exact bin look.
        im = ax.imshow(
            data,
            origin='lower',
            interpolation='nearest',
            aspect='auto',
            extent=(xmin, xmax, ymin, ymax),
            cmap=new_cmap,
            vmin=0, vmax=1
        )
        im.set_rasterized(True)  # harmless for PNG; ensures raster in vector backends
        #im.set_antialiased(antialiased)

        # ---- Draw gridlines in one artist via LineCollection (fast) ----
        # Vertical lines at each x-edge:
        # (only draw a reasonable subset if x_edges is huge; but here we keep them all.)
        v_segments = [((xe, ymin), (xe, 0)) for xe in x_edges]
        vcoll = LineCollection(v_segments, colors='grey', linewidths=0.2, antialiased=antialiased)
        ax.add_collection(vcoll)

        # Horizontal dashed lines every hline_step:
        y_positions = np.arange(0, data.shape[0]-1, hline_step[i])
        h_segments = [((xmin, y), (xmax, y)) for y in y_positions]
        hcoll = LineCollection(h_segments, colors='grey', linewidths=0.2, linestyles='--', antialiased=antialiased)
        ax.add_collection(hcoll)

        # Baseline at y=0:
        ax.axhline(0, color='k', linestyle='-', linewidth=0.5, antialiased=antialiased)

        # Text labels at mid-groups (keep minimal; text is expensive)
        for y in y_positions:
            ax.text(xmin, y + hline_step[i]/2, getlabel(i, y)+' ',
                    color='black', va='center', ha='right',
                    fontsize=6 if i < 2 else 5)

        ax.set_title(plot_titles[i])
        ax.set_xlabel(x_label)
        ax.set_yticklabels([])
        ax.tick_params(axis='y', which='both', size=0)

        # Tight data limits (imshow with extent already sets them, but be explicit):
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

    # ----- Colorbar (reuse the image mappable to avoid extra ScalarMappable) -----
    cbar_ax = fig.add_subplot(gs[1, 1])
    cbar_ax.tick_params(size=0)
    # Use the last 'im' (all share same cmap/norm). If you want exact control, keep a ref list.
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_ticks([0.0, 1.0])
    cbar.set_ticklabels([f"{t:.2f}" for t in [0.0, 1.0]])
    cbar.set_ticks(np.linspace(0, 1, len(colors)))  # if you actually want every color tick

    fig.suptitle(f'Run {run}{spec0}', fontsize=13, y=0.985)

    # Manual layout (fast; doesn’t trigger an extra “tight” pass)
    fig.subplots_adjust(left=0.04, right=0.98, top=0.925, bottom=0.035, wspace=0.1)

    # ----- Save ONCE at the same DPI (fastest path) -----
    name_spec = f'_{spec2}' if spec2 else ''
    outpath = os.path.join(output_dir, f'full_canvas2{name_spec}.png')

    # Faster PNG writing: lower compression, no transparency
    fig.savefig(
        outpath,
        dpi=final_dpi,               # same as figure dpi → avoids extra renderer
        transparent=False,
        pil_kwargs={"compress_level": png_compress_level, "optimize": False}
    )

    plt.close(fig)
    return outpath
         
         
     
         
         

if __name__ == '__main__':
    print('Making canvas5')
    make_canvas5()

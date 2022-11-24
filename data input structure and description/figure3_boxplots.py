#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
variables needed:
data_dict = dictionary containing all the dataframes of the datasets we want to make predictions on 
            (e.g.: {"data X":data_X_df, "data Y":data_y_df, ...})

signature_dict = dictionary with list of genes of the signatures and subtypes,
                 e.g. {"Signature X": (list_genes_signature_X, subtypes_related_to_signature_X),...}
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon 
import matplotlib.lines as mlines
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy import stats

###############################################################################
# Figure 3
###############################################################################

from statsmodels.sandbox.stats.multicomp import multipletests

all_pval_list = []
for s_genes in signature_dict.keys(): # genes, subtypes
    # first take signatures 
    g = signature_dict[s_genes][0]
    print("\n signatures: ", s_genes)
    
    for s_data in signature_dict.keys():
        # now, for each data:
        print("data: ", s_data)
        df = data_dict[s_data].copy()
        sub = signature_dict[s_data][1]
        
        df = df[set(df.columns).intersection(g)]     
        df.index  = sub
        
        if s_data=="Moffitt et al.":
            pval = []
            for i in df.columns:
                F, p = stats.ranksums(df.loc["Basal", i], df.loc["Classical", i])
                pval.append(p)

        elif s_data=="Collisson et al.":
            pval = []
            for i in df.columns:
                F, p = stats.kruskal(df.loc["Exocrine-like PDA", i], df.loc["QM-PDA", i], df.loc["Classical PDA", i])
                pval.append(p)

        elif s_data=="Bailey et al.":
            pval = []
            for i in df.columns:
                F, p = stats.kruskal(df.loc["ADEX", i], df.loc["Pancreatic Progenitor", i], df.loc["Immunogenic", i], df.loc["Squamous", i])
                pval.append(p)
        elif s_data=="Puleo et al.":
            pval = []
            for i in df.columns:
                F, p = stats.kruskal(df.loc["Desmoplastic", i], df.loc["ImmuneClassical", i], df.loc["PureBasal-like", i], df.loc["PureClassical", i], df.loc["StromaActivated", i])
                pval.append(p)
        
        a, correct_pval, b, c = multipletests(pvals=pval, alpha=0.05, method="bonferroni")
        all_pval_list.append(-np.log10(correct_pval))   
         
            
all_boxes = all_pval_list
fig, ax1 = plt.subplots(figsize=(10, 5))
fig.canvas.set_window_title('A Boxplot Example')
fig.subplots_adjust(left=0.075, right=0.95, top=0.9, bottom=0.25)

bp = ax1.boxplot(all_boxes, notch=0, sym='+', vert=1, whis=1.5)

plt.setp(bp['boxes'], color='grey', linewidth=1)
plt.setp(bp['whiskers'], color='grey', linewidth=1)
plt.setp(bp['medians'], color='grey', linewidth=1)
plt.setp(bp['fliers'], color='gainsboro', marker='o', linewidth=1, markersize=4, alpha=0.5)

ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

ax1.set_axisbelow(True)
ax1.set_title("Comparison of signatures expression level across subtypes", fontsize=16)
ax1.set_xlabel("")
ax1.set_ylabel('-log10(adjusted pvalue)', fontsize=14)
ax1.set_xticks(list(range(1,10)))

box_colors = ['orange', 'darkcyan', 'deeppink', 'blue']
num_boxes = 12
medians = np.empty(num_boxes)
for i in range(num_boxes):
    box = bp['boxes'][i]
    boxX = []
    boxY = []
    for j in range(5):
        boxX.append(box.get_xdata()[j])
        boxY.append(box.get_ydata()[j])
    box_coords = np.column_stack([boxX, boxY])
    # Alternate between Dark Khaki and Royal Blue
    ax1.add_patch(Polygon(box_coords, facecolor=box_colors[i % 4], alpha=0.5))
    # Now draw the median lines back over what we just filled in

ax1.set_xticklabels([])

import matplotlib.lines as mlines
gold = mlines.Line2D([], [], color='orange', marker='s', markersize=11, label='Moffitt et al.', alpha=0.5)
cyan = mlines.Line2D([], [], color='darkcyan', marker='s', markersize=11, label='Collisson et al.', alpha=0.5)
pink = mlines.Line2D([], [], color='deeppink', marker='s', markersize=11, label='Bailey et al.', alpha=0.5)
blue = mlines.Line2D([], [], color='blue', marker='s', markersize=11, label='Puleo et al.', alpha=0.5)
ax1.legend(handles=[gold, cyan, pink, blue], title="Subtypes", loc=1, fontsize=12, title_fontsize=12)

ax1.text(1.6, -1.3, "Moffitt et al.", fontsize=13)  
ax1.text(6.6, -1.3, "Collisson et al.", fontsize=13) 
ax1.text(9.6, -1.3, "Bailey et al.", fontsize=13) 
ax1.text(12.6, -1.3, "Puleo et al.", fontsize=13) 

ax1.axhline(-np.log10(0.05), color="darkred", linestyle='--', lw=1.2)

ax1.text(4.6, -2.1, "Signatures", fontsize=14) 
                     
for i in [4.5, 8.5, 12.5]:            
    ax1.axvline(i, color="grey", linestyle='--', lw=0.5)

ax1.set_yticks([0, 2, 4, 6, 8, 10])
ax1.set_yticklabels([0, 2, 4, 6, 8, 10], fontsize=12)
ax1.grid(color='grey', axis='y', linestyle='--', linewidth=0.3, alpha=0.7)

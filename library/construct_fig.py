#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##                                Library
##
##########################################################################################
##########################################################################################

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

from Bio import Phylo

import seaborn as sns ; sns.set_style("ticks") 

plt.rcParams['text.color'] = 'black'
plt.rcParams['svg.fonttype'] = 'none'  # Editable SVG text
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.weight'] = 'light'

##########################################################################################
##########################################################################################
##
##                                Functions
##
##########################################################################################
##########################################################################################

def proportion_proteobacteria(df_found, LIST_SYSTEMS, ax, column_wanted, size='xx-small'):

    """
    Function that plot the proportion of each systems found

    :param df_found: The dataframe with the information about the detected systems
    :type: pandas.Dataframe
    :param LIST_SYSTEMS: list of the name of all the systems of the analysis (ex : T2SS, T4P ...)
    :type: list of str
    :param ax: The axe on which to plot the sub figure
    :type: matplotlib.axes.Axes
    :param column_wanted: Name of the columns used to do the presence absence table
    :type: str
    :param size: Font size in points or as a string (e.g., 'large')
    :type: float or str
    :return: Nothing
    """

    # XXX Creation du dataframe intermediare
    df_figure = pd.DataFrame(0, index=['Proteobacteria (𝛼, 𝛽, 𝛾)','Rest'] ,columns=LIST_SYSTEMS)
    df_figure.loc['Proteobacteria (𝛼, 𝛽, 𝛾)'] = df_found[df_found \
                                                              .Phylum.isin(["Gammaproteobacteria","Alphaproteobacteria","Betaproteobacteria"])] \
                                                              .loc[:,column_wanted].value_counts()
    
    df_figure.loc['Rest'] = df_found[~df_found.Phylum.isin(["Gammaproteobacteria","Alphaproteobacteria","Betaproteobacteria"])] \
                                                .loc[:,column_wanted].value_counts()
    df_figure.fillna(0, inplace=True)

    # XXX Le plot avec les couleurs choisies
    colors = [(0.65098041296005249, 0.80784314870834351, 0.89019608497619629, 1.0), 
              (0.3997693305214246, 0.6478123867044262, 0.80273742044673246, 1.0)]
    
    plot = df_figure.T.plot(kind="bar", stacked=True, color=colors, rot=0, legend=False, ax=ax)

    # XXX Ajout des valeurs en haut du plot et on ne prend que les valeurs finale
    index_true_patches = int(len(plot.patches)/2)
    
    # Patches of the dark blue bars
    true_patches = plot.patches[index_true_patches:]
    # Patches of the light blue bars    
    false_patches = plot.patches[:index_true_patches]
    
    for p in true_patches:
        x=p.get_bbox().get_points()[:,0]
        y=p.get_bbox().get_points()[1,1] 
        
        # Change the value of y if y == 0
        y = y if (y != 0) else false_patches[true_patches.index(p)].get_bbox().get_points()[1,1]

        plot.annotate('{}'.format(int(y)), (x.mean(), y),
                ha='center', va='bottom', size=size)


    plot.tick_params(axis='x', which='both', bottom=False,top=False, labelbottom=False)
    plot.xaxis.set_label_text("")
    
    
    plot.tick_params(axis = 'y', 
                     which = 'both', 
                     labelsize = size, 
                     )
    
    plot.yaxis.set_label_text(label = "Detected {}".format(column_wanted),
                             fontdict = {'size': size})
    
    sns.despine(ax=ax, 
                bottom=True,
                trim=True)
    
    return

##########################################################################################
##########################################################################################

def heatmap_df(report_pivot, count_df, list_wanted, ax, size='xx-small', rotation=0, color_gradient='red') :
    
    '''
    Function to plot the heatmap on the figure based on the presence absence data

    :param report_pivot: Pivot table of the report that represent the presence absence table of the gene/systems 
    :type: pandas.DataFrame
    :param count_df: Number of the genomes for each phyla order in the same order as the heatmap
    :type: pandas.DataFrame
    :param list_wanted: List of all the genes/systems wanted to apear on the figure
    :type: list of str
    :param ax: The axe on which to plot the sub figure
    :type: matplotlib.axes.Axes
    :param size: Font size in points or as a string (e.g., 'large')
    :type: float or str
    :param rotation: The angle to which the label of the heatmap to rotate
    :type: int
    :param color_gradient:
    :type:
    :return: Nothing
    '''
    
    cmap = sns.light_palette(color_gradient, as_cmap=True)
    
    try_missing =  list(set(list_wanted) - set(report_pivot.index))
    
    if try_missing :
        for missing in try_missing :
            report_pivot.loc[missing] = 0
            
    df_annot = report_pivot.reindex(list_wanted)
        
    df_heatmap = df_annot.div(count_df.Count, axis=0)
    
    sns.heatmap(df_heatmap,
                cmap = cmap,
                linewidths = 1,
                linecolor = (0.3997693305214246, 0.6478123867044262, 0.80273742044673246, 1.0),
                annot = df_annot,
                annot_kws = {'color':'black', 'fontsize':size},
                cbar = False,
                ax = ax,
                fmt = "d",
                yticklabels = False,
               )
    
    ax.set_xticklabels(ax.get_xticklabels(), rotation=rotation)
    
    # The mesh is the figure itself here, so to change the facecolor of the cell in the heatmap
    # we need to parse the mesh as in the seaborn instance code
    # So now the 0 are white
    
    mesh = ax.collections[0]
    all_values_mesh = mesh.get_array()
    all_color_mesh = mesh.get_facecolors()
    new_color = []
    len_values = len(all_values_mesh)
    
    ax.collections[0].set_facecolor('none')
    for i in range(len_values) :
        if all_values_mesh[i] == 0 :
            new_color.append('white')
        else :
            new_color.append(all_color_mesh[i])

    mesh.set_facecolor(new_color)

    # Modify axis and ticks
    ax.xaxis.set_ticks_position('top')
    ax.tick_params(axis='x', 
                   which='both', 
                   labelsize=size,
                   length = 0,
                  )
    
    ax.xaxis.set_label_text("")
    ax.yaxis.set_label_text("")

    return

##########################################################################################
##########################################################################################

def number_genomes(ax, count_df, size='xx-small') :
    
    '''
    Function to plot the number of genome in the database by phyla

    :param ax: The axe on which to plot the sub figure
    :type: matplotlib.axes.Axes    
    :param count_df: Number of the genomes for each phyla order in the same order as the heatmap
    :type: pandas.DataFrame
    :param size: Font size in points or as a string (e.g., 'large')
    :type: float or str
    :return: Nothing    
    '''
    
    cmap = ListedColormap(['white'])
        
    df_count = count_df[['Count']].rename(columns={'Count':'Number\nof\ngenomes'})
    
    sns.heatmap(df_count,
                cmap = cmap,
                linewidths = 1,
                linecolor = (0.3997693305214246, 0.6478123867044262, 0.80273742044673246, 1.0),
                annot = df_count,
                annot_kws = {'color':'black', 'fontsize':size},
                cbar = False,
                ax = ax,
                fmt = ".0f",
                yticklabels = False,
               )
    
    # Modify axis and ticks
    ax.xaxis.set_ticks_position('top')
    ax.tick_params(axis='x', 
                   which='both', 
                   labelsize=size,
                   length = 0,
                  )
    
    ax.xaxis.set_label_text("")
    ax.yaxis.set_label_text("")
    
    #ax.set_aspect(0.3)

    return

##########################################################################################
##########################################################################################

def draw_tree_prok(tree_file, ax) :
    
    '''
    Function to plot the tree of the phyla in the side of the heatmap from a tree file

    :param tree_file: Path to the rooted tree file in newick format : ((A,B),(C, D), E)
    :type: str
    :param ax: The axe on which to plot the sub figure
    :type: matplotlib.axes.Axes    
    :return: Nothing
    '''
    

    with plt.rc_context(rc={'lines.linewidth': 1.5, "font.size":10}):
        
        tree_prok = Phylo.read(tree_file, "newick")

        tree_prok.name = ""

        Phylo.draw(tree_prok, 
                   axes=ax, 
                   do_show=False, 
                   show_confidence=False,
                  )


        xmax = np.max([t.get_position()[0] for t in ax.texts])
        biggest_leaves = np.max([len(t.get_text()) for t in ax.texts])
        leaves = []

        for t in ax.texts:
            X = t.get_position()[0]
            Y = t.get_position()[1]

            # handle leaves              

            ax.hlines(Y, X, xmax, "lightgrey", linestyle="dotted")
            # right align label
            t.set_x(xmax)
            # optional
            # t.set_text("{:>{width}}".format(t.get_text().strip(), width=biggest_leaves))
            # t.set_fontproperties("monospace")
            t.set_zorder(-1)
            leaves.append(t.get_text().strip().replace('_', ' ').replace('T.', 'Thermodesulfobium')) # get leaves as they are when reading the tree from top to bottom

            
    ax.tick_params(axis='both', 
                   which='both', 
                   bottom = False,
                   left = False, 
                   labelbottom = False,
                   labelleft = False
                  )
    
    ax.xaxis.set_label_text("")
    
    ax.yaxis.set_label_text("")
    
    sns.despine(ax=ax, 
                bottom=True,
                left=True)  

    ax.set_ylim(len(leaves)+0.5, 0.5)

    return leaves

##########################################################################################
##########################################################################################

def set_legend(count_df, ax, color_species=[], color_gradient='red', size='xx-small'):
    
    '''
    Function that plot the legend of the headmap and barplot on the bottom left

    :param count_df: Number of the genomes for each phyla order in the same order as the heatmap
    :type: pandas.DataFrame
    :param ax: The axe on which to plot the sub figure
    :type: matplotlib.axes.Axes  
    :param color_species: List of the color for the group of phyla used in the barplot
    :type: list of rgb color
    :param color_gradient: color for the gradient used in the heatmap
    :type: color in str, rgb or hex
    '''
    
    axins = inset_axes(ax,
                    width="50%",  # width = 50% of parent_bbox width
                    height="15%",  # height : 15%
                    loc='upper center')

    axins.set_title(label = 'Colour key (% of genomes with at least one genes)',
                    fontdict = {'fontsize':size})

    # do the gradient legend oan the first ax
    cmap = sns.light_palette(color_gradient, as_cmap=True)

    cbar = plt.colorbar(ScalarMappable(cmap=cmap), 
                 cax = axins,
                 orientation = 'horizontal',
                  )

    cbar.ax.tick_params(labelsize=size)

    cbar.set_ticks([0,0.25,0.5,0.75,1])
    cbar.set_ticklabels(["0", "25", "50", "75", '100'])


    ax.tick_params(axis='both', 
                      which='both', 
                      left = False, 
                      bottom = False, 
                      labelleft = False,
                      labelbottom = False,
                      )

    sns.despine(ax=ax, 
                left=True,
                bottom = True
               )  

    # ax.set_title('Colour key (% of\ngenomes with at least one genes)', size = 'x-small')

    
    # Do the square on the second
    mini_tab = pd.DataFrame(0, index=['Proteobacteria (𝛼, 𝛽, 𝛾)','Rest'] ,columns=["Count"])
    mini_tab.loc['Proteobacteria (𝛼, 𝛽, 𝛾)'] = count_df.loc['Gammaproteobacteria'] + count_df.loc['Betaproteobacteria'] + count_df.loc['Alphaproteobacteria']
    mini_tab.loc['Rest'] = count_df.sum() - mini_tab.loc['Proteobacteria (𝛼, 𝛽, 𝛾)']
    
    if color_species == [] :
        color = [(0.65098041296005249, 0.80784314870834351, 0.89019608497619629, 1.0), 
                 (0.3997693305214246, 0.6478123867044262, 0.80273742044673246, 1.0)]
    else :
        color = color_species
        
    legend = [r'Proteobacteria ($\alpha$, $\beta$, $\gamma$) ({} genomes)'.format(int(mini_tab.loc['Proteobacteria (𝛼, 𝛽, 𝛾)'].values)), 
                  'Rest of the dataset ({} genomes)'.format(int(mini_tab.loc['Rest'].values))]


    handles = [mpatches.Patch(color=color[i], label=legend[i]) for i in range(2)]

    ax.legend(handles = handles, 
              frameon = False,
              fontsize = size,
              loc = 'lower center'
              )

    return 

##########################################################################################
##########################################################################################

def Phylum_wanted_search(x, *wanted_phylum) :

    """
    Function to search the phylum in a lineage and return the name of the phylum found

    :param x: row of the dataframe
    :type: pandas.DataFrame
    :param wanted_phylum: list of all the phylum present in the tree
    :type: list
    :return: The phylum of the species
    :rtype: str
    """
    a = list(set(x.Lineage.split(';')) & set(wanted_phylum))

    if not a :
        a = list(set(x.Lineage.replace(' ', ';').split(';')) & set(wanted_phylum))

    return a[0] if a else 'phylum_not_in_list'

##########################################################################################
##########################################################################################

def make_pivot(report_like, columns_pivot, wanted_phylum, order_columns) :
    
    '''
    Function to create the presence absence from a table of detected gene/systemes that contains at least the columns needed

    Columns needed = 'Lineage', 'Species_Id' and column you want to count (e.g. Reference_system, Gene...)

    The Species_Id column is here to identify the detected genes/systems that belong to the same species

    :param report_like: Dataframe that contains information about the lineage, the species_id
    :type: pandas.DataFrame
    :param columns_pivot: Name of the columns to do the pivot table on
    :type: str
    :param wanted_phylum: Name of the phylum present in the tree
    :type: list of str
    :param order_columns: ordered list with the columns in the order wanted in the final figure
    :type: list of str
    :return: A presence/absence dataframe ordered as wanted by the user and as the tree is
    :rtype: pandas.DataFrame
    '''
    
    report_like = report_like[report_like.Phylum.isin(wanted_phylum)].reset_index(drop=True)
    
    report_like = report_like.drop_duplicates(['Species_Id', columns_pivot]).reset_index(drop=True)
    
    report_like.loc[:,'one'] = 1
    
    report_pivot = pd.pivot_table(report_like, values='one', index=['Phylum'], columns=[columns_pivot], aggfunc=np.sum, fill_value=0)
    
    if order_columns :
        report_pivot = report_pivot[order_columns]    
    
    return report_pivot

##########################################################################################
##########################################################################################

def create_countdf(INFO_TAB, wanted_phylum) :
    
    '''
    Function that will count the number of genomes in each phylum. It use a table that 
    have at least an Replicon_name present in the database, an Id or name for each species
    and the lineage for each genome

    Table with at least columns Replicon_name, Species_Id, Lineage

    :param INFO_TAB: Table with the information about the genome, the species and the lineage
    :type: pandas.DataFrame
    :param wanted_phylum: Name of the phylum present in the tree
    :type: list of str
    :return: The dataframe with the count of the phylum in the database that you used
    :rtype: pandas.DataFrame
    '''
    
    INFO_TAB.loc[:,'Phylum'] = INFO_TAB.apply(Phylum_wanted_search, args=(wanted_phylum), axis=1)
    INFO_TAB.loc[:,'Kingdom'] = INFO_TAB.Lineage.apply(lambda x : x.split(';')[0])
    
    count_df = INFO_TAB[INFO_TAB.Phylum.isin(wanted_phylum)] \
                            .drop_duplicates('Species_Id').groupby(['Kingdom', 'Phylum']) \
                            .Phylum.count().reset_index(name='Count')
    
    return count_df.set_index('Phylum')

##########################################################################################
##########################################################################################

def make_figure_distribution(file_output, report_like, column_pivot, order_columns, INFO_TAB, file_tree, size='xx-small', rotation=0, 
                            main_rowspan = 5, tree_colspan = 2, heatmap_colspan = 8, barplot_rowspan = 1, base_span = 1,
                            width=15, height=10, color='red', add_span=0) :
    
    '''
    report at least columns : Replicon_name, column you want to count (e.g. Reference_system, Gene...)
    INFO_TAB at least : Replicon_name, Species_Id, Lineage

    :param fill_output: Path to the output
    :type: str
    :param report_like: Dataframe that contains information about the lineage, the species_id
    :type: pandas.DataFrame
    :param columns_pivot: Name of the columns to do the pivot table on
    :type: str
    :param order_columns: ordered list with the columns in the order wanted in the final figure
    :type: list of str
    :param INFO_TAB: Table with the information about the genome, the species and the lineage
    :type: pandas.DataFrame    
    :param file_tree: Path to the rooted tree file in newick format : ((A,B),(C, D), E)
    :type: str    
    :param size: Font size in points or as a string (e.g., 'large')
    :type: float or str    
    :param rotation: The angle to which the label of the heatmap to rotate
    :type: int
    :param main_rowspan: Number of rows you wanted the tree and heatmap to be span
    :type: int
    :param tree_colspan: Number of columns you want the tree to be span
    :type: int
    :param heatmap_colspan: Number of columns you want the heatmap to be span
    :type: int
    :param barplot_rowspan: Number of rows you want the barplot to be span bellow the heatmap
    :type: int
    :param base_span: Number of columns or row a figure is span by default (mostly for the legend and the number of genomes)
    :type: int
    :param width: Width total of the figure in point
    :type: int 
    :param height: Height total of the figure in point
    :type: int
    :params color: Color of the heatmap
    :type: color in str, rgb or hex
    :params add_span: blank columns to add after the tree if needed
    :type: int
    :return: Nothing
    '''
    
    fig = plt.figure(figsize=(width,height))
    fig.set_tight_layout(True)

    nrow = barplot_rowspan + main_rowspan
    ncol = base_span + tree_colspan + heatmap_colspan

    ax1 = plt.subplot2grid(shape = (nrow, ncol), 
                           loc = (0,0), 
                           rowspan = main_rowspan, 
                           colspan = tree_colspan)

    ax2 = plt.subplot2grid(shape=(nrow, ncol), 
                           loc=(0,tree_colspan + add_span), 
                           rowspan = main_rowspan, 
                           colspan = base_span)

    ax3 = plt.subplot2grid(shape = (nrow, ncol), 
                           loc = (0,tree_colspan + add_span + base_span), 
                           rowspan = main_rowspan, 
                           colspan = heatmap_colspan)

    ax4 = plt.subplot2grid(shape = (nrow, ncol), 
                             loc = (main_rowspan,0), 
                             rowspan = base_span, 
                             colspan = tree_colspan + add_span + base_span)

    ax5 = plt.subplot2grid(shape=(nrow, ncol), 
                           loc=(main_rowspan,tree_colspan + add_span + base_span), 
                           rowspan = barplot_rowspan, 
                           colspan = heatmap_colspan)


    # Draw tree on ax1
    leaves_order = draw_tree_prok(file_tree, ax1)


    # Read INFO_TAB
    info_df = pd.read_table(INFO_TAB)

    info_df = info_df.replace({'Candidatus Gracilibacteria':'CPR',
                               'Candidatus Saccharibacteria':'CPR',
                               'Candidatus Bipolaricaulota':'Bipolaricaulota',
                               'Candidatus Omnitrophica':'Omnitrophica',
                               'Candidatus Cloacimonetes':'Cloacimonetes', 
                               'Candidatus Dependentiae':'Dependentiae'}, regex=True)

    # Read report_like
    report_df = pd.read_table(report_like)

    # Reorder count_df based on leaves order
    count_df = create_countdf(info_df, leaves_order)
    count_df = count_df.reindex(leaves_order)
    count_df = count_df[['Count']]
    
    
    if 'Phylum' not in report_df.columns.tolist() :
        set_columns = set(info_df.columns.tolist()) - set(report_df.columns.tolist())
        infoDF_columns = list(set_columns) + ['Replicon_name']
        report_df = report_df.merge(info_df.loc[:,infoDF_columns], on='Replicon_name')

    # Convert name to name in tree
    report_df.loc[:,'Phylum'] = report_df.apply(Phylum_wanted_search, args=(leaves_order),  axis=1)
    
    report_pivot = make_pivot(report_df, column_pivot, leaves_order, order_columns)

    # Draw number of genomes in dataset
    number_genomes(ax2, count_df, size=size)

    # draw big heatmap
    heatmap_df(report_pivot, count_df, leaves_order, ax3, size=size, rotation=rotation, color_gradient=color)

    # Draw legends

    set_legend(count_df, ax4, color_gradient=color, size=size)
    
    # Draw proportion
    proportion_proteobacteria(report_df, report_pivot.columns, ax5, column_pivot, size=size)

    
    #plt.tight_layout()
    plt.savefig(file_output)

    return

##########################################################################################
##########################################################################################

def make_figure_distribution_gridspec(file_output, report_like, column_pivot, order_columns, INFO_TAB, file_tree, size='xx-small', rotation=0, 
                            main_rowspan = 5, tree_colspan = 2, heatmap_colspan = 8, barplot_rowspan = 1, base_span = 1,
                            width=15, height=10, color='red', add_span=0) :
    
    '''
    report at least columns : Replicon_name, column you want to count (e.g. Reference_system, Gene...)
    INFO_TAB at least : Replicon_name, Species_Id, Lineage

    :param fill_output: Path to the output
    :type: str
    :param report_like: Dataframe that contains information about the lineage, the species_id
    :type: pandas.DataFrame
    :param columns_pivot: Name of the columns to do the pivot table on
    :type: str
    :param order_columns: ordered list with the columns in the order wanted in the final figure
    :type: list of str
    :param INFO_TAB: Table with the information about the genome, the species and the lineage
    :type: pandas.DataFrame    
    :param file_tree: Path to the rooted tree file in newick format : ((A,B),(C, D), E)
    :type: str    
    :param size: Font size in points or as a string (e.g., 'large')
    :type: float or str    
    :param rotation: The angle to which the label of the heatmap to rotate
    :type: int
    :param main_rowspan: Number of rows you wanted the tree and heatmap to be span
    :type: int
    :param tree_colspan: Number of columns you want the tree to be span
    :type: int
    :param heatmap_colspan: Number of columns you want the heatmap to be span
    :type: int
    :param barplot_rowspan: Number of rows you want the barplot to be span bellow the heatmap
    :type: int
    :param base_span: Number of columns or row a figure is span by default (mostly for the legend and the number of genomes)
    :type: int
    :param width: Width total of the figure in point
    :type: int 
    :param height: Height total of the figure in point
    :type: int
    :params color: Color of the heatmap
    :type: color in str, rgb or hex
    :params add_span: blank columns to add after the tree if needed
    :type: int
    :return: Nothing
    '''
    
    fig = plt.figure(figsize=(width,height))
    # fig.set_tight_layout(True)

    nrow = barplot_rowspan + main_rowspan
    ncol = base_span + tree_colspan + heatmap_colspan

    gs = gridspec.GridSpec(nrow, ncol)

    ax1 = fig.add_subplot(gs[:main_rowspan, :tree_colspan])
    ax2 = fig.add_subplot(gs[:main_rowspan, tree_colspan:(tree_colspan+base_span)])
    ax3 = fig.add_subplot(gs[:main_rowspan, -heatmap_colspan:])
    ax4 = fig.add_subplot(gs[-1, :(tree_colspan+base_span)])
    ax5 = fig.add_subplot(gs[-1, -heatmap_colspan:])

    # Draw tree on ax1
    leaves_order = draw_tree_prok(file_tree, ax1)

    # Read INFO_TAB
    info_df = pd.read_table(INFO_TAB)

    info_df = info_df.replace({'Candidatus Gracilibacteria':'CPR',
                               'Candidatus Saccharibacteria':'CPR',
                               'Microgenomates group incertae sedis':'CPR',
                               'Candidatus Bipolaricaulota':'Bipolaricaulota',
                               'Candidatus Cloacimonetes':'Cloacimonetes',
                               'Candidatus Dependentiae':'Dependentiae'}, regex=True)

    # Read report_like
    report_df = pd.read_table(report_like)

    # Reorder count_df based on leaves order
    count_df = create_countdf(info_df, leaves_order)
    count_df = count_df.reindex(leaves_order)
    count_df = count_df[['Count']]
    
    
    if 'Phylum' not in report_df.columns.tolist() :
        set_columns = set(info_df.columns.tolist()) - set(report_df.columns.tolist())
        infoDF_columns = list(set_columns) + ['Replicon_name']
        report_df = report_df.merge(info_df.loc[:,infoDF_columns], on='Replicon_name')
    else :
        report_df.loc[:,'Phylum'] = report_df.apply(Phylum_wanted_search, args=(leaves_order),  axis=1)
    
    report_pivot = make_pivot(report_df, column_pivot, leaves_order, order_columns)

    # Draw number of genomes in dataset
    number_genomes(ax2, count_df, size=size)

    # draw big heatmap
    heatmap_df(report_pivot, count_df, leaves_order, ax3, size=size, rotation=rotation, color_gradient=color)
    
    # Draw legends
    set_legend(count_df, ax4, color_gradient=color, size=size)

    # Draw proportion
    proportion_proteobacteria(report_df, report_pivot.columns, ax5, column_pivot, size=size)

    
    #plt.tight_layout()
    gs.tight_layout(fig)
    plt.savefig(file_output)

    return

##########################################################################################
##########################################################################################

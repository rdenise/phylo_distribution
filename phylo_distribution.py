#!/usr/bin/env python3
# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##                                Library
##
##########################################################################################
##########################################################################################

import argparse
from textwrap import dedent
import sys, os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'library'))

from utils import create_folder
from construct_fig import make_figure_distribution

##########################################################################################
##########################################################################################
##
##                                Main
##
##########################################################################################
##########################################################################################

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     description=dedent("""Phylegentic Dictribution Figure Creation""") )

general_option = parser.add_argument_group(title = "General input dataset options")
general_option.add_argument("-r",'--report',
                            metavar="<REPORT>",
                            dest="report",
                            help="Tabulated File that contain the detection of the genes/systems with the columns Replicon_name \
                                  (name of the replicon [for exemple each chromosome of a species is a different replicon]) \
                                  and the columns that you want to see the distribution",
                            required=True)
general_option.add_argument("-i",'--infotab',
                            metavar="<INFO_TAB>",
                            dest="infotab",
                            help="Name of the tabulated file that contain the lineage information for each replicon. Need to have a columns\
                                 name Replicon_name (same as in the report), Species_Id (could be the name of the species, but need \
                                 to be unique for each genome), Lineage (full lineage of the genome)",
                            required=True)
general_option.add_argument("-c",'--columns_distribution',
                            metavar="<NAME OF THE COLUMNS OF INTEREST>",
                            dest="columns_distribution",
                            type=str,
                            help="The name of the column in the report that will be use to do the distribution. The column contains the annotation for the sequence",
                            required=True)
general_option.add_argument("-t",'--tree_file',
                            metavar="<FILE_TREE>",
                            dest="tree_file",
                            help="The name of tree file of the phyla in newick format",
                            required=True)
general_option.add_argument("-o",'--output',
                            default=None,
                            dest="output",
                            metavar='<FOLDER>',
                            help="Name of the figure where to put the database (default: $PWD/phylo_distribution.pdf)")
general_option.add_argument("-order",'--order_columns',
                            metavar=("<GENE_NAME/SYSTEM_NAME>"),
                            nargs='+',
                            dest="order_columns",
                            help="Ordered list of the genes/systems you want in the heatmap (default: The order of appearance in the report)",
                            default='')



figure_option = parser.add_argument_group(title = "Figure options")
figure_option.add_argument("-size",'--size_font',
                            metavar="<SIZE>",
                            dest="size_font",
                            default=8,
                            help="Font size in points or as a string (e.g., 'large') (default: 8pt)")
figure_option.add_argument("-main_rowspan",'--main_rowspan',
                            metavar="<INT>",
                            dest="main_rowspan",
                            type=int,
                            default=5,
                            help="Number of row to span the main figure (tree and heatmap) (default: 5)")
figure_option.add_argument("-tree_colspan",'--tree_colspan',
                            metavar="<INT>",
                            dest="tree_colspan",
                            type=int,
                            default=2,
                            help="Number of columns to span the tree (default: 2)")
figure_option.add_argument("-heatmap_colspan",'--heatmap_colspan',
                            metavar="<INT>",
                            dest="heatmap_colspan",
                            type=int,
                            default=8,
                            help="Number of columns to span the heatmap (default: 8)")
figure_option.add_argument("-barplot_rowspan",'--barplot_rowspan',
                            metavar="<INT>",
                            dest="barplot_rowspan",
                            type=int,
                            default=1,
                            help="Number of rows to span the barplot under the heatmap (default: 1)")
figure_option.add_argument("-base_span",'--base_span',
                            metavar="<INT>",
                            dest="base_span",
                            type=int,
                            default=1,
                            help="Number of rows/columns to span that are not cited above (defaul: 1)")
figure_option.add_argument("-height",'--height',
                            metavar="<INT>",
                            dest="height",
                            type=int,
                            default=10,
                            help="Total height of the figure in point (default: 10)")
figure_option.add_argument("-width",'--width',
                            metavar="<INT>",
                            dest="width",
                            type=int,
                            default=12,
                            help="Total width of the figure in point (default: 12)")
figure_option.add_argument("-p",'--prettify',
                            action='store_true',
                            dest="prettify",
                            help='Allow the use of a try to prettify the figure in case the leaf tree are in the number of genome',
                            default=False)
figure_option.add_argument("-color",'--color_heatmap',
                            metavar=("<COLOR>"),
                            dest="color_heatmap",
                            default='red',
                            help="Color of the gradient in the heatmap (default: red)")

##########################################################################################

args = parser.parse_args()

##########################################################################################

if args.output :
    OUTPUT = args.output
else :
    OUTPUT = os.path.join(os.getcwd(), 'phylo_distribution.pdf')

##########################################################################################

create_folder(os.path.dirname(OUTPUT))

##########################################################################################

REPORT = args.report
COLUMN = args.columns_distribution
INFO_TAB = args.infotab
TREE_FILE = args.tree_file
ORDER_COLUMNS = args.order_columns

##########################################################################################

SIZE = args.size_font
MAIN_ROWSPAN = args.main_rowspan
TREE_COLSPAN = args.tree_colspan
HEATMAP_COLSPAN = args.heatmap_colspan
BARPLOT_ROWSPAN = args.barplot_rowspan
BASE_SPAN = args.base_span
HEIGHT = args.height
WIDTH = args.width
COLOR = args.color_heatmap
PRETTIFY = args.prettify

##########################################################################################
# Needed in case of bugs

def main():

    make_figure_distribution(file_output = OUTPUT,
                         report_like = REPORT, 
                         column_pivot = COLUMN, 
                         order_columns = ORDER_COLUMNS,
                         INFO_TAB = INFO_TAB,
                         file_tree = TREE_FILE,
                         size=SIZE,
                         main_rowspan = MAIN_ROWSPAN, tree_colspan = TREE_COLSPAN, 
                         heatmap_colspan = HEATMAP_COLSPAN, barplot_rowspan = BARPLOT_ROWSPAN, 
                         base_span = BASE_SPAN,
                         prettify=PRETTIFY, height=HEIGHT, width=WIDTH
                        )

##########################################################################################

if __name__ == "__main__":
    main()

#!/usr/local/bin/python
"""
***********************************************
- PROGRAM: draw_features_heatmap.py
- CONTACT: Gaurav Jain(gaurav.jain@dzne.edu)
***********************************************
"""
print (__doc__)

import argparse
import os.path
import sys
import textwrap
import string
import re
import numpy  as np
import matplotlib as mpl
#mpl.use('Agg') # to use matplotlib without X11
mpl.use('pdf') # to use matplotlib without X11
import matplotlib.pyplot as plt
import math
import time
import subprocess
from gjainPyLib import *     # import all the functions from the Gaurav`s python scripts/library	
from numpy import loadtxt, dtype,float32
from termcolor import colored
from collections import *
from scipy import stats
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import numpy.ma as ma
################ USER CONFIGURATION ###################
# 1)  Save matplotlib figure with text as vectorized text and not as a path
#     The matplotlib documentation in http://matplotlib.org/users/customizing.html 
#     states that the default value for the parameter svg.fonttype is 'path', 
#     which means that characters will be converted to paths when exporting to svg format.
matplotlib.rcParams['svg.fonttype']    = 'none'
matplotlib.rcParams['pdf.fonttype']    =  42
matplotlib.rcParams['patch.edgecolor'] = 'white'
plt.style.use('seaborn-white')
#######################################################

def main():
    # Get input options
    args         = check_options()
    input_file   = args.input_file
    output_file  = args.output_file
    l2fc         = args.l2fc
    padj         = args.padj
    base_mean    = args.base_mean
    xclustering  = args.xclustering
    yclustering  = args.yclustering
    rnaseq       = args.rnaseq
    simpleDE     = args.simpleDE
    column_order = args.column_order
    if column_order:
        column_order = column_order.split(",")
        print("Column Order: ")
        for c in column_order:
            print("\t- {0}".format(c))
    
    # Get the differentially expressed features heatmaps
    get_feature_enrichment(input_file, get_file_info(output_file)[3], l2fc, padj, base_mean, xclustering, yclustering, rnaseq, column_order, simpleDE)

################ USER DEFINED FUNCTIONS ###################
def get_feature_enrichment(input_file, output_file, l2fc_param, padj_param, basm_param, xc, yc, rnaseq, column_order, simpleDE):
    ''' Parse the input file and get the log2foldchange for each miRNA '''

    # Get the information in a dictionary
    feature_enrichment = defaultdict(float)
    with open(input_file,'rU') as fg:
        # head input/DE/Blood_DifferentialExpression.csv | cut -f1-5
        # feature	treatmentMean	controlMean	baseMean	log2FoldChange	lfcSE	stat	pvalue	padj	X2wHCCA3_mm_smallrna_sr_Cemil
        # mmu-miR-301a-3p	1245.33333333	1242.28571429	8977.0584113	-0.000465075976182	0.0648336184922	-0.00717337682205	0.994276522471	0.998959090355	8584.89149252
        # mmu-miR-412-5p	18.8333333333	18.5714285714	131.6066911	-0.000127955807554	0.098081456794	-0.00130458714355	0.998959090355	0.998959090355	123.861151712
        # ....

        # Get the header and save it
        header = fg.readline().strip()
        
        # Get the output file handle
        #fo = open(output_file + "_differentially_expressed.txt" , 'w')
        fo = open("{0}_differentially_expressed.txt".format(get_file_info(output_file)[3]), 'w')

        # Write the filtered smallRNAs to the output file
        fo.write(header + "\n")

        # Initialize vars 
        total_smallrna = 0
        skip_iv = 0

        # get the zmatrix
        zmatrix = list()
        row_header=list()
        if simpleDE:
            # 0:feature	1:baseMean	2:log2FoldChange	3:lfcSE	4:stat	5:pvalue	6:padj  7:sample_9097_mm_rna_sr_Raaz_Y_1_S3_L001_R1_001
            column_header = header.split("\t")[7:]
            if column_order:
                column_header = column_order
        elif rnaseq:
            # 0:GeneName	1:feature	2:baseMean	3:log2FoldChange	4:lfcSE	5:stat	6:pvalue	7:padj	8:Sample_9097_mm_rna_sr_Raaz_Y_1_S3_L001_R1_001
            column_header = header.split("\t")[8:]
            if column_order:
                column_header = column_order
        else:
            #  0:feature	1:treatmentMean	2:controlMean	3:baseMean	4:log2FoldChange	5:lfcSE	6:stat	7:pvalue	8:padj   9:Sample_*
            column_header = header.split("\t")[9:]

        # Loop through rest of the file
        for line in useful_lines(fg):
            featureid        = str(line.split("\t")[0])
            if simpleDE:
                # 0:feature	1:baseMean	2:log2FoldChange	3:lfcSE	4:stat	5:pvalue	6:padj  7:sample_9097_mm_rna_sr_Raaz_Y_1_S3_L001_R1_001
                baseMean         = line.split("\t")[1]
                log2FoldChange   = line.split("\t")[2]
                padj             = line.split("\t")[6]
                de_between       = line.split("\t")[7:]
            elif rnaseq:
                # 0:GeneName	1:feature	2:baseMean	3:log2FoldChange	4:lfcSE	5:stat	6:pvalue	7:padj	8:Sample_9097_mm_rna_sr_Raaz_Y_1_S3_L001_R1_001
                baseMean         = line.split("\t")[2]
                log2FoldChange   = line.split("\t")[3]
                padj             = line.split("\t")[7]
                de_between       = line.split("\t")[8:]
            else:
                # 0:feature	1:treatmentMean	2:controlMean	3:baseMean	4:log2FoldChange	5:lfcSE	6:stat	7:pvalue	8:padj	9:X1743ACC_mm_smallrna_sr_magda_E_4
                baseMean         = line.split("\t")[3]
                log2FoldChange   = line.split("\t")[4]
                padj             = line.split("\t")[8]
                de_between       = line.split("\t")[9:]

            zde_between      = stats.zscore(np.array(de_between, dtype=float), ddof=1, axis=0) # axis=0 columns, axis=1 rows
            
            # Skip mirs with invalid log fold change i.e. "NA"
            if log2FoldChange == 'NA' or baseMean == 'NA' or padj == 'NA' or log2FoldChange == '' or padj == '':
                skip_iv += 1
                continue
            log2FoldChange = np.float(log2FoldChange)

            # Skip mirs with log2 fold change that are between (-0.5, 0.5)
            if -1*l2fc_param < log2FoldChange < l2fc_param:
                skip_iv += 1
                continue

            # Filter for baseMean < 100 
            if float(baseMean) < basm_param:
                skip_iv += 1
                continue

            # Filter for padjusted > 0.05
            if padj != '':
                if float(padj) > padj_param:
                    skip_iv += 1
                    continue
            else:
                skip_iv += 1
                continue

            # Get the mirs enrichment in the datastructure 
            feature_enrichment[featureid] = de_between
            zmatrix.append(zde_between)
            row_header.append(featureid)
            total_smallrna += 1
            
            # Write the filtered smallRNAs to the output file
            fo.write(line + "\n")

        # Run clustering if the zmatrix is not empty
        if np.array(zmatrix).size:
            print("\t- Drawing the heatmap...", file=sys.stderr)
            if xc:  
                row_method    = 'average'
            else:
                row_method    = None

            if yc:
                column_method = 'single'
            else:
                column_method    = None

            hierarchical_cluster_heatmap(np.array(zmatrix, dtype=float), row_header, column_header, output_file, row_method, column_method)

            # Print the filtered data summary
            print("\n- Conditions used:\n\t- log2FC    <> {0}\n\t- basemean  >= {1}\n\t- padjusted <= {2}\n\t- xClustering: {3}\n\t- yClustering: {4}".format(l2fc_param, basm_param, padj_param, xc, yc))
            print("\n- Total features : {0}\n". format(total_smallrna))
        
            # Close the output file handles
            fo.close()
        else:
            print("\n *********** NOTE **********")
            print("- Threshold conditions used:\n\t- log2FC    <> {0}\n\t- basemean  >= {1}\n\t- padjusted <= {2}".format(l2fc_param, basm_param, padj_param))
            print("\n- There is no data left after applying the filtering conditions. Your cut off values are too stringent. Please check threshold.")
            print(" ***************************")

def hierarchical_cluster_heatmap(x, row_header, column_header, filename="test", row_method='average', column_method='single', row_metric='cityblock', column_metric='euclidean'):    
    """
    This below code is based in large part on the protype methods:
    http://old.nabble.com/How-to-plot-heatmap-with-matplotlib--td32534593.html
    http://stackoverflow.com/questions/7664826/how-to-get-flat-clustering-corresponding-to-color-clusters-in-the-dendrogram-cre

    x is an m by n ndarray, m observations, n genes
    """

    ### Scale the Matplotlib window size
    default_window_hight = 40
    default_window_width = 40
    fig = plt.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
    color_bar_w = 0.0001 ### Sufficient size to show
        
    ## calculate positions for all elements
    # ax1, placement of dendrogram 1, on the left of the heatmap
    #if row_method != None: w1 = 
    [ax1_x, ax1_y, ax1_w, ax1_h] = [0.05,0.22,0.2,0.6]   ### The second value controls the position of the matrix relative to the bottom of the view
    width_between_ax1_axr = 0.004
    height_between_ax1_axc = 0.004 ### distance between the top color bar axis and the matrix
    
    # axr, placement of row side colorbar
    [axr_x, axr_y, axr_w, axr_h] = [0.31,0.1,color_bar_w,0.6] ### second to last controls the width of the side color bar - 0.015 when showing
    axr_x = ax1_x + ax1_w + width_between_ax1_axr
    axr_y = ax1_y; axr_h = ax1_h
    width_between_axr_axm = 0.004

    # axc, placement of column side colorbar
    [axc_x, axc_y, axc_w, axc_h] = [0.4,0.63,0.5,color_bar_w] ### last one controls the hight of the top color bar - 0.015 when showing
    axc_x = axr_x + axr_w + width_between_axr_axm
    axc_y = ax1_y + ax1_h + height_between_ax1_axc
    height_between_axc_ax2 = 0.004

    # axm, placement of heatmap for the data matrix
    [axm_x, axm_y, axm_w, axm_h] = [0.4,0.9,2.5,0.5]
    axm_x = axr_x + axr_w + width_between_axr_axm
    axm_y = ax1_y; axm_h = ax1_h
    axm_w = axc_w

    # ax2, placement of dendrogram 2, on the top of the heatmap
    [ax2_x, ax2_y, ax2_w, ax2_h] = [0.3,0.72,0.6,0.15] ### last one controls hight of the dendrogram
    ax2_x = axr_x + axr_w + width_between_axr_axm
    ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
    ax2_w = axc_w

    # axcb - placement of the color legend
    [axcb_x, axcb_y, axcb_w, axcb_h] = [0.07,0.88,0.18,0.015]

    # Compute and plot top dendrogram
    if column_method != None:
        print("\n- Performing hiearchical clustering for columns using:\n\t- method: {0}\n\t- metric: {1}".format(column_method, column_metric))
        start_time = time.time()
        #d2 = dist.pdist(x.T)
        d2 = dist.pdist(np.transpose(x))
        D2 = dist.squareform(d2)
        ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=False)
        Y2 = sch.linkage(D2, method=column_method, metric=column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
        Z2 = sch.dendrogram(Y2)
        ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
        ax2.set_xticks([]) ### Hides ticks
        ax2.set_yticks([])
        time_diff = str(round(time.time()-start_time,1))
        print("\n\t- Column clustering completed in %s seconds" % time_diff)
    else:
        ind2 = ['NA']*len(column_header) ### Used for exporting the flat cluster data
        
    # Compute and plot left dendrogram.
    if row_method != None:
        print("\n- Performing hiearchical clustering for rows using:\n\t- method: {0}\n\t- metric: {1}".format(row_method, row_metric))
        start_time = time.time()
        d1 = dist.pdist(x)
        D1 = dist.squareform(d1)  # full matrix
        ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=False) # frame_on may be False
        Y1 = sch.linkage(D1, method=row_method, metric=row_metric) ### gene-clustering metric - 'average', 'single', 'centroid', 'complete'
        #Z1 = sch.dendrogram(Y1, orientation='right')
        Z1 = sch.dendrogram(Y1, orientation='left')
        ind1 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
        ax1.set_xticks([]) ### Hides ticks
        ax1.set_yticks([])
        time_diff = str(round(time.time()-start_time,1))
        print("\n\t- Row clustering completed in %s seconds" % time_diff)
    else:
        ind1 = ['NA']*len(row_header) ### Used for exporting the flat cluster data
        
    # Plot distance matrix.
    axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h],aspect="auto")  # axes for the data matrix
    xt = x
    if column_method != None:
        idx2 = Z2['leaves'] ### apply the clustering for the array-dendrograms to the actual matrix data
        xt = xt[:,idx2]
        ind2 = ind2[idx2] ### reorder the flat cluster to match the order of the leaves the dendrogram
    if row_method != None:
        idx1 = Z1['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
        xt = xt[idx1,:]   # xt is transformed x
        ind1 = ind1[idx1] ### reorder the flat cluster to match the order of the leaves the dendrogram

    ### taken from http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894
    #im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cmap, norm=norm) ### norm=norm added to scale coloring of expression with zero = white or black

    buworb_cmap=mpl.colors.LinearSegmentedColormap.from_list('buworb',colors=['darkblue','blue','lightblue','white','orange','red','darkred'])
    #buworb_cmap.set_bad(color = '#7cfc00', alpha = 1.) # 0.82=grey
    buworb_cmap.set_bad(color = 'grey', alpha = 1.) # 0.82=grey
    plt.register_cmap(cmap=buworb_cmap)
    
    # Set x-limit and y-limit
    plt.xlim(xmax=np.shape(xt)[1])
    plt.ylim(ymax=np.shape(xt)[0])
    
    maxmin = max(abs(np.nanmin(xt)), abs(np.nanmax(xt)))
    cmin   = -1 * maxmin
    cmax   = maxmin

    #---------------------------------------------------------------#
    '''
    # for specific purpose when masked value of -1000 is used.
    maxmin = abs(np.nanmax(xt))
    cmin   = -1 * maxmin
    cmax   = maxmin
    # Set masked elements in different color. So mask it and plot it Lawn Green color
    xt = np.ma.masked_invalid(ma.masked_values(xt,-1000))
    '''
    #---------------------------------------------------------------#

    pclm = axm.pcolormesh(xt,cmap='buworb', vmin=float(cmin), vmax=float(cmax))

    ### Orient the ticks outsite
    axm.tick_params(axis='y', direction='out')
    axm.tick_params(axis='x', direction='out')

    ### Hides left and top ticks
    axm.yaxis.tick_right()
    axm.xaxis.tick_bottom()

    if row_method != None:
        # For clustered rows
        if len(row_header) > 300: ### Don't show gene names when more than 300 rows
            print("\t- Skipping row names since they are greater than 300: {0}".format(len(row_header)))
    else:
        # When not clustering rows
        if len(row_header) > 100: ### Don't show gene names when more than 100 rows
            print("\t- Skipping row names since they are greater than 100: {0}".format(len(row_header)))

        
    # Get the new row headers
    new_row_header=[]
    for i in range(x.shape[0]):
        if row_method != None:
            # For clustered rows
            if len(row_header) <= 300: ### Don't show gene names when more than 300 rows
                new_row_header.append(row_header[idx1[i]])
        else:
            # When not clustering rows
            if len(row_header) <= 100: ### Don't show gene names when more than 100 rows
                new_row_header.append(row_header[i])

    # Get the new columns headers
    new_column_header=[]
    for i in range(x.shape[1]):
        if column_method != None:
            # For clustered columns
            new_column_header.append(column_header[idx2[i]])
        else: 
            # When not clustering columns
            new_column_header.append(column_header[i])
    

    ##############################################
    # put the major ticks at the middle of each cell
    axm.set_yticks(np.arange(x.shape[0])+0.5, minor=False)
    axm.set_xticks(np.arange(x.shape[1])+0.5, minor=False)
    
    # Want a more natural, table-like display
    #axm.invert_yaxis()
    
    axm.set_yticklabels(new_row_header   , minor=False)
    axm.set_xticklabels(new_column_header, minor=False, rotation=90, verticalalignment='top')

    ##############################################
    # Plot color legend
    norm = mpl.colors.Normalize(vmin=float(cmin), vmax=float(cmax))
    axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
    cb = mpl.colorbar.ColorbarBase(axcb, cmap=buworb_cmap, orientation='horizontal', norm=norm)
    axcb.tick_params(labelsize=10) 
    axcb.set_title("ColorKey", size=12)
    cb.set_label("Z-Score (normalized Counts)", size=12)
    mtitle = get_file_info(filename)[1]
    plt.suptitle(mtitle, fontsize=20)  

    #exportFlatClusterData(filename + '.pdf', new_row_header,new_column_header,xt,ind1,ind2)

    ### Render the graphic
    if len(row_header)>100 or len(column_header)>100:
        plt.rcParams['font.size'] = 3
    else:
        plt.rcParams['font.size'] = 10

    # Get output pdf dir
    pdfdir = "{0}/pdf".format(get_file_info(filename)[0])
    create_dir(pdfdir)

    # Save the plots
    plt.savefig("{0}.png"    .format(filename)        ,  bbox_inches='tight')
    plt.savefig("{0}/{1}.pdf".format(pdfdir, get_file_info(filename)[1]),  bbox_inches='tight')
    plt.close()

def getColorRange(x):
    """ Determines the range of colors, centered at zero, for normalizing cmap """
    vmax=x.max()
    vmin=x.min()
    if vmax<0 and vmin<0: direction = 'negative'
    elif vmax>0 and vmin>0: direction = 'positive'
    else: direction = 'both'
    if direction == 'both':
        vmax = max([vmax,abs(vmin)])
        vmin = -1*vmax
        return vmax,vmin
    else:
        return vmax,vmin
    

def exportFlatClusterData(filename, new_row_header,new_column_header,xt,ind1,ind2):
    """ Export the clustered results as a text file, only indicating the flat-clusters rather than the tree """
    
    filename = string.replace(filename,'.pdf','.cluster')
    export_text = open(filename,'w')
    column_header = string.join(['UID','row_clusters-flat']+new_column_header,'\t')+'\n' ### format column-names for export
    export_text.write(column_header)
    column_clusters = string.join(['column_clusters-flat','']+ list(map(str, ind2)),'\t')+'\n' ### format column-flat-clusters for export
    export_text.write(column_clusters)
    
    ### The clusters, dendrogram and flat clusters are drawn bottom-up, so we need to reverse the order to match
    new_row_header = new_row_header[::-1]
    xt = xt[::-1]
    
    ### Export each row in the clustered data matrix xt
    i=0
    for row in xt:
        export_text.write(string.join([new_row_header[i],str(ind1[i])]+list(map(str, row)),'\t')+'\n')
        i+=1
    export_text.close()
    
    ### Export as CDT file
    filename = string.replace(filename,'.cluster','.cdt')
    export_cdt = open(filename,'w')
    column_header = string.join(['UNIQID','NAME','GWEIGHT']+new_column_header,'\t')+'\n' ### format column-names for export
    export_cdt.write(column_header)
    eweight = string.join(['EWEIGHT','','']+ ['1']*len(new_column_header),'\t')+'\n' ### format column-flat-clusters for export
    export_cdt.write(eweight)
    
    ### Export each row in the clustered data matrix xt
    i=0
    for row in xt:
        export_cdt.write(string.join([new_row_header[i]]*2+['1']+list(map(str, row)),'\t')+'\n')
        i+=1
    export_cdt.close()

def useful_lines(f):
    ''' Filter out useless lines from the blast output file '''
    for line in f:
        line = line.strip()
        if line.startswith('#'):
            continue
        if not line:
            continue
        yield line

def check_options():
    ''' Checks the options to the program '''

    # Create parser object
    parser = argparse.ArgumentParser(add_help=True,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog=textwrap.dedent('''\
        ----------------- SAMPLE USAGE ------------------
        - python scripts/draw_features_heatmap.py -if=output/filtered_data/deseq2/mirna/results_DEseq2/ACC18mwt_over_ACC4mwt_DE_RESULTS.txt -of=output/filtered_data/deseq2/mirna/heatmaps_DEseq2/ACC18mwt_over_ACC4mwt_DE_heatmaps.txt -xc -yc
        - python scripts/draw_features_heatmap.py -if=output/filtered_data/deseq2/mirna/results_DEseq2/Blood18mwt_over_Blood4mwt_DE_RESULTS.txt -of=output/filtered_data/deseq2/mirna/heatmaps_DEseq2/Blood18mwt_over_Blood4mwt_DE_heatmaps.txt -lf=0.25 -pj=0.05 -bm=25
        -------------------------------------------------
        CONTACT: 
        	Gaurav Jain
        	gaurav.jain@dzne.de
        -------------------------------------------------
        '''))

    # Add arguments 
    parser.add_argument("-if", metavar="--ifile" , help="*SmallRNA DEseq2 file", dest="input_file"  , type=str, required=True)
    parser.add_argument("-of", metavar="--ofile" , help="*Output file name"    , dest="output_file" , type=str, required=True)
    parser.add_argument("-lf", metavar="--l2fc"  , help=" Log2 Fold change (Default <> 0.5 )", dest="l2fc", type=float, default=0.5)
    parser.add_argument("-pj", metavar="--padj"  , help=" Adjusted pvalue  (Default <= 0.05)", dest="padj", type=float, default=0.05)
    parser.add_argument("-bm", metavar="--bmean" , help=" Base Mean        (Default >= 100 )", dest="base_mean", type=int, default=100)
    parser.add_argument("-co", metavar="--colod" , help=" Column order. Ex: A1,A2,A4,A3,C3,C2,C1", dest="column_order", type=str)
    parser.add_argument('-xc', "--xclustering"   , help="if set, cluster x-rows", action='store_true', default=False)
    parser.add_argument('-yc', "--yclustering"   , help="if set, cluster y-rows", action='store_true', default=False)
    parser.add_argument('-rq', "--rnaseq"        , help="if set, draw heatmap for rnaseq data", action='store_true', default=False)
    parser.add_argument('-sm', "--simpleDE"      , help="if set, draw heatmap for simple DE analysis", dest="simpleDE", action='store_true', default=False)

    # Print the help message only if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)

    # Save the STDOUT output in a log file
    if parser.parse_args().output_file:
        create_dir(get_file_info(parser.parse_args().output_file)[0])
        logfile = get_file_info(parser.parse_args().output_file)[3] + ".log"
    else:
        cwd = os.getcwd()
        logfile = "{0}/{1}.log".format(cwd,get_file_info(sys.argv[0])[1])
    logf = open(logfile, 'w')
    sys.stdout = Log(logf, sys.stdout)

    # Parse command line with parse_args and store it in an object
    args = parser.parse_args()
    print_initial_arguments(parser)
    return args

# main function
if __name__=="__main__":
      main()

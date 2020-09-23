# pwd
cd /home/rad/users/gaurav/projects/workflows/nfatacseq
# Perform tissue specific analysis for Colon and Duo

# mkdir -p tissueSpecific/{Colon,Duo}/input
# cd tissueSpecific/Colon/input/
# ln -s /media/rad/HDD1/atacseq/miguel/fignColCancer/results/bwa/mergedLibrary/MT1033-Colon-EtOH_R1.mLb.clN.sorted.bam MT1033-Colon-EtOH_R1.mLb.clN.sorted.bam
# ln -s /media/rad/HDD1/atacseq/miguel/fignColCancer/results/bwa/mergedLibrary/MT1033-Colon-TAM_R1.mLb.clN.sorted.bam MT1033-Colon-TAM_R1.mLb.clN.sorted.bam
# ln -s /media/rad/HDD1/atacseq/miguel/fignColCancer/results/bwa/mergedLibrary/MT1033-Colon-EtOH_R1.mLb.clN.sorted.bam.bai MT1033-Colon-EtOH_R1.mLb.clN.sorted.bam.bai
# ln -s /media/rad/HDD1/atacseq/miguel/fignColCancer/results/bwa/mergedLibrary/MT1033-Colon-TAM_R1.mLb.clN.sorted.bam.bai MT1033-Colon-TAM_R1.mLb.clN.sorted.bam.bai
# ln -s /media/rad/HDD1/atacseq/miguel/fignColCancer/results/bwa/mergedLibrary/macs/broadPeak/MT1033-Colon-EtOH_R1.mLb.clN_peaks.broadPeak MT1033-Colon-EtOH_R1.mLb.clN_peaks.broadPeak
# ln -s /media/rad/HDD1/atacseq/miguel/fignColCancer/results/bwa/mergedLibrary/macs/broadPeak/MT1033-Colon-TAM_R1.mLb.clN_peaks.broadPeak MT1033-Colon-TAM_R1.mLb.clN_peaks.broadPeak
# cd -

# cd tissueSpecific/Duo/input/
# ln -s /media/rad/HDD1/atacseq/miguel/fignColCancer/results/bwa/mergedLibrary/MT1033-Duo-EtOH_R1.mLb.clN.sorted.bam MT1033-Duo-EtOH_R1.mLb.clN.sorted.bam
# ln -s /media/rad/HDD1/atacseq/miguel/fignColCancer/results/bwa/mergedLibrary/MT1033-Duo-TAM_R1.mLb.clN.sorted.bam MT1033-Duo-TAM_R1.mLb.clN.sorted.bam
# ln -s /media/rad/HDD1/atacseq/miguel/fignColCancer/results/bwa/mergedLibrary/MT1033-Duo-EtOH_R1.mLb.clN.sorted.bam.bai MT1033-Duo-EtOH_R1.mLb.clN.sorted.bam.bai
# ln -s /media/rad/HDD1/atacseq/miguel/fignColCancer/results/bwa/mergedLibrary/MT1033-Duo-TAM_R1.mLb.clN.sorted.bam.bai MT1033-Duo-TAM_R1.mLb.clN.sorted.bam.bai
# ln -s /media/rad/HDD1/atacseq/miguel/fignColCancer/results/bwa/mergedLibrary/macs/broadPeak/MT1033-Duo-EtOH_R1.mLb.clN_peaks.broadPeak MT1033-Duo-EtOH_R1.mLb.clN_peaks.broadPeak
# ln -s /media/rad/HDD1/atacseq/miguel/fignColCancer/results/bwa/mergedLibrary/macs/broadPeak/MT1033-Duo-TAM_R1.mLb.clN_peaks.broadPeak MT1033-Duo-TAM_R1.mLb.clN_peaks.broadPeak
# cd -


# Parameters for the script
species="mouse"
user="miguel"
projName="fignColCancer"
outdir="/media/rad/HDD1/atacseq"
jobdir="/home/rad/users/gaurav/projects/workflows/nfatacseq"
projDir="${outdir}/${user}/${projName}"

# 1) Generate tissue specific consensus peaks region file and raw count matrix
for tissue in Colon Duo;
do
  # Generate consensus peaks region file and raw count matrix
  peaksDir="/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/${tissue}/input"
  bamDir="/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/${tissue}/input"
  consensusDir="/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/${tissue}/consensus"
  bash /home/rad/users/gaurav/projects/workflows/nfchipseq/scripts/generate_rawCount_mergedPeak_files.sh ${projDir} broadPeak ${peaksDir} ${consensusDir} ${bamDir}

  # Parse concensus raw matrix and boolean matrix to get annoation files
  bname="${projName}_consensus_peaks"
  rawCountsTxtFile="${consensusDir}/${bname}_rawCounts.txt"
  peaksAnnTxtFile="${consensusDir}/${bname}_annotation.txt"
  origAnnFile="${consensusDir}/interimFiles/${projName}_consensus_peaks.mLb.clN.boolean.txt"
  echo "bash scripts/parse_nfatacseq_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${consensusDir} ${bamDir} ${origAnnFile} ${jobdir}"
  bash scripts/parse_nfatacseq_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${consensusDir} ${bamDir} ${origAnnFile} ${jobdir}
  echo "- Consensus bed file: ${consensusPeaksBed}"
  echo "- Raw peaks count   : ${rawCountsTxtFile}"
  echo "- Peaks annotation  : ${peaksAnnTxtFile}"
done

# 3.3) Normalize raw matrix with vst
R
suppressPackageStartupMessages(library("DESeq2", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("data.table", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("splitstackshape", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("dplyr", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("ggplot2", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("reshape2", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("pander", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library("Hmisc", warn.conflicts=FALSE, quietly=TRUE))
suppressPackageStartupMessages(library(RColorBrewer, warn.conflicts=FALSE, quietly=TRUE))

# inputFile     <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/Duo/consensus/fignColCancer_consensus_peaks_rawCounts.txt"
# annotFile     <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/Duo/consensus/fignColCancer_consensus_peaks_annotation.txt"
# rldoutputFile <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/Duo/consensus/fignColCancer_consensus_peaks_rldCounts.txt"

inputFile     <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/Colon/consensus/fignColCancer_consensus_peaks_rawCounts.txt"
annotFile     <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/Colon/consensus/fignColCancer_consensus_peaks_annotation.txt"
rldoutputFile <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/Colon/consensus/fignColCancer_consensus_peaks_rldCounts.txt"

# Import data from featureCounts
countdata   <- fread(inputFile, header=TRUE)
# sampleTable <- fread(annotFile, header=TRUE)

# Filter out peaks where not a single sample has more than 50 reads.
countdata = countdata[apply(countdata[,5:ncol(countdata)], 1, max) > 50,]

# Subset genomic ranges from the counts DT
grangesDT  <- countdata[, c('PeakChrom','PeakStart','PeakEnd','PeakID')]
# Filter out peaks where not a single sample has more than 50 reads.
grangesDT <- grangesDT[apply(countdata[,4:ncol(countdata)], 1, max) > 50,]

# Remove PeakChrom PeakStart  PeakEnd columns
countdata[, c('PeakChrom','PeakStart','PeakEnd','PeakID'):=NULL]

# Create the sampletable from the column names
targetsDT <- data.table(colnames(countdata))
SampleTable <- cSplit(targetsDT, 'V1','-', drop=F)
#             V1               V1_1    V1_2     V1_3    
#  ------------------------- -------- ------- --------- 
#   1: MT1033-Colon-EtOH_R1   MT1033   Colon   EtOH_R1  
#   2: MT1033-Colon-TAM_R1    MT1033   Colon   TAM_R1   
#   3: MT1033-Duo-EtOH_R1     MT1033   Duo     EtOH_R1  
#   4: MT1033-Duo-TAM_R1      MT1033   Duo     TAM_R1   
SampleTable$V1_3 <-gsub('_R1', '', SampleTable$V1_3)

# Rename columns
names(SampleTable)[names(SampleTable) == "V1"  ] <- "SampleID"
names(SampleTable)[names(SampleTable) == "V1_1"] <- "MouseID"
names(SampleTable)[names(SampleTable) == "V1_2"] <- "TissueID"
names(SampleTable)[names(SampleTable) == "V1_3"] <- "Treatment"
#          SampleID           MouseID  TissueID Treatment
#  ------------------------- -------- ------- --------- 
#   1: MT1033-Colon-EtOH_R1   MT1033   Colon    EtOH_R1  
#   2: MT1033-Colon-TAM_R1    MT1033   Colon    TAM_R1   
#   3: MT1033-Duo-EtOH_R1     MT1033   Duo      EtOH_R1  
#   4: MT1033-Duo-TAM_R1      MT1033   Duo      TAM_R1 

# # Drop useless columns
# SampleTable <- subset(SampleTable, select = -c(V1_1, V1_3, V1_4))

# Get the dds object
dds <- DESeqDataSetFromMatrix(colData  = SampleTable, countData=countdata, design = ~Treatment)

# Relative Log Transformation
# Transform data to log space and visualize samples
rld <- assay(rlogTransformation(dds, blind = TRUE, fitType='local'))

# Add the genomic ranges in base R
rldBindDT <- cbind(grangesDT, rld)

# Round all normalized counts to 4 decimal places
rldOutDT  <- data.table(rldBindDT %>% mutate_at(vars(-PeakChrom,-PeakStart,-PeakEnd,-PeakID), funs(round(., 4))))

# Save the counts file to output csv file
fwrite(rldOutDT, file = rldoutputFile , sep = '\t', quote = FALSE)

################################
################################
ipython
#****************************************************************************************************
import datatable as dt
from datatable import *
import scanpy as sc
import time
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.lines import Line2D
filled_markers = itertools.cycle(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X'])
fillstyles     = ('full', 'left', 'right', 'bottom', 'top', 'none')


pd.set_option('display.max_rows', 20)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

####################### USER DEFINIED FUNCTIONS ###########################
# Generate PCA with 
def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.05, point['y']+0.5, str(int(point['val'])))

def show_values_on_bars(axs, h_v="v", space=0.4):
        '''https://stackoverflow.com/questions/43214978/seaborn-barplot-displaying-values'''
        def _show_on_single_plot(ax):
                if h_v == "v":
                        for p in ax.patches:
                                _x = p.get_x() + p.get_width() / 2
                                _y = p.get_y() + p.get_height()
                                value = int(p.get_height())
                                ax.text(_x, _y, value, ha="center") 
                elif h_v == "h":
                        for p in ax.patches:
                                _x = p.get_x() + p.get_width() + float(space)
                                _y = p.get_y() + p.get_height()
                                value = int(p.get_width())
                                ax.text(_x, _y, value, ha="left")

                if isinstance(axs, np.ndarray):
                        for idx, ax in np.ndenumerate(axs):
                                _show_on_single_plot(ax)
                else:
                        _show_on_single_plot(axs)

################################################################################

# Input and output files
input_file = "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/Duo/consensus/fignColCancer_consensus_peaks_rldCounts.txt"
rankoutput_file = "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/Duo/consensus/fignColCancer_consensus_peaks_annotation.txt"
tissue="Duo"
# input_file = "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/Colon/consensus/fignColCancer_consensus_peaks_rldCounts.txt"
# annot_file = "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/tissueSpecific/Colon/consensus/fignColCancer_consensus_peaks_annotation.txt"
# tissue="Colon"
projName   = "fignColCancer"
bname      = "{0}_{1}_mLcP_rld".format(projName,tissue) # mLcP = merged_library_consensus_peaks
outputDir  = "{0}/final".format(get_file_info(input_file)[0]); create_dir(outputDir)
plotsDir   = "{0}/plots".format(outputDir); create_dir(plotsDir)
dataDir    = "{0}/data".format(outputDir); create_dir(dataDir)
sample_pc  = 0.90

# Import data into the datatable
start_time = time.time()
peaksDF    = pd.read_csv(input_file, sep="\t", index_col=3)
annotDF    = pd.read_csv(annot_file, sep="\t", index_col=3)
print("%s seconds" % (time.time() - start_time))

# Drop the genomic ranges columns from the peaksDT
peaksDF.drop(columns = ['PeakChrom','PeakStart','PeakEnd'], inplace=True)

# Get the total sum for each peak
librarySizeDF = peaksDF.sum(axis=0)

# Filter peaks sum of 20 in all the samples
origPeaksDF = peaksDF.copy()

# Filter peaks with sum of peaks less than 75 for all samples
peaksDF = peaksDF[peaksDF.sum(axis=1) >= sample_pc]

# Get the variance for each peak
peaksDF['variance'] = peaksDF.var(axis=1)

# Rank the variance (expressed as percentile rank)
peaksDF['pctRank']  = peaksDF['variance'].rank(pct=True)

# Save the variance and ranks genes
peaksDF[['variance','pctRank']].to_csv("{0}/{1}.ranks".format(dataDir,bname), index=True, index_label = 'PeakID', header=True, sep="\t", float_format='%.4g')

# Get the filtered dataframe by taking top 10% of ranked peaks
hvPeaksDF = peaksDF[peaksDF['pctRank'] >= sample_pc]

# Get the dimensions of original and top peaks df
peaksDF.shape    # peaksDF.shape
hvPeaksDF.shape # (2259, 28)

# Drop the variance pctRnak columns from the top ranked genes
peaksDF.drop  (columns=['variance','pctRank'], inplace=True)
hvPeaksDF.drop(columns=['variance','pctRank'], inplace=True)

# Save the highly variable peaks to a file
hvPeaksDF.to_csv("{0}/{1}_hvPeaks{2}pc.txt".format(dataDir,bname, np.int(sample_pc*100)), index=True, index_label = 'PeakID', header=True, sep="\t", float_format='%.4g')

######################## PLOTS ###########################################
# 2) Principal component analysis
# Standardize the Data: Since PCA yields a feature subspace that maximizes the variance along the axes, 
# it makes sense to standardize the data, especially, if it was measured on different scales. 
# x = StandardScaler().fit_transform(x)
x       = hvPeaksDF.T
pca     = PCA(n_components=2)
pcs     = pca.fit_transform(x)
pcaDF   = pd.DataFrame(data = pcs, columns = ['PC1', 'PC2'])

# Add index to pcaDF
pcaDF['pcaIndex'] = hvPeaksDF.columns.tolist()

# # Add groups information
# # https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
# pcaDF.set_index('pcaIndex', inplace=True)
# pcaDF = pd.concat([pcaDF,librarySizeDF], axis=1)
# pcaDF.rename(columns={0:'LibrarySize'}, inplace=True)
# pcaDF = pd.concat([pcaDF,annotDF], axis=1)
# pcaDF['sno'] = np.arange(len(pcaDF)) + 1
# n = pcaDF['sno'].tolist()


# 2.1) Perform a Scree Plot of the Principal Components
# A scree plot is like a bar chart showing the size of each of the principal components. 
# It helps us to visualize the percentage of variation captured by each of the principal components
percent_variance = np.round(pca.explained_variance_ratio_* 100, decimals =2)
columns = ['PC1', 'PC2']
screeDF = pd.DataFrame(percent_variance, index=['PC1', 'PC2'])
screeDF.reset_index(inplace=True)
screeDF.rename(columns={ screeDF.columns[1]: "variance" }, inplace = True)
sns_t = sns.barplot(x='variance',y='index', data=screeDF, palette="Blues_d", orient='h')
show_values_on_bars(sns_t, "h", 0.3)
plt.xlabel('Percentate of Variance Explained')
plt.ylabel('Principal Components')
plt.title('PCA Scree Plot')
# Hide the right and top spines
sns.despine(top=True, right=True, left=False, bottom=False)
# or
ax=plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
screePlotPdf = "{0}/01_{1}_PCA_Scree_plot.pdf".format(plotsDir, bname)
plt.savefig(screePlotPdf,bbox_inches = 'tight')
plt.close('all')

# 2.2) Visualize 2D Projection
filled_markers = itertools.cycle(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X'])
fig = plt.figure()
ax = fig.add_subplot(1,1,1) 
# pHue = 'Treatment'; pSize = 20; pStyle = 'Treatment'
# g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue=pHue, markers=filled_markers)
g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF)
box = g.axes.get_position() # get position of figure
g.axes.set_position([box.x0, box.y0, box.width , box.height* 0.85]) # resize position
g.axes.legend(loc='center right', bbox_to_anchor=(1.10,1.10), ncol=4, prop={'size': 6})# Put a legend at the top
ax.set_xlabel('PC1 ({0:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize = 15)
ax.set_ylabel('PC2 ({0:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize = 15)
ax.set_title('VST normalized peaks PCA', fontsize = 20)
pcaPlotPdf = "{0}/02_{1}_PCA.pdf".format(plotsDir, bname)
plt.savefig(pcaPlotPdf,bbox_inches = 'tight')
plt.close('all')


# 2.5)Generate the sample correlation heatmap from top 10% varied peaks
h = sns.clustermap(peaksDF.corr(),figsize=(10, 10),cmap='gist_heat_r')
# plt.show()
heatmapPlotPdf = "{0}/05_{1}_sampleCorrelation_heatmap.pdf".format(plotsDir, bname)
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# 2.6) Generate the heatmap of top 5% varied peaks
sns.set(font_scale=0.5)
h = sns.clustermap(hvPeaksDF,z_score=0,cmap=sns.diverging_palette(220, 20, n=7),figsize=(10, 20)); 
h.ax_heatmap.set_xticklabels(h.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
# plt.show()
heatmapPlotPdf = "{0}/06_{1}_HighlyVariablePeaks{2}pc_heatmap.pdf".format(plotsDir, bname, np.int(sample_pc*100))
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

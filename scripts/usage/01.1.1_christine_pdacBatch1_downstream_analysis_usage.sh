# pwd
cd /home/rad/users/gaurav/projects/workflows/nfatacseq

# Parameters for the script
species="mouse"
user="christine"
projName="AGRad_ATACseq_MUC001"
outdir="/media/rad/HDD1/atacseq"
jobdir="/home/rad/users/gaurav/projects/workflows/nfatacseq"
projDir="${outdir}/${user}/${projName}"

#########################################################################################
# 2) PERFORM CUSTOM PEAKCALLING
#########################################################################################
# 2.1) Perform cutoff analysis to get the threshold for the peaks
# Output paramerts for downstream analysis
projDir="${outdir}/${user}/${projName}"
bamDir=${projDir}/results/bwa/mergedLibrary
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
qvalueDir="0_1" # qvalue = 0.1
qvalue="${qvalueDir//_/.}" 
customPeaksDir="${analysisDir}/customPeaks/mergedLibrary"; mkdir -p ${customPeaksDir}
peaksOutDir=${customPeaksDir}/q${qvalueDir}

# Custom peakcalling with calculated threshold using the --cutoff-analysis parameter
# The output columns are: 
  # 1) p-score -log10(pvalue) cutoff, 
  # 2) q-score -log10(qvalue) cutoff, 
  # 3) number of peaks called, 
  # 4) total length in basepairs of peaks called, 
  # 5) average length of peak in basepair
# --nolambda in the absence of a control is only valid for broad peaks, so diffuse histone modification ChIPs. If you apply it on sharp peak ChIPs like transcription factors, you deactivate the local background estimation that aims to decide if a peak is true or false based on the signal landscape surrounding the putative peak.

# for b in ${bamDir}/*.bam;
# do
#   bname=$(basename ${b} .mLb.clN.sorted.bam)
#   /home/rad/miniconda/bin/macs2 callpeak -f BAM --seed=39751 -g mm  --keep-dup all  --cutoff-analysis --nolambda --broad --broad-cutoff  ${qvalue} --outdir ${peaksOutDir} -t ${b} -n ${bname}
# done

# From the cutoff analysis, we identified a cutoff of qvalue=0.01
qvalueDir="0_01" ; # qvalue=0.01 or qscore ~ 1.96
qvalue="${qvalueDir//_/.}" 
customPeaksDir="${analysisDir}/customPeaks/mergedLibrary"; mkdir -p ${customPeaksDir}
peaksOutDir=${customPeaksDir}/q${qvalueDir}
for b in ${bamDir}/*.bam;
do
  bname=$(basename ${b} .mLb.clN.sorted.bam)
  macs2 callpeak -f BAM --seed=39751 -g mm  --keep-dup all  --nolambda --broad -q ${qvalue} --broad-cutoff  ${qvalue} --outdir ${peaksOutDir} -t ${b} -n ${bname}
done

# 2.2) Generate consensus peaks region file and raw count matrix
consensusDir="${peaksOutDir}/consensus"
echo -e "- ${projDir}\n- ${peaksOutDir}\n- ${consensusDir}"
bash /home/rad/users/gaurav/projects/workflows/nfchipseq/scripts/generate_rawCount_mergedPeak_files.sh ${projDir} broadPeak ${peaksOutDir} ${consensusDir}

# 2.3) 
#########################################################################################
# 3) DOWNSTREAM ANALYSIS
#########################################################################################
# 3.1) Output paramerts for downstream analysis
projDir="${outdir}/${user}/${projName}"
bamDir=${projDir}/results/bwa/mergedLibrary
analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
qvalueDir="0_01" ; # qvalue=0.01 or qscore ~ 1.96
qvalue="${qvalueDir//_/.}" 
customPeaksDir="${analysisDir}/customPeaks/mergedLibrary"; mkdir -p ${customPeaksDir}
peaksOutDir=${customPeaksDir}/q${qvalueDir}
consensusDir="${peaksOutDir}/consensus"

bname="${projName}_consensus_peaks"
rawCountsTxtFile="${consensusDir}/${bname}_rawCounts.txt"
peaksAnnTxtFile="${consensusDir}/${bname}_annotation.txt"
origAnnFile="${consensusDir}/interimFiles/${projName}_consensus_peaks.mLb.clN.boolean.txt"

# 3.2) Parse concensus raw matrix and boolean matrix to get annoation files
echo "bash scripts/parse_nfatacseq_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${consensusDir} ${bamDir} ${origAnnFile} ${jobdir}"
bash scripts/parse_nfatacseq_consensus_peaks_annotation.sh ${species} ${user} ${projName} ${consensusDir} ${bamDir} ${origAnnFile} ${jobdir}

# Output files are:
echo "- Consensus bed file: ${consensusPeaksBed}"
echo "- Raw peaks count   : ${rawCountsTxtFile}"
echo "- Peaks annotation  : ${peaksAnnTxtFile}"


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

inputFile     <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/customPeaks/mergedLibrary/q0_01/consensus/fignColCancer_consensus_peaks_rawCounts.txt"
annotFile     <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/customPeaks/mergedLibrary/q0_01/consensus/fignColCancer_consensus_peaks_annotation.txt"
rldoutputFile <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/customPeaks/mergedLibrary/q0_01/consensus/fignColCancer_consensus_peaks_rldCounts.txt"
vstoutputFile <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/customPeaks/mergedLibrary/q0_01/consensus/fignColCancer_consensus_peaks_vstCounts.txt"
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
# Save the counts file to output csv file
write.table(rld   , file = rldoutputFile , sep = '\t', quote = F)

# VST 
vst      <- varianceStabilizingTransformation(dds, blind=TRUE)  
vsd      <- assay(vst)
# Add the genomic ranges in base R
vstBindDT <- cbind(grangesDT, vsd)
# Round all normalized counts to 4 decimal places
vstOutDT  <- data.table(vstBindDT %>% mutate_at(vars(-PeakChrom,-PeakStart,-PeakEnd,-PeakID), funs(round(., 4))))
# Save the counts file to output csv file
fwrite(vstOutDT   , file = vstoutputFile , sep = '\t', quote = F)

# PCA
pca = prcomp(t(rld))
print(summary(pca))
pcaData = as.data.frame(pca$x)
pcaData$sample=rownames(pcaData)
si <- as.data.frame(SampleTable)
rownames(si) <- si$SampleID; si$SampleID <- NULL
pcaData=merge(pcaData, si, by.x=0, by.y=0)
percentVar = round(100 * (pca$sdev^2 / sum( pca$sdev^2 ) ))
p=ggplot(data=pcaData, aes(x = PC1, y = PC2, group=Treatment)) + geom_point(aes(size=3, shape=TissueID, color=Treatment))
p=p+xlab(paste0("PC1: ", percentVar[1], "% variance"))
p=p+ylab(paste0("PC2: ", percentVar[2], "% variance"))
normPlotPDF  <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/analysis_q0_01/fignColCancer_rld_PCA.pdf"
ggsave(normPlotPDF, width=15)

# Data visualization of normalized counts
lf = melt(rld, id.vars=c())
pander(head(lf))
pg1 <- ggplot(data=lf, aes(x=Var2, y=value)) + geom_boxplot(aes(group=Var2)) + xlab("SampleID") + ylab("Normalized Count") + coord_flip()
normPlotPDF  <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/analysis_q0_01/fignColCancer_rld_boxplot.pdf"
ggsave(normPlotPDF)
pg2 <- ggplot(data=lf, aes(x=value, after_stat(density))) + geom_freqpoly(aes(group=Var2, color=Var2), bins=30) + xlab("Normalized Count") + ylab("density")
normPlotPDF  <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/analysis_q0_01/fignColCancer_rld_normalization.pdf"
ggsave(normPlotPDF, width=15)

# Sample cluster dendrogram
hclustPlotPDF  <- "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/analysis_q0_01/fignColCancer_rld_hclust_dendrogram.pdf"
dists          <- dist(t(rld))
pdf(hclustPlotPDF); par(mar=c(20,4,1,1));
plot(hclust(dists))
dev.off()

# Heatmap of sample-to-sample distances using the Poisson Distance
suppressMessages(library("PoiClaClu"))
suppressMessages(library("pheatmap"))
poisd <- PoissonDistance(t(vsd))

# In order to plot the sample distance matrix with the rows/columns arranged by the distances in our distance matrix, we manually provide sampleDists 
# to the clustering_distance argument of the pheatmap function. Otherwise the pheatmap function would assume that the matrix contains the data values 
# themselves, and would calculate distances between the rows/columns of the distance matrix, which is not desired.
samplePoisDistMatrix <- as.matrix( poisd$dd)

# Change the row names of the distance matrix to contain Treatment and Control instead of sample ID, so that we have all this information in view when looking at the heatmap.
rownames(samplePoisDistMatrix) <- paste( attr(vsd,"colData")$condition, rownames(attr(vsd,"colData")),sep="-")
colnames(samplePoisDistMatrix) <- NULL
# colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
colors <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
# colors <- colorRampPalette(rev(c('gold','darkorange','darkred')))(256)
# colors <- colorRampPalette(c('darkblue','blue','lightblue','white','orange','red','darkred'))(1024)
pheatmap(samplePoisDistMatrix,
        clustering_distance_rows=poisd$dd,
        clustering_distance_cols=poisd$dd,
        col=colors,
        fontsize=8)

#par(new=TRUE)
# Turn off device driver (to flush output to PNG/PDF file)
dev.off()

################################
# Get the top 1% most varying genes and plot the pca
ipython
#****************************************************************************************************
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

pd.set_option('display.max_rows', 5)
pd.set_option('display.max_columns', 8)
pd.set_option('display.width', 1000)

# Input and output files
input_file = "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/customPeaks/mergedLibrary/q0_01/consensus/fignColCancer_consensus_peaks_rldCounts.txt"
rankoutput_file = "/media/rad/HDD1/atacseq/miguel/fignColCancer/analysis/customPeaks/mergedLibrary/q0_01/consensus/fignColCancer_consensus_peaks_rldCounts.ranks"

# Import data into the dataframe
peaksDF    = pd.read_csv(input_file, sep="\t", index_col=0)

# Get the total sum for each peak
librarySizeDF = peaksDF.sum(axis=0)

# Get the variance for each peak
peaksDF['variance'] = peaksDF.var(axis=1)

# Rank the variance (expressed as percentile rank)
peaksDF['pctRank']  = peaksDF['variance'].rank(pct=True)

# Save the variance and ranks genes
peaksDF[['variance','pctRank']].to_csv(rankoutput_file, index=True, header=True, sep="\t", float_format='%.2f')

# Get the filtered dataframe by taking top 10% of ranked peaks
topPeaksDF = peaksDF[peaksDF['pctRank'] >= 0.90]

# Get the dimensions of original and top peaks df
peaksDF.shape    # peaksDF.shape
topPeaksDF.shape # (1600, 20)

# Drop the variance pctRnak columns from the top ranked genes
topPeaksDF.drop(columns=['variance','pctRank'], inplace=True)

# Principal component analysis
features = topPeaksDF.columns.tolist()
x        = topPeaksDF.T
# Standardize the Data: Since PCA yields a feature subspace that maximizes the variance along the axes, 
# it makes sense to standardize the data, especially, if it was measured on different scales. 
x = StandardScaler().fit_transform(x)

# PCA Projection to 2D
pca   = PCA(n_components=3)
pcs   = pca.fit_transform(x)
pcaDF = pd.DataFrame(data = pcs, columns = ['PC1', 'PC2','PC3'])


# Add the cellLines information
pcaDF = pd.concat([pcaDF, pd.DataFrame(data=features, columns=['SampleID'])], axis = 1)

# Add groups information
# https://www.geeksforgeeks.org/python-pandas-split-strings-into-two-list-columns-using-str-split/
pcaDF['TissueID']  = pcaDF["SampleID"].str.split("-", n = 0, expand = True)[1] # Colon, Duo
pcaDF['Treatment'] = pcaDF["SampleID"].str.split("-", n = 0, expand = True)[2].str.replace("_R1","") # EtOH, TAM

# Perform a Scree Plot of the Principal Components
# A scree plot is like a bar chart showing the size of each of the principal components. 
# It helps us to visualize the percentage of variation captured by each of the principal components
percent_variance = np.round(pca.explained_variance_ratio_* 100, decimals =2)
columns = ['PC1', 'PC2', 'PC3']
screeDF = pd.DataFrame(percent_variance, index=['PC1', 'PC2', 'PC3'])
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
# plt.show()
screePlotPdf = "{0}_PCA_Scree_plot.pdf".format(get_file_info(input_file)[3])
# screePlotPdf = "{0}_PCA_Scree_plot_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(screePlotPdf,bbox_inches = 'tight')
plt.close('all')

# Visualize 2D Projection
fig = plt.figure()
ax = fig.add_subplot(1,1,1) 
g = sns.scatterplot(x='PC1', y='PC2', data=pcaDF, hue='Treatment', style='TissueID', size='Treatment')
box = g.axes.get_position() # get position of figure
g.axes.set_position([box.x0, box.y0, box.width , box.height* 0.85]) # resize position
g.axes.legend(loc='center right', bbox_to_anchor=(1.10,1.10), ncol=4, prop={'size': 6})# Put a legend at the top
ax.set_xlabel('PC1 ({0:.2f}%)'.format(pca.explained_variance_ratio_[0]*100), fontsize = 15)
ax.set_ylabel('PC2 ({0:.2f}%)'.format(pca.explained_variance_ratio_[1]*100), fontsize = 15)
ax.set_title('Top 1% variance ranked peaks PCA', fontsize = 20)
pcaPlotPdf = "{0}_PCA_plot.pdf".format(get_file_info(input_file)[3])
# pcaPlotPdf = "{0}_PCA_plot_ohne_004_samples.pdf".format(get_file_info(input_file)[3])
plt.savefig(pcaPlotPdf,bbox_inches = 'tight')
# plt.show()
plt.close('all')

# Generate the heatmap of top 1% varied peaks
sns.set(font_scale=0.5)
h = sns.clustermap(topPeaksDF,z_score=0,cmap=sns.diverging_palette(220, 20, n=7),figsize=(10, 20)); 
h.ax_heatmap.set_xticklabels(h.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
# plt.show()
heatmapPlotPdf = "{0}_top10pc_heatmap.pdf".format(get_file_info(input_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')

# Generate the sample correlation heatmap from top 1% varied peaks
# sns.set(font_scale=0.5)
h = sns.clustermap(topPeaksDF.corr(),figsize=(10, 10),cmap='gist_heat_r', vmax=1.1, vmin=-0.1, annot=True)
# plt.show()
heatmapPlotPdf = "{0}_sampleCorrelation_top10pc_heatmap.pdf".format(get_file_info(input_file)[3])
plt.savefig(heatmapPlotPdf,bbox_inches = 'tight')
plt.close('all')




# Generate PCA ohne 004_* samples
# Remove all columns containing the string 004_
originalPeaksDF=peaksDF.copy()
peaksDF = peaksDF.loc[:,~peaksDF.columns.str.contains('004_', case=False)] 
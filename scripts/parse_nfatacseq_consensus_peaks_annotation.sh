#!/bin/bash

# USAGE: bash scripts/parse_nfchip_consensus_peaks_annotation.sh mouse christine pdacBatch1 /media/rad/HDD1/nfchip /media/rad/HDD1/nfchip/christine/pdacBatch1/results/bwa/mergedLibrary /media/rad/HDD1/nfchip/christine/pdacBatch1/results/bwa/mergedLibrary/macs/narrowPeak/consensus/abcam/abcam.consensus_peaks.boolean.txt /home/rad/users/gaurav/projects/nfPipelines/nfchipseq

species=${1:-"mouse"}
user=${2:-"anja"}
projName=${3:-"rep_tALLcellLineMm"}
outdir=${4:-"/media/rad/HDD1/atacseq"}
bamdir=${5:-""}
origAnnFile=${6:-""}  # /media/rad/HDD1/atacseq/anja/nfTALLMm/results/bwa/mergedLibrary/macs/broadPeak/consensus/consensus_peaks.mLb.clN.boolean.txt
jobdir=${7:-"/home/rad/users/gaurav/projects/nfPipelines/nfchipseq"}

# Get relevant directories
bname="${projName}_consensus_peaks"

# Input argument and relevant directories
echo "######################################################################"
echo "# Input argument and relevant directories "
echo ""
echo "species=${species}"
echo "user=${user}"
echo "projName=${projName}"
echo "outdir=${outdir}"
echo "bamdir=${bamdir}"
echo "origAnnFile=${origAnnFile}"
echo "jobdir=${jobdir}"
echo ""
echo "######################################################################"
echo ""

# Re-index bam files 
echo "- Re-index bam files "
ls ${bamdir}/*.bam | parallel --progress --eta -j 16 "samtools index {}"

# Generate consensus file with PeakID to analysis folder
echo "- Generate consensus file with PeakID to analysis folder"
consensusPeaksBed="${outdir}/${bname}.bed"
cut -f1-4 ${origAnnFile} | egrep -v "start|end" > ${consensusPeaksBed}

# Get the tab seaprated raw counts matrix for consensus peaks on all samples using deeptools
echo "- Getting the tab seaprated raw counts matrix for consensus peaks on all samples"
echo "    - This may take some time as it is using multiBamSummary from deeptools"
rawCountsBname="${outdir}/${bname}_rawCounts"
rawCountsTxtFile="${rawCountsBname}.txt"
multiBamSummary BED-file --BED ${consensusPeaksBed} --bamfiles ${bamdir}/*.bam --smartLabels -out ${rawCountsBname}.npz --outRawCounts ${rawCountsBname}.txt -p 64

# Sort the input file with the header line intact
echo "- Sorting the input file with the header line intact"
(head -n 1 ${rawCountsTxtFile} && tail -n +2 ${rawCountsTxtFile} | sort -k1,1V -k2,2g -k3,3g) > ${rawCountsTxtFile}.tmp && mv ${rawCountsTxtFile}.tmp ${rawCountsTxtFile}

# Remove single quotes from the header
echo "- Remove single quotes from the header"
sed -i "s/[']//g" ${rawCountsTxtFile}

# Rename columns. \b is word boundary to search for full word only
echo "- Rename columns and add chr to the chromosome column"
sed -i "s/#chr\b/PeakChrom/" ${rawCountsTxtFile}
sed -i "s/start\b/PeakStart/" ${rawCountsTxtFile}
sed -i "s/end\b/PeakEnd/" ${rawCountsTxtFile}
sed -i 's/^/chr/' ${rawCountsTxtFile}
sed -i "s/chrPeakChrom\b/PeakChrom/" ${rawCountsTxtFile}
sed -i 's/^/chr/' ${consensusPeaksBed}

# Change float to integer in the python one liner
python -c "import sys; import pandas as pd; import datatable as dt; input_file=sys.argv[1]; input_bed=sys.argv[2]; outtxt_file=sys.argv[3]; peaksDT = dt.fread(input_file, sep='\t', header=True, nthreads=16); bedDT = dt.fread(input_bed, sep='\t', header=False, nthreads=16); peaksDF = peaksDT.to_pandas(); bedDF = bedDT.to_pandas(); peaksDF.insert (3, 'PeakID', peaksDF['PeakChrom'].str.cat(peaksDF['PeakStart'].apply(str), sep='_').str.cat(peaksDF['PeakEnd'].apply(str), sep='_')); bedDF.insert (3, 'PeakID', bedDF['C0'].str.cat(bedDF['C1'].apply(str), sep='_').str.cat(bedDF['C2'].apply(str), sep='_')); peaksDF.set_index('PeakID', inplace=True); bedDF.set_index('PeakID', inplace=True); bedDF.drop(columns=['C0','C1','C2'], inplace=True); peaksDF = pd.merge(peaksDF, bedDF, left_index=True, right_index=True, how='outer'); col = peaksDF.pop('C3'); peaksDF.insert(3, col.name, col); peaksDF.rename(columns={ peaksDF.columns[3]: 'PeakID' }, inplace = True); peaksDF.to_csv(outtxt_file, header=True, index=False, sep='\t', float_format='%.0f'); " ${rawCountsTxtFile} ${consensusPeaksBed} ${rawCountsTxtFile}.tmp

mv ${rawCountsTxtFile}.tmp ${rawCountsTxtFile}

# Remove additional characters from the header
sed -i "s/.mLb.clN.bam//g" ${rawCountsTxtFile}
sed -i "s/.mLb.clN.sorted.bam//g" ${rawCountsTxtFile}
sed -i "s/.mLb.clN.sorted//g" ${rawCountsTxtFile}

# Get parsed boolean annotation file
peaksAnnTabFile="${outdir}/${bname}_annotation.tab"
# Rscript -e 'library(data.table); inputfile  <- commandArgs(TRUE)[1]; outputfile <- commandArgs(TRUE)[2]; peaksDT    <- fread(inputfile, header=TRUE, sep="\t"); dropCols   <- grep("(fc|qval|pval)$", colnames(peaksDT)); peaksDT[,(dropCols) :=NULL]; to.replace <- names(which(sapply(peaksDT, is.logical))); for (var in to.replace) peaksDT[, (var):= as.numeric(get(var))]; peaksDT[,PeakID:=paste(interval_id,chr,start,end, sep="_")]; peaksDT[,"interval_id":=NULL]; setnames(peaksDT, old = c("chr","start", "end", "num_peaks", "num_samples"), new = c("PeakChrom", "PeakStart", "PeakEnd","NumPeaks","NumSamples")); newColOrder <-  c(colnames(peaksDT)[0:3],"PeakID", colnames(peaksDT)[5:length(colnames(peaksDT))-1]); setcolorder(peaksDT, newColOrder);names(peaksDT) = gsub(pattern = ".mLb.clN.", replacement = "_", x = names(peaksDT));peaksDT[,PeakChrom := paste0("chr",PeakChrom)];fwrite(peaksDT, outputfile, sep = "\t");' ${origAnnFile} ${peaksAnnTabFile}
Rscript -e 'library(data.table); inputfile  <- commandArgs(TRUE)[1]; outputfile <- commandArgs(TRUE)[2]; peaksDT    <- fread(inputfile, header=TRUE, sep="\t"); dropCols   <- grep("(fc|qval|pval)$", colnames(peaksDT)); peaksDT[,(dropCols) :=NULL]; to.replace <- names(which(sapply(peaksDT, is.logical))); for (var in to.replace) peaksDT[, (var):= as.numeric(get(var))]; peaksDT[,PeakID:=paste(interval_id)]; peaksDT[,"interval_id":=NULL]; setnames(peaksDT, old = c("chr","start", "end", "num_peaks", "num_samples"), new = c("PeakChrom", "PeakStart", "PeakEnd","NumPeaks","NumSamples")); newColOrder <-  c(colnames(peaksDT)[0:3],"PeakID", colnames(peaksDT)[5:length(colnames(peaksDT))-1]); setcolorder(peaksDT, newColOrder);names(peaksDT) = gsub(pattern = ".mLb.clN.", replacement = "_", x = names(peaksDT));peaksDT[,PeakChrom := paste0("chr",PeakChrom)];fwrite(peaksDT, outputfile, sep = "\t");' ${origAnnFile} ${peaksAnnTabFile}

# Add genomic annotation using chipseeker
peaksAnnTxtFile="${outdir}/${bname}_annotation.txt"

# Create a log file
analysisLogDir="${outdir}/logs/analysisLogs"; mkdir -p ${analysisLogDir}
mkdir -p ${analysisLogDir}
annotationLogFile="${analysisLogDir}/${bname}_annotation.log"

# Annotate consensus peaks
echo "species  = ${species}"  2>&1 | tee    ${annotationLogFile}
echo "user     = ${user}"     2>&1 | tee -a ${annotationLogFile}
echo "projName = ${projName}" 2>&1 | tee -a ${annotationLogFile}
echo "outdir   = ${outdir}"   2>&1 | tee -a ${annotationLogFile}
echo "jobdir   = ${jobdir}"   2>&1 | tee -a ${annotationLogFile}

echo "-------------------------------------------------------------------"
echo "Rscript ${jobdir}/scripts/R_annotate_peaks_chipseeker.R -if=${peaksAnnTabFile} -of=${peaksAnnTxtFile} -sp=${species} 2>&1 | tee -a ${annotationLogFile}"
echo "-------------------------------------------------------------------"
Rscript ${jobdir}/scripts/R_annotate_peaks_chipseeker.R -if=${peaksAnnTabFile} -of=${peaksAnnTxtFile} -sp=${species} 2>&1 | tee -a ${annotationLogFile}

# Remove the peaksAnnTabFile file as the data is there in the peaksAnnTxtFile
rm ${peaksAnnTabFile}

# Output files are:
echo "######################################################################" 2>&1 | tee    ${annotationLogFile}
echo "" 2>&1 | tee    ${annotationLogFile}
echo "- Consensus bed file: ${consensusPeaksBed}" 2>&1 | tee    ${annotationLogFile}
echo "- Raw peaks count   : ${rawCountsTxtFile}"  2>&1 | tee    ${annotationLogFile}
echo "- Peaks annotation  : ${peaksAnnTxtFile}"   2>&1 | tee    ${annotationLogFile}
echo "" 2>&1 | tee    ${annotationLogFile}
echo "######################################################################" 2>&1 | tee    ${annotationLogFile}

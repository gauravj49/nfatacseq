cd /home/rad/users/gaurav/projects/workflows/nfatacseq

# Parameters for the script
species="mouse"
user="helena"
projName="cnsTreg"
outdir="/media/rad/HDD1/atacseq"
jobdir="/home/rad/users/gaurav/projects/workflows/nfatacseq"

#########################################################################################
# 1) RUN THE PIPELINE
#########################################################################################

# Parameters to run the pipeline
projDir="${outdir}/${user}/${projName}"
fastqDir="${projDir}/fastq"; mkdir -p ${fastqDir}

# # 1.1) Copy the fastq files
# Provided by Helena on a external hard drive

# 1.2) Run the pipeline
cd ${projDir}
nextflow run /home/rad/users/gaurav/projects/nfPipelines/nfatacseq --input ${projDir}/nfCnsTreg_design.csv --genome GRCm38 --single_end -name ${projName}1

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

for b in ${bamDir}/*.bam;
do
  bname=$(basename ${b} .mLb.clN.sorted.bam)
  /home/rad/miniconda/bin/macs2 callpeak -f BAM --seed=39751 -g mm  --keep-dup all  --cutoff-analysis --nolambda --broad --broad-cutoff  ${qvalue} --outdir ${peaksOutDir} -t ${b} -n ${bname}
done


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




















# 3.3) Custom peakcalling with calculated threshold using the --cutoff-analysis parameter
# The output columns are: 
  # 1) p-score (log10(-p)) cutoff, 
  # 2) q-score (log10(-q)) cutoff, 
  # 3) number of peaks called, 
  # 4) total length in basepairs of peaks called, 
  # 5) average length of peak in basepair

customPeaksDir="${analysisDir}/customPeaks"; mkdir -p ${customPeaksDir}

# Cat the following commands in a temporary file and run using gnu parallel
# # modify manually for H3k4me1 marks to run peaks as broad peaks
for pd in 0_0001 0_00005 0_000001
do
  pvalue="${pd//_/.}"
  echo ${pvalue}
  peaksOutDir=${customPeaksDir}/p${pd}
  macs2 callpeak -f BAM --seed=39751 -g mm  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t  /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/results/bwa/mergedLibrary/foxp1ovko_200519062115_A0000504_FKO372Foxp1_mmu_chipseq_se_R1.mLb.clN.sorted.bam -c /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/results/bwa/mergedLibrary/foxp1ovko_200519062115_A0000504_FKO372Input_mmu_chipseq_se_R1.mLb.clN.sorted.bam -n foxp1ovko_200519062115_A0000504_FKO372Foxp1_mmu_chipseq_se_R1_${pd} --outdir ${peaksOutDir}
  macs2 callpeak -f BAM --seed=39751 -g mm  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t  /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/results/bwa/mergedLibrary/foxp1ovko_200519062115_A0000504_FKO387Foxp1_mmu_chipseq_se_R1.mLb.clN.sorted.bam -c /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/results/bwa/mergedLibrary/foxp1ovko_200519062115_A0000504_FKO387Input_mmu_chipseq_se_R1.mLb.clN.sorted.bam -n foxp1ovko_200519062115_A0000504_FKO387Foxp1_mmu_chipseq_se_R1_${pd} --outdir ${peaksOutDir}
  macs2 callpeak -f BAM --seed=39751 -g mm  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t  /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/results/bwa/mergedLibrary/foxp1ovko_200519062115_A0000504_PKF1OE057Foxp1_mmu_chipseq_se_R1.mLb.clN.sorted.bam -c /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/results/bwa/mergedLibrary/foxp1ovko_200519062115_A0000504_PKF1OE057Input_mmu_chipseq_se_R1.mLb.clN.sorted.bam -n foxp1ovko_200519062115_A0000504_PKF1OE057Foxp1_mmu_chipseq_se_R1_${pd} --outdir ${peaksOutDir}
  macs2 callpeak -f BAM --seed=39751 -g mm  -p ${pvalue} --keep-dup auto  --cutoff-analysis -t  /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/results/bwa/mergedLibrary/foxp1ovko_200519062115_A0000504_PKF1OE062Foxp1_mmu_chipseq_se_R1.mLb.clN.sorted.bam -c /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe/results/bwa/mergedLibrary/foxp1ovko_200519062115_A0000504_PKF1OE062Input_mmu_chipseq_se_R1.mLb.clN.sorted.bam -n foxp1ovko_200519062115_A0000504_PKF1OE062Foxp1_mmu_chipseq_se_R1_${pd} --outdir ${peaksOutDir}
done
# cat > ${customPeaksDir}/${projName}_macs2_callpeak_p005.sh
# # Cltr+C
# chmod 775 ${customPeaksDir}/${projName}_macs2_callpeak_p005.sh
# parallel :::: ${customPeaksDir}/${projName}_macs2_callpeak_p005.sh

# Annotate using chipseeker
for pd in 0_005 0_0001 0_00005 0_000001 
do
 peaksOutDir=${customPeaksDir}/p${pd}
 for peaksInputFile in $(ls ${peaksOutDir}/*.narrowPeak)
 do
  peaksAnnTxtFile=${peaksOutDir}/$(basename ${peaksInputFile} .narrowPeak)_annotation.txt
  Rscript ${jobdir}/scripts/R_annotate_peaks_chipseeker.R -if=${peaksInputFile} -of=${peaksAnnTxtFile} -sp=${species} -mf -ac
 done
done


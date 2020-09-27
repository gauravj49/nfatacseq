#!/bin/bash
# USAGE: bash scripts/generate_rawCount_mergedPeak_files.sh /media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe narrowPeak

# Input commandline parameters
projDir=${1:-"/media/rad/HDD1/nfchip/thorsten/foxp1FkoPkf1oe"}
peaksShape=${2:-"narrowPeak"} # or broadPeak
peaksDir=${3:-"${projDir}/results/bwa/mergedLibrary/macs/${peaksShape}"}
analysisDir=${4:-"${projDir}/analysis"}; mkdir -p ${analysisDir}

# Parameters to run the pipeline
projName=$(basename ${projDir})
bamFiles=${projDir}/results/bwa/mergedLibrary/*.bam
# peaksDir=${projDir}/results/bwa/mergedLibrary/macs/${peaksShape}
# Get relevant directories
# analysisDir="${projDir}/analysis"; mkdir -p ${analysisDir}
interimFilesDir="${analysisDir}/interimFiles"; mkdir -p ${interimFilesDir}

# 1.1) Get the concensus peaks from all groups
## MergedIntervalTxtFile is file created using commands below:
## 1) broadPeak
## sort -k1,1 -k2,2n <MACS_BROADPEAK_FILES_LIST> | mergeBed -c 2,3,4,5,6,7,8,9 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > merged_peaks.txt
## 2) narrowPeak
## sort -k1,1 -k2,2n <MACS_NARROWPEAK_FILE_LIST> | mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > merged_peaks.txt
# Get list of all peak files separated by comma
sampleNameList=""
macsPeaksFileList=""
for f in $(ls ${peaksDir}/*.${peaksShape});
do
  echo ${f}
  bname=$(basename ${f} _peaks.${peaksShape})
  sampleNameList="${bname},${sampleNameList}"
  macsPeaksFileList="${f} ${macsPeaksFileList}"
done
# Remove the last comma
sampleNameList=$(echo "${sampleNameList%,*} ${sampleNameList##*,}")

# Get Merged Interval Txt File 
echo "- 1.1) Get the concensus peaks from all groups"
MergedPeaksTxt="${interimFilesDir}/${projName}_merged_peaks_all_samples.txt"
if [ ${peaksShape} = "narrowPeak" ]; then
  # For narrow peaks 
  sort -k1,1V -k2,2g ${macsPeaksFileList} | mergeBed -c 2,3,4,5,6,7,8,9,10 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > ${MergedPeaksTxt}
else
  # For broad peaks
  sort -k1,1V -k2,2g ${macsPeaksFileList} | mergeBed -c 2,3,4,5,6,7,8,9 -o collapse,collapse,collapse,collapse,collapse,collapse,collapse,collapse > ${MergedPeaksTxt}
fi
echo -e "\t- ${MergedPeaksTxt}"; echo ""

# 1.2) Get the consensus peaks annotation and boolean matrix
echo "- 1.2) Get the consensus peaks annotation and boolean matrix"
boolAnnPeaksFile="${interimFilesDir}/${projName}_consensus_peaks.mLb.clN.boolean.txt"
if [ ${peaksShape} = "narrowPeak" ]; then
  python /home/rad/users/gaurav/projects/workflows/nfchipseq/bin/macs2_merged_expand.py ${MergedPeaksTxt} ${sampleNameList} ${boolAnnPeaksFile} --is_narrow_peak
else
  python /home/rad/users/gaurav/projects/workflows/nfchipseq/bin/macs2_merged_expand.py ${MergedPeaksTxt} ${sampleNameList} ${boolAnnPeaksFile}
fi
echo -e "\t- ${boolAnnPeaksFile}"; echo ""

python -c "import sys; import pandas as pd; import datatable as dt; input_file=sys.argv[1]; outtxt_file=sys.argv[2]; peaksDT = dt.fread(input_file, sep='\t', header=True, nthreads=16); peaksDF = peaksDT.to_pandas(); peaksDF.insert (3, 'interval_id_new', peaksDF['interval_id'].str.cat(peaksDF['chr'], sep='_').str.cat(peaksDF['start'].apply(str), sep='_').str.cat(peaksDF['end'].apply(str), sep='_')); peaksDF.drop(columns=['interval_id'], inplace=True); peaksDF.rename(columns={ peaksDF.columns[3]: 'interval_id' }, inplace = True); peaksDF.to_csv(outtxt_file, header=True, index=False, sep='\t', float_format='%.0f')" ${boolAnnPeaksFile} ${boolAnnPeaksFile}.tmp
mv ${boolAnnPeaksFile}.tmp ${boolAnnPeaksFile}

# 1.3) Get the concensus peaks from all groups
echo "- 1.3) Get the concensus peaks from all groups (bed and saf format)"
# Get the output consensus bed file
#  === ============== ============== 
#   1       16228066       16228163  
#   1       44951210       44951269  
#   1       85495021       85495094  
#   1      106171472      106171616  
#  === ============== ============== 
consensusPeaksBed="${interimFilesDir}/${projName}_consensus_peaks.bed"
cut -f1-3 ${boolAnnPeaksFile}| sort -k1,1V -k2,2g | egrep -v "start|end"> ${consensusPeaksBed}

# Get the output consensus SAF (Simplified Annotation Format) file
# ============ ===== ============== ============== ========= 
#   GeneID      Chr        Start         End         Strand  
# ============ ===== ============== ============== ========= 
# Interval_1     1       16228066       16228163      +     
# Interval_2     1       44951210       44951269      +     
# Interval_3     1      106171472      106171616      +     
# Interval_4     1      156615557      156615645      +     
# ============ ===== ============== ============== ========= 
consensusPeaksSaf="${interimFilesDir}/${projName}_consensus_peaks.saf"
python -c "import sys; import pandas as pd; import datatable as dt; input_file=sys.argv[1]; outtxt_file=sys.argv[2]; 
peaksDT = dt.fread(input_file, sep='\t', header=False, nthreads=16);peaksDF = peaksDT.to_pandas(); 
peaksDF.rename(columns={ peaksDF.columns[0]: 'Chr' }, inplace = True); peaksDF.rename(columns={ peaksDF.columns[1]: 'Start' }, inplace = True); peaksDF.rename(columns={ peaksDF.columns[2]: 'End'   }, inplace = True); peaksDF.insert (0, 'GeneID', ['Interval_{0}'.format(x) for x in range(peaksDF.shape[0])]); peaksDF.insert (4, 'Strand', '+'); peaksDF.to_csv(outtxt_file, header=True, index=False, sep='\t', float_format='%.0f')" ${consensusPeaksBed} ${consensusPeaksSaf}
echo -e "\t- ${consensusPeaksBed}"; 
echo -e "\t- ${consensusPeaksSaf}"; 
echo ""

# 1.4) Get the tab seaprated raw counts of the merged peaks for all the samples
echo "- 1.4) Get the tab seaprated raw counts of the merged peaks for all the samples"
origConsFile="${interimFilesDir}/${projName}_consensus_peaks.mLb.clN.featureCounts.txt"
echo "featureCounts -F SAF -O --fracOverlap 0.2 -T 6 -p --donotsort -a ${consensusPeaksSaf} -o ${origConsFile} $(echo ${bamFiles})"
featureCounts "-F" "SAF" "-O" "--fracOverlap" "0.2" "-T" "6" "-p" "--donotsort" "-a" ${consensusPeaksSaf} -o ${origConsFile} $(echo ${bamFiles})
echo -e "\t- ${origConsFile}"; echo ""
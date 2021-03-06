# pwd
cd /home/rad/users/gaurav/projects/workflows/nfatacseq/

# NOTE: This will be depricated from 2021_01Jan_01 and the pipeline usage will be done in a separate usage file
# rm -rf work/ results/ .nextflow*

##########################################################
# 1) Anja
##########################################################
# 1.1) nfTALLmm
# ls /media/nas/temporary/PUB_CRCs/tallCRCs/mouse/fastq/*.gz | parallel --progress --eta -j 32 "rsync --ignore-existing -arzRP {} /media/rad/HDD1/atacseq/anja/nfTALLmm/fastq"
cd /media/rad/HDD1/atacseq/anja/nfTALLmm
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/anja/nfTALLmm/nfTALLmm_design.csv --genome GRCm38 --single_end --narrow_peak -resume
cd -

# 1.2) nfTALLhs
# ls /media/nas/temporary/PUB_CRCs/tallCRCs/human/fastq/*.gz | parallel --progress --eta -j 32 "rsync --ignore-existing -arzRP {} /media/rad/HDD1/atacseq/anja/nfTALLhs/fastq"
cd /media/rad/HDD1/atacseq/anja/nfTALLhs
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/anja/nfTALLhs/nfTALLhs_design.csv --genome GRCh38 --single_end --narrow_peak -resume
cd -

##########################################################
# 2) Christine
##########################################################
# 2.1) CK-ATACseq samples
cd /media/rad/HDD1/atacseq/christine/ckatac
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/christine/ckatac/nfCKatac_design.csv --genome GRCm38 --single_end -name ckatac1 
cd -

# 2.2) AGRad_ATACseq_MUC001 samples
cd /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/christine/AGRad_ATACseq_MUC001/nfAGRad_ATACseq_MUC001_design.csv --genome GRCm38 --single_end -name AGRad_ATACseq_MUC001 
cd - 

# 2.2)AGBuchholz samples
cd /media/rad/HDD1/atacseq/christine/AGBuchholz
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/christine/AGBuchholz/AGBuchholz_design.csv --genome GRCm38 -name AGBuchholz 
cd - 

# 2.3)Test PE 150 public samples from GEO GSE145705
cd /media/rad/HDD1/atacseq/christine/GSE145705
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/christine/GSE145705/GSE145705_design.csv --genome GRCm38 -name GSE145705 
cd - 


##########################################################
# 3) helena
##########################################################
# 3.1)CNS Treg samples
cd /media/rad/HDD1/atacseq/helena/cnsTreg
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/helena/cnsTreg/nfCnsTreg_design.csv --genome GRCm38 --single_end -name cnsTreg1

##########################################################
# 4) Miguel
##########################################################
# 4.1)CRC colorectal cancer samples
cd /media/rad/HDD1/atacseq/miguel/fignColCancer
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/miguel/fignColCancer/nffignColCancer_design.csv --genome GRCm38 --single_end -name fignColCancer

# 4.2)TM colorectal cancer samples
cd /media/rad/HDD1/atacseq/miguel/TMColCancer
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/miguel/TMColCancer/nfTMColCancer_design.csv --genome GRCm38 -name TMColCancer

##########################################################
# 5) Roman
##########################################################
# 5.1)Fign induction samples
cd /media/rad/HDD1/atacseq/roman/inductionFign
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/roman/inductionFign/nfinductionFign_design.csv --genome GRCm38 --single_end -name inductionFign

# 5.2)Fign organoids samples
cd /media/rad/HDD1/atacseq/roman/organoidsFign
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/roman/organoidsFign/nforganoidsFign_design.csv --genome GRCm38 --single_end -name organoidsFign

##########################################################
# 6) Gaurav
##########################################################
# 6.1) Test paired-end samples
cd /media/rad/HDD1/atacseq/gaurav/petest
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/gaurav/petest/nfpetest_design.csv --genome GRCm38 -name TMColCancer

##########################################################
# 7) Thorsten
##########################################################
# 7.1)Thorsten FKO/PKF1OE
# # Copy the fastq files
# cd /media/nas/raw/TUM_Nextseq/200924_NB501802_0283_AHKW3YBGXG/Data/Intensities/BaseCalls/A0000559_fastq_ATAC
# ls A0000*{FKO,PKF}*.fastq.gz | parallel --progress --eta -j 16 "rsync --ignore-existing -arzPR {} /media/rad/HDD1/atacseq/thorsten/atacfkopkf1oe/fastq"
# cd -

# Renamed PKF1O057E to PKF1OE057
# mv /media/rad/HDD1/atacseq/thorsten/atacfkopkf1oe/fastq/A0000559-021-PKF1O057E-P14-atac_S3_R1_001.fastq.gz /media/rad/HDD1/atacseq/thorsten/atacfkopkf1oe/fastq/A0000559-021-PKF1OE057-P14-atac_S3_R1_001.fastq.gz

# Run the pipeline
cd /media/rad/HDD1/atacseq/thorsten/atacfkopkf1oe
nextflow run /home/rad/users/gaurav/projects/workflows/nfatacseq --input /media/rad/HDD1/atacseq/thorsten/atacfkopkf1oe/nfatacfkopkf1oe_design.csv --genome GRCm38 --single_end -name atacfkopkf1oe


---
title: "Deep-amplicon sequencing of the complete beta tubulin gene in T. trichiura"
author: "Javier Gandaseggui"
date: '08/08/2024'
raw data: SRA BioProject ID PRJNA1145892
output file: md document
---

# Code used to perform the analysis
## Creating working environment

```bash
#Creating working directory
mkdir /lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI
cd /lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI
WORKING_DIR=/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI

#And working enviroment
mkdir 01_REF 02_RAW 03_MAP

#gather the raw files in 02_RAW/ diretory
#bring the T. trichuria reference genome to 01_REF/ directory (trichuris_trichiura.renamed.fa)

#Loading enviroments
module load PaM/environment
/data/pam/installs/scripts/setup_farm.sh
```

## Mapping fo the raw files to the T. trichiura reference genome

```bash
#Mapping using Nextflow
cd 03_MAP
module load mapping-helminth/v1.0.9

#We have to generate a manifest file with the location of the paired reads and the sample name, which is:
"
ID,R1,R2
4-10-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-10-HEL_S4_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-10-HEL_S4_L001_R2_001.fastq
4-12-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-12-HEL_S5_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-12-HEL_S5_L001_R2_001.fastq
4-14-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-14-HEL_S6_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-14-HEL_S6_L001_R2_001.fastq
4-15-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-15-HEL_S7_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-15-HEL_S7_L001_R2_001.fastq
4-19-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-19-HEL_S9_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-19-HEL_S9_L001_R2_001.fastq
4-21-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-21-HEL_S8_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-21-HEL_S8_L001_R2_001.fastq
4-25-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-25-HEL_S10_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-25-HEL_S10_L001_R2_001.fastq
4-26-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-26-HEL_S11_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-26-HEL_S11_L001_R2_001.fastq
4-27-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-27-HEL_S12_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-27-HEL_S12_L001_R2_001.fastq
4-28-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-28-HEL_S13_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-28-HEL_S13_L001_R2_001.fastq
4-2-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-2-HEL_S1_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-2-HEL_S1_L001_R2_001.fastq
4-31-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-31-HEL_S14_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-31-HEL_S14_L001_R2_001.fastq
4-32-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-32-HEL_S15_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-32-HEL_S15_L001_R2_001.fastq
4-33-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-33-HEL_S16_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-33-HEL_S16_L001_R2_001.fastq
4-35-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-35-HEL_S17_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-35-HEL_S17_L001_R2_001.fastq
4-40-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-40-HEL_S18_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-40-HEL_S18_L001_R2_001.fastq
4-45-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-45-HEL_S19_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-45-HEL_S19_L001_R2_001.fastq
4-46-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-46-HEL_S20_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-46-HEL_S20_L001_R2_001.fastq
4-4-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-4-HEL_S2_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-4-HEL_S2_L001_R2_001.fastq
4-51-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-51-HEL_S21_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-51-HEL_S21_L001_R2_001.fastq
4-53-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-53-HEL_S22_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-53-HEL_S22_L001_R2_001.fastq
4-54-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-54-HEL_S23_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-54-HEL_S23_L001_R2_001.fastq
4-55-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-55-HEL_S24_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-55-HEL_S24_L001_R2_001.fastq
4-59-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-59-HEL_S25_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-59-HEL_S25_L001_R2_001.fastq
4-61-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-61-HEL_S26_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-61-HEL_S26_L001_R2_001.fastq
4-62-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-62-HEL_S27_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-62-HEL_S27_L001_R2_001.fastq
4-7-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-7-HEL_S3_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/4-7-HEL_S3_L001_R2_001.fastq
6-10-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-10-HEL_S31_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-10-HEL_S31_L001_R2_001.fastq
6-11-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-11-HEL_S32_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-11-HEL_S32_L001_R2_001.fastq
6-12-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-12-HEL_S33_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-12-HEL_S33_L001_R2_001.fastq
6-13-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-13-HEL_S34_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-13-HEL_S34_L001_R2_001.fastq
6-16-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-16-HEL_S35_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-16-HEL_S35_L001_R2_001.fastq
6-20-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-20-HEL_S36_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-20-HEL_S36_L001_R2_001.fastq
6-23-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-23-HEL_S37_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-23-HEL_S37_L001_R2_001.fastq
6-26-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-26-HEL_S38_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-26-HEL_S38_L001_R2_001.fastq
6-27-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-27-HEL_S39_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-27-HEL_S39_L001_R2_001.fastq
6-29-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-29-HEL_S40_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-29-HEL_S40_L001_R2_001.fastq
6-2-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-2-HEL_S28_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-2-HEL_S28_L001_R2_001.fastq
6-39-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-39-HEL_S41_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-39-HEL_S41_L001_R2_001.fastq
6-3-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-3-HEL_S29_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-3-HEL_S29_L001_R2_001.fastq
6-40-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-40-HEL_S42_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-40-HEL_S42_L001_R2_001.fastq
6-45-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-45-HEL_S43_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-45-HEL_S43_L001_R2_001.fastq
6-46-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-46-HEL_S44_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-46-HEL_S44_L001_R2_001.fastq
6-48-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-48-HEL_S45_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-48-HEL_S45_L001_R2_001.fastq
6-49-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-49-HEL_S46_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-49-HEL_S46_L001_R2_001.fastq
6-9-HEL,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-9-HEL_S30_L001_R1_001.fastq,/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/02_RAW/ungziped/6-9-HEL_S30_L001_R2_001.fastq
"

#And mapping using the command
mapping-helminth --input deep_amp_manifest.csv --reference /lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI/01_REF/trichuris_trichiura.renamed.fa
#bsub -q oversubscribed -J jg_pool_sub -o jg_pool_sub.out -e jg_pool_sub.err -R "select[mem>4000] rusage[mem=4000]" -M4000 'mapping-helminth --input /lustre/scratch125/pam/teams/team333/jg34/STOP_pooled_samples/03_MAP/pooled_samples_manifest.csv --reference /lustre/scratch125/pam/teams/team333/jg34/STOP_pooled_samples/01_REF/STH_combined.renamed.fa'

#According to multiqc results, we have to remove sicne those samlpes has very low percentaje of mapped reads:
#(looks like the ampicon was not specific)
4-10-HEL
4-2-HEL
4-26-HEL
4-32-HEL
4-4-HEL
6-16-HEL
6-9-HEL
```

## Let's check the position of the amplicon of the beta tubuline in the new genome

```bash
cd ${WORKING_DIR}/01_REF/
module load mummer4/4.0.0rc1
module load bedtools/2.31.0--hf5e1c6e_3
module load samtools

#From mapping, I know that the reads went to the 95_trichuris_trichiura chr
#Let's narrow the annealong to that chr in order to avoid unespecific matches
samtools faidx ${REF_DIR}/trichuris_trichiura.renamed.fa 95_trichuris_trichiura > 95_trichuris_trichiura.fa

#Let's see waht we found
#Download the gen from NCBI Gen BAnl Accs no.: AF034219.1 - dowload fasta and renamed to Trichuris-beta-tub.fa
nucmer -c 100 -p nucmer ${WORKING_DIR}/01_REF/95_trichuris_trichiura.fa ${WORKING_DIR}/01_REF/Trichuris-beta-tub.fa
show-coords -r -c -l nucmer.delta > nucmer.coords
cat nucmer.coords

"NUCMER

    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]
===============================================================================================================================
10684111 10686714  |     2608        1  |     2604     2608  |    99.31  | 11299416     2611  |     0.02    99.89  | 95_trichuris_trichiura     Trichuris-beta-tub"

#Let's generate a bed tool with the coordinates from the exons and introns from a nucmer coordinates
#First we get the sequences of the exons from ncbi and generate a fasta file named beta_tub_exons.fa
#We have prepared a single line fasta with the ionformation we have in NCBI, let's preapre a multilne fasta
fold -b -w 60 beta_tub_exons.fa > beta_tub_exons_folded.fa
nucmer -c 100 -p beta_tub_exons ${WORKING_DIR}/01_REF/95_trichuris_trichiura.fa ${WORKING_DIR}/01_REF/beta_tub_exons_folded.fa
show-coords -r -c -l beta_tub_exons.delta > beta_tub_exons.coords

"NUCMER

    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]
===============================================================================================================================
10684167 10684696  |      530        1  |      530      530  |    99.25  | 11299416      530  |     0.00   100.00  | 95_trichuris_trichiura     EXON6
10684751 10684937  |      187        1  |      187      187  |   100.00  | 11299416      187  |     0.00   100.00  | 95_trichuris_trichiura     EXON5
10684987 10685436  |      450        1  |      450      450  |    99.78  | 11299416      450  |     0.00   100.00  | 95_trichuris_trichiura     EXON4
10685483 10685727  |      245        1  |      245      245  |    99.18  | 11299416      245  |     0.00   100.00  | 95_trichuris_trichiura     EXON3
10685780 10685999  |      220        1  |      220      220  |   100.00  | 11299416      220  |     0.00   100.00  | 95_trichuris_trichiura     EXON2
10686294 10686647  |      355        1  |      354      355  |    99.44  | 11299416      355  |     0.00   100.00  | 95_trichuris_trichiura     EXON1"

#Now, we conver this to bed format
awk -F '|' '!/^=/{gsub(/^[ \t]+|[ \t]+$/, "", $0); split($1, coords, " +"); split($15, tags, " +"); print tags[1] "\t" coords[1] "\t" coords[2] "\t" tags[2] "\t0\t+"}' beta_tub_exons.coords > beta_tub_exons_folded.bed
#And some brief modification in nano
#Now we use that file to obtain the introns
samtools faidx 95_trichuris_trichiura.fa
cut -f 1,2 95_trichuris_trichiura.fa.fai > 95_trichuris_trichiura.genome.bed
bedtools complement -i beta_tub_exons_folded.bed -g 95_trichuris_trichiura.genome.bed > beta_tub_introns.bed
#And some nano modifications to get only the amplicon based on Trichuris-beta-tub.fa sequence
```


## Estimating coverage and variant frequency using grenedalf
Grenedalf is specifically designed to analyse pools of individuals

```bash
#Let's created the folder
mkdir ${WORKING_DIR}/05_DIV/
mkdir ${WORKING_DIR}/05_DIV/grenedalf
cd ${WORKING_DIR}/05_DIV/grenedalf

#Let's generate the enviromentla variables
BAM_LIST=${WORKING_DIR}/04_VARIANTS/bam.list #this include the mapped files
REFERENCE=${WORKING_DIR}/01_REF/trichuris_trichiura.renamed.fa


#Loading the modules
module load grenedalf/0.3.0
module load samtools/1.9--h91753b0_8

#Let's select those SNPs in the chr of interest (95_trichuris_trichiura) - some reads went to other chr, so they were not specific
mkdir bam_files_chr

while read BAM; do \
SAMPLE=$( echo ${BAM} | awk -F '/' '{print $NF}' | sed -e 's/.bam//g' )
samtools index ${BAM}
samtools view -b -h ${BAM} "95_trichuris_trichiura" > bam_files_chr/${SAMPLE}.chr.bam;
done < ${BAM_LIST}

#Let's estimate the variant frequency
mkdir freq
grenedalf frequency \
--sam-path bam_files_chr/ --write-sample-counts --write-sample-coverage --write-sample-ref-freq --out-dir freq/ \
--file-prefix per_sample

#Now per sample and per site diversity
mkdir div
grenedalf diversity \
--sam-path bam_files_chr/ \
--filter-sample-min-count 100 --window-type single --out-dir div/ --pool-sizes 10 \
--measure theta-pi --file-prefix per_sample

grenedalf diversity \
--sam-path bam_files_chr/ \
--filter-sample-min-count 10 --window-type genome --out-dir div/ --pool-sizes 10 \
--filter-region-bed beta_tub.bed --measure theta-pi --file-prefix per_sample_total
```


################################################################################################################################################################












################################################################################################################################################################

#Let's call the variants using genotypes as previously done
#working environment
WORKING_DIR=/lustre/scratch125/pam/teams/team333/jg34/DEEP_AMP_TRI
mkdir ${WORKING_DIR}/04_VARIANTS
cd ${WORKING_DIR}/04_VARIANTS

#Loeading modules
module load gatk/4.1.4.1
module load vcftools/0.1.16-c4
module load samtools/1.14--hb421002_0

#indexing the ref
REFERENCE=${WORKING_DIR}/01_REF/trichuris_trichiura.renamed.fa
REF_DIR=${WORKING_DIR}/01_REF
# CreateSequenceDictionary
samtools faidx ${REF_DIR}/trichuris_trichiura.renamed.fa
gatk --java-options "-Xmx8g -Xms4g -Djava.io.tmpdir=${TMP_DIR}" CreateSequenceDictionary \
     --REFERENCE ${REF_DIR}/trichuris_trichiura.renamed.fa \
     --OUTPUT ${REF_DIR}/trichuris_trichiura.renamed.dict \
     --spark-runner LOCAL

#Make GVCFs per sample
#new folder
mkdir ${WORKING_DIR}/04_VARIANTS/GVCFS
cd ${WORKING_DIR}/04_VARIANTS/GVCFS
# create bam list using full path to bams - this allows bams to be anywhere
ls ${WORKING_DIR}/03_MAP/results/*/*.bam > ${WORKING_DIR}/04_VARIANTS/bam.list
BAM_LIST=${WORKING_DIR}/04_VARIANTS/bam.list

#Now we I edit the bam ist mannually to remove those samples that have no mmaped reads in the multiqc:
#4-10-HEL, 4-2-HEL, 4-26-HEL, 4-32-HEL, 4-4-HEL, 6-16-HEL, 6-9-HEL

#Makin the GVCF per sample
while read BAM; do \
	SAMPLE=$( echo ${BAM} | awk -F '/' '{print $NF}' | sed -e 's/.bam//g' )
	mkdir ${SAMPLE}_GATK_HC_GVCF
	mkdir ${SAMPLE}_GATK_HC_GVCF/LOGFILES
	gatk HaplotypeCaller \
          --input ${BAM} \
          --output ${SAMPLE}_GATK_HC_GVCF/${SAMPLE}.tmp.gvcf.gz \
          --reference ${REFERENCE} \
          --heterozygosity 0.015 \
          --indel-heterozygosity 0.01 \
          --annotation DepthPerAlleleBySample --annotation Coverage --annotation ExcessHet --annotation FisherStrand --annotation MappingQualityRankSumTest --annotation StrandOddsRatio --annotation RMSMappingQuality --annotation ReadPosRankSumTest --annotation DepthPerSampleHC --annotation QualByDepth \
          --min-base-quality-score 20 --minimum-mapping-quality 30 --standard-min-confidence-threshold-for-calling 30 \
          --emit-ref-confidence GVCF > ${SAMPLE}_GATK_HC_GVCF/run_hc_${SAMPLE}.tmp.job;
done < ${BAM_LIST}

#Combining all gvcf in a single file
ls -1 ${WORKING_DIR}/04_VARIANTS/GVCFS/*/*gz > ${WORKING_DIR}/04_VARIANTS/GVCFS/gvcf.list
GVCF_LIST=${WORKING_DIR}/04_VARIANTS/GVCFS/gvcf.list
cd ${WORKING_DIR}/04_VARIANTS/

echo -e "gatk CombineGVCFs -R ${REFERENCE} \\" > run_merge_vcfs.tmp.sh
     while read GVCF; do
          echo -e "--variant ${GVCF} \\" >> run_merge_vcfs.tmp.sh;
     done < ${GVCF_LIST}
echo -e "--output merged.gvcf.gz" >> run_merge_vcfs.tmp.sh

bash run_merge_vcfs.tmp.sh

#And genotype
gatk GenotypeGVCFs \
--reference ${REFERENCE} \
--V merged.gvcf.gz \
--heterozygosity 0.015 \
--indel-heterozygosity 0.01 \
--annotation DepthPerAlleleBySample --annotation Coverage --annotation ExcessHet --annotation FisherStrand --annotation MappingQualityRankSumTest --annotation StrandOddsRatio --annotation RMSMappingQuality --annotation ReadPosRankSumTest --annotation DepthPerSampleHC --annotation QualByDepth \
--output deep_amp_tri_genome.vcf

vcftools --vcf deep_amp_tri_genome.vcf --out miss_indv --missing-indv


#Now we have to filter - are we sure? Let's filter anyway
vcftools \
--vcf deep_amp_tri_genome.vcf \
--remove-filtered-geno-all \
--remove-filtered-all \
--min-alleles 1 \
--max-alleles 4 \
--maf 0.001 \
--recode \
--recode-INFO-all \
--out deep_amp_tri_genome.filtered

################################################################################################################################################################



#############

#Let's generate a plot with the beta tubuline variations obserned and per-site diversity
mkdir ${WORKING_DIR}/05_DIV
cd ${WORKING_DIR}/05_DIV
VCF=${WORKING_DIR}/04_VARIANTS/deep_amp_tri_genome.filtered.noInDels.recode.vcf
module load vcftools/0.1.16-c4

################################################################################################################################################################



#Let's see what is going on with the second beta tub
mkdir ${WORKING_DIR}/06_BETA_2
cd ${WORKING_DIR}/06_BETA_2

#Copying Steve nivlear vcf of trichuris
cp /lustre/scratch125/pam/teams/team333/sd21/trichuris_trichiura/04_VARIANTS/GATK_HC_MERGED/Trichuris_trichiura.cohort.nuclear_variants.final.recode.vcf.gz .
cp /lustre/scratch125/pam/teams/team333/sd21/trichuris_trichiura/04_VARIANTS/GATK_HC_MERGED/mod_human_samples.list .

#Now, we have to figure out the position of the other beta tub, as well as the exosn and introns
module load mummer4/4.0.0rc1
module load samtools
module load bedtools/2.31.0--hf5e1c6e_3

#A TBLASTN in worm base parasite showed the other isotype named TTRE8258
#Download it and generate the file gen.fa

nucmer -c 100 -p beta_2 ${WORKING_DIR}/01_REF/trichuris_trichiura.renamed.fa gen.fa
show-coords -r -c -l beta_2.delta > beta_2.coords
cat beta_2.coords

"NUCMER

    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]
===============================================================================================================================
27478327 27481371  |        1     3045  |     3045     3045  |   100.00  | 29164577     3045  |     0.01   100.00  | 102_trichuris_trichiura    TTRE_chr2"

#It is in scaffold 102
samtools faidx ${WORKING_DIR}/01_REF/trichuris_trichiura.renamed.fa 102_trichuris_trichiura > 102_trichuris_trichiura.fa

#Form WBP we also download the information of the exons and name the file as exon.fa
nucmer -c 100 -p beta_tub_2_exons 102_trichuris_trichiura.fa exon.fa
show-coords -r -c -l beta_tub_2_exons.delta > beta_tub_2_exons.coords
cat beta_tub_2_exons.coords

"NUCMER

    [S1]     [E1]  |     [S2]     [E2]  |  [LEN 1]  [LEN 2]  |  [% IDY]  |  [LEN R]  [LEN Q]  |  [COV R]  [COV Q]  | [TAGS]
===============================================================================================================================
27478377 27478487  |        1      111  |      111      111  |   100.00  | 29164577      111  |     0.00   100.00  | 102_trichuris_trichiura    TTRE8258_t1
27478650 27478869  |        1      220  |      220      220  |   100.00  | 29164577      220  |     0.00   100.00  | 102_trichuris_trichiura    TTRE8258_t1
27478923 27479167  |        1      245  |      245      245  |   100.00  | 29164577      245  |     0.00   100.00  | 102_trichuris_trichiura    TTRE8258_t1
27479219 27479525  |        1      307  |      307      307  |   100.00  | 29164577      307  |     0.00   100.00  | 102_trichuris_trichiura    TTRE8258_t1
27479572 27479717  |        1      146  |      146      146  |   100.00  | 29164577      146  |     0.00   100.00  | 102_trichuris_trichiura    TTRE8258_t1
27480871 27481057  |        1      187  |      187      187  |   100.00  | 29164577      187  |     0.00   100.00  | 102_trichuris_trichiura    TTRE8258_t1
27481108 27481321  |        1      214  |      214      214  |   100.00  | 29164577      214  |     0.00   100.00  | 102_trichuris_trichiura    TTRE8258_t1"

#Now, we conver this to bed format
awk -F '|' '!/^=/{gsub(/^[ \t]+|[ \t]+$/, "", $0); split($1, coords, " +"); split($15, tags, " +"); print tags[1] "\t" coords[1] "\t" coords[2] "\t" tags[2] "\t0\t+"}' beta_tub_2_exons.coords > beta_tub_2_exons.bed
#And some brief modification in nano
#Now we use that file to obtain the introns
samtools faidx 102_trichuris_trichiura.fa
cut -f 1,2 102_trichuris_trichiura.fa.fai > 102_trichuris_trichiura.genome.bed
bedtools complement -i beta_tub_2_exons.bed -g 102_trichuris_trichiura.genome.bed > beta_tub_2_introns.bed

#Let's select the infor using Steve scripts
module load vcftools/0.1.16-c4

#I have to change the chr names
#in the vcf file is: Trichuris_trichiura_2_001,length=29164577,assembly=trichuris_trichiura.fa
sed '/^102_trichuris_trichiura/s/^/Trichuris_trichiura_2_001,length=29164577,assembly=trichuris_trichiura.fa/' beta_tub_2_gen.bed > beta_tub_2_gen.renamed.bed

vcftools \
      --gzvcf Trichuris_trichiura.cohort.nuclear_variants.final.recode.vcf \
      --bed beta_tub_2_gen.renamed.bed \
      --site-pi \
      --keep mod_human_samples.list \
      --maf 0.01 \
      --out BZ_nuc_div

#> After filtering, kept 30 out of 61 Individuals
#> After filtering, kept 205 out of a possible 6933531 Sites

vcftools \
      --gzvcf Trichuris_trichiura.cohort.nuclear_variants.final.recode.vcf \
      --bed beta_tub_2_gen.renamed.bed \
      --freq \
      --keep mod_human_samples.list \
      --maf 0.01 \
      --out BZ_allele_freq
#> After filtering, kept 30 out of 61 Individuals
#> After filtering, kept 205 out of a possible 6933531 Sites


################################################################################################################################################################

#SNPeff
mkdir ${WORKING_DIR}/07_SNPeff
module load vcftools/0.1.16-c4
module load snpeff/5.1d--hdfd78af_0
module load bcftools/1.14--h88f3f91_0

#Let's select the SNPs that are higher to 10% in grenedalf - that data come from R
#I did this bed file in R
"chrom     start     end
95_trichuris_trichiura	10684175	10684175
95_trichuris_trichiura	10684243	10684243
95_trichuris_trichiura	10684295	10684295
95_trichuris_trichiura	10684296	10684296
95_trichuris_trichiura	10684309	10684309
95_trichuris_trichiura	10684354	10684354
95_trichuris_trichiura	10684360	10684360
95_trichuris_trichiura	10684369	10684369
95_trichuris_trichiura	10684370	10684370
95_trichuris_trichiura	10684449	10684449
95_trichuris_trichiura	10684471	10684471
95_trichuris_trichiura	10684520	10684520
95_trichuris_trichiura	10684635	10684635
95_trichuris_trichiura	10684845	10684845
95_trichuris_trichiura	10684857	10684857
95_trichuris_trichiura	10684929	10684929
95_trichuris_trichiura	10685444	10685444
95_trichuris_trichiura	10685606	10685606
95_trichuris_trichiura	10685714	10685714
95_trichuris_trichiura	10685740	10685740
95_trichuris_trichiura	10685817	10685817
95_trichuris_trichiura	10686038	10686038
95_trichuris_trichiura	10686121	10686121
95_trichuris_trichiura	10686196	10686196
95_trichuris_trichiura	10686201	10686201
95_trichuris_trichiura	10686256	10686256
95_trichuris_trichiura	10686375	10686375
95_trichuris_trichiura	10686478	10686478
95_trichuris_trichiura	10686492	10686492
95_trichuris_trichiura	10686502	10686502
95_trichuris_trichiura	10686505	10686505
95_trichuris_trichiura	10686528	10686528
95_trichuris_trichiura	10686559	10686559
95_trichuris_trichiura	10686568	10686568
95_trichuris_trichiura	10686580	10686580
95_trichuris_trichiura	10686637	10686637
95_trichuris_trichiura	10686656	10686656
95_trichuris_trichiura	10686661	10686661
95_trichuris_trichiura	10686668	10686668"

#Filter the varants
vcftools \
--gzvcf ../04_VARIANTS/deep_amp_tri_genome.vcf \
--bed grenedalf_snps_10.bed \
--out grenedalf_snps_10 \
--recode --keep-INFO-all

#Let's firue out the chr name based on the lech of our ref genome
REFERENCE=${WORKING_DIR}/01_REF/trichuris_trichiura.renamed.fa
awk '/^>/ { if (seqlen) {
              print seqlen
              }
            print

            seqtotal+=seqlen
            seqlen=0
            seq+=1
            next
            }
    {
    seqlen += length($0)
    }     
    END{print seqlen
        print seq" sequences, total length " seqtotal+seqlen
    }' ${REFERENCE}

#The length of the chromosome where the beta tub is in our samples it is: >95_trichuris_trichiura = 11299416

#Get the annotation
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS19/species/trichuris_trichiura/PRJEB535/trichuris_trichiura.PRJEB535.WBPS19.annotations.gff3.gz
gunzip trichuris_trichiura.PRJEB535.WBPS19.annotations.gff3.gz
#In the annotation we have:
##sequence-region TTRE_chr1_scaffold1 1 11299416

#We have to renme the vcf file using TTRE_chr1_scaffold1 instead of 95_trichuris_trichiura
echo "95_trichuris_trichiura TTRE_chr1_scaffold1" >> chr_name_conv.txt
bcftools annotate --rename-chrs chr_name_conv.txt grenedalf_snps_10.recode.vcf > grenedalf_snps_10.renamed.vcf

#Let's reun snpEff
#Make data directory and a trichuris directory inside
mkdir data
mkdir data/trichuris
#cp the annotitaion to trichuris file and renamed it as genes.gff
#nake he following file using nano
"# Haemonchus contortus chromosomes V4
trichuris.genome : trichuris"
#mk a genomes diretory, get the trichuris ref genom and renamed it as trichuris.fa
#build the database
snpEff build -noCheckCds -gff3 -v trichuris #I have to use that flag, not sure why
#run snpEff
#All annotations
snpEff -v trichuris grenedalf_snps_10.renamed.vcf > grenedalf_snps_10.ann.vcf
snpEff -v -onlyProtein trichuris grenedalf_snps_10.renamed.vcf > grenedalf_snps_10.ann.protein.vcf
#Let's do the file "my_transcripts.txt" as:
"transcript:TTRE11895_t1"
#And run
snpEff -v -onlyTr my_transcripts.txt trichuris grenedalf_snps_10.renamed.vcf > grenedalf_snps_10.ann.trans.vcf
#same result, still don't know the meaning









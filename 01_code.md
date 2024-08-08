---
title: "Deep-amplicon sequencing of the complete beta tubulin gene in T. trichiura"
author: "Javier Gandasegui"
date: '08/08/2024'
raw data: SRA BioProject ID PRJNA1145892
output file: md
---

# Main code

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
## Now, we use the output files to generate some plots and perform some analysys in R
### First, the plots for coverage and variant freq are prapred

```R
library(tidyverse)
library(ggsci)
library(patchwork)
library(viridis)
library(ggrepel)
library(gcookbook)
library(janitor)
library(tidyverse)
library(ggpubr)
library(ggridges)
library(scales)
library(readxl)

#setting WD
setwd("~/R/DEEP_amp_seq")

#Let's plot the beta tubuline scheme alongised the coverage and variant freq
#Let's generate a scheme of the beta tubuline
exons <- read.table("beta_tub_var/beta_tub_exons_folded.bed", header=F)
colnames(exons) <- c("chrom", "start", "end")

introns <- read.table("beta_tub_var/beta_tub_introns.bed", header=F)
colnames(introns) <- c("chrom", "start", "end")

resistant_snps <- read.table("beta_tub_var/btubulin.canonicalresistantSNPs.bed",header=T)

#Let's explore COV and SNP frequencies
per_samplefrequency <- read_csv("grenedalf/per_samplefrequency.csv")
#Let's write a csv for supp material
freq_to_csv <- per_samplefrequency
colnames(freq_to_csv) <- str_remove(colnames(freq_to_csv), '.chr.1')
colnames(freq_to_csv) <- str_replace(colnames(freq_to_csv), 'REF_CNT', 'ref_read_count')
colnames(freq_to_csv) <- str_replace(colnames(freq_to_csv), 'ALT_CNT', 'alt_read_count')
colnames(freq_to_csv) <- str_replace(colnames(freq_to_csv), 'COV', 'coverage')
colnames(freq_to_csv) <- str_replace(colnames(freq_to_csv), 'FREQ', 'ref_freq')
write_csv(freq_to_csv, 'Table_S2.csv')

#Let's estimate pre and post coverage per position
cov_pre <- per_samplefrequency %>%
  select(POS, contains('COV')) %>%
  mutate(., pre_cov_mean = rowMeans(select(., contains("4-")), na.rm = TRUE))

cov_pre_post <- cov_pre %>%
  select(contains('COV')) %>%
  mutate(., post_cov_mean = rowMeans(select(., !contains("4-")), na.rm = TRUE)) %>%
  select(post_cov_mean) %>%
  cbind(cov_pre, .)

cov_pre_post_all <- cov_pre %>%
  select(contains('COV')) %>%
  mutate(., all_cov_mean = rowMeans(.,na.rm = TRUE)) %>%
  select(all_cov_mean) %>%
  cbind(cov_pre_post, .)

#And plot
cov_to_plot <- cov_pre_post_all %>%
  select(POS, pre_cov_mean, post_cov_mean) %>%
  tidyr::gather(., point, coverage, -POS)

cov_to_plot$point <- str_replace(cov_to_plot$point, "pre_cov_mean", "Pre-treatment")
cov_to_plot$point <- str_replace(cov_to_plot$point, "post_cov_mean", "Post-treatment")
cov_to_plot$point <- factor(cov_to_plot$point, levels=c("Pre-treatment", "Post-treatment"))

plot_cov  <- ggplot() +
  geom_rect(data=introns, aes(xmin=start,ymin=-220, xmax=end, ymax=220), fill="grey50") +
  geom_rect(data=exons, aes(xmin=start,ymin=-1450, xmax=end, ymax=1450), fill="grey90") +
  #geom_segment(data=resistant_snps, aes(x=pos, xend=pos, y=-1450, yend=1450),col="red", size=1) +
  #geom_text_repel(data=resistant_snps, aes(x=pos, y=0, label=name), col="red", box.padding = 0.5, max.overlaps = Inf, nudge_y = 5000) +
  geom_line(data =cov_to_plot, aes(x=POS, y = coverage, col = point), size=1, alpha=0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlim(10684100, 10686800) +
  scale_color_discrete(name = " ") +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  ylab("Covegare") + xlab("Position")

#Let's estimate pre and post SNP freq
freq_pre <- per_samplefrequency %>%
  select(POS, contains('FREQ')) %>%
  mutate(., pre_freq_mean = rowMeans(select(., contains("4-")), na.rm = TRUE))

freq_pre_post <- freq_pre %>%
  select(contains('FREQ')) %>%
  mutate(., post_freq_mean = rowMeans(select(., !contains("4-")), na.rm = TRUE)) %>%
  select(post_freq_mean) %>%
  cbind(freq_pre, .)

freq_pre_post_all <- freq_pre %>%
  select(contains('FREQ')) %>%
  mutate(., all_freq_mean = rowMeans(.,na.rm = TRUE)) %>%
  select(all_freq_mean) %>%
  cbind(freq_pre_post, .)

#And plot
freq_to_plot <- freq_pre_post_all %>%
  select(POS, pre_freq_mean, post_freq_mean) %>%
  tidyr::gather(., point, frequency, -POS)

freq_to_plot_freq <- as_tibble(sapply(freq_to_plot$frequency, function(x) if(is.numeric(x)) 1-x else ""))
freq_to_plot <- cbind(select(freq_to_plot, POS, point), freq_to_plot_freq)

freq_to_plot$point <- str_replace(cov_to_plot$point, "pre_freq_mean", "Pre-treatment")
freq_to_plot$point <- str_replace(cov_to_plot$point, "post_freq_mean", "Post-treatment")
freq_to_plot$point <- factor(freq_to_plot$point, levels=c("Pre-treatment", "Post-treatment"))

plot_freq <- ggplot() +
  geom_rect(data=introns, aes(xmin=start,ymin=-0.015, xmax=end, ymax=0.015), fill="grey50") +
  geom_rect(data=exons, aes(xmin=start,ymin=-0.10, xmax=end, ymax=0.10), fill="grey90") +
  geom_segment(data=resistant_snps, aes(x=pos, xend=pos, y=-0.10, yend=0.10),col="red", size=1) +
  geom_text_repel(data=resistant_snps, aes(x=pos, y=0, label=name), col="red", box.padding = 0.5, max.overlaps = Inf, nudge_y = 0.5) +
  geom_line(data=freq_to_plot, aes(x=POS, y = value, col = point), size=1, alpha=0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlim(10684100, 10686800) +
  scale_color_discrete(name = " ") +
  ylab("Variant frequency") + xlab("Position")  +
  scale_y_continuous(labels = label_number(accuracy = 0.001))
```

### And now, we continue with analysing the nuc. diversity

```R
#Let's add metadata to explore diversity, EPG and Ct values
irsi_amp <- read_excel("databases/irsi_amp_mod.xlsx") #this metadata contains partitipant information and it is not public

#let's add the egg count data per sample to get the plots of pi and hapo vs eggs
stool_data_1 <- read_excel("databases/MARS_bases_de_dados_resultados_V0V1V2_Urina_e_Fezes_V4V6_Fezes_08.04.2020.xlsx",
                           sheet = "MARS-1_v0v1v2v4v6") #same here

stool_data_2 <- read_excel("databases/MARS_bases_de_dados_resultados_V0V1V2_Urina_e_Fezes_V4V6_Fezes_08.04.2020.xlsx",
                           sheet = "MARS-2_v0v1v2v4v6") #and here

stool_data <- rbind(stool_data_1, stool_data_2)

egg_count_v4 <- irsi_amp %>% filter(pre_post == 'pre') %>%
  left_join(., stool_data, by = c('nida' = 'NIDA_V4')) %>%
  select("n_estudo","nida","code","pre_post",
         "triplex_tri","micro_tri","PCR_pyro","167",             
         "198","200", "tri_kk1_v4_ovos", "tri_kk2_v4_ovos",
         "tri_kk3_v4_ovos", "tri_kk4_v4_ovos")
colnames(egg_count_v4) <- str_remove(colnames(egg_count_v4), 'v4_')

egg_count_v6 <- irsi_amp %>% filter(pre_post == 'post') %>%
  left_join(., stool_data, by = c('nida' = 'NIDA_V6')) %>%
  select("n_estudo","nida","code","pre_post",
         "triplex_tri","micro_tri","PCR_pyro","167",             
         "198","200", "tri_kk1_v6_ovos", "tri_kk2_v6_ovos",
         "tri_kk3_v6_ovos", "tri_kk4_v6_ovos")
colnames(egg_count_v6) <- str_remove(colnames(egg_count_v4), 'v6_')

egg_count <- rbind(egg_count_v4, egg_count_v6)

egg_count$tri_kk1_ovos <- as.numeric(egg_count$tri_kk1_ovos)
egg_count$tri_kk2_ovos <- as.numeric(egg_count$tri_kk2_ovos)
egg_count$tri_kk3_ovos <- as.numeric(egg_count$tri_kk3_ovos)
egg_count$tri_kk4_ovos <- as.numeric(egg_count$tri_kk4_ovos)

egg_count <- egg_count %>%
  mutate(epg = (tri_kk1_ovos + tri_kk2_ovos + tri_kk3_ovos + tri_kk4_ovos)/4*24)

qPCR_egg_count <- read_excel("databases/MARS_PCR_all_1190_Leiden_V210104.xlsx") %>%
  select(pin, ct_tt, pcr_tt) %>%
  left_join(egg_count, .,  by = c('nida' = 'pin'))

#Let's see diversity and intensity of infection
#Let's load per sample diversity obtained in Grenedalf
per_sample_total_div <- read_csv("grenedalf/per_sample_totaldiversity.csv")

df2 <- data.frame(t(select(per_sample_total_div, contains('pi_rel'))[]))
colnames(df2) <- per_sample_total_div[, 1]
per_sample_div <- rownames_to_column(df2, var = 'sample')
per_sample_div$point <- ifelse(grepl("4-", per_sample_div$sample), "Pre-treatment", "Post-treatment")
colnames(per_sample_div) <- c('sample', 'div', 'point')
per_sample_div$point <- factor(per_sample_div$point, levels=c("Pre-treatment", "Post-treatment"))


per_sample_div$sample <- str_remove(per_sample_div$sample, '-HEL.chr.1.theta_pi_rel')
qPCR_egg_count$ct_tt <- as.numeric(qPCR_egg_count$ct_tt)

qPCR_egg_count$pre_post <- str_replace(qPCR_egg_count$pre_post, "pre", "Pre-treatment")
qPCR_egg_count$pre_post <- str_replace(qPCR_egg_count$pre_post, "post", "Post-treatment")
qPCR_egg_count$pre_post <- factor(qPCR_egg_count$pre_post, levels=c("Pre-treatment", "Post-treatment"))

coor <- qPCR_egg_count %>%
  select(code, epg, ct_tt, pre_post) %>%
  filter(epg > 0) %>%
  left_join(per_sample_div, ., by = c('sample' = 'code')) %>%
  select(div, epg, ct_tt)

cor.test(coor$div, coor$epg, method=c('spearman'))
cor.test(coor$div, coor$ct_tt, method=c('spearman'))

div_vs_epg <- qPCR_egg_count %>%
  select(code, epg, ct_tt, pre_post) %>%
  left_join(per_sample_div, ., by = c('sample' = 'code')) %>%
  ggplot(aes(x=log(epg), y=div, col = pre_post)) +
  geom_point(size = 3) +
  #geom_smooth(method = lm, se = F, col= 'grey30', alpha = 0.6) +
  annotate("text", x = 6.9, y = 0.0048, col = 'grey32',size = 3,
           label = "rho = 0.39, p = 0.036") +
  scale_color_discrete(name = " ") +
  scale_color_discrete(name = " ") +
  labs(x = "Log(epg)" , y = "Nucleotide diversity (Pi)") +
  theme_bw()

div_vs_ct <- qPCR_egg_count %>%
  select(code, epg, ct_tt, pre_post) %>%
  left_join(per_sample_div, ., by = c('sample' = 'code')) %>%
  ggplot(aes(x=ct_tt, y=div, col = pre_post)) +
  geom_point(size = 3) +
  #geom_smooth(method = lm, se = F, col= 'grey30', alpha = 0.6) +
  annotate("text", x = 31, y = 0.0048, col = 'grey32', size = 3,
           label = "rho = -0.26, p = 0.172") +
  scale_color_discrete(name = " ") +
  labs(x = "Ct value" , y = "Nucleotide diversity (Pi)") +
  theme_bw()

#Let's find the paired samples
qPCR_egg_count %>%
  filter(pre_post == 'Pre-treatment') %>%
  inner_join(., filter(qPCR_egg_count, pre_post == 'Post-treatment'), by = 'n_estudo') %>%
  select(code.x, code.y) %>% print()

#Pairs are:
#4-12   6-20  
#4-15   6-10  
#4-19   6-12  
#4-46   6-26  
#4-54   6-39  
#4-59   6-45

per_sample_div$pair[per_sample_div$sample == '4-12'] <- 'a'
per_sample_div$pair[per_sample_div$sample == '6-20'] <- 'a'

per_sample_div$pair[per_sample_div$sample == '4-15'] <- 'b'
per_sample_div$pair[per_sample_div$sample == '6-10'] <- 'b'

per_sample_div$pair[per_sample_div$sample == '4-19'] <- 'c'
per_sample_div$pair[per_sample_div$sample == '6-12'] <- 'c'

per_sample_div$pair[per_sample_div$sample == '4-46'] <- 'd'
per_sample_div$pair[per_sample_div$sample == '6-26'] <- 'd'

per_sample_div$pair[per_sample_div$sample == '4-54'] <- 'e'
per_sample_div$pair[per_sample_div$sample == '6-39'] <- 'e'

per_sample_div$pair[per_sample_div$sample == '4-59'] <- 'f'
per_sample_div$pair[per_sample_div$sample == '6-45'] <- 'f'

#Let's see differences in per-sample diversity pre and post treatment
boxplot_div <- ggplot(as_tibble(per_sample_div), aes(x=point, y=div, fill = point)) +
  geom_boxplot(alpha=0.2) +
  geom_point(aes(y=div, col = point)) +
  geom_line(aes(group=pair), color = 'grey32') +
  labs(y = "Nucleotide diversity (pi)", x = '', col = '') +
  #ylim(0.002, 0.0125) +
  annotate("text", x = 1.5, y = 0.01, col = 'grey32',size = 3,
           label = "p = 0.914") +
  theme_bw() + theme(legend.position = "none")
```

### Amd now I get some partitipant data and infection intensity
This infroamtion is used to create table 1

```R
library(haven)
part_data <- read_dta("databases/data_for_Javier_with_ages_OM_VN_01JUL2024.dta") #this mnetadata is not public

#Now get uniq study numbers from irsiamp database
part_data_table <- irsi_amp[!duplicated(irsi_amp$n_estudo), ] %>%
  select(n_estudo, pre_post) %>%
  left_join(., part_data, by = c('n_estudo' = 'study_numer')) %>%
  select(n_estudo, pre_post, age, SEXO, Posto_admistrativo) 

part_data_table$pre_post[part_data_table$n_estudo == "0459" | part_data_table$n_estudo == "0481" |
                          part_data_table$n_estudo == "0261" | part_data_table$n_estudo == "0603" |
                          part_data_table$n_estudo == "0699" | part_data_table$n_estudo == "0295"] <- 'paired'

part_data_table$SEXO <- tolower(part_data_table$SEXO)
part_data_table$SEXO <- str_replace(part_data_table$SEXO, 'homem', 'Male')
part_data_table$SEXO <- str_replace(part_data_table$SEXO, 'mulher', 'Female')

part_data_table$Posto_admistrativo <- str_replace_all(str_to_title(tolower(part_data_table$Posto_admistrativo)), ' ', '_')
part_data_table$Posto_admistrativo <- str_replace(part_data_table$Posto_admistrativo, '3_De_Fevreiro', '3_De_Fevereiro')

library(tableone)
colnames(part_data_table)
listVars<-c("age","SEXO","Posto_admistrativo")
catVars<-c("SEXO","Posto_admistrativo")
table1<-CreateTableOne(data=part_data_table,vars=listVars,factorVars=catVars,
                       addOverall=T,test=F,strata=c('pre_post'))

#Now we have to estimate the KK mean and Ct-median for each value
#KK
overall <- qPCR_egg_count %>%
  filter(epg > 0)
nrow(overall)
mean(overall$epg)
range(overall$epg)

pre_epg <- part_data_table %>%
  filter(pre_post == 'pre') %>%
  left_join(., qPCR_egg_count, by = 'n_estudo') %>%
  filter(epg > 0)
nrow(pre_epg)
mean(pre_epg$epg)
range(pre_epg$epg)

post_epg <- part_data_table %>%
  filter(pre_post == 'post') %>%
  left_join(., qPCR_egg_count, by = 'n_estudo') %>%
  filter(epg > 0)
nrow(post_epg)
mean(post_epg$epg)
range(post_epg$epg)

pair_pre_epg <- part_data_table %>%
  filter(pre_post == 'paired') %>%
  left_join(., qPCR_egg_count, by = 'n_estudo') %>%
  filter(epg > 0) %>%
  filter(pre_post.y == 'Pre-treatment')
nrow(pair_pre_epg)
mean(pair_pre_epg$epg)
range(pair_pre_epg$epg)

pair_post_epg <- part_data_table %>%
  filter(pre_post == 'paired') %>%
  left_join(., qPCR_egg_count, by = 'n_estudo') %>%
  filter(epg > 0) %>%
  filter(pre_post.y == 'Post-treatment')
nrow(pair_post_epg)
mean(pair_post_epg$epg)
range(pair_post_epg$epg)

#qPCR
median(qPCR_egg_count$ct_tt)
range(qPCR_egg_count$ct_tt)

pre_ct <- part_data_table %>%
  filter(pre_post == 'pre') %>%
  left_join(., qPCR_egg_count, by = 'n_estudo')
nrow(pre_ct)
mean(pre_ct$ct_tt)
range(pre_ct$ct_tt)

post_ct <- part_data_table %>%
  filter(pre_post == 'post') %>%
  left_join(., qPCR_egg_count, by = 'n_estudo')
nrow(post_ct)
mean(post_ct$ct_tt)
range(post_ct$ct_tt)

pair_pre_ct <- part_data_table %>%
  filter(pre_post == 'paired') %>%
  left_join(., qPCR_egg_count, by = 'n_estudo') %>%
  filter(pre_post.y == 'Pre-treatment')
nrow(pair_pre_ct)
mean(pair_pre_ct$ct_tt)
range(pair_pre_ct$ct_tt)

pair_post_ct <- part_data_table %>%
  filter(pre_post == 'paired') %>%
  left_join(., qPCR_egg_count, by = 'n_estudo') %>%
  filter(pre_post.y == 'Post-treatment')
nrow(pair_post_ct)
mean(pair_post_ct$ct_tt)
range(pair_post_ct$ct_tt)
```

## Now, the impact of the variants in the proteis is assessed
### Let's get first the vcf file

```bash
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
```

## And now we run SnpEff

```bash
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
```
### Now, Table S3 is generatet 

```
#Let's generate table S3 with the SNPeff info
#Let's add intron exon info for all SNPs, just in case
df_freqs <- select(freq_pre_post_all, pre_freq_mean, post_freq_mean)

df_freqs <- as_tibble(sapply(df_freqs, function(x) if(is.numeric(x)) 1-x else ""))
df_freqs <- as_tibble(sapply(df_freqs, function(x) if(is.numeric(x)) round(x, digits = 3) else ""))
df_freqs <- cbind(select(per_samplefrequency, POS, REF, ALT), df_freqs)
df_freqs$POS <- as.numeric(df_freqs$POS)

df_freqs$gen[df_freqs$POS <= 10684167] <- 'Intron'
df_freqs$gen[df_freqs$POS >= 10684168 & df_freqs$POS <= 10684696] <- "Exon"
df_freqs$gen[df_freqs$POS >= 10684697 & df_freqs$POS <= 10684751] <- "Intron"  
df_freqs$gen[df_freqs$POS >= 10684752 & df_freqs$POS <= 10684937] <- "Exon"  
df_freqs$gen[df_freqs$POS >= 10684938 & df_freqs$POS <= 10684987] <- "Intron"  
df_freqs$gen[df_freqs$POS >= 10684988 & df_freqs$POS <= 10685436] <- "Exon"  
df_freqs$gen[df_freqs$POS >= 10685437 & df_freqs$POS <= 10685483] <- "Intron"  
df_freqs$gen[df_freqs$POS >= 10685484 & df_freqs$POS <= 10685727] <- "Exon"  
df_freqs$gen[df_freqs$POS >= 10685728 & df_freqs$POS <= 10685780] <- "Intron"
df_freqs$gen[df_freqs$POS >= 10685781 & df_freqs$POS <= 10685999] <- "Exon" 
df_freqs$gen[df_freqs$POS >= 10686000 & df_freqs$POS <= 10686294] <- "Intron" 
df_freqs$gen[df_freqs$POS >= 10686295 & df_freqs$POS <= 10686647] <- "Exon"
df_freqs$gen[df_freqs$POS >= 10686648] <- 'Intron'

#Let's do a bed file with gernedalf variants
vcf_grenedalf <- read_table('grenedalf_snps_10.txt')
vcf_grenedalf$POS
vcf_grenedalf <- cbind(vcf_grenedalf, vcf_grenedalf$POS)
colnames(vcf_grenedalf) <- c('CHR', 'start', 'end')
write_delim(vcf_grenedalf, 'grenedalf_snps_10.bed', delim = '\tab')

vcf <- read.vcfR("grenedalf_snps_10_prot.ann.vcf", verbose = FALSE)
Rscript ./Parse_SnpEff.r grenedalf_snps_10.ann.trans.vcf snpeff.trans.csv 
snpeff <- read_csv("snpeff.trans.csv")

right_join(select(per_samplefrequency, POS, REF, ALT), vcf_grenedalf, by = c('POS' = 'start')) %>%
  select(POS, REF, ALT) %>%
  left_join(., select(df_freqs, POS, gen), by = 'POS') %>%
  left_join(., select(snpeff, POS, REF, ALT, effect, impact, gene, mut_cdna, mut_aa), by = 'POS') %>%
write.csv('Table_S3.csv') # then is sligly modify for publication
```

## Let's analyse the second beta tubulin
### First, we fin dthe location of the second isotype in the reference genome

```bash
#Let's see what is going on with the second beta tub
mkdir ${WORKING_DIR}/06_BETA_2
cd ${WORKING_DIR}/06_BETA_2

#Copying Steve nuclear vcf of trichuris (see https://github.com/stephenrdoyle/ancient_trichuris/)
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
```

### And now is R the the informaiton is visuallised

```R
# the location and gene data is loaded
exons <- read.table("beta_tub_2/beta_tub_2_exons.bed", header=T)
introns <- read.table("beta_tub_2/beta_tub_2_introns.bed", header=T)
nuc <- read.table("beta_tub_2/BZ_nuc_div.sites.pi", header=T)

resistant_snps <- read.table("beta_tub_2/res_snps_location.txt",header=T)

allele_freq <- read.table("beta_tub_2/BZ_allele_freq.frq", header=F, sep="\t", skip=1)
colnames(allele_freq) <- c("chrom", "POS", "alleles", "total_alleles", "REF", "VAR")
allele_freq[c('ref', 'ref_freq')] <- str_split_fixed(allele_freq$REF, ':', 2)
allele_freq[c('var', 'var_freq')] <- str_split_fixed(allele_freq$VAR, ':', 2)
allele_freq <- allele_freq[c("chrom", "POS", "alleles", "total_alleles", "ref", "ref_freq", "var", "var_freq")]
allele_freq$var_freq <- as.numeric(allele_freq$var_freq)

nuc_freq <- left_join(nuc, allele_freq, by = 'POS')

#Plotting
beta_tub_2_shceme <- ggplot() +
  geom_rect(data=introns, aes(xmin=start,ymin=-0.05, xmax=end, ymax=0.05), fill="grey50") +
  geom_rect(data=exons, aes(xmin=start,ymin=-0.15, xmax=end, ymax=0.15), fill="grey90") +
  ylim(-0.2,1) +
  labs(x="Genomic position (bp)", y="Nucleotide diversity (pi)", color="Variant frequency") +
  geom_segment(data=resistant_snps, aes(x=pos, xend=pos, y=-0.15, yend=0.15),col="red", size=1) +
  geom_text_repel(data=resistant_snps, aes(x=pos, y=0, label=name), col="red", box.padding = 0.5, max.overlaps = Inf, nudge_y = 0.75) +
  geom_segment(data=nuc_freq, aes(x=POS, xend=POS, y=0, yend=PI, col=var_freq), size=1) +
  geom_point(data=nuc_freq, aes(x=POS, y=PI, col=var_freq), size=3) +
  theme_bw() + scale_colour_viridis()  + scale_fill_viridis()


#Genome wide diversity
data <- read.table("beta_tub_2/trichuris_allsites_pi.txt", header=T)

chr <- filter(data,chromosome=="Trichuris_trichiura_2_001")
chr <- mutate(chr, colour = ifelse( window_pos_1>=27460001 & window_pos_1<27480001, "1", "0.5"))
chr <- chr %>%
  group_by(pop) %>%
  dplyr::mutate(position = 1:n())
#And plot
plot_1 <- ggplot(chr, aes(position*20000, avg_pi, col=colour, size=colour)) +
  geom_point() +
  facet_grid(pop~.) +
  geom_vline(xintercept=c(27478327,27481371), linetype="dashed", size=0.5) +
  labs(x = "Genomic position (bp)" , y = "Nucleotide diversity") +
  scale_colour_manual(values = c("grey", "brown1")) +
  scale_size_manual(values=c("0.5" = 0.5, "1" = 2)) +
  theme_bw() + theme(legend.position="none")
# rearrange Pi data to determine where btubulin sits in the distribution of pi
chr_sort <- arrange(chr,avg_pi)
chr_sort <-
  chr_sort %>%
  group_by(pop) %>%
  dplyr::mutate(position = 1:n())
chr_sort_quantile <-
  chr_sort %>%
  group_by(pop) %>%
  dplyr::summarise(enframe(quantile(avg_pi, c(0.05,0.95)), "quantile", "avg_pi"))
#And plot
plot_2 <- ggplot(chr_sort) +
  geom_hline(data=chr_sort_quantile,aes(yintercept=c(avg_pi)), linetype="dashed", size=0.5) +
  geom_point(aes(position, avg_pi, col=colour, size=colour)) + facet_grid(pop~.) + theme_bw() +
  scale_colour_manual(values = c("grey", "brown1")) +
  labs(x = "Position sorted by Pi" , y = "Nucleotide diversity") +
  scale_size_manual(values=c("0.5" = 0.5, "1" = 2)) +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

gw_beta_tub_2 <- plot_1 + plot_2
```





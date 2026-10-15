# Table of Contents
- [Data analyses](#data-analyses)
  - [List of tools/software](#list-of-toolssoftware)
  - [Adaptor Trimming](#adaptor-trimming)
  - [Splitting](#splitting)
  - [Genome download and indexing](#genome-download-and-indexing)
  - [Alignment and re-merging of Hi-C data](#alignment-and-re-merging-of-hi-c-data)
  - [Parsing](#parsing)
  - [Generating contact maps](generating-contact-maps)
- [Downstream analyses](#downstream-analyses)
  - [Recommended tools/software](#recommended-toolssoftware)
  - [Compartment calling](#compartment-calling)
  - [Loop calling](#loop-calling)
    - [Differential loop calling](#differential-loop-calling)
  - [Intersection with ChIP datasets](#intersection-with-chip-datasets)
 
# Data analyses
This GitHub page is a collection of scripts and general code used for the analysis of Hi-C datasets. Differential analyses were given with a comparison between CS1AN CSB+ and CSB- cell lines as an example.

## List of tools/software
Make sure these dependencies are installed in your system:
- anaconda3
- python (version 3.9)
- Parallel v20240322
- fastp v0.23.4
- Tools
- FastQC v0.11.9
- seqkit v2.8.2
- BWA-MEM 0.7.18
- SAMtools v1.22.1
- BEDTools v2.31.1
- Pairtools v1.1.0
- Java v21.0.8
- Juicer_tools v1.22.01

In anaconda3, make separate environments for the following Python suites used for the data analyses:
- deepTools v3.5.6
- MACS2 v2.2.9.1

For other analyses, install or clone Github repositories:
- R
- HOMER v5.1

For visualization:
- IGV v2.19.8

## Adaptor Trimming
Raw fastq files were trimmed to remove adaptor reads before alignment.
```bash
echo "FILENAME_1_R1 FILENAME_1_R2 FILENAME_2_R1 FILENAME_2_R2" > FILENAME.txt

#PBS -lwalltime=10:00:00
#PBS -lselect=1:ncpus=10:mem=128gb

module load parallel/20240322-GCCcore-13.2.0
eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)"

fastq_dir=/PATH/TO/TRIMMED-FILES
tmp_dir=$EPHEMERAL/tmp_trimming
mkdir -p $tmp_dir

sort -u <(sed 's/_R[12]$//' $PBS_O_WORKDIR/FILENAME.txt) | \
parallel -j 10 \
"/rds/general/user/nfa23/home/anaconda3/bin/fastp -i $fastq_dir/{}_R1.fq.gz -I $fastq_dir/{}_R2.fq.gz -o $tmp_dir/trimmed_{}_R1.fq.gz -O $tmp_dir/trimmed_{}_R2.fq.gz --detect_adapter_for_pe -l 20 -w 1 --compression 1 -j $tmp_dir/{}.fastp.json -h $tmp_dir/{}.fastp.html"

cp $tmp_dir/trimmed_*.fq.gz $tmp_dir/*.fastp.json $tmp_dir/*.fastp.html $PBS_O_WORKDIR
```

To check the quality of trimming, quality control was performed with FastQC
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=20:mem=128gb

module load tools/prod
module load FastQC/0.11.9-Java-11

fastq_dir=/PATH/TO/TRIMMED-FILES
tmp_dir=$EPHEMERAL/tmp_fastqc_post
mkdir -p $tmp_dir

fastqc -t 20 -d $TMPDIR -o $tmp_dir $fastq_dir/*.fq.gz

cp $tmp_dir/*.html $tmp_dir/*.zip "$PBS_O_WORKDIR"
```
## Splitting
Trimmed fastq files were split to allow parallellization of the alignment workflow.
```bash
#PBS -lwalltime=16:00:00
#PBS -lselect=1:ncpus=10:mem=128gb

eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)"
conda activate seqkit_env
module load parallel/20240322-GCCcore-13.2.0

fastq_dir=/PATH/TO/TRIMMED-FILES
tmp_dir=$EPHEMERAL/tmp_split
mkdir -p $tmp_dir

ls $fastq_dir/trimmed_*_R1.fq.gz | xargs -n1 basename | sed 's/^trimmed_//; s/_R1\.fq\.gz$//' | \
parallel -j 10 \
"mkdir -p $tmp_dir/{} && seqkit split2 -1 $fastq_dir/trimmed_{}_R1.fq.gz -2 $fastq_dir/trimmed_{}_R2.fq.gz -s 4000000 -O $tmp_dir/{} -f -e .gz -j 1"

cp -r $tmp_dir/* $PBS_O_WORKDIR
```
## Genome download and indexing
All files were aligned to the human reference genome hg38 (UCSC). Here, the link to the genome of reference is given.
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=1:mem=900gb

wget --no-verbose --directory-prefix hg38_genome "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
gunzip hg38.fa.gz

cp * $PBS_O_WORKDIR
```

For Hi-C datasets, NGS data was aligned to BWA-MEM. Below, the indexing of hg38 genome using BWA-MEM is shown. 
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=1:mem=900gb

eval "$(~/anaconda3/bin/conda shell.bash hook)"
conda activate bwa_env

cd $PBS_O_WORKDIR

bwa index -a bwtsw -p hg38 $PBS_O_WORKDIR/hg38.fa

cp * $PBS_O_WORKDIR
```

## Alignment and re-merging of Hi-C data
To align multiple files in a parallel manner, a <b>submission loop</b> was created on top of the usual PBS job. The submission loop can be started using the bash command on Terminal. Alignment submission loop can be submitted while splitting is still running, given the conditional command that the alignment could only start after all the splitting job has been completed. Upon alignment, the files were re-merged per biological replicate.

### Alignment Submission Loop
```bash
SPLIT_JOBID=xxxx.pbs-7

out_base=/PATH/TO/ALIGNED-FILES

for sample in C1_T1 C1_T2 C2_T1 C2_T2 W1_T1 W1_T2 W2_T1 W2_T2 W3_T1 W3_T2; do
  mkdir -p $out_base/$sample
  qsub -W depend=afterok:$SPLIT_JOBID -v SAMPLE=$sample Alignment_HiC.pbs
done
```

### Alignment
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=128:mem=128gb

module load parallel/20240322-GCCcore-13.2.0
eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)"
conda activate bwa_env

REF=/rds/general/project/diantonio/live/tools/genomes/hg38_ucsc_bwa/hg38
split_dir=/PATH/TO/SPLIT-FILES
tmp_dir=$EPHEMERAL/tmp_align/$SAMPLE
mkdir -p $tmp_dir

ls $split_dir/$SAMPLE/*R1.part*.fq.gz | sed 's/.*R1\.\(part_[0-9]*\)\.fq\.gz/\1/' | \
parallel -j 7 \
"bwa mem -SP5M -t 18 $REF $split_dir/$SAMPLE/trimmed_${SAMPLE}_R1.{}.fq.gz $split_dir/$SAMPLE/trimmed_${SAMPLE}_R2.{}.fq.gz > $tmp_dir/${SAMPLE}.{}.sam"

cp $tmp_dir/*.sam /PATH/TO/ALIGNED-FILES/$SAMPLE/
```
### Re-merging Submission Loop
```bash
for biorep in C1 C2 W1 W2 W3; do
  qsub -v BIOREP=$biorep MergingBioRep.pbs
done
```

### Re-merging
```bash
#PBS -lwalltime=04:00:00
#PBS -lselect=1:ncpus=4:mem=16gb

module load SAMtools/1.18-GCC-12.3.0

align_dir=/PATH/TO/ALIGNED-FILES
out_dir=/PATH/TO/MERGED-FILES/$BIOREP
mkdir -p $out_dir

samtools merge -O bam -@ 4 $out_dir/${BIOREP}.bam \
  $align_dir/${BIOREP}_T1/${BIOREP}_T1.merged.bam \
  $align_dir/${BIOREP}_T2/${BIOREP}_T2.merged.bam
```

## Parsing
Aligned, merged files were filtered, sorted, deduplicated, and parsed to determine valid contact pairs. This is one of the most time-consuming steps of the analysis, and as such the process were done only on compressed files. MAPQ value is the quality cutoff for filtering; higher value corresponds to more stringent filtering.

### Parsing Submission Loop
```bash
for biorep in C1 C2 W1 W2 W3; do
  qsub -v BIOREP=$biorep Parsing_HiC_Compressed.pbs
done
```
### Parsing
```bash
#PBS -lwalltime=48:00:00
#PBS -lselect=1:ncpus=8:mem=128gb

eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)"
conda activate hic_env

biorep_dir=/PATH/TO/MERGED-FILES/$BIOREP
final_out_dir=/PATH/TO/PARSED-FILES/$BIOREP
tmp_dir=$EPHEMERAL/tmp_parsing/$BIOREP
mkdir -p $tmp_dir $final_out_dir

chromsizes=/rds/general/project/diantonio/live/tools/genomes/hg38_ucsc_bwa/hg38.chr.genome
mapq=30
out=$tmp_dir/$BIOREP

samtools view -h -@ 8 $biorep_dir/${BIOREP}.bam | pairtools parse -c $chromsizes --min-mapq $mapq --assembly hg38 --add-columns mapq --nproc-out 8 -o ${out}.MAPQ30.parsed.pairsam.gz

pairtools sort --nproc 8 --nproc-out 8 -o ${out}.MAPQ30.sorted.pairsam.gz ${out}.MAPQ30.parsed.pairsam.gz
rm ${out}.MAPQ30.parsed.pairsam.gz

pairtools dedup --mark-dups --nproc-in 8 --nproc-out 8 -o ${out}.MAPQ30.deduped.pairsam.gz ${out}.MAPQ30.sorted.pairsam.gz --output-stats ${out}.MAPQ30.stats
rm ${out}.MAPQ30.sorted.pairsam.gz

pairtools select '(pair_type == "UU" or pair_type == "UR" or pair_type == "RU") and (mapq1>=30) and (mapq2>=30)' --nproc-in 8 --nproc-out 8 --output ${out}.MAPQ30.filtered.pairsam.gz ${out}.MAPQ30.deduped.pairsam.gz
rm ${out}.MAPQ30.deduped.pairsam.gz

pairtools split --nproc-in 8 --output-pairs ${out}.MAPQ30.nodups.pairs.gz ${out}.MAPQ30.filtered.pairsam.gz
rm ${out}.MAPQ30.filtered.pairsam.gz

zcat ${out}.MAPQ30.nodups.pairs.gz | grep -v '#' | perl -F'\t' -lane 'print join "\t", @F[0,1,2,5,3,4,6]' > ${out}.MAPQ30.Valid.pairs
gzip ${out}.MAPQ30.Valid.pairs

cp ${out}.MAPQ30.nodups.pairs.gz ${out}.MAPQ30.Valid.pairs.gz ${out}.MAPQ30.stats $final_out_dir/
```

## Generating contact maps
Contact maps can be generated using either Juicer tools (Aiden Lab) or Cooler tools. This depends on the format compatibility with the downstream analyses. Here, the code for contact map generated using Juicer is shown for resolution from (10000 to 10000000 bp)/

### Juicer-Cooler Submission Loop
```bash

parsing_dir=/PATH/TO/PARSED-FILES
condmerge_dir=/PATH/TO/PARSED-FILES-MERGED-PER-CELL-CONDITION
out_dir=/PATH/TO/CONTACT-MAP

for sample in C1 C2 W1 W2 W3; do
  qsub -v PAIRS=$parsing_dir/$sample/$sample.nodups.pairs.gz,OUT=$out_dir/$sample/$sample.hic Juicer_HiC.pbs
  qsub -v PAIRS=$parsing_dir/$sample/$sample.nodups.pairs.gz,OUT=$out_dir/$sample/$sample.mcool Cooler_HiC.pbs
done

for cond in CSB WT; do
  qsub -v PAIRS=$condmerge_dir/${cond}_combined.nodups.pairs.gz,OUT=$out_dir/CondMerge/${cond}_combined.hic Juicer_HiC.pbs
  qsub -v PAIRS=$condmerge_dir/${cond}_combined.nodups.pairs.gz,OUT=$out_dir/CondMerge/${cond}_combined.mcool Cooler_HiC.pbs
done
```

### Juicer
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=16:mem=64gb

module load Java/21.0.8

chromsizes=/rds/general/project/diantonio/live/tools/genomes/hg38_ucsc_bwa/hg38_ordered.chr.genome
mkdir -p $(dirname $OUT)

java -Xmx56G -jar /rds/general/user/nfa23/projects/diantonio/live/Naura/ExternalBioinformatics/hESC_BG4HiC/Parsing/juicer_tools_1.22.01.jar pre \
  -q 5 --threads 16 \
  -r 10000,20000,25000,50000,100000,250000,500000,1000000,10000000 \
  $PAIRS $OUT $chromsizes
```

### Cooler
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=16:mem=64gb

eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)" 
conda activate hic_env

pairs_dir=/rds/general/project/diantonio/live/Naura/HiC/CS1AN_WTvCSB_Novogene/Parsing/ParsingMAPQ5/CondMerge/
chromsizes=/rds/general/project/diantonio/live/tools/genomes/hg38_ucsc_bwa/hg38.chr.genome
out_dir=/PATH/TO/COOLER-MAPS
mkdir -p $out_dir

cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
  ${chromsizes}:10000 $PAIRS ${BASE}.10000.cool

cooler zoomify --nproc 16 --resolutions 10000,20000,50000,100000,250000,500000,1000000,10000000 \
  -o $OUT ${BASE}.10000.cool
```

# Downstream analyses
After generating the Hi-C contact map, several downstream analyses could be performed; compartment calling (CALDER2), loop calling and differential calling (FitHiChIP, FitHiChIP), intersections with ChIP-seq data (BEDTools), etc.

## Recommended tools/software
Depending on the type of analysis, different tools could be utilized. Several softwares/tools used to analyze the data of CSB+ vs CSB- cells ware shown below:
- CALDER2 v0.7
- DiffComp v0.0.3.post1.dev11+g2878167a9
- FitHiChIP v
- DiffHiChIP v
  
## Compartment calling
Compartment calling could be performed using CALDER2, whereas differential calling between two cell line conditions could then be performed using DiffComp (CSO Group, EPFL). This will show the distribution of A/B compartments, which denote active/inactive regions of the chromatin, ranked and clustered per-subcompartment. Differential analysis will allow the analysis of differences between two cell lines. The code below shows an example pipeline for the comparison between CSB+ (CSB) and CSB- (WT) cells.

### CALDER2 compartment calling
```bash
#PBS -lwalltime=08:00:00
#PBS -lselect=1:ncpus=8:mem=64gb

eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)"
conda activate calder_env

Rscript /rds/general/user/nfa23/projects/diantonio/live/Naura/HiC/CS1AN_WTvCSB_Novogene/CALDER_Compartment/Calder_HiC.R
```
```R
library(CALDER)

samples <- list(
  CSB = "/rds/general/project/diantonio/live/Naura/HiC/CS1AN_WTvCSB_Novogene/ContactMap/CondMerge/CSB_combined.hic",
  WT = "/rds/general/project/diantonio/live/Naura/HiC/CS1AN_WTvCSB_Novogene/ContactMap/CondMerge/WT_combined.hic"
)

for (s in names(samples)) {
  CALDER(
    contact_file_hic = samples[[s]],
    chrs = c(1:22, "X"),
    bin_size = 50000,
    genome = "hg38",
    save_dir = paste0("/rds/general/project/diantonio/live/Naura/HiC/CS1AN_WTvCSB_Novogene/CALDER_Compartment/", s),
    single_binsize_only = TRUE,
    save_intermediate_data = FALSE,
    n_cores = 8,
    sub_domains = FALSE
  )
}
```

### DiffComp for differential analysis
```bash
#PBS -lwalltime=04:00:00
#PBS -lselect=1:ncpus=4:mem=16gb

eval "$(/rds/general/user/nfa23/home/anaconda3/bin/conda shell.bash hook)"
conda activate diffComp

out_dir=/rds/general/project/diantonio/live/Naura/HiC/CS1AN_WTvCSB_Novogene/CALDER_Compartment/diffComp
mkdir -p $out_dir

calder_dir=/rds/general/project/diantonio/live/Naura/HiC/CS1AN_WTvCSB_Novogene/CALDER_Compartment
calder_ind_dir=/rds/general/project/diantonio/live/Naura/HiC/CS1AN_WTvCSB_Novogene/CALDER_Compartment/CALDER_Individual

cores \
  $calder_dir/CSB/sub_compartments/CSB_all_sub_compartments.bed \
  $calder_dir/WT/sub_compartments/WT_all_sub_compartments.bed \
  50000 \
  $out_dir/CSBvWT_output \
  --control1_path $calder_ind_dir/C1/sub_compartments/all_sub_compartments.bed,$calder_ind_dir/C2/sub_compartments/all_sub_compartments.bed \
  --control2_path $calder_ind_dir/W2/sub_compartments/all_sub_compartments.bed,$calder_ind_dir/W3/sub_compartments/all_sub_compartments.bed \
  --genome hg38 \
  --bed_path $out_dir/CSBvWT_CoREs.bed \
  --signal_path $out_dir/CSBvWT_signal \
  --verbose
```

## Loop calling
Loop callings can be performed with several softwares (e.g. FitHiChIP, Chromosight), and for several modes (e.g. ALL2ALL, P2P). Below, peak calling using FitHiChIP, and subsequent differential analysis using DiffHiChIP, are shown.

Different configuration files for various modes of loop calling can be found in the Configuration File folder. ALL2ALL refers to all detected loop interactions within the 3D contact map, P2P refers to detected loop interactions with ChIP peaks designated as anchors. While not show, P2ALL mode could also be performed, where the ChIP peaks are designated as a hub for looping interactions that are not necessarily between two peaks.

```bash
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=8:mem=64gb

cd $PBS_O_WORKDIR

module load SAMtools/1.18-GCC-12.3.0
module load BEDTools/2.31.0-GCC-12.3.0

eval "$(~/anaconda3/bin/conda shell.bash hook)"
conda activate FitHiChIP_env

export PATH=/rds/general/user/nfa23/home/anaconda3/envs/HiC-Pro_v3.1.0/bin:$PATH
export PATH=/rds/general/user/nfa23/home/InsHiC-Pro_v3.1.0/HiC-Pro_3.1.0/bin:$PATH
export PATH=/rds/general/user/nfa23/home/anaconda3/envs/FitHiChIP_env/bin:$PATH

FITHICHIP=/rds/general/user/nfa23/home/FitHiChIP/FitHiChIP_HiCPro.sh

mkdir -p $OUTDIR

echo "Starting FitHiChIP $CONDITION: $(date)"
bash $FITHICHIP -C $CONFIG
echo "Finished FitHiChIP $CONDITION: $(date)"
```

### Differential loop calling
```bash
#!/bin/bash
#PBS -lwalltime=24:00:00
#PBS -lselect=1:ncpus=8:mem=64gb

cd /rds/general/user/nfa23/home/DiffHiChIP/src
eval "$(~/anaconda3/bin/conda shell.bash hook)"
conda activate FitHiChIP_env

Rscript /rds/general/user/nfa23/home/DiffHiChIP/src/DiffLoop_Code_FitHiChIP.R \
  --InpTable=/rds/general/project/diantonio/live/Naura/HiC/CS1AN_WTvCSB_Novogene/FitHiChIP/DiffHiChIP_config.txt \
  --OutDir=/rds/general/project/diantonio/live/Naura/HiC/CS1AN_WTvCSB_Novogene/FitHiChIP/DiffHiChIP_CSBvWT_All2All \
  --ChrSizeFile=/rds/general/project/diantonio/live/tools/genomes/hg38_ucsc_bwa/hg38.chr.genome \
  --FDRThr=0.01 \
  --Model=1 \
  --DistStrat=1 \
  --DiffFDRThr=0.05 \
  --LFCThr=1
```
## Intersection with ChIP datasets

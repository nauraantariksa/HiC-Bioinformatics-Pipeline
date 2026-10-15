# Table of Contents
- [Data analyses](#data-analyses)
  - [List of tools/software](#list-of-toolssoftware)
  - [Adaptor Trimming](#adaptor-trimming)
  - [Splitting](#splitting)
  - [Genome download and indexing](#genome-download-and-indexing)
  - [Alignment and re-merging of Hi-C data](#alignment-and-re-merging-of-hi-c-data)
  - [Parsing](#parsing)
 
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
- Picard v2.25.1
- BEDTools v2.31.1

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
SPLIT_JOBID=xxxx.pbs-7

out_base=/PATH/TO/ALIGNED-FILES

for sample in C1_T1 C1_T2 C2_T1 C2_T2 W1_T1 W1_T2 W2_T1 W2_T2 W3_T1 W3_T2; do
  mkdir -p $out_base/$sample
  qsub -W depend=afterok:$SPLIT_JOBID -v SAMPLE=$sample Alignment_HiC.pbs
done

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

biorep_dir=/rds/general/project/diantonio/live/Naura/HiC/CS1AN_WTvCSB_Novogene/Merging/$BIOREP
final_out_dir=/rds/general/project/diantonio/live/Naura/HiC/CS1AN_WTvCSB_Novogene/ParsingMAPQ30/$BIOREP
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

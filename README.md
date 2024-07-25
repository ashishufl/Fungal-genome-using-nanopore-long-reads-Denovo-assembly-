# Fungal-genome-using-nanopore-long-reads-Denovo-assembly-
Minion long reads denovo assembly polish using illumina reads

This is a pipeline to do genome assemble, starting from basecalling of raw reads/signals, to processing and filtering of reads and assembly of a fungal genome. The assembly was aslo polished with illumina short reads for gap filling.

This was done in cloud computer using bash scripts.

cat *.fastq > Combined_fastq.fastq


#!/bin/sh
#SBATCH --account=plantpath
#SBATCH --qos=plantpath-b
#SBATCH --job-name=110407_3_1_1_Race_0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pcvgt@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=85gb
#SBATCH --time=96:00:00

pwd; hostname; date

module load gcc/5.2.0
module load porechop

porechop -i /ufrc/loria/share/pcvgt/110407_3_1_1_Race_0/20200803_1417_MN33357_FAM92903_9d1ad5ea/fast5/completefastq/TotalComplete.fastq -o PoreChopout


#!/bin/sh
#SBATCH --account=plantpath
#SBATCH --qos=plantpath-b
#SBATCH --job-name=110407_3_1_1_Race_0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pcvgt@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=85gb
#SBATCH --time=96:00:00

pwd; hostname; date

module load gcc/5.2.0
module load porechop

porechop -i /ufrc/loria/share/pcvgt/110407_3_1_1_Race_0/20200803_1417_MN33357_FAM92903_9d1ad5ea/fast5/completefastq/TotalComplete.fastq -o PoreChopout
From James Fulton to Everyone:  09:36 AM
#!/bin/sh
#SBATCH --account=plantpath
#SBATCH --qos=plantpath-b
#SBATCH --job-name=110407_3_1_1_Race_0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pcvgt@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=70gb
#SBATCH --time=72:00:00
#SBATCH --output=110407_3_1_1_Race_0_nanoplot_%j.out
pwd; hostname; date

module load nanoplot
NanoPlot --fastq /ufrc/loria/share/pcvgt/110407_3_1_1_Race_0/20200803_1417_MN33357_FAM92903_9d1ad5ea/fast5/completefastq/PoreChop/110407_3_1_1_RACE_0_filter_HAC_mod_90quality_10kblength.fastq -o NanoplotOUT --verbose


#!/bin/sh

#SBATCH --account=plantpath

#SBATCH --qos=plantpath-b

#SBATCH --job-name=150524_Race1

#SBATCH --mail-type=END,FAIL

#SBATCH --mail-user=pcvgt@ufl.edu

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=4

#SBATCH --mem=10gb

#SBATCH --time=72:00:00

#SBATCH --output=150524_Race1_filtlong_%j.out

pwd; hostname; date

 

module load  filtlong

filtlong --min_mean_q 90 --min_length 10000 150524_RACE_1_PoreChopout.fastq > 150524_RACE_1_filter_HAC_mod_90quality_10kblength.fastq

#!/bin/sh
#SBATCH --account=plantpath
#SBATCH --qos=plantpath-b
#SBATCH --job-name=150524_RACE_1_smartdenovo
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pcvgt@ufl.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=80gb
#SBATCH --time=96:00:00
#SBATCH --output=smartdenovo_%j.out
pwd; hostname; date

module load smartdenovo


smartdenovo.pl -p 150524_RACE_1 -c 1 /ufrc/loria/share/pcvgt/150524_Race1/20200729_1632_MN33357_FAM97078_3f67d641/fast5/Complete/PoreChop/150524_RACE_1_filter_HAC_mod_90quality_10kblength.fastq > 150524_RACE_1.mak
make -f 150524_RACE_1.mak




http://cab.cc.spbu.ru/quast/


Quest..
#!/bin/sh
#SBATCH --account=plantpath
#SBATCH --qos=plantpath-b
#SBATCH --job-name=Quast
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=adhikariashish@ufl.edu
#SBATCH --ntasks=1
#SBATCH --mem=4g
#SBATCH --time=4:00:00
#SBATCH --output=Quast_%j.out
pwd; hostname; date


module load quast
 
quast.py -o Quast_racon_2k -l "quast_output" -t 2 /blue/goss/adhikariashish/2Kb_data_assembly/racon_2kb/racon_2kb.fasta

# minimap
minimap2 -t 8 -ax map-ont --secondary=no /blue/... trimmedont-reads.fq > aln.sam      # for Oxford Nanopore reads




#synteny mapping
what reference to use?


****
fastx_trimmer
gunzip and trimmed by 
-m 140 


#!/bin/sh

#SBATCH --account=plantpath

#SBATCH --qos=plantpath-b

#SBATCH --job-name=110407_3_1_1_Race_0

#SBATCH --mail-type=END,FAIL

#SBATCH --mail-user=pcvgt@ufl.edu

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --mem=45gb

#SBATCH --time=96:00:00

#SBATCH --output=110407_3_1_1_Race_0_pilon_%j.out

pwd; hostname; date

 

 

module load bwa/0.7.17

module load samtools

module load pilon

 

bwa index /ufrc/loria/share/pcvgt/110407_3_1_1_Race_0/20200803_1417_MN33357_FAM92903_9d1ad5ea/pilon/smartdenovo/round1/pilon_round1.fasta

 

 

bwa mem -t 14 /ufrc/loria/share/pcvgt/110407_3_1_1_Race_0/20200803_1417_MN33357_FAM92903_9d1ad5ea/pilon/smartdenovo/round1/pilon_round1.fasta /ufrc/loria/share/pcvgt/Illumina_Data_2Batch/NS1977-J8HMN-Pool/110407_3_1_1_Race_0_7-1_5-1_L001_ds.9b27c20a1a6d4e2a9cbd90b323aea422/fastx_trimmer/repeat_run/FASTxOUT_R1.fastq /ufrc/loria/share/pcvgt/Illumina_Data_2Batch/NS1977-J8HMN-Pool/110407_3_1_1_Race_0_7-1_5-1_L001_ds.9b27c20a1a6d4e2a9cbd90b323aea422/fastx_trimmer/repeat_run/FASTxOUT_R2.fastq | samtools view - -Sb | samtools sort - -@14 -o /ufrc/loria/share/pcvgt/110407_3_1_1_Race_0/20200803_1417_MN33357_FAM92903_9d1ad5ea/pilon/smartdenovo/illumina_mapping/mapping_pilon2.sorted.bam

 

samtools index /ufrc/loria/share/pcvgt/110407_3_1_1_Race_0/20200803_1417_MN33357_FAM92903_9d1ad5ea/pilon/smartdenovo/illumina_mapping/mapping_pilon2.sorted.bam

 

 

export _JAVA_OPTIONS="-Xmx45g"

 

pilon --genome /ufrc/loria/share/pcvgt/110407_3_1_1_Race_0/20200803_1417_MN33357_FAM92903_9d1ad5ea/pilon/smartdenovo/round1/pilon_round1.fasta --fix all --changes --frags /ufrc/loria/share/pcvgt/110407_3_1_1_Race_0/20200803_1417_MN33357_FAM92903_9d1ad5ea/pilon/smartdenovo/illumina_mapping/mapping_pilon2.sorted.bam --threads 8 --output /ufrc/loria/share/pcvgt/110407_3_1_1_Race_0/20200803_1417_MN33357_FAM92903_9d1ad5ea/pilon/smartdenovo/round2/pilon_round2 | tee /ufrc/loria/share/pcvgt/110407_3_1_1_Race_0/20200803_1417_MN33357_FAM92903_9d1ad5ea/pilon/smartdenovo/round2/pilon_round2




#!/bin/sh

#SBATCH --account=plantpath

#SBATCH --qos=plantpath-b

#SBATCH --job-name=RepeatMod

#SBATCH --mail-type=END,FAIL

#SBATCH --mail-user=

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=5

#SBATCH --mem=15gb

#SBATCH --time=96:00:00

#SBATCH --output=RepeatMod_%j.out

 

module load repeatmodeler

REF=location_of_genome

REF_BASE="desired_name"

BuildDatabase -name $REF_BASE $REF

RepeatModeler -pa 5 -database $REF_BASE

 

 

#!/bin/sh

#SBATCH --account=plantpath

#SBATCH --qos=plantpath-b

#SBATCH --job-name=RepeatMasker

#SBATCH --mail-type=END,FAIL

#SBATCH --mail-user=

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=8

#SBATCH --mem=25gb

#SBATCH --time=96:00:00

#SBATCH --output=RepeatMasker_%j.out

 

pwd; hostname; date

 

module load repeatmasker/4.0.5

 

 

RepeatMasker -gff -lib

 

name/consensi.fa.classified genome_file

https://www.animalgenome.org/bioinfo/resources/manuals/RepeatMasker.html




""maker""
#after running maker for first round; run these codes

ml snap maker

gff3_merge -d *index.log

maker2zff *all.gff

fathom genome.ann genome.dna -validate > snap_validate_output.txt

fathom genome.ann genome.dna -categorize 1000

fathom uni.ann uni.dna -export 1000 -plus

forge export.ann export.dna

hmm-assembler.pl /orange/goss/adhikariashish/MIGS/Pilon2/Maker/Maker_again/pilon_round2.fasta.masked . > pilon_round2.hmm


# again edit maker.opt file 
line 34, 41 and 42

then add ml snap to maker script and run it.


*maker augustus
ml maker
gff3_merge -d *index.log
maker2zff *all.gff
fathom genome.ann genome.dna -validate > snap_validate_output.txt

fathom genome.ann genome.dna -categorize 1000

fathom uni.ann uni.dna -export 1000 -plus

nano zff2augustus_gbk.pl
https://github.com/hyphaltip/genome-scripts/blob/master/gene_prediction/zff2augustus_gbk.pl (copy paste these codes)
ml perl
perl zffaugustus_gbk.pl > augustus.gbk

ml augustus 
randomSplit.pl augustus.gbk 100
new_species.pl --species=bipolaris


nano augustus.sh
ml augustus
export config file
new_species.pl --species=bipolaris

nano augustus_train.sh
....
ml maker augustus
export Augustus_cONFIG_PATH=/home/pcvgt/augutus_conf_dir/config

etraining --species=.bipolaris augustus.gbk.train
augustus --species=bipolaris augustus.gbk.test | tee first_training.out




nano opyimize_augustus.sh
....
increase memory and time (up to 7 days)

ml maker augustus
export Augustus_cONFIG_PATH=/home/pcvgt/augutus_conf_dir/config

optimize_augustus.pl --species=bipolaris augustus.gbk.train


*Second training after optimizing


*Then run maker again
change opts file: change augustus species=bi_speices
		   keep preds=1


*last run with genemark
#make a script for genemark
for genemark, we use the assembly file


rm -f........
ml genemark_es
gmes_petap.pl -ES -fungus -cores 10 -sequence /orange/goss/adhikariashish/assemby.fasta








*circos ordered pic
CM017956.1_RagTag,CM017956.1,CM017957.1_RagTag,CM017957.1,CM017958.1_RagTag,CM017958.1,CM017959.1_RagTag,CM017959.1,CM017960.1_RagTag,CM017960.1,CM017961.1_RagTag,CM017961.1,CM017962.1_RagTag,CM017962.1,CM017963.1_RagTag,CM017963.1,CM017964.1_RagTag,CM017964.1,CM017966.1_RagTag,CM017966.1,CM017967.1_RagTag,CM017967.1,CM017968.1_RagTag,CM017968.1,CM017969.1_RagTag,CM017969.1,CM017971.1_RagTag,CM017971.1,SRZH01000017.1_RagTag,SRZH01000017.1,utg2847_pilon_pilon,utg538_pilon_pilon,utg65_pilon_pilon,CM017970.1,SRZH01000018.1,SRZH01000019.1,SRZH01000020.1,SRZH01000021.1,SRZH01000022.1



CM017956.1_RagTag,CM017956.1,CM017957.1_RagTag,CM017957.1,CM017958.1_RagTag,CM017958.1,CM017959.1_RagTag,CM017959.1,CM017960.1_RagTag,CM017960.1,CM017961.1_RagTag,CM017961.1,CM017962.1_RagTag,CM017962.1,CM017963.1_RagTag,CM017963.1,CM017964.1_RagTag,CM017964.1,CM017966.1_RagTag,CM017965.1,CM017967.1_RagTag,CM017966.1,CM017968.1_RagTag,CM017967.1,CM017969.1_RagTag,CM017968.1,CM017971.1_RagTag,CM017969.1,SRZH01000017.1_RagTag,CM017970.1,utg2847_pilon_pilon,CM017971.1,utg538_pilon_pilon,SRZH01000017.1,utg65_pilon_pilon,SRZH01000018.1,SRZH01000019.1,SRZH01000020.1,SRZH01000021.1,SRZH01000022.1





**Blast Script
ml ncbi_blast
makeblastdb –in your_sequences_db –dbtype nucl or prot

blastn -query your_sequences –db your_database –out your_results –outfmt 1 (full format),  6 (tab format) 

*if you want to blast proteins change blastn for blastp or use the other programs blastx, tblastn , tblastx, depending of your query and/or db



##### SNP filter website https://speciationgenomics.github.io/filtering_vcfs/




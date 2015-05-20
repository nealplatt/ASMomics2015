# ASMomics 2015
# Mitogenome assembly in vespertilionid and phylostomid bats
# RNPlatt
# 20 May 2015
# v 0.1


#tmp variables
FASTX_DIR=/lustre/work/apps/fastx_toolkit-0.0.14/bin/


#The first thing we nee to do is set up a directory structure.  There is no 
#  hard and fast rule, but I try to stick with the Noble (2009) structure.
#  http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424


#Here we are setting variable names for major directories.  This will streamline
#  downstreal processing and reduce the chance for errors.  It is easier to get 
#  it right once rather than a few dozen times.  
HOME_DIR="/lustre/scratch/roplatt/asm/mitogenomes"
DATA_DIR="$HOME_DIR/data"
RESULTS_DIR="$HOME_DIR/results"
BIN_DIR="$HOME_DIR/bin"

mkdir $HOME_DIR $DATA_DIR $RESULTS_DIR $BIN_DIR


#Now we will move into the data directory to begin processing
cd $DATA_DIR

#First we need to download the data from an FTP site.
wget <address>..batUCE_batch1_2015-01.tgz

#To save space and increase transfer speeds, the sequence data is compressed and
#  archived.  Expand the data into a useable format.
tar -xvzf batUCE_batch1_2015-01.tgz

#Lets take an opportunity to look at the data...
#First, lets look at the all the files, their sizes, and permissions
ls -l *

#The file permission associated with these files give the owner (you) permission
#  to write onto/over them.  Whenever you are working from RAW data, it is best
#  to preserve them in their original state.  Remove the write permissions.
chmod a-w *.fq

#Its always nice to know the size of your files, but the previous command 
#  outputs this data in bytes.  It is more intuitive to look at this in a "human"
#  readable number.
ls -lh *

#You should notice that your file permission have changed, the "w" is gone, and
#  file sizes make more sense (G = gigabytes, M = megabytes, K = kilobytes).

#Now lets look at the actual data...
head braCav_17_R1.fq

# description of fastq slides here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#We know our files are large, and they are in the FastQ format, but how big are
#  they
wc *.fq
#This gives us the TOTAL number of lines, words, characters, and the filename.
#  but we know that the .FQ format has 4 lines per entry.  Lets get a better idea

wc -l *.fq | awk '{print $1/4"\t"$2}'       #<- first use of the pipe , explain !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#or by millions
wc -l *.fq | awk '{print $1/4000000"\t"$2}' 

#Also, notice that the number of reads are identical between libraries
#  R1 vs. R2.  This expected from paired end data.

################################################################################
################################################################################
##                                                                            ##
##         Quality control via FastQC and FastX                               ##
##                                                                            ##
################################################################################
################################################################################
#
#Here we will use a series of programs to investigate the quality of the raw 
#  reads, filter dirty data, then see how that changed the overall quality of
#  the remaining (filtered) reads. To do this we will use two different packages
#  FastX_toolkit and FastQC.

# *** REMEMBER - we will be manipulating and creating new data.  This will be done
#       in the results/seqQC directory ***

#Pick a sample that you will use for processing.  Designate the R1 and R2 reads
#  here.  Do not include a file extension (ex. .fq)
BAT_R1=braCav_17_R1
BAT_R2=braCav_17_R2

#create the directory and set a variable
SEQ_DIR="$RESULTS_DIR/seqQC"

mkdir $SEQ_DIR
cd $SEQ_DIR

fastqc                          \
        --outdir .              \
        --threads 10            \
        $DATA_DIR/$BAT_R1.fq    \
        $DATA_DIR/$BAT_R2.fq
#Open the FastQC .html output on your local computer.  This can be done by 
#  transfering the file through FileZilla or copying it to DropBox 
cp $BAT_R1*.fastq.html ~/Dropbox
cp $BAT_R2*.fastq.html ~/Dropbox

#discuss fastqc output and the meaning behind each of the figures. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Even though the data looks to be high quality, it is good practice to filter
#  out low quality reads and trim adapter seuqences.  This will be done in a
#  series of steps using the FastX toolkit.  There is more than one way to
#  skin a cat...this is my preferred way.


#Here we clip adapters and remove reads with Q20 over half the read...for R1
ADAPTER_SEQ=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

$FASTX_DIR/fastx_clipper                                \
        -n                                              \
        -a $ADAPTER_SEQ                                 \
        -i $DATA_DIR/$BAT_R1.fq                         \
        | $FASTX_DIR/fastq_quality_filter               \
                -Q33                                    \
                -q 20                                   \
                -p 50                                   \
                -o $SEQ_DIR/$BAT_R1"_filtered.fq"  

#Here we clip adapters and remove reads with Q20 over half the read...for R1
$FASTX_DIR/fastx_clipper                                \
        -n                                              \
        -a $ADAPTER_SEQ                                 \
        -i $DATA_DIR/$BAT_R2.fq                         \
        | $FASTX_DIR/fastq_quality_filter               \
                -Q33                                    \
                -q 20                                   \
                -p 50                                   \
                -o $SEQ_DIR/$BAT_R2"_filtered.fq"  




#Adapters TruSeq Universal Adapter:
#  i7:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
#  i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT*GTGTAGATCTCGGTGGTCGCCGTATCATT







#assemble mito genomes using de novo
NAME[0]=Index_Bat_9
NAME[1]=Index_Bat_10
NAME[2]=Index_Bat_11
NAME[3]=Index_Bat_12
NAME[4]=Index_Bat_13
NAME[5]=Index_Bat_14
NAME[6]=Index_Bat_15
NAME[7]=Index_Bat_16
NAME[8]=Index_Bat_18
NAME[9]=Index_Bat_19
NAME[10]=Index_Bat_21
NAME[11]=Index_Bat_28

for (( i = 0; i < ${#NAME[@]}; i++))
do
    /lustre/work/apps/trinityrnaseq_r20140717/Trinity \
        --seqType fa \
        --JM 30G \
        --left  /lustre/scratch/roplatt/asm/mitoGenomes/data/${NAME[$i]}_1.fq \
        --right /lustre/scratch/roplatt/asm/mitoGenomes/data/${NAME[$i]}_2.fq \
        --CPU 20 \
        --output "/lustre/scratch/roplatt/asm/mitoGenomes/results/${NAME[$i]}.trinity" \
        --full_cleanup 
done



NAME[0]=Index_Bat_9
NAME[1]=Index_Bat_21
NAME[2]=Index_Bat_28
NAME[3]=Index_Bat_12
NAME[4]=Index_Bat_13
NAME[5]=Index_Bat_14
NAME[6]=Index_Bat_15
NAME[7]=Index_Bat_16
NAME[8]=Index_Bat_18
NAME[9]=Index_Bat_19
NAME[10]=Index_Bat_10
NAME[11]=Index_Bat_11

for (( i = 0; i < ${#NAME[@]}; i++))
do
/lustre/work/apps/trinityrnaseq_r20140717/Trinity \
        --seqType fq \
        --JM 20G \
        --left "/lustre/scratch/roplatt/asm/mitoGenomes/data/"${NAME[$i]}"_1.fq" \
        --right "/lustre/scratch/roplatt/asm/mitoGenomes/data/"${NAME[$i]}"_2.fq" \
        --CPU 20 \
        --output ${NAME[$i]}
done



cat bat10_1.fa bat10_2.fa >bat10_X.fa

 /lustre/work/apps/blast/bin/makeblastdb -in bat10_mtGenome.fa -dbtype nucl


 /lustre/work/apps/blast/bin/blastn -query bat10_X.fa -db bat10_mtGenome.fa -out bat10_blastn.out -outfmt 6 -evalue 1e-3


 time /lustre/work/apps/trinityrnaseq_r20140717/Trinity \
	--seqType fq \
	--JM 20G \
	--left /lustre/scratch/roplatt/asm/mitoGenomes/data/Index_Bat_11_1.fq \
	--right /lustre/scratch/roplatt/asm/mitoGenomes/data/Index_Bat_11_2.fq \
	--CPU 20 \
	--output bat11


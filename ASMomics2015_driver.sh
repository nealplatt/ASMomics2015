# ASMomics 2015
# Mitogenome assembly in vespertilionid and phylostomid bats
# RNPlatt
# 20 May 2015
# v 0.1


#tmp variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
##         Quality control via FastQC and Trimmomatic                         ##
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
BAT=braCav_17
BAT_R1=$BAT"_R1"
BAT_R2=$BAT"_R2"

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

#First lets get the Trimmomatic binary (ths is sent to your bin directory
wget \
        --directory-prefix=$BIN_DIR \
        http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip

#uncompress the file
unzip $BIN_DIR/Trimmomatic-0.33.zip

#and set a variable for the path
TRIM_BIN=$BIN_DIR/Trimmomatic-0.33


#First get the adapter sequeancs into a file for Trimmomatic file
echo ">i7">TruSeq4-PE.fa
echo "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG">>TruSeq4-PE.fa
echo ">i5">>TruSeq4-PE.fa
echo "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT" >>TruSeq4-PE.fa

#Now we can run Trimmomatic.  This is a fairly complicated filtering protocol
#  and I would suggest reading the manual to have a better understanding.
#  http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
#  In short, we are removing poritons of reads that match to the sequencing
#  adapter, (2) removing reads below Q20 at the beginning and end of the reads
#  (3) clipping reads once the quality score drops below Q20 in any 4bp window
#  and (4) removing any read that is less than 33 bp after these filters have
#  been applied.  More importantly, Trimmomatic is "pair aware" and when a read
#  from R1 is culled (placed in the UNPaired file) the R2 is as well. [10 min]

java -jar $TRIM_BIN/trimmomatic-0.33.jar \
        PE \
        -threads 10 \
        -phred64 \
        $DATA_DIR/$BAT_R1.fq \
        $DATA_DIR/$BAT_R2.fq \
        "$BAT_R1"_filterPaired.fq \
        "$BAT_R1"_filterUNPaired.fq \
        "$BAT_R2"_filterPaired.fq \
        "$BAT_R2"_filterUNPaired.fq \
        ILLUMINACLIP:TruSeq4-PE.fa:2:30:10 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:4:20 \
        MINLEN:33
       
#After filtering, this is a good time to look at the data.  Use "less" to scroll
#  through the data.  Notice that read lengths and file sizes vary, but the 
#  number of sequences per group matches

#Take this opportunity to re-run FastQC and look at the differences in the .html
#  output.

fastqc \
        --outdir . \
        --threads 10 \
        "$BAT_R1"_filterPaired.fq \
        "$BAT_R1"_filterUNPaired.fq \
        "$BAT_R2"_filterPaired.fq \
        "$BAT_R2"_filterUNPaired.fq

#Is there anything we can do to increase the quality of the sequences?  Look
#  through the Trimmomatic manual and and FastQC results.  Can we re-run
#  Trimmomatic to improve our overall quality?

#
java -jar $TRIM_BIN/trimmomatic-0.33.jar \
        PE \
        -threads 10 \
        -phred64 \
        $DATA_DIR/$BAT_R1.fq \
        $DATA_DIR/$BAT_R2.fq \
        "$BAT_R1"_filterPaired.fq \
        "$BAT_R1"_filterUNPaired.fq \
        "$BAT_R2"_filterPaired.fq \
        "$BAT_R2"_filterUNPaired.fq \
        ILLUMINACLIP:TruSeq4-PE.fa:2:30:10 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:4:20 \
        MINLEN:33 \
        HEADCROP:6

#Lets go ahead and combine our unpaired reads into a singletons file
cat "$BAT_R1"_filterUNPaired.fq "$BAT_R2"_filterUNPaired.fq >$BAT"_RX_UNpaired.fq"

fastqc \
        --outdir . \
        --threads 10 \
        "$BAT_R1"_filterPaired.fq \
        "$BAT_R1"_filterUNPaired.fq \
        "$BAT_R2"_filterPaired.fq \
        "$BAT_R2"_filterUNPaired.fq


################################################################################
################################################################################
####                                                                        ####
####                De Novo Asembly w/ Trinity                              ####
####                                                                        ####
################################################################################
################################################################################
#Now for the heavy lifting.  We are going to use our quality filtered data along
#  with a Trinity (http://trinityrnaseq.github.io/) to assemble our data into 
#  contigs.  Since we are not using a refernece genome, these are "de novo" 
#  assemblies.  There is no magic bullet when it comes to assembly.  While the 
#  parameters used here are "good enough" a significant amount optimization could
#  occur here.  

#Lets organize our assemblies into a seperate directory (in Results)
ASSEMBLY_DIR="$RESULTS_DIR/assemblies"
mkdir $ASSEMBLY_DIR
cd $ASSEMBLY_DIR

#start the De Novo Assembly
 /lustre/work/apps/trinityrnaseq_r20140717/Trinity \
        --seqType fq \
        --JM 15G \
        --left  $SEQ_DIR/$BAT_R1"_filterPaired.fq" \
        --right $SEQ_DIR/$BAT_R2"_filterPaired.fq" \
        --single $SEQ_DIR/$BAT"_RX_UNpaired.fq" \
        --CPU 10 \
        --output "$ASSEMBLY_DIR/$BAT" \
        --full_cleanup 

#on average these assemblies take ~an hour to complete.

#lets use some preliminary tools to investigate the assembly.  For me, it is
#  easier to analyze the data in

#there are several ways to find potential mitochondrial sequences in your data.
#  Method 1a & 1b: Use NCBI Blast
#  Method 2: search based on length (prior knowledge).


#Method 1a: Blasting against NCBI mitogenomes.  It is possible to take your
#  asesemblies and upload them to NCBI for the BLAST search.  This is not a 
#  feasible option.  Try it, but expect it to fail

#Method 1b: Download mitochondiral genomes from other bat species, then run a
#  local blast on your EC2 instance.  

#Go to NCBI and download several of the bat mitochondrial genomes.  Upload
#  them to your EC2 instance for a local blast.

ANNOTATION_DIR="$RESULTS_DIR/annotations"
mkdir $ANNOTATION_DIR
cd $ANNOTATION_DIR

#I named my mitochondiral genome file <chiropteraMitoGenomes_2015-06-05.fasta>.
#  Since this is a file that will NOT be altered and was generated by others, I
#  store it in the "data" directory.

#The first thing you need to do is to create a blast index.

MITOGENOMES=$DATA_DIR/chiropteraMitoGenomes_2015-06-05.fasta
makeblastdb -in $MITOGENOMES -dbtype nucl

#Now run your blast search (as an example I am using the Eptesicus furnalis
#  assembly).
/lustre/work/apps/blast/bin/blastn \
    -db $MITOGENOMES \
    -query ../assemblies/eptFur_21.Trinity.fasta \
    -out eptFur_21_vsChiropterMitoGenomes.blastn.out \
    -outfmt 6 \
    -num_threads 8
 

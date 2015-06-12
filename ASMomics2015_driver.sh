# ASMomics 2015
# Mitogenome assembly in vespertilionid and phylostomid bats
# RNPlatt
# 12 June 2015
# v 1.0

################################################################################
################################################################################
##                                                                            ##
##                      PUTTING OUR HOUSE IN ORDER                            ##
##                                                                            ##
################################################################################
################################################################################
#The first thing we need to do is make sure our Biolinux installation and all of
#  the relevant programs are up-to-date.  

#These commands update all of the repositories, then upgrade any programs with
#  recent releases.  This is very important. For example, we will be assembling
#  with Trinity...which relies on Bowtie (v1), but the version included in the
#  Biolinux AMI is Bowtie (v0).  Updating Bowtie (and everything else) now, will 
#  make our lives MUCH easier.
apt-get update
apt-get upgrade


#Organizing the data for even small projects (like this one) is a challenge. 
#  There is no hard and fast rule, for directory structure/organization, but I 
#  try to stick with the Noble (2009) reccomendations.
#  http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424


#Here we are setting variable names for major directories.  This will streamline
#  downstreal processing and reduce the chance for errors.  It is easier to get 
#  it right once rather than a few dozen times.  
HOME_DIR="/data/ASMomics2015" #<- ENTER THE NAME OF YOUR HOME DIR HERE
DATA_DIR="$HOME_DIR/data"
RESULTS_DIR="$HOME_DIR/results"
BIN_DIR="$HOME_DIR/bin"

mkdir $DATA_DIR $RESULTS_DIR $BIN_DIR

#Lets set the number of available processors for consistent use in downstream
#  analyses.
PROCS=$(nproc)

#Now we will move into the data directory to begin processing
cd $DATA_DIR


                        #~~~~~~~~~~~~~~~~~~~~#
                        #  GETTING THE DATA  #
                        #~~~~~~~~~~~~~~~~~~~~#
#First we need to download the data from my dropbox account.
# md5 9d3120b7e2a532d448b276e055eaafe2
wget -O batUCE_batch1_2015-02.tgz \
    https://www.dropbox.com/s/uhoi3y0xzvczobi/batUCE_batch1_2015-02.tgz?dl=0
#as a back up.  You can also get the data here:
# wget ftp://crocgenomes.org/pub/ASM/batUCE_batch1_2015-02.tgz

#To save space and increase transfer speeds, the sequence data is compressed and
#  archived.  Expand the data into a useable format.
tar -xvzf batUCE_batch1_2015-02.tgz

#Lets take an opportunity to look at the data...
#First, lets look at the all the files, their sizes, and permissions
ls -l *

#The file permission associated with these files give the owner (you) permission
#  to write onto/over them.  Whenever you are working from RAW data, it is best
#  to preserve them in their original state.  Remove the write permissions.
chmod a-w *.tgz

#Its always nice to know the size of your files, but the previous command 
#  outputs this data in bytes.  It is more intuitive to look at this in a "human"
#  readable number.
ls -lh *

#You should notice that your file permission have changed, the "w" is gone, and
#  file sizes make more sense (G = gigabytes, M = megabytes, K = kilobytes).

#Now lets look at the actual data...
head traCir_R*.fq

#We know our files are large, and they are in the FastQ format, but how big are
#  they
wc *.fq
#This gives us the TOTAL number of lines, words, characters, and the filename.
#  but we know that the .FQ format has 4 lines per entry.  Lets get a better idea

wc -l *.fq | awk '{print $1/4"\t"$2}'       
#This is the first use of the "|" (pipe).  This is a very valuable tool and
#  allows us to direct output from one progmram/command as input into another.

#or by millions
wc -l *.fq | awk '{print $1/4000000"\t"$2}' 

#Also, notice that the number of reads are identical between libraries
#  R1 vs. R2.  This expected from paired end data.


################################################################################
################################################################################
##                                                                            ##
##              QUALITY CONTROL OF RAW SEQUENCE DATA                          ##
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Pick a sample that you will use for processing.  Designate the R1 and R2 reads
#  here.  Do not include a file extension (ex. .fq)

#braCav = Brachyphylla cavernarum
#chrAur = Chotopterus auritus
#eptFur = Eptisicus furinalis
#macCal = Macrotus califronicus
#micMeg = Micronycteris megalotis
#nycHum = Nycticeus humeralis
#perSub = Perimyotis subflavus
#phyDis = Phylostomus discolor
#rhiAle = Rhinophylla alethina
#tonSau = Tonatia saurophila
#traCir = Trachops cirrhosus
#uroBil = Uroderma bilobatum

BAT=traCir                      #<---------------enter the bat abbreviation here
BAT_R1=$BAT"_R1"
BAT_R2=$BAT"_R2"

#create the directory and set a variable
SEQ_DIR="$RESULTS_DIR/seqQC"

mkdir $SEQ_DIR
cd $SEQ_DIR

fastqc \
        --outdir . \
        --threads $PROCS \
        $DATA_DIR/$BAT_R1.fq \
        $DATA_DIR/$BAT_R2.fq
#Open the FastQC .html output on your local computer.  This can be done by 
#  transfering the file through FileZilla or copying it to DropBox 
#  cp $BAT_R1*.fastq.html <location of Dropbox>
#  cp $BAT_R2*.fastq.html <location of Dropbox>

#Even though the data looks to be high quality, it is good practice to filter
#  out low quality reads and trim adapter seuqences.  This will be done in a
#  series of steps using the FastX toolkit.  There is more than one way to
#  skin a cat...this is my preferred way.

#First lets get the Trimmomatic binary. Since this is a program, lets put it in
#  the /bin directory
cd $BIN_DIR
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip

#uncompress the file
unzip Trimmomatic-0.33.zip

#and set a variable for the path
TRIM_BIN=$BIN_DIR/Trimmomatic-0.33

#Now lets go back to our SEQ_DIR
cd $SEQ_DIR

#First get the adapter sequeancs into a file for Trimmomatic file
echo ">i7">TruSeq4-PE.fa
echo "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG">>TruSeq4-PE.fa
echo ">i5">>TruSeq4-PE.fa
echo "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT" >>TruSeq4-PE.fa


                        #~~~~~~~~~~~~~~~~~~~~#
                        # FILTERING THE DATA #
                        #~~~~~~~~~~~~~~~~~~~~#
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
    -threads $PROCS \
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
        --threads $PROCS \
        "$BAT_R1"_filterPaired.fq \
        "$BAT_R1"_filterUNPaired.fq \
        "$BAT_R2"_filterPaired.fq \
        "$BAT_R2"_filterUNPaired.fq

#Is there anything we can do to increase the quality of the sequences?  Look
#  through the Trimmomatic manual and and FastQC results.  Can we re-run
#  Trimmomatic to improve our overall quality?
java -jar $TRIM_BIN/trimmomatic-0.33.jar \
        PE \
        -threads $PROCS \
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
        "$BAT_R2"_filterPaired.fq \
        $BAT"_RX_UNpaired.fq"


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


                        #~~~~~~~~~~~~~~~~~~~~#
                        #  INSTALL TRINITY   #
                        #~~~~~~~~~~~~~~~~~~~~#
#But oh crap.  We need to install trinity.  Again, this is a program.  It belongs
#  in the /bin folder
cd $BIN_DIR

#Now dowload and uncompress the source code.
wget https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.0.6.tar.gz
tar -xvzf v2.0.6.tar.gz

#and move on to compiling trinity
cd trinityrnaseq-2.0.6

#take a deep breath
make
#a bunch of stuff will flash across the screen and it is going to be very scary...

#So.  Great.  Now Trinity is installed.  Lets add it to our $PATH so that we can
#  call it directly from any subdirectory.
export PATH=$PATH:$(pwd)


                        #~~~~~~~~~~~~~~~~~~~~#
                        #     ASSEMBLE       #
                        #~~~~~~~~~~~~~~~~~~~~#
#Ok.  Now we can get back to our assembly directory and start working
cd $ASSEMBLY_DIR

#start the De Novo Assembly
Trinity \
        --seqType fq \
        --max_memory 20G \
        --left  $SEQ_DIR/$BAT_R1"_filterPaired.fq" \
        --right $SEQ_DIR/$BAT_R2"_filterPaired.fq" \
        --single $SEQ_DIR/$BAT"_RX_UNpaired.fq" \
        --CPU $PROCS \
        --output $ASSEMBLY_DIR/$BAT"_trinity" \
        --full_cleanup 

#

#on average these assemblies take ~an hour to complete.


################################################################################
################################################################################
####                                                                        ####
####                 Exploring the Assembly                                 ####
####                                                                        ####
################################################################################
################################################################################
#lets use some preliminary tools to investigate the assembly. Our main goal is
#  to find assembled mitochondrial genomes, and to calculate the amount of
#  coverage we are getting.

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

ASSEMBLED_SEQS=$ASSEMBLY_DIR/$BAT"_trinity.Trinity.fasta"

ANNOTATION_DIR="$RESULTS_DIR/annotations"
mkdir $ANNOTATION_DIR
cd $ANNOTATION_DIR


                        #~~~~~~~~~~~~~~~~~~~~#
                        #      BLAST         #
                        #~~~~~~~~~~~~~~~~~~~~#
#Lets extract some mitoGenomes from NCBI to run a local blast.  Do that now.
#  I named my mitochondiral genome file chiropteraMitoGenomes_2015-06-10.fasta.
#  Since this is a file that will NOT be altered and was generated by others, I
#  store it in the "data" directory. You MAY have to give yourself permision.

MITOGENOMES=$DATA_DIR/chiropteraMitoGenomes_2015-06-10.fasta #<--- name of NCBI file

chmod a+w /data/ASMomics2015/data/

#The first thing you need to do is to create a blast index.
makeblastdb -in $MITOGENOMES -dbtype nucl

#Now run your blast search 
blastn \
    -db $MITOGENOMES \
    -query $ASSEMBLED_SEQS \
    -out $ANNOTATION_DIR/$BAT"_assemVsMito_blastn.out" \
    -outfmt 6 \
    -num_threads $PROCS
 
#To check out the length distributions of your blast hits.  Use the field info
#  below.  

# 1 Query
# 2 Target
# 3 Percent identity
# 4 Alginment length
# 5 Number of mismatches
# 6 Number of gaps
# 7 Query start position
# 8 Query end position
# 9 Subject start position
# 10 Subject end position
# 11 Evalue
# 12 Bit score

#Since the file is pretty small browse through it using "less" to find the
#  assemble contig that is the best hits to the mitogenome.  (There are 
#  other/better ways to do this: see below).

less $ANNOTATION_DIR/$BAT"_assemVsMito_blastn.out"

#Instead of using "less" you can use unix to "sort" the file based on specific
#  columns.
sort -n $ANNOTATION_DIR/$BAT"_assemVsMito_blastn.out" -nr -k12 | head

#Find the assembled contig with the longest/best hit to the mitogenomes and
#  extract it.  This can be done manually or with Samtools.
#  Which of the sequences is the best hit to the FULL mitoGenome.  

MITO_CONTIG="TR2755|c0_g1_i1"     #<------ Enter the mito contig here.


                        #~~~~~~~~~~~~~~~~~~~~#
                        #      MAPPING       #
                        #~~~~~~~~~~~~~~~~~~~~#
#Now extract it into it's own seperate file.
samtools faidx \
    $ASSEMBLED_SEQS \
    $MITO_CONTIG \
    >$BAT"_mitoGenome.fas"


#Some of the mitogenomes will vary in quality. Lets look at sequence coverage.
#  First we need to map all of our reads back to the mitochondrial contig
bowtie2-build \
    $BAT"_mitoGenome.fas" \
    $BAT"_mitoGenome.fas"

#Now lets map our cleaned data sequnce data to the miotchondrial contig.
bowtie2 \
    --threads $PROCS \
    -x $BAT"_mitoGenome.fas" \
    -1 $SEQ_DIR/$BAT_R1"_filterPaired.fq" \
    -2 $SEQ_DIR/$BAT_R2"_filterPaired.fq" \
    -U $SEQ_DIR/$BAT"_RX_UNpaired.fq" \
    -S $BAT"_mitoGenome.SAM"


#Convert the SAM file to BAM
samtools view \
    -Sb \
    $BAT"_mitoGenome.SAM" \
    | samtools sort \
        -f \
        - \
        $BAT"_rawMapped_sorted.BAM"
 
#... and then convert BAM to BED.  Welcome to bioinformatics.
 bedtools bamtobed \
    -i $BAT"_rawMapped_sorted.BAM" \
    >$BAT"_rawMapped_sorted.BED" 

#A Bed and Sam files are ASCII files while BAM is binary.  Here we will take the
#  opportunity to discuss the file structure of BED and SAM.

#And since nothing is easy we have to create a "special" genome file.  Check out
#  the bedtools doccumentation to see the file type.  
fasta_formatter \
    -i $BAT"_mitoGenome.fas" \
    -t \
    | awk '{print $1"\t"length($2)}' \
     >$BAT"_mitoGenome.tab"


#Now we can FINALLY calculate the overall coverage of our genome....

#at each position...
bedtools genomecov \
    -i $BAT"_rawMapped_sorted.BED" \
    -g $BAT"_mitoGenome.tab" \
    -d

#at each position (with reads displayed)...
samtools \
    mpileup \
    $BAT"_rawMapped_sorted.BAM"


#Or across the entire genome (with AWK).
bedtools genomecov \
    -i $BAT"_rawMapped_sorted.BED" \
    -g $BAT"_mitoGenome.tab" \
    -d \
    | awk '{sum += $3; n++ } END {print sum / n;}' 
#The result is the average coverage across the ENTIRE ASSEMBLED mitogenome.

#Want to look at it more closely (time permitting).  Lets use IGV.
#  1) Install the java based program based on your operating system based on 
#     those available here. https://www.broadinstitute.org/software/igv/download
#  2) Download the assembled mitogenome (through Filezilla/FTP) and load it 
#     as a genome
#  3) Download the *_rawMapped_sorted.BED then File->load from file.
#  4) Now expand the track.  How does the read distribution look?
                       

                        #~~~~~~~~~~~~~~~~~~~~#
                        #      THE END?      #
                        #~~~~~~~~~~~~~~~~~~~~#
#At this point there are a million other things that can be done.  From
#  annotating individual genes to various phylogenetic analyses.  Since we don't
#  have an entire semester we will stop here.  If there is time, consider
#  RERUNNING these analyses with a different species...Or the same species and 
#  see if you get a similar result.  It is possible by only modifying a few of
#  the given commands, but you will have to pay attention.  You could download
#  DNA reads from NCBI's SRA and reconstroct mitochondrial genomes.  If you are
#  interested, here some example syntax to get you started with the white tiger
#  mitogenome.

#Download the data from the SRA.  This will take forever, so dont try and do this
#  today.  Skip this step.
# cd $DATA_DIR
# wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR171/SRR1712667/SRR1712667.sra

#Data stored on the SRA is saved in a special format (again with the formats).
#  It needs to be converted from .sra to .fastq.  This will generate two HUGE
#  files (>110 Gb).  This is overkill for what we want to do and beyond our 
#  capabilites for the day.  Skip.
# fastq-dump --split-files SRR1712667.sra

#Since we only need a very small portion of all the DNA reads to get resonable
#  coverage of the mitogenome, lets subsamle only 10M reads from the R1 and R2
#  ends of the fragment.  Again, for today...skip.
# head -40000000 $DATA_DIR/SRR1712667_1.fastq >$SEQ_DIR/SRR1712667_10M_1.fastq
# head -40000000 $DATA_DIR/SRR1712667_2.fastq >$SEQ_DIR/SRR1712667_10M_2.fastq

cd $SEQ_DIR
#You can donwnload the data from the previous steps here:

#Quick note:  I will remove this link ~9/1/2015, so be aware.
wget -O SRR1712667_10M.tgz \
    https://www.dropbox.com/s/4vkpk3bbcv6rk6m/SRR1712667_10M.tgz?dl=0

# Uncompress it.
tar -xvzf SRR1712667_10M.tgz

#Now you can use this data, and the by making slight modifications to the script
#  given above, you could assemble the mitogenome from the tiger.  All the
#  the relevant assembly information can be found here: 
#  http://www.ncbi.nlm.nih.gov/sra/SRR1712667/.  Make sure to note that the 
#  fragments have an insert size of 400 bp.
                    

                        #~~~~~~~~~~~~~~~~~~~~~#
                        #   TO BE CONTINUED   #
                        #~~~~~~~~~~~~~~~~~~~~~#















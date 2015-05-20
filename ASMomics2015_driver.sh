# ASMomics 2015
# Mitogenome assembly in vespertilionid and phylostomid bats
# RNPlatt
# 20 May 2015
# v 0.1

#The first thing we nee to do is set up a directory structure.  There is no 
#  hard and fast rule, but I try to stick with the Noble (2009) structure.
#  http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424

mkdir mitogenomes
mkdir mitogenomes/data mitogenomes/results

#Here we are setting variable names for major directories.  This will streamline
#  downstreal processing and reduce the chance for errors.  It is easier to get 
#  it right once rather than a few dozen times.  
DATA_DIR="mitogenomes/data"
RESULTS_DIR="mitogenomes/results"

#Now we will move into the data directory to begin processing
cd $DATA_DIR

#First we need to download the data from an FTP site.
wget <address>..initialBatUCE_2015-01.tgz

#To save space and increase transfer speeds, the sequence data is compressed and
#  archived.  Expand the data into a useable format.
tar -xvzf initialBatUCE_2015-01.tgz

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
head Index_Bat_9_1.fq

# description of fastq slides here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#We know our files are large, and they are in the FastQ format, but how big are
#  they
wc *.fq
#This gives us the TOTAL number of lines, words, characters, and the filename.
#  but we know that the .FQ format has 4 lines per entry.  Lets get a better idea

wc -l *.fq | awk '{print $1/4"\t"$2}'                                           <- first use of the pipe , explain

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
#       in the results directory ***

cd $RESULTS_DIR


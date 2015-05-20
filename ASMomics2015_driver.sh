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
#  file sizes make more sense.



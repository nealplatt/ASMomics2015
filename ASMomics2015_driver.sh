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

DATA_DIR="mitogenomes/data"
RESULTS_DIR="mitogenomes/results"

cd $DATA_DIR
wget <address>..initialBatUCE_2015-01.tgz



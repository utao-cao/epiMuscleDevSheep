#/bin/bash -e
###
 # @Authour: utao.cao
 # @input: 1. accept homer output log file; 2. outfile prefix [string]
### 
log_file=$1  summaryReportPrefix=$2
# cat logbar_motif.txt | grep  "similar to" | grep -E "1e-[1-9][0-9]+"  >> summary_${2}_motifResults
if [[ -f summary_${summaryReportPrefix}_motifResults ]];then
	ls $log_file >> summary_${summaryReportPrefix}_motifResults.txt
else
	ls $log_file >> summary_${summaryReportPrefix}_motifResults.txt
fi
cat $log_file | grep  "similar to" | grep -E "1e-[1-9][0-9]+"  >> summary_${summaryReportPrefix}_motifResults.txt

# finally arrange
# cat summary_atacIDR_motifResults.txt | sed -e "s@/.*@@g" | cut -d " " -f 1,7 | sed  's/([-\?\.A-Za-z0-9]*)//g'
# cat summary_atacIDR_motifResults.txt | sed -e "s@/.*@@g" | cut -d " " -f 1,7 | sed -e 's/([-\?\.A-Za-z0-9]*)//g' -e 's@[[:space:]][0-9]*@@g'
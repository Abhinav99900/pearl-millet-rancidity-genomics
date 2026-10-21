#!/bin/bash
# Extract VCF records lying within the candidate gene models (mRNA +/- 100 bp).
# Usage: bash 01_extract_gene_regions.sh <input.vcf.gz> <output.vcf> uni/universe_regions.txt
# Run for the filtered and the raw (unfiltered) VCF of each variety, e.g.:
#   bash 01_extract_gene_regions.sh SO_14026_A1.filtered.vcf.gz vcf/u732.filt.vcf uni/universe_regions.txt
#   bash 01_extract_gene_regions.sh SO_14026_A1.raw.vcf.gz      vcf/u732.raw.vcf  uni/universe_regions.txt
#   (and likewise u538.filt.vcf / u538.raw.vcf for GHB-538)
zcat "$1" | awk -v R="$3" 'BEGIN{while((getline l < R)>0){split(l,a," ");for(b=int(a[2]/10000);b<=int(a[3]/10000);b++){k=a[1]":"b; L[k]=L[k]" "a[2]","a[3]}}}
/^#/{print;next}{k=$1":"int($2/10000); if(k in L){n=split(L[k],iv," "); for(i=1;i<=n;i++){split(iv[i],p,","); if($2>=p[1]&&$2<=p[2]){print;break}}}}' > "$2"

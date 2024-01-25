# Sequentially batch process all fastq.gz files in the data folder
find ./data -type f -name '*fastq.gz*' | sed 's|./data/||' | while read -r file; do
    Rscript --vanilla CrispRscreenAnalysis_V1_EI.R "$file"
done


# test a single file
# Rscript --vanilla CrispRscreenAnalysis_V1_EI.R Sample_Het-1/Het-1_S7_L003_R1_001.fastq.gz

# parallel
# ??? parallel does improve the speed maybe the alignment file is shared and can be access by only one process at the same time
find ./data -type f -name '*fastq.gz*' | sed 's|./data/||' | xargs -I {} -P 5 -n 1 Rscript --vanilla CrispRscreenAnalysis_V1_EI.R {}

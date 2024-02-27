## ----setup, include=FALSE-------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, error = FALSE, message = FALSE)


## ----eval=FALSE, include=FALSE--------------------------------------------------------------------------------------------------------------------
## rm(list = ls())
## # OneDrive_dir <- "~/OneDrive - Inside MD Anderson/Box Sync/Genomic"
## # change main directory
## # OneDrive_dir <- "/Users/cgu3/Documents/Crisp"
## # Functions_dir <- paste0(OneDrive_dir,"/My_R_functions")
## # project_dir <- paste0(OneDrive_dir,"/FarhadDanesh/CRISPR/CRISPR_screen_KidneyGene_Jan2024")
## # result_dir <- paste0(project_dir,"/output")
## # data_dir <- paste0(project_dir,"/data")
## # main_dir <- paste0(project_dir,"/script/Subread")
## # setwd(main_dir)


## -------------------------------------------------------------------------------------------------------------------------------------------------
# install.packages("here")
library(here)
setwd(here())
data_dir <- paste0(getwd(), '/data')
project_dir <- getwd()
result_dir <- paste0(getwd(),"/output")


## ----eval=FALSE, include=FALSE--------------------------------------------------------------------------------------------------------------------
## # this chuck is intented to used for batch processing after this rmarkdown is converted to a R script
args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  # expect args[1] to be relative path to the fastq.gz input file
  args[2] = gsub(".fastq.gz", '.csv', args[1])
}

## -------------------------------------------------------------------------------------------------------------------------------------------------
library(ShortRead, verbose = FALSE, quietly = TRUE)
library(Rsubread, verbose = FALSE, quietly = TRUE)

## ----eval=FALSE, include=FALSE--------------------------------------------------------------------------------------------------------------------
## library(Rsubread)
myFQs <- paste0(data_dir, '/', args[1])
myMapped <- align("Addgene",myFQs,output_file = gsub(".fastq.gz",".bam",myFQs),
                  nthreads=4,unique=TRUE,nBestLocations=1,type = "DNA",TH1 = 1,
                  maxMismatches = 0,indels = 0)



## -------------------------------------------------------------------------------------------------------------------------------------------------
library(GenomicAlignments, verbose = FALSE, quietly = TRUE)
temp <- readGAlignments(gsub(".fastq.gz",".bam",myFQs))
temp


## -------------------------------------------------------------------------------------------------------------------------------------------------
counts <- data.frame(table(seqnames(temp)),row.names = "Var1")
counts[1:2,,drop=FALSE]


output_path <- paste0(result_dir, "/", args[2])
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(counts, file = output_path, row.names = TRUE, col.names = TRUE)


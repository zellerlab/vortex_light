#!/usr/bin/env Rscript


library(optparse)
library(tidyverse)

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

option_list = list(
  make_option(c("--kraken2_res_path"), type="character", default=NULL, 
              help="Path to folder with kraken2 result tables", metavar="character"),
  make_option(c("--PathSeq_res_path"), type="character", default=NULL, 
              help="Path to folder with PathSeq result tables", metavar="character"),
  make_option(c("--libsize_res_path"), type="character", default=NULL, 
              help="Path to folder with libsize results", metavar="character"),
  make_option(c("--lib_layout_res_path"), type="character", default=NULL, 
              help="Path to folder with library_layout results", metavar="character"),
  make_option(c("--flagstats_res_path"), type="character", default=NULL, 
              help="Path to folder with flagstats", metavar="character"),
  make_option(c("--N_raw_counts_path"), type="character", default=NULL, 
              help="Path to folder with number of raw_counts", metavar="character"),
  
  
  make_option(c("--out_folder"), type="character", default=NULL, 
              help="output folder path", metavar="character"),
  make_option(c("--libdir"), type="character", default=".", help="path to functions", metavar="character"),
  make_option(c("--gtdb_markers"), type="character", default=NULL, help="path to gtdb markers", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

source(file.path(opt$libdir, "functions_read_in_profiled_data_210702.R"))

out.folder <- opt$out_folder
message("kraken2")
#kraken2
if(!(is.null(opt$kraken2_res_path))){
  if(length(list.files(opt$kraken2_res_path))>0){
    res.kraken2 <- .f_read_in_files_kraken2(path_to_folder = opt$kraken2_res_path,
                                            tax.level = "genus")
    saveRDS(res.kraken2,paste0(out.folder,"/res_kraken2.rds"))  
  }else{
    message("kraken2 path empty")
  }
}

#pathseq
message("pathseq")
if(!(is.null(opt$PathSeq_res_path))){
  if(length(list.files(opt$PathSeq_res_path))>0){
  res.pathseq <- .f_read_in_files_PathSeq(path_to_folder = opt$PathSeq_res_path,
                                                    tax.level = "genus")
  saveRDS(res.pathseq,paste0(out.folder,"/res_pathseq.rds"))
  } else{
    message("PathSeq path empty")
  }
}

#libsize
if(!(is.null(opt$libsize_res_path))){
  if(length(list.files(opt$libsize_res_path))>0){
    res.libsize <- .f_read_in_libsize(path_to_folder = opt$libsize_res_path)
    saveRDS(res.libsize,paste0(out.folder,"/res_libsize.rds"))
  } else{
    message("libsize path empty")
  }
}

#liblayout
if(!(is.null(opt$lib_layout_res_path))){
  if(length(list.files(opt$lib_layout_res_path))>0){
    res.lib_layout <- .f_read_in_lib_layout(path_to_folder = opt$lib_layout_res_path)
    saveRDS(res.lib_layout,paste0(out.folder,"/res_lib_layout.rds"))
  } else{
    message("liblayout path empty")
  }
}

#flagstats --> total numbers of reads before QC
#No needed anymore in vknight build >= 827290457b
if(!(is.null(opt$flagstats_res_path))){
  if(length(list.files(opt$flagstats_res_path))>0){
    res.flagstats <- .f_read_in_flagstats(path_to_folder = opt$flagstats_res_path)
    saveRDS(res.flagstats,paste0(out.folder,"/res_flagstat_total_reads.rds"))
  } else{
    message("flagstats path empty")
  }
}

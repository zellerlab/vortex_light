###
#Functions to collate profiler outputs and library size information into .rds files
#
###
library(stringr)
library(progress)

.f_read_in_files_kraken2 <- function(path_to_folder,tax.level){
  ### Read in kraken2 result files and return matrix with counts per bacteria and sample
  
  #convert tax.level to kraken2 compatible taxonoic symbol in order to subset the DF
  tax.sym <- toupper(str_extract(tax.level,pattern = "^[A-z]"))
  
  var.names <- c("pct.total","counts.sum","counts.only.here","tax.symbol","tax.ID","tax.name")
  file_list <- list.files(path=path_to_folder)  
  ### initialize output df
  counts.df <- tibble(tax.name = character(0))
  ### iterate over every file and select counts at the selected tax level
  pb <- progress_bar$new(total=length(file_list))
  total.counts.mapped <- tibble(Sample_ID = !!gsub(x = file_list,pattern = ".txt",replacement = ""),tot.counts.mapped = double(length(file_list)))
  i <- 1
  for(i in seq(1,length(file_list))){
    #message(i)
    ### read in file
    #chec if file is empty
    if(file.info(paste0(path_to_folder,file_list[i]))$size==0){
      message(paste0("\nSample ",file_list[i]," is empty - skipping"))
      pb$tick()
      next
    }
    #c.f <- read_tsv(paste0(path_to_folder,file_list[i]),col_names = var.names,col_types = cols())
    c.f <-read.csv(paste0(path_to_folder,file_list[i]),header = FALSE,sep = "\t",comment.char = "",check.names = FALSE,col.names = var.names)
    c.f$tax.name <- str_replace_all(c.f$tax.name,pattern = " ",replacement = "")
    c.total.counts.mapped <- c.f %>% filter(tax.name == "root") %>% pull(counts.sum)
    if(!(is_empty(c.total.counts.mapped))){
      total.counts.mapped[i,2] <- c.total.counts.mapped
    }  
    ### extract all bacterial counts
    bac.start <- which(c.f[,6] == "Bacteria")
    if(is_empty(bac.start)){
      message(paste0("\nNo bacteria profiled in sample ",file_list[i], " - skipping"))
      pb$tick()
      next
    }
    
    virus.start <- which(c.f[,6] == "Viruses") #get row in table at which "viruses" start
    archea.start <- bac.end <- which(c.f[,6] == "Archaea") #get row in table at which "archea" start
    
    bac.end <- min(virus.start,archea.start)-1 #get last row with entries for bacteria
    archea.end <- archea.end <- min(virus.start)-1
    
    if(!(is.finite(bac.end))){
      bac.end <- nrow(c.f)
    }
    if(!is.finite(archea.end)){
      archea.end <- bac.end
    }
    
    #subset to keep only the entries mapped to bacteria AND archea
    bac.reads <- c.f[(bac.start:archea.end),]  
    
    ### select counts and tax names
    c.counts.df <- 
      bind_rows(
        bac.reads %>% filter(tax.symbol == tax.sym) %>% select(tax.name,counts.sum), #select counts at genus level
        bac.reads %>% filter(tax.symbol == "D") %>% select(tax.name,counts.sum) %>% #select "Bacterial" and "Archeal" counts, sum them togehter as "Bacteria" and merge with counts data
          mutate(tax.name = "Bacteria") %>% 
          group_by(tax.name) %>% 
          summarise(counts.sum = sum(counts.sum))) %>% 
      arrange(-counts.sum) %>% 
      rename(!!gsub(x = file_list[i],pattern = ".txt",replacement = "") := counts.sum)
    
    # Add to data from other samples
    counts.df <- suppressMessages(full_join(counts.df,c.counts.df,by="tax.name"))
    pb$tick()
    #print(i)
  }
  ### convert to output matrix
  counts.df <- counts.df %>% mutate_at(vars(-tax.name),as.numeric)
  counts.mat <- as.matrix(counts.df[,-1])
  rownames(counts.mat) <- (counts.df$tax.name)
  counts.mat[is.na(counts.mat)] <- 0
  #remove everything after the first point in the samplenames
  colnames(counts.mat) <- sub(colnames(counts.mat),pattern = ".kraken2_report",replacement = "")
  
  res.list <- list(counts.mat,total.counts.mapped)
  names(res.list) <- c("counts.mat","total.counts.mapped")
  #return(counts.mat)
  return(res.list)
}


.f_read_in_files_PathSeq <- function(path_to_folder,tax.level){
  ### Read PathSeq output files, select tax level (eg genus) and return 2 matrices:
  #1) matrix with score_normalized values
  #2) with count values of !unambiguously! mapped reads
  
  
  file_list <- list.files(path=path_to_folder,pattern = "\\.txt$")  
  if(is_empty(file_list)){
    file_list <- list.files(path=path_to_folder,pattern = "\\.scores$")  
  }
  
  tax.sym <- tax.level
  
  ### initialize output df
  score.df <- tibble(tax.name = character(0),tax_id=double(0))
  counts.df <- tibble(tax.name = character(0),tax_id=double(0))
  
  ### iterate over every file and select counts at the selected tax level
  empty.counter <- 0
  pb <- progress_bar$new(total=length(file_list))
  for(i in seq(1,length(file_list))){
    ## read in file
    #c.f <-read_tsv(paste0(path_to_folder,file_list[i]),col_types = cols())
    c.f <-read.csv(paste0(path_to_folder,file_list[i]),header = TRUE,sep = "\t",comment.char = "",check.names = FALSE)
    
    #check if file is empty
    if(nrow(c.f) == 0){
      message(paste0("\nSample ",file_list[i]," is empty - skipping"))
      pb$tick()
      empty.counter <- empty.counter+1
      next
    }
    
    ### select counts and tax names
    #Edit 22-09-01: When Archaeal reads also should be considered: 
    #Select "Archaea" and "Bacteria", sum up "score" and "unambigous" of the "superkingdom" and then normalize every tax.lvl counts against this number
    
    c.combined.df <- 
      c.f %>% 
      filter(type == "superkingdom",
             kingdom %in% c("Bacteria","Archaea")) %>% 
      mutate(name = "Bacteria") %>% 
      group_by(name) %>% 
      summarise(score = sum(score),
                unambiguous = sum(unambiguous)) %>% 
      add_column(tax_id = 0000) %>%  #just to have a tax_id
      bind_rows(.,c.f %>% filter(kingdom %in% c("Bacteria","Archaea"),
                                 type == tax.sym)) %>% 
      mutate(score_normalized = score/score[1]*100)
    
    c.score.df <-
      c.combined.df %>% 
      select(name,tax_id,score_normalized) %>%
      rename(!!gsub(x = file_list[i],pattern = ".pathseq.txt",replacement = "") := score_normalized,
             tax.name = name)
    
    # also unambiguous counts
    c.counts.df <-
      c.combined.df %>% 
      select(name,tax_id,unambiguous) %>%
      rename(!!gsub(x = file_list[i],pattern = ".pathseq.txt",replacement = "") := unambiguous,
             tax.name = name)
    
    
    #Join with data from other samples
    score.df <- full_join(score.df,c.score.df,by=c("tax.name","tax_id"))
    counts.df <- full_join(counts.df,c.counts.df,by=c("tax.name","tax_id"))
    
    pb$tick()
  }
  
  #'In the PathSeq DB of Dohlman et al there seems to be a naming issue:
  #'Both,tax_ID 1937007 and tax_ID 1980680 correspond to the genus "Ileibacterium". 
  #'According to NCBI taxonomy, tax_ID 1980680 refers to "Alterileibacterium
  #'--> correct the naming manually
  score.df <- score.df %>% mutate(tax.name = case_when(tax_id == 1980680 ~ "Alterileibacterium",
                                                       TRUE ~ tax.name)) %>% 
    select(-tax_id)
  counts.df <- counts.df %>% mutate(tax.name = case_when(tax_id == 1980680 ~ "Alterileibacterium",
                                                         TRUE ~ tax.name)) %>% 
    select(-tax_id)
  
  
  ### convert to output matrices and return list
  score.df <- score.df %>% mutate_at(vars(-tax.name),as.numeric)
  score.mat <- as.matrix(score.df[,-1])
  rownames(score.mat) <- (score.df$tax.name)
  score.mat[is.na(score.mat)] <- 0
  #remove all-zero rows
  score.mat <- score.mat[rowSums(score.mat)>0,,drop=F]
  #remove everything after the first point in the samplenames
  colnames(score.mat) <- sub(colnames(score.mat),pattern = ".pathseq.scores",replacement = "")
  
  counts.df <- counts.df %>% mutate_at(vars(-tax.name),as.numeric)
  counts.mat <- as.matrix(counts.df[,-1])
  rownames(counts.mat) <- (counts.df$tax.name)
  counts.mat[is.na(counts.mat)] <- 0
  counts.mat <- counts.mat[rowSums(counts.mat)>0,,drop=F]
  #remove everything after the first point in the samplenames
  colnames(counts.mat) <- sub(colnames(counts.mat),pattern = ".pathseq.scores",replacement = "")
  
  
  res.list <- list(score.mat,counts.mat)
  names(res.list) <- c("score_normalized.mat","counts_unambiguous.mat")
  
  return(res.list)
}

.f_read_in_libsize <- function(path_to_folder){
  file_list <- list.files(path_to_folder)
  
  df.libsize <- tibble()
  for(i in seq(1,length(file_list))){
    tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t",comment.char = "",header = F)
    tmp.file$ID <- sub(file_list[i],pattern = ".libsize.txt",replacement = "")
    colnames(tmp.file)[1] <- "libsize"
    df.libsize <- bind_rows(df.libsize,tmp.file)
  }
  df.libsize <- df.libsize %>% relocate(ID)
  return(df.libsize)
}

.f_read_in_lib_layout <- function(path_to_folder){
  file_list <- list.files(path_to_folder)
  df.lib_layout <- tibble()
  for(i in seq(1,length(file_list))){
    tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t",comment.char = "",header = F)
    tmp.file$ID <- sub(file_list[i],pattern = ".is_paired.txt",replacement = "")
    colnames(tmp.file)[1] <- "lib_layout"
    df.lib_layout <- bind_rows(df.lib_layout,tmp.file)
  }
  df.lib_layout <- df.lib_layout %>% relocate(ID)
  return(df.lib_layout)
}

.f_read_in_flagstats <- function(path_to_folder){
  #No needed anymore in vknight build >= 827290457b
  #read in the total number of reads from the flagstats output -> Corresponds to ALL reads in the FastQ input file (before QC)
  file_list <- list.files(path_to_folder,pattern = "*.flagstats.txt")
  df.flagstats <- tibble()
  for(i in seq(1,length(file_list))){
    tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t",comment.char = "",header = F)
    #Extract total number of reads (1st number in 1st row of the flagstats output)
    tmp.total.reads <- as.numeric(gsub("\\ +.*", "",tmp.file[1,1]))
    names(tmp.total.reads) <- sub(file_list[i],pattern = ".flagstats.txt",replacement = "")
    
    tmp.res <- enframe(tmp.total.reads,name = "ID",value = "total.libsize")
    
    df.flagstats <- bind_rows(df.flagstats,tmp.res)
  }
  return(df.flagstats)
}

.f_read_in_raw_counts_number <- function(path_to_folder){
  #'take raw_counts computed by vknight (>build 827290457b) when preprocessing the files with multiQC
  #'exports a tibble with Sample_ID and raw_read_count
  file_list <- list.files(path_to_folder)
  
  df.raw_counts <- tibble(Sample_ID=character(),raw_counts=double())
  for(i in seq(1,length(file_list))){
    tmp.file <- read.table(paste0(path_to_folder,file_list[i]),sep = "\t") %>% 
      separate(V1,sep = " ",into = c("paired_info","V1")) %>% 
      mutate(paired_info = as.numeric(paired_info))
    tmp.raw_counts <- tibble(Sample_ID = str_remove(file_list[i],pattern = ".txt"),
                             raw_counts = tmp.file$paired_info * tmp.file$V2)
    
    df.raw_counts <- bind_rows(df.raw_counts,tmp.raw_counts)
  }
  return(df.raw_counts)
}



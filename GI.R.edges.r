#!/usr/bin/env Rscript
suppressMessages(library("optparse"))


option_list = list(
  make_option(c("-i", "--interactions"), type="character", default=NULL, 
              help="Looping bedpe8", metavar="character"),
  make_option(c("-b", "--bed"), type="character", default=NULL, 
              help="BED regions will be used to generate the network", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);





suppressMessages(library(GenomicInteractions))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

suppressMessages(library(reticulate))

options(repr.matrix.max.rows=10)

bed_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)

  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }

  if(length(df)<3){
    stop("File has less than 3 columns")
  }

  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]

  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }

  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end, names=id)))
  }
    
    return(gr)

}




interactions = opt$interactions

nodes = opt$bed

out = opt$out



pet = makeGenomicInteractionsFromFile(
    interactions,
    type="bedpe", 
    experiment_name="LNCaP", 
    description="LNCaP")


cre = bed_to_granges(nodes)

annotation.features = list(cre = cre)


annotateInteractions(pet, annotation.features)


print(
  categoriseInteractions(pet)
)




df_all = as.data.frame(pet[isInteractionType(pet, c("cre"), c("cre"))])


pd = import("pandas")

df = r_to_py(df_all)
df$to_csv(out, sep="\t", index=FALSE)

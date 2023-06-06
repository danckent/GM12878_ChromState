library( GenomicRanges)
library( tibble)


load( "~/Data/GM12878/Chromatin State/Chrom_state_per_gene_for_PCA.RData")

DNase   <- read.table( "/home/dankent/Data/GM12878/ENCFF097LEF_DNase.bed", header = F)
DNase   <- rowid_to_column( DNase, "ID")

# Function to find overlaps of individual chromatin states and % coverage on the fragments and join them in a df
MetaMerge <- function( Genes, DNase) {
  
  test <- makeGRangesFromDataFrame( Genes, keep.extra.columns = TRUE, start.field = "Gene.start..bp.", end.field = "Gene.end..bp.", seqnames.field = "Chromosome.scaffold.name")
  test1 <- makeGRangesFromDataFrame( DNase, keep.extra.columns = TRUE, start.field = "V2", end.field = "V3", seqnames.field = "V1")
  
  hits <- findOverlaps( test, test1) # Finds all overlaps
  match_hit <- data.frame( test[ queryHits( hits)] , data.frame(test1[subjectHits(hits)] ),stringsAsFactors=T)
  bait_match_hit <- data.frame( test[ queryHits( hits)] ,stringsAsFactors=T)
  CS_match_hit <- data.frame( test1[ subjectHits( hits)]  ,stringsAsFactors=T)
  
  Bait_Score_CS <- makeGRangesFromDataFrame( match_hit, keep.extra.columns = T)
  x <- makeGRangesFromDataFrame( bait_match_hit, keep.extra.columns = T)
  y <- makeGRangesFromDataFrame( CS_match_hit, keep.extra.columns = T)
  overlaps <- pintersect( x,y)
  percentageOverlap <- width( overlaps) / width( x)
  match_hit <- data.frame( test[ queryHits( hits)] , data.frame( mcols( test1[ subjectHits( hits)] ) ), data.frame( percentageOverlap), stringsAsFactors=T)
  
  df = as( match_hit, "data.frame")
  return( df)
  
}
DNase_test <- MetaMerge( Genes, DNase)

colnames( DNase_test) <-  c( "seqnames", "start", "end", "width", "strand","ID", "TSS", "HGNC.symbol", "Gene.type", "Gene.stable.ID", "GO.accession", "GO.name", "ID_DNase",
                            "string name", "uint score", "strand", "float signalValue", "float pValue", "float qValue", "int peak", "percentageOverlapDNase")
#Sum all DNase percentages across the each gene by ID 
Gene_makeup_DNase <-  aggregate(x = DNase_test$percentageOverlapDNase, by = list(DNase_test$Gene.stable.ID), FUN = sum)
Gene_makeup_DNase <- data.frame(Gene_makeup_DNase %>% 
                             arrange(Group.1))
#filter the merged Gene & DNAse df for one gene per row  
DNase_unique <- data.frame(DNase_test %>% 
  select(c( "seqnames", "start", "end", "width", "TSS", "strand","ID", "Gene.stable.ID", "HGNC.symbol", "Gene.type", "GO.name",
            "ID_DNase", "string name", "uint score", "strand", "float signalValue", "float pValue", "float qValue", "int peak", "percentageOverlapDNase")) %>% 
  group_by(Gene.stable.ID) %>% 
  distinct(Gene.stable.ID,.keep_all = TRUE) %>% 
  arrange(Gene.stable.ID))

# add the combined percentage overlap for DNase data to the dataframe of single genes per row 
Gene_makeup_DNase <- cbind(DNase_unique, Gene_makeup_DNase)
#clean up
save(Gene_makeup_DNase, file = "~/Data/GM12878/Gene_with_DNase")
Gene_makeup_DNase <- Gene_makeup_DNase %>% 
                             select(c( "seqnames", "start", "end", "width", "ID", "Gene.stable.ID", "HGNC.symbol", "Gene.type", "x"))

# What about all the genes that don't have any percentage open? 

nrow(Genes) #22811 includes all chromosome scaffold 
nrow(Gene_makeup_DNase) #13453
mean(Gene_makeup_DNase$x) #0.01551847 (1.5% of genes are open)

#### overlap the DNase genes with Chromatin state ####

# Function to find overlaps of individual chromatin states and % coverage on the fragments and join them in a df
MetaMerge <- function( Gene_makeup_DNase, ChromState) {
  
  test <- makeGRangesFromDataFrame( Gene_makeup_DNase, keep.extra.columns = TRUE)
  test1 <- makeGRangesFromDataFrame( ChromState, keep.extra.columns = TRUE, start.field = "V2", end.field = "V3", seqnames.field = "V1")
  
  hits <- findOverlaps( test, test1) # Finds all overlaps
  match_hit <- data.frame( test[ queryHits( hits)] , data.frame(test1[subjectHits(hits)] ),stringsAsFactors=T)
  bait_match_hit <- data.frame( test[ queryHits( hits)] ,stringsAsFactors=T)
  CS_match_hit <- data.frame( test1[ subjectHits( hits)]  ,stringsAsFactors=T)
  
  Bait_Score_CS <- makeGRangesFromDataFrame( match_hit, keep.extra.columns = TRUE)
  x <- makeGRangesFromDataFrame( bait_match_hit, keep.extra.columns = TRUE)
  y <- makeGRangesFromDataFrame( CS_match_hit, keep.extra.columns = TRUE)
  overlaps <- pintersect( x,y)
  percentageOverlap <- width( overlaps) / width( x)
  match_hit <- data.frame( test[ queryHits( hits)] , data.frame( mcols( test1[ subjectHits( hits)] ) ), data.frame( percentageOverlap), stringsAsFactors=T)
  
  df = as( match_hit, "data.frame")
  return( df)
  
}
DNase_Chrom <- MetaMerge( Gene_makeup_DNase, ChromState)

# Convert chromatin states to columns with percentageOverlaps as the values. 
x=as.character( paste( "seqnames+start+end+width+ID+Gene.stable.ID+HGNC.symbol+DNase"))
f                         <- as.formula( paste( paste( x, collapse = " + "), "~ CS"))
colnames( DNase_Chrom) <-  c( "seqnames", "start", "end", "width", "strand","ID", "Gene.stable.ID", "HGNC.symbol", "Gene.type", "DNase", "CS", "percentageOverlap")
DNase_Chrom %>% reshape2::dcast( f, value.var = "percentageOverlap", fun.aggregate = sum) %>% select(c("seqnames", "start", "end", "width", "ID", "Gene.stable.ID", 
  "HGNC.symbol", "DNase", "E1", "E2", "E3",  "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11")) -> DNase_correct
DNase_correct <- DNase_correct %>% 
  arrange(seqnames)

Gene_DNase_Chrom <- DNase_correct

save( Gene_DNase_Chrom, file = "/home/dankent/Data/GM12878/Chromatin State/Chrom_DNase_per_gene_for_PCA.RData")

library( dplyr)
library( tibble)
library( data.table)
library( GenomicRanges)

Genes       <- read.table( "/home/dankent/Data/mart_export (8).txt", header = TRUE, sep = "\t")
Genes$Chromosome.scaffold.name <- paste0( 'chr',Genes$Chromosome.scaffold.name)
#remove duplicate of ENSG00000109072 (appears as different ID)
Genes <- Genes[-23323,]
Genes       <- subset( Genes, Genes$Gene.type=="protein_coding")
Genes       <- rowid_to_column( Genes, "ID")
Genes$width <- Genes$Gene.end..bp.- Genes$Gene.start..bp.
ChromState  <- read.table( "/home/dankent/Data/GM12878/Chromatin State/Segmentation_GM12878/GM12878_11_segments.bed", header = F)

# Function to find overlaps of individual chromatin states and % coverage on the fragments and join them in a df
MetaMerge <- function( Genes, ChromState) {
  
  test <- makeGRangesFromDataFrame( Genes, keep.extra.columns = T, start.field = "Gene.start..bp.", end.field = "Gene.end..bp.", seqnames.field = "Chromosome.scaffold.name")
  test1 <- makeGRangesFromDataFrame( ChromState, keep.extra.columns = T, start.field = "V2", end.field = "V3", seqnames.field = "V1")
  
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
Genes_Chrom <- MetaMerge( Genes, ChromState)

# Convert chromatin states to columns with percentageOverlaps as the values. 
x=as.character( paste( "seqnames+start+end+width+TSS+ID+Gene.stable.ID+HGNC.symbol+GO.term"))
f                         <- as.formula( paste( paste( x, collapse = " + "), "~ CS"))
colnames( Genes_Chrom) <-  c( "seqnames", "start", "end", "width", "strand","ID", "TSS", "HGNC.symbol", "Gene type", "Gene.stable.ID", "GO.accenssion","GO.term", "CS", "percentageOverlap")
Genes_Chrom %>% reshape2::dcast( f, value.var = "percentageOverlap", fun.aggregate = sum) -> Gene_makeup_CS
Gene_makeup_CS <- arrange(Gene_makeup_CS, seqnames, .by_group = TRUE)

save( Gene_makeup_CS, Genes, ChromState, file = "/home/dankent/Data/GM12878/Chromatin State/Chrom_state_per_gene_for_PCA.RData")

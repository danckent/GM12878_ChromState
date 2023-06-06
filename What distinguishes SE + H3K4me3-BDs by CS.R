library(dplyr)
library(factoextra)

ChromState  <- read.table( "/home/dankent/Data/GM12878/Chromatin State/Segmentation_GM12878/GM12878_11_segments.bed", header = F)

# Function to find overlaps of individual chromatin states and % coverage on the fragments and join them in a df
MetaMerge <- function( test, ChromState) {
  
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

Hybrid_Chrom <- MetaMerge( test = Hybrid_region_GR, ChromState)

# Convert chromatin states to columns with percentageOverlaps as the values. 
x=as.character( paste( "seqnames+start+end+width"))
f                         <- as.formula( paste( paste( x, collapse = " + "), "~ CS"))
colnames( Hybrid_Chrom) <-  c( "seqnames", "start", "end", "width", "strand","width", "ID", "CS", "percentageOverlap")
Hybrid_Chrom %>% reshape2::dcast( f, value.var = "percentageOverlap", fun.aggregate = sum) -> Hybrid_Chrom
Hybrid_Chrom <- arrange(Hybrid_Chrom, seqnames, .by_group = TRUE)
Hybrid_Chrom <- Hybrid_Chrom %>% 
  mutate(Category = "Hybrid")  %>% 
  mutate(E5 = 0) %>% 
  select("seqnames", "start", "end", "width", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "Category")

#### SE 

SE_Chrom <- MetaMerge( test = SE_GM12878[,1], ChromState)

# Convert chromatin states to columns with percentageOverlaps as the values. 
x=as.character( paste( "seqnames+start+end+width"))
f                         <- as.formula( paste( paste( x, collapse = " + "), "~ CS"))
colnames( SE_Chrom) <-  c( "seqnames", "start", "end", "width", "strand", "ID", "CS", "percentageOverlap")
SE_Chrom %>% reshape2::dcast( f, value.var = "percentageOverlap", fun.aggregate = sum) -> SE_Chrom
SE_Chrom <- arrange(SE_Chrom, seqnames, .by_group = TRUE)
SE_Chrom <- SE_Chrom %>% 
  mutate(Category = "SE") %>% 
  select("seqnames", "start", "end", "width", "E1", "E2", "E3", "E4", "E6", "E7", "E8", "E9", "E10", "E11", "Category")
SE.pca <- prcomp(SE_Chrom[ ,c( 5:14)], scale. = TRUE)
fviz_eig(SE.pca, addlabels = TRUE, ylim = c(0, 25), linecolor = "orange", ncp = 11)

fviz_pca_biplot(SE.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                col.var = "orange",
                pointshape = 21,
                pointsize = 3,
                col.ind = SE_Chrom$width,
                repel = TRUE)

#### BD 

BD_Chrom <- MetaMerge( test = BD_GM12878_GR[,1], ChromState)

# Convert chromatin states to columns with percentageOverlaps as the values. 
x=as.character( paste( "seqnames+start+end+width"))
f                         <- as.formula( paste( paste( x, collapse = " + "), "~ CS"))
colnames( BD_Chrom) <-  c( "seqnames", "start", "end", "width", "strand", "ID", "CS", "percentageOverlap")
BD_Chrom %>% reshape2::dcast( f, value.var = "percentageOverlap", fun.aggregate = sum) -> BD_Chrom
BD_Chrom <- arrange(BD_Chrom, seqnames, .by_group = TRUE)
BD_Chrom <- BD_Chrom %>% 
  mutate(Category = "BD") %>% 
  select("seqnames", "start", "end", "width", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "Category")
BD.pca <- prcomp(data.frame(BD_Chrom[ ,c( 5:15)]), scale. = TRUE) 
fviz_eig(BD.pca, addlabels = TRUE, ylim = c(0, 25), linecolor = "orange", ncp = 11)

fviz_pca_biplot(BD.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                col.var = "orange",
                pointshape = 21,
                pointsize = 3,
                col.ind = BD_Chrom$width,
                repel = TRUE)

get_eigenvalue(BD.pca)


CA_BD <- FactoMineR::CA(data.frame(BD_Chrom[ , 5:15])) 
fviz_ca_biplot(CA_BD,
               select.row = list(contrib = 20),
               repel = TRUE) 

### combine SE + BD ####

SE_Chrom <- SE_Chrom %>% 
  mutate(Category = "SE") %>% 
  mutate(E5 = 0) %>% 
  select("seqnames", "start", "end", "width", "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "Category")

comb <- rbind(SE_Chrom, BD_Chrom, Hybrid_Chrom)

ggplot(comb, aes(x=E5+E6, y=E11+E10, colour = Category)) +
  geom_point() 

comb.pca <- prcomp(data.frame(comb[ ,c( 5:15)]), scale. = TRUE)
fviz_pca_biplot(comb.pca, 
                # Fill individuals by groups
                geom.ind = "point",
                col.var = "orange",
                pointshape = 21,
                pointsize = 3,
                col.ind = comb$Category,
                repel = TRUE)

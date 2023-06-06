library(dplyr)
library(regioneR)

Chrom_state_GM12878_BD <- read.table("Segmentation_GM12878/E10_E11_chromHMM_BD_GM12878_size_filtered.bed", header = TRUE)
Chrom_state_GM12878_BD <- makeGRangesFromDataFrame(Chrom_state_GM12878_BD, start.field = "V2", end.field = "V3", seqnames.field = "V1")
seqlevelsStyle(Chrom_state_GM12878_BD)      <- "UCSC"
# Only ranges in autosomes
Chrom_state_GM12878_BD <- filterChromosomes(Chrom_state_GM12878_BD, chr.type = "autosomal")

Chrom_state_GM12878_SE <- read.table("Segmentation_GM12878/E9_chromHMM_SE_GM12878_size_filtered.bed", header = TRUE)
Chrom_state_GM12878_SE <- makeGRangesFromDataFrame(Chrom_state_GM12878_SE, start.field = "V2", end.field = "V3", seqnames.field = "V1")
seqlevelsStyle(Chrom_state_GM12878_SE)      <- "UCSC"
# Only ranges in autosomes
Chrom_state_GM12878_SE <- filterChromosomes(Chrom_state_GM12878_SE, chr.type = "autosomal")

#### SE distance to BD #####
chrom_SE_Distance_to_BD <- distanceToNearest(Chrom_state_GM12878_BD, Chrom_state_GM12878_SE)
chrom_SE_Distance_to_BD <- data.frame(chrom_SE_Distance_to_BD)
chrom_SE_Distance_to_BD_10kb <- chrom_SE_Distance_to_BD[chrom_SE_Distance_to_BD$distance<=10000,]

#combine df to show which match which 
chrom_Combined_nearest_BD_and_SE <- cbind(data.frame(Chrom_state_GM12878_BD[chrom_SE_Distance_to_BD_10kb$queryHits]),data.frame(Chrom_state_GM12878_SE[chrom_SE_Distance_to_BD_10kb$subjectHits]), chrom_SE_Distance_to_BD_10kb$distance)

#clean df 
chrom_Combined_nearest_BD_and_SE <- chrom_Combined_nearest_BD_and_SE[c(1,2,3,4,7,8,11)]
colnames(chrom_Combined_nearest_BD_and_SE) <- c("chr", "startBD", "endBD", "widthBD", "startSE", "endSE", "distance")
chrom_Combined_nearest_BD_and_SE_GR <- makeGRangesFromDataFrame(chrom_Combined_nearest_BD_and_SE, keep.extra.columns = TRUE, start.field = "startBD", end.field = "endBD")



#Open the Ensebl Gene list in RStudio
Gene_List_Ensembl <- read.table( "/home/dankent/Data/HG_Esembl98_Genes.tab", sep = "\t", header = TRUE)
#Add 'chr' into the column containing the chromosome notation
Gene_List_Ensembl$Chromosome.scaffold.name <- paste0( 'chr',Gene_List_Ensembl$Chromosome.scaffold.name) 

#custom order
Gene_List_Ensembl$Chromosome.scaffold.name <- factor( Gene_List_Ensembl$Chromosome.scaffold.name, ordered = TRUE, levels =  c( "chr1","chr2", "chr3", 
                                                                                                                               "chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22", 
                                                                                                                               "chrX", "chrY"))
Gene_List_Ensembl <- Gene_List_Ensembl[ order( Gene_List_Ensembl$Chromosome.scaffold.name), ]

#check for NAs and omit based on genome patch builds
Gene_List_Ensembl[ !complete.cases( Gene_List_Ensembl),]
Gene_List_Ensembl <- na.omit( Gene_List_Ensembl)

#duplicate start and end columns for inclusion into metadata 
Gene_List_Ensembl$startgene = Gene_List_Ensembl$Gene.start..bp.
Gene_List_Ensembl$endgene = Gene_List_Ensembl$Gene.end..bp.

#Add unique ID to each row
Gene_List_Ensembl <- tibble::rowid_to_column(Gene_List_Ensembl, "ID")

colnames(Gene_List_Ensembl) <- c("ID", "Gene Stable ID", "HGNC Symbol", "Chromosome", "Start", "End", "Type")

#Make genomic ranges objects for the data frame
GR_Genes <- makeGRangesFromDataFrame( df = Gene_List_Ensembl, keep.extra.columns = TRUE)


chrom_Genes <- find_overlaps(chrom_Combined_nearest_BD_and_SE_GR, GR_Genes)
chrom_Genes <- data.frame(chrom_Genes)
match2 <- chrom_Genes$`Gene Stable ID`

find_overlaps(chrom_Combined_nearest_BD_and_SE_GR, Combined_nearest_BD_and_SE_GR)
find_overlaps(Combined_nearest_BD_and_SE_GR, chrom_Combined_nearest_BD_and_SE_GR)
#only those that overlap 
hgnc <- chrom_Genes[c(2,3,5,7,8,9),]


join_overlap_inner(chrom_Combined_nearest_BD_and_SE_GR, Combined_nearest_BD_and_SE_GR)


#Dell server

#sort
bedtools sort -chrThenSizeA -i /home/dankent/Data/GM12878/Segmentation_GM12878/GM12878_11_segments.bed

#subset data for BD and SE 
awk '/E9/ {print}' /home/dankent/Data/GM12878/Segmentation_GM12878/GM12878_11_segments.bed > Segmentation_GM12878/E9_chromHMM_SE.bed
awk '/c/ && $4 == "E10" || $4 == "E11" {print}' /home/dankent/Data/GM12878/Segmentation_GM12878/GM12878_11_segments.bed > Segmentation_GM12878/E10_E11_chromHMM_BD.bed

chmod u+xrw Segmentation_GM12878/E10_E11_chromHMM_BD.bed
chmod u+xrw Segmentation_GM12878/E9_chromHMM_SE.bed

#merge all regions within 1000bp of one another 
bedtools merge -d 1000 -i /home/dankent/Data/GM12878/Segmentation_GM12878/E10_E11_chromHMM_BD.bed > Segmentation_GM12878/E10_E11_Merged.bed
bedtools merge -d 7500 -i /home/dankent/Data/GM12878/Segmentation_GM12878/E9_chromHMM_SE.bed > Segmentation_GM12878/E9_Merged.bed

Rscript --vanilla /home/dankent/Data/GM12878/Segmentation_GM12878/size_filter.R



#HP server
#make new temp folder outside my home. This is the only way I could get it to work. 

mkdir /dankent_tmp

java -jar /progs/ChromHMM-1.18/ChromHMM.jar BinarizeBam -center /progs/ChromHMM-1.18/CHROMSIZES/hg19.txt /dankent_tmp /dankent_tmp/myDataFile_GM12878-c.txt /dankent_tmp

# create the folder for segmentaiton file
mkdir /dankent_tmp/Segmentation_GM12878

# do segmentation 
java -jar /progs/ChromHMM-1.18/ChromHMM.jar MakeSegmentation -printposterior -printstatesbyline /dankent_tmp/model_11_All_cell_types.txt /dankent_tmp /dankent_tmp/Segmentation_GM12878

#make browser files
java -jar /progs/ChromHMM-1.18/ChromHMM.jar MakeBrowserFiles /dankent_tmp/Segmentation_GM12878 GM12878 /dankent_tmp






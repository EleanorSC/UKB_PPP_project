#16-11-23

setwd("/Scratch_Area/Sarah/UKBB_proteomics_3K")

#read in proteomics data
d1=read.table ("olink_data.txt", header = TRUE)
dim (d1)

#keep ins_index=0 (baseline data)
d2=d1[d1$ins_index=="0",]
dim (d2)

#remove ins_index

d3 = subset(d2, select = -c(ins_index) )
dim (d3)

#convert long format to wide format
d4=reshape (d3, idvar="eid", timevar="protein_id", direction="wide")

dim (d4)
#[1] 53029  2924

write.table (d4, file="UKBB_proteomics_baseline", row.names=F)


#30-1-24

setwd("/Scratch_Area/Sarah/UKBB_proteomics_3K")

#read in proteomics baseline data
d4=read.table ("UKBB_proteomics_baseline", header =TRUE)
dim (d4)
#[1] 53029  2924

#read in meta data (No. of proteins measured, Plate, Well, UKB-PPP consortium selected participant)

d5=read.table ("/Scratch_Area/Sarah/UKBB_proteomics_3K/UKB_PROTEOMICS_FOR_SARAH_JAN2024_GD.txt", header = TRUE)
dim (d5)
#[1] 502356     12

#keep baseline metadata
d6=subset(d5, select = c(eid, X30900.0.0, X30901.0.0, X30902.0.0, X30903.0.0))
dim (d6)
#[1] 502356     5 

#replace NAs with 0 to indicate non-consortium participants
d7= replace (d6, is.na(d6),0)

#merge proteomics data with meta data 
d8=merge(d7, d4, by.x="eid", by.y="eid")
dim (d8)
#[1] 53029  2928

#read in batch data
d9=read.table ("olink_batch_number.dat.txt", header = TRUE)

#keep batches 1-6
d10=d9[d9$Batch >0 & d9$Batch < 7,]
dim (d10)
#[1] 536   2

#merge batch number with proteomics data keeping only batches 1-6 
d11=merge (d10, d8, by.x="PlateID", by.y="X30901.0.0")
dim (d11)
#[1] 45441  2929

#remove consortium selected individuals
d12=d11[d11$X30903.0.0 == 0,]
dim (d12)
#[1] 40441  2929

#keep individuals with measures for more than or equal to 90% of proteins
d13=d12[d12$X30900.0.0 >2631,]
dim (d13)
#[1] 36787  2929

#keep proteins with measures for more than or equal to 90% of individuals
miss = c()
for(i in 1:ncol(d13)) {
  if(length(which(is.na(d13[,i])))>0.1*nrow(d13)) miss = append(miss,i)
}
d14=d13[,-miss,]
dim (d14)
#[1] 36787  2923



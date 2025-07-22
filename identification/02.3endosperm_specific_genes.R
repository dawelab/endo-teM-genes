specific <- read.delim("~/Documents/lab/b73.all.6_38.txt", comment.char="#")

#calculate the max TPM in all the endosperm tissues
specific$maxendo<-apply(specific[,25:41],1,max)

#calculte the max TPM in all the sporophyte tissues.
specific$max<-apply(specific[,c(14,16:23)],1,max)

#calcilate the expression FC between endosperm and sporophyte tissue.
specific$FC<-(specific$maxendo+1)/(specific$max+1)

#select the endosperm specific genes using FC >=5 and endosperm TPM >=20
specific<-specific[which(specific$maxendo>=20&specific$FC>=5),]



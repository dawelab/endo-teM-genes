

specific<-b73.all.6_38
specific$maxendo<-apply(specific[,25:41],1,max)
specific$max<-apply(specific[,c(14,16:23)],1,max)
specific$FC<-(specific$maxendo+1)/(specific$max+1)
specific<-specific[which(specific$maxendo>=20&specific$FC>=5),]

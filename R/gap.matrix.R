gap.matrix <-
function(partition, data){
		nseq<-length(data$alignments[[partition]])
		nbp<-length(data$alignments[[partition]][[1]])
		gap.mat<-matrix(nrow=nseq, ncol=nbp)
		rownames(gap.mat)<-names(data$alignments[[partition]])
		for(i in 1:nseq){
			gap.mat[i,]<-!(data$alignments[[partition]][[i]]=="A"|data$alignments[[partition]][[i]]=="a"|data$alignments[[partition]][[i]]=="C"|data$alignments[[partition]][[i]]=="c"|data$alignments[[partition]][[i]]=="G"|data$alignments[[partition]][[i]]=="g"|data$alignments[[partition]][[i]]=="T"|data$alignments[[partition]][[i]]=="t")
			}		
		return(gap.mat)
		}

commandlinemaker <-
function(model, partition, data){
		nsamples<-length(data$log[,1])
		commandline<-rep("-or", times=nsamples)
		seqlen<-length(data$alignment[[partition]][[1]])
		commandline<-cbind(commandline, paste("-l", seqlen, sep=""))
		commandline<-cbind(commandline, "-mGTR")
		if(model$alpha){
			if(!(paste(partition, "alpha", sep=".")%in%names(data$log))){
				cat("WARNING:GENE TREE NAMES DO NOT CORRESPOND TO LOG FILE COLUMN NAMES\n")
				cat("the following names from the log file:\n", names(data$log)[grep(partition, names(data$log))], "\n", "do not match these names:\n", data$genes, "\nfrom the input data. change the log file names to correct the issue.\n")
				}
			commandline<-cbind(commandline, paste("-a", data$log[,paste(partition, "alpha", sep=".")], sep=""))
			p.alpha<-data$log[,paste(partition, "alpha", sep=".")]
			kg<-rep(4, times=nsamples)
			}else{
				p.alpha<-rep(1, times=nsamples)
				kg<-rep(1, times=nsamples)
				}
			
		if(model$pInv){
			commandline<-cbind(commandline, paste("-i", data$log[,paste(partition, "pInv", sep=".")], sep=""))
			p.inv<-data$log[,paste(partition, "pInv", sep=".")]
			}else{
				p.inv<-rep(0, times=nsamples)
				}

		if(model$basefreq=="estimated"){
			freqs<-cbind(data$log[,paste(partition, "frequencies1", sep=".")], data$log[,paste(partition, "frequencies2", sep=".")], data$log[,paste(partition, "frequencies3", sep=".")], data$log[,paste(partition, "frequencies4", sep=".")])
			p.freqs<-freqs
			freqs<-paste(freqs[,1], freqs[,2], freqs[,3], freqs[,4], sep=",")
			freqs2<-paste("-f", freqs, sep="")
			commandline<-cbind(commandline, freqs2)
			}
		if(model$basefreq=="fixed"){
			cat("for partition ", partition, "setting base frequencies as equal. support for empirical base frequencies not yet added.\n")
			freqs<-rep(0.25, times=nsamples)
			p.freqs<-cbind(freqs, freqs, freqs, freqs)
			freqs<-paste(freqs, freqs, freqs, freqs, sep=",")
			freqs2<-paste("-f", freqs, sep="")
			commandline<-cbind(commandline, freqs2)
			}


		if(model$rmat=="gtr"){
			rmatrix<-cbind(data$log[,paste(partition, "ac", sep=".")], data$log[,paste(partition, "ag", sep=".")], data$log[,paste(partition, "at", sep=".")], data$log[,paste(partition, "cg", sep=".")], 1, data$log[,paste(partition, "gt", sep=".")])
			p.rmat<-rmatrix
			rmatrix<-paste(rmatrix[,1],rmatrix[,2],rmatrix[,3],rmatrix[,4],rmatrix[,5],rmatrix[,6], sep=",")
			rmatrix2<-paste("-r", rmatrix, sep="")
			commandline<-cbind(commandline, rmatrix2)
			}
		if(model$rmat=="hky"){
			rmatrix<-cbind(1, data$log[,paste(partition, "kappa", sep=".")], 1, 1, data$log[,paste(partition, "kappa", sep=".")], 1)
			p.rmat<-rmatrix
			rmatrix<-paste(rmatrix[,1],rmatrix[,2],rmatrix[,3],rmatrix[,4],rmatrix[,5],rmatrix[,6], sep=",")
			rmatrix2<-paste("-r", rmatrix, sep="")
			commandline<-cbind(commandline, rmatrix2)			
			}
		if(model$rmat=="trn"){
			rmatrix<-cbind(1, data$log[,paste(partition, "kappa1", sep=".")], 1, 1, data$log[,paste(partition, "kappa2", sep=".")], 1)
			p.rmat<-rmatrix
			rmatrix<-paste(rmatrix[,1],rmatrix[,2],rmatrix[,3],rmatrix[,4],rmatrix[,5],rmatrix[,6], sep=",")
			rmatrix2<-paste("-r", rmatrix, sep="")
			commandline<-cbind(commandline, rmatrix2)
			}
	colnames(commandline)<-NULL
	out<-list("seqgen"=commandline, "phangorn"=list("p.rmat"=p.rmat, "p.freqs"=p.freqs, "p.inv"=p.inv, "p.alpha"=p.alpha, "kg"=kg))
	return(out)
	}

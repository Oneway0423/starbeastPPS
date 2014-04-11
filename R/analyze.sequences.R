analyze.sequences <-
function(data, seqgendir="/directory/containing/seqgen"){
		nsamples<-length(data$log[,1])
		nparts<-length(data$genes.to.partitions[,1])
		loci<-data$genes
		parts<-data$genes.to.partitions
		
		rescaler<-function(tree, rate=1){
			tree$edge.length<-tree$edge.length*tree$rate[,rate]
			return(tree)
			}

		###create gap/ambiguity matrices
		gaps<-list()
		for(i in 1:nparts){
			gaps[[parts[i,2]]]<-gap.matrix(parts[i,2], data)
			}
		###get gap/ambiguity columns for exclusion. 
		gap.cols<-list()
		for(i in 1:nparts){
			gap.cols[[parts[i,2]]]<-gap.columns(gaps[[parts[i,2]]])
			}
		
		cat("gaps done\n")
		
		###get alignment order
		order<-list()
		for(i in 1:nparts){
			order[[i]]<-names(data$alignment[[i]])
			}
		
		###generate seqgen and phangorn inputs
		models<-list()
		commands<-list()
		for(i in 1:nparts){
			models[[i]] <- get.model(data, parts[i,2])
			commands[[parts[i,2]]] <- commandlinemaker(models[[i]], partition=parts[i,2], data)
			}
		cat("commands done\n")


		###rescale time trees to substitution trees		
		partnumbers<-1
		for(i in 2:nparts){
			if(parts[(i-1),1]==parts[i,1]){
				partnumbers<-c(partnumbers, (partnumbers[i-1]+1))
				}else{partnumbers<-c(partnumbers, 1)}
			}
		parts<-cbind(parts, partnumbers)
		empirical.trees<-list()
		for(i in 1:nparts){
			empirical.trees[[parts[i,2]]]<-data$gene.trees[[parts[i,1]]]
			empirical.trees[[parts[i,2]]]<-lapply(empirical.trees[[parts[i,2]]], rescaler, rate=as.numeric(parts[i,3]))
			cat("trees for partition", parts[i,1:2],"rescaled\n")
			}

		###calculate empirical likelihoods
		empirical.lik<-matrix(nrow=nsamples, ncol=nparts, dimnames=list(c(), parts[,2]))
		empirical.unconst<-matrix(nrow=nsamples, ncol=nparts, dimnames=list(c(), parts[,2]))
		empirical.sites<-vector(length=nparts)
		for(i in 1:nparts){
			cat("Excluding ambiguous sites and using ", sum(gap.cols[[parts[i,2]]])/length(gap.cols[[parts[i,2]]]), " of alignment for calculations for locus ", parts[i,2], ".\n")
			empirical.dat<-phyDat(data$alignments[[parts[i,2]]])
			empirical.dat<-as.character(empirical.dat)
			empirical.dat<-empirical.dat[,gap.cols[[parts[i,2]]]]
			empirical.dat<-phyDat(empirical.dat)
			for(j in 1:nsamples){
				pml.out<-pml(tree=empirical.trees[[parts[i,2]]][[j]], data=empirical.dat, bf=commands[[parts[i,2]]]$phangorn$p.freqs[j,], Q=commands[[parts[i,2]]]$phangorn$p.rmat[j,], inv=commands[[parts[i,2]]]$phangorn$p.inv[j], k=commands[[parts[i,2]]]$phangorn$kg[j], shape=commands[[parts[i,2]]]$phangorn$p.alpha[j])
				empirical.lik[j,i]<-pml.out$logLik
				w<-pml.out$weight
				w<-w[w > 0]
				ll0 = sum(w * log(w/sum(w)))
				empirical.unconst[j,i]<-ll0
				}
			empirical.sites[i]<-var.sites(empirical.dat)
			cat("empirical probabilities for ", parts[i,2], " are done\n")
			}
		difference.emp<-(empirical.unconst-empirical.lik)
		cat("empirical probabilities are done\n")

		###simulate new sequences, reorder them, recreate gaps, and calculate likelihoods
		simulated.lik<-matrix(nrow=nsamples, ncol=nparts, dimnames=list(c(), parts[,2]))
		simulated.unconst<-matrix(nrow=nsamples, ncol=nparts, dimnames=list(c(), parts[,2]))
		simulated.sites<-matrix(nrow=nsamples, ncol=nparts, dimnames=list(c(), parts[,2]))
		for(i in 1:nparts){
			for(j in 1:nsamples){
				seq.1<-seqgenrunner(commands[[i]]$seqgen[j,], empirical.trees[[parts[i,2]]][[j]], seqgendir)
				seq.1<-read.phylip.seqgen(seq.1)
				seq<-seq.1$org.code
				rownames(seq)<-seq.1$seqname
				seq<-seq[order[[i]],]
				seq[gaps[[i]]]<-"-"
				seq<-seq[,gap.cols[[parts[i,2]]]]
				simulated.dat<-phyDat(seq)
				pml.out<-pml(tree=empirical.trees[[parts[i,2]]][[j]], data=simulated.dat, bf=commands[[parts[i,2]]]$phangorn$p.freqs[j,], Q=commands[[parts[i,2]]]$phangorn$p.rmat[j,], inv=commands[[parts[i,2]]]$phangorn$p.inv[j], k=commands[[parts[i,2]]]$phangorn$kg[j], shape=commands[[parts[i,2]]]$phangorn$p.alpha[j])
				simulated.lik[j,i]<-pml.out$logLik
				w<-pml.out$weight
				w<-w[w > 0]
				ll0 = sum(w * log(w/sum(w)))
				simulated.unconst[j,i]<-ll0
				simulated.sites[j,i]<-var.sites(simulated.dat)
				}
			cat("simulated probabilities for ", parts[i,2], " are done\n")
			}
		difference.sim<-(simulated.unconst-simulated.lik)
		cat("simulated probabilities are done\n")
		
		sitesdist<-simulated.sites
		for(i in 1:length(sitesdist[,1])){sitesdist[i,]<-sitesdist[i,]-empirical.sites}

		empirical<-list("phy.likelihood"=empirical.lik, "multi.likelihood"=empirical.unconst, "gcstat"=difference.emp, "v.sites"=empirical.sites)
		simulated<-list("phy.likelihood"=simulated.lik, "multi.likelihood"=simulated.unconst, "gcstat"=difference.sim, "v.sites"=simulated.sites)
		difference<-list("phy.likelihood"=simulated.lik-empirical.lik, "multi.likelihood"=simulated.unconst-empirical.unconst, "gcstat"=difference.sim-difference.emp, "v.sites"=sitesdist)
		output<-list("empirical"=empirical, "simulated"=simulated, "test.stat"=difference)
		class(output)<-"sequenceteststats"
		return(output)
		}

read.sptree.samples <-
function(file){
	#####parts of this function were adapted from a function in package 'phyloch', written by Christoph Heibl
	#####read in trees without stats
	tree<-read.nexus(file)
	many<-class(tree)=="multiPhylo"
	
	#####get tree strings
    X <- scan(file = file, what = "", sep = "\n", quiet = TRUE)	
    X <- X[grep("tree STATE", X)]					  
    X <- gsub("tree STATE_[[:digit:]]+[[:space:]]\\[[^]]*\\] = \\[&R\\] ", "", X)		

	######a function to produce branch width vectors and order them appropriately

	extract.order.stats<-function(treestring){
		
		ntax<-length(unlist(strsplit(treestring, "\\[")))/2
		nodes<-(ntax+1):((2*ntax)-1)
		Y<-treestring

		####adds internal node labels
		for(i in 1:length(nodes)){
			repl<-paste(")", nodes[i], "[", sep="")
			Y<-sub(")\\[", repl, Y)
			}	
		
		meta<-unlist(strsplit(Y, "\\[|\\]"))[grep("dmv", unlist(strsplit(Y, "\\[|\\]")))]
		metacols<-length(unlist(strsplit(meta[1], ",")))
		meta<-gsub("&|dm.=|\\{|\\}", "", meta)

		Ysub<-gsub("\\[[^]]*\\]", "\\[\\]", Y)
		Ysub<-unlist(strsplit(Ysub, ",|)"))
		Ysub<-gsub("\\(|\\)|;|\\[|\\]", "", Ysub)
		Ysub[length(Ysub)]<-paste(Ysub[length(Ysub)], NA, sep=":")

		branchdata<-array(dim=c(length(meta), 2+metacols))

		for(i in 1:length(meta)){
			branchdata[i,]<-c(unlist(strsplit(Ysub[i], ":")), unlist(strsplit(meta[i], ",")))
			}

		if(metacols==3){colnames(branchdata)<-c("br", "length", "dmt", "dmv_start", "dmv_end")}
		if(metacols==1){colnames(branchdata)<-c("br", "length", "dmv")}
		rownames(branchdata)<-branchdata[,1]

		string <- gsub("\\[[^]]*\\]", "", Y)
		stree<-read.tree(text=string)
		translate<-cbind(stree$node.label, (ntax+1):nodes[length(nodes)])
		translate<-rbind(translate, cbind(stree$tip.label, 1:ntax))
		rownames(translate)<-translate[,2]
		translate2<-translate[as.character(stree$edge[,2]),]
		branchdata2<-branchdata[translate2[,1],]

		branchdata2<-rbind(branchdata2, branchdata[length(branchdata[,2]),])
		
		rownames(branchdata2)<-NULL
		return(branchdata2)
		
		
		}

	branchdata<-lapply(X, extract.order.stats)

	#######append branch width stats to trees.  vector "check" is included as a verification 
	#######and should exactly match vector "edge.length"
	#######this is the LONGEST part of the function.  can I do this faster???!?
	if(length(branchdata[[1]][1,])==5){	
		if(class(tree)=="multiPhylo"){
		for(i in 1:length(tree)){
#			tree[[i]][["check"]]<-as.numeric(branchdata[[i]][,2])
			tree[[i]][["dmt"]]<-as.numeric(branchdata[[i]][,3])
			tree[[i]][["dmv_start"]]<-as.numeric(branchdata[[i]][,4])
			tree[[i]][["dmv_end"]]<-as.numeric(branchdata[[i]][,5])
			}
			}
		if(class(tree)=="phylo"){
#			tree[["check"]]<-as.numeric(branchdata[[1]][,2])
			tree[["dmt"]]<-as.numeric(branchdata[[1]][,3])
			tree[["dmv_start"]]<-as.numeric(branchdata[[1]][,4])
			tree[["dmv_end"]]<-as.numeric(branchdata[[1]][,5])
			}
		}
	if(length(branchdata[[1]][1,])==3){	
		if(class(tree)=="multiPhylo"){
		for(i in 1:length(tree)){
#			tree[[i]][["check"]]<-as.numeric(branchdata[[i]][,2])
			tree[[i]][["dmv"]]<-as.numeric(branchdata[[i]][,3])
			}
			}
		if(class(tree)=="phylo"){
#			tree[["check"]]<-as.numeric(branchdata[[1]][,2])
			tree[["dmv"]]<-as.numeric(branchdata[[1]][,3])
			}
		
		}
	return(tree)
	}

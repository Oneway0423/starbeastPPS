mstree <-
function(phy, msdir, nseq, nreps, samplescheme, ploidy=1){
	
	#####phy is the population tree in which gene trees are to be simulated. IN THIS SCRIPT, ALL TREES MUST BE SPECIES TREES WITH ASSOCIATED 'DMV' BRANCH WIDTH VALUES
	#####msdir is the path to enclosing directory of ms.  ms should be modified to increase the number of digits in the branch length calculations.  
	#####nseq is the total number of tips to be simulated
	#####nreps is the number of trees to generate given the tree and sampling scheme
	#####samplescheme is a dataframe with two columns and n rows.  The first column contains the tip labels of the population tree and the second contains the number of alleles to be sampled from that population
	#####ploidy refers to *BEAST's "ploidy" in the xml files and modifies the DMV values. When all loci have the same ploidy, it should be left as 1. When ploidy varies, it should be 0.5 for mitochondrial and 2 for diploid nuclear. 
	
	#####gets tip labels for descendents of node in phy
	node.tips<-function (phy, node){
 		n<- length(phy$tip.label)
 		if (node<= n){node}else{
 			l<- numeric()
 			d<- phy$edge[which(phy$edge[,1]==node),2]
 			for(j in d){if(j<= n){l<- c(l, j)}else{l<-c(l, node.tips(phy,j))}}
         	l}}
	
	nspec<-length(phy$tip.label)

	values<-phy$dmv*2*ploidy	###dmv values are now converted to theta=4Nu. dividing branches by theta yields species tree branch lengths in 4N generations, which ms uses for simulation. 
	values<-values  ####as on the next line, "dmv" is Nu haploid alleles, or 2Nu individuals.  Therefore, multiply by 2 to get the standard 4Nu population scaled mutation rate.  

	branchscalar<-values[length(values)]
	phy$edge.length<-phy$edge.length/branchscalar ###transforms branch lengths by theta.  ##BEAST's dmv is equal to Nu where N is the pop size in alleles.  in other words 2Nu in diploid individuals.  so multiply dmv by 2 to get branch lengths in 4N generations for ms.  
	values<-values/values[length(values)]	####dmvs should are now a fraction of the root dmv.
	sort(branching.times(phy))->bt;
	scaleindex<-c(phy$edge[,2], (nspec+1))

	####need to set initial pop sizes -I npops samplescheme[,2] <-n .... -n...> comlinetree etc...
	#### -n pop scale<size=(scale*No)>
	###-en time pop scale<newsize=(scale*No)>
	
	initialpops<-c()
	for(i in 1:nspec){
		initialpops<-c(initialpops, "-n", i, values[which(phy$edge[,2]==which(phy$tip.label==samplescheme[i,1]))])
		
		}
	
	
	comlinetree<-c();
	for(i in 1:length(bt)){
		child <- sort(phy$edge [phy$edge[,1] == names(bt[i]), 2])
		tips<-sort( c( sort(node.tips(phy, child[1]))[1], sort(node.tips(phy, child[2]))[1] ))
		popi<-which(samplescheme[,1]==phy$tip.label[tips[2]])
		popj<-which(samplescheme[,1]==phy$tip.label[tips[1]])
		scale<-values[which(scaleindex==names(bt[i]))]
		
		comlinetree<-c(comlinetree, "-ej", bt[i], popi, popj, "-en", bt[i], popj, scale)	
		}
	names(comlinetree)<-NULL
	
	commandline<-c("cd", msdir, ";", "./ms", nseq, nreps, "-T", "-I", length(phy$tip.label), samplescheme[,2], initialpops, comlinetree, "| grep \\;")
	junk<-system(paste(unlist(commandline), collapse = " "), intern=TRUE)
	read.tree(text=junk)->trees
	
	if(nreps>1){
		for(i in 1:nreps){
			trees[[i]]$edge.length<-trees[[i]]$edge.length*branchscalar
			}
		}else{trees$edge.length<-trees$edge.length*branchscalar}
	return(trees)
#	return(paste(unlist(commandline), collapse = " "))	
	}

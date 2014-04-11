deep.coal <-
function(sptree, gtree, association){
	
	species<-sptree$tip.label
	nspec<-length(species)
	nbranch<-(2*nspec)-1

	tiplist<-list()
	for(i in 1:nspec){
		tiplist[[species[i]]]<-association[which(association[,1]==species[i]),2]
		}

	####this subfunction returns the tip numbers for a given node
	node.tips<-function (phy, node){
 		n<- length(phy$tip.label)
 		if (node<= n){node}else{
 			l<- numeric()
 			d<- phy$edge[which(phy$edge[,1]==node),2]
 			for(j in d){if(j<= n){l<- c(l, j)}else{l<-c(l, node.tips(phy,j))}}
       		l}}
	#####end subfunction
	####subfunction that gets a tree containing descendants of all gene lineages that pass through a branch 
	######do I really need that node conditional??? or should I just deal with it later?
	get.lin<-function(sptree, gtree, association, node){	
		species<-sptree$tip.label
		if((length(species)+1)!=node){
			spec.desc<-species[node.tips(sptree, node)]
			gtree.desc<-c()
			for(i in 1:length(spec.desc)){
				gtree.desc<-c(gtree.desc, association[which(association[,1]==spec.desc[i]),2])
				}
			match(gtree.desc, gtree$tip.label)->index
			pruned<-drop.tip(gtree, gtree$tip.label[-index])
			}else{pruned<-gtree}		
		return(pruned)
		}
	#####end subfunction
	#####prepare branch demographic info...
	demo.maker<-function(sptree, nspec){
		demo<-cbind(sptree$edge[,2], sptree$edge.length)
		demo<-rbind(demo, c((nspec+1), Inf))
		demo<-cbind(demo, sptree$dmv)
		demo<-demo[order(demo[,1]),]
		sbt<-branching.times(sptree)
		sbt<-sbt[order(as.numeric(names(sbt)))]
		sbt<-c(rep(0, nspec), sbt)
		demo<-cbind(demo, sbt)
		colnames(demo)<-c("node", "length", "dmv", "sbt")
		rownames(demo)<-c(1:length(sbt))
		return(demo)
		}
	#####

	#####now how to calculate likelihoods!??!?!
		####for a single branch...
	
	b.prob<-function(sptree, gtree, demo, node){
		pruned<-get.lin(sptree, gtree, association, node)
		gbt<-sort(branching.times(pruned))
		gbt<-c(0, gbt)
		gbt<-cbind(gbt, length(gbt):1)
		start<-demo[node,"sbt"]
		end<-(demo[node, "sbt"]+demo[node, "length"])
		enter<-gbt[gbt[,1]==max(gbt[gbt[,1]<=start,1]),2]
		exit<-gbt[gbt[,1]==max(gbt[gbt[,1]<=end,1]),2]
		
		exit<-exit-1
		#cat(exit, "\n")###test
		return(exit)
		}

	######now actually calculate the whole gene tree probability
		####remember to exclude species tip branches with only one allele sampled
	
	demo<-demo.maker(sptree, nspec)
	numtips<-lapply(tiplist, length)
	if(any(numtips==1)){
		demo2<-demo[-which(numtips==1),]}	
		else{demo2<-demo}
	
	
	deep.coalescences<-c()
	nodes<-demo2[,"node"]
	names(nodes)<-NULL
	for(i in 1:length(nodes)){
		deep.coalescences<-c(deep.coalescences, b.prob(sptree, gtree, demo, nodes[i]))
		}
	
	return((sum(deep.coalescences)))
	}

gap.columns <-
function(gapmatrix){
		dims<-dim(gapmatrix)
		rowSums(gapmatrix)->nonmissing.rows
		nonmissing.rows<-nonmissing.rows<dims[2]
		gap.mat<-gapmatrix[nonmissing.rows,]
		colSums(gap.mat)==0->gap.cols
		return(gap.cols)
		}

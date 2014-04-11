var.sites <-
function(phydat){
		mat<-as.character(phydat)
		dims<-dim(mat)
		mat<-mat[rowSums(mat=="-"|mat=="?")!=dims[2],]
		uniqEle<-function(vec){
			uni<-sum(!duplicated(vec))
			uni
			}
		out<-apply(mat, MARGIN=2, uniqEle)
		var.s<-length(out[out>1])
		var.s
		}

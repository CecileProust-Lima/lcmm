#ex: X 2 modalites
#type=0 -----> factor(X)=X
#type=1 ------>factor(X)=X.1  
#
factor.names <- function(x,type){
	if(type==1){
		x <- matrix(x,nrow=1)
		Names <- apply(x,MARGIN=2,FUN=function(x){
			if(length(grep("factor",x))!= 0){
				pos1 <- grep("\\(",unlist(strsplit(x,split="")))+1
				pos2 <- grep("\\)",unlist(strsplit(x,split="")))-1
				compris.factor <- substr(x,start=pos1,stop=pos2)
				after.factor <- substr(x,start=(pos2+2),stop=length(unlist(strsplit(x,split=""))))
				paste(compris.factor,after.factor,sep=".")
			}else{
				x
			}
		}
		)
		
	}else{
		x <- matrix(x,nrow=1)
		Names <- apply(x,MARGIN=2,FUN=function(x){
			if(length(grep("factor",x))!= 0){
				pos1 <- grep("\\(",unlist(strsplit(x,split="")))+1
				pos2 <- grep("\\)",unlist(strsplit(x,split="")))-1
				compris.factor <- substr(x,start=pos1,stop=pos2)
				compris.factor
			}else{
				x
			}
		}
		)		
	}

#	Names <- interraction(Names,1)
#	Names <- interraction(Names,2)
	
	return(Names)
}
	
#ex: interraction A:B
#type=1 -----> :
#type=2 ------> *  
# #
# interraction <- function(x,type){
# 	if(type==1){
# 		# interraction .
# 		Names <- NULL
# 		for(i in 1:length(x)){
# 			pos <- grep(":",unlist(strsplit(x[i],split="")))
# 			if(length(pos)>0){
# 				xgauche <- substr(x[i],start=1,stop=(pos-1))
# 				xdroite <- substr(x[i],start=(pos+1),stop=length(unlist(strsplit(x[i],split=""))))
# 				out <- c(xgauche,xdroite)
# 			}else{
# 				out <- x[i]
# 			}
# 			Names <- c(Names,out)
# 		}
# 		return(Names)
# 	}
# 	if(type==2){
# 		# interraction .
# 		Names <- NULL
# 		for(i in 1:length(x)){
# 			pos <- grep("\\*",unlist(strsplit(x[i],split="")))
# 			if(length(pos)>0){
# 				xgauche <- substr(x[i],start=1,stop=(pos-1))
# 				xdroite <- substr(x[i],start=(pos+1),stop=length(unlist(strsplit(x[i],split=""))))
# 				out <- c(xgauche,xdroite)
# 			}else{
# 				out <- x[i]
# 			}
# 			Names <- c(Names,out)
# 		}
# 		return(Names)
# 	}
# 
# }





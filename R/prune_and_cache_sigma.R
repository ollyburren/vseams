#!/usr/bin/RScript


library(GenomicRanges)
library(corpcor)
library(mvtnorm)
library(Matrix)

#############
###FUNCTIONS#
#############



prune.snps<-function(r2,gwas.gr,thr){
	clusters<-'1'
	if(length(r2)<2){
		names(clusters)<-colnames(r2)
	}else{
		D<-as.dist(1-r2)
		hc <- hclust(D, method="complete")
		clusters <- cutree(hc, h=1-thr)
	}
	names(clusters)<-gsub("^.*::","",names(clusters))
	clusters<-clusters[which(names(clusters) %in% start(gwas.gr))]
	
	gwas.gr[start(gwas.gr) %in% names(clusters[!duplicated(clusters)])]
	
}

attempt.pos.def<-function(mat,diag.val=1.0001){
  print(paste("diag.val",diag.val))
	if(!is(mat,"Matrix"))
		stop("mat is not a Matrix!")
	if(diag.val >= 1.1){
	  print("Matrix is not positive definite. Finding closest approximation..")
		diag(mat)<-1
		return(as(make.positive.definite(mat),"Matrix"))
	}
	diag(mat)<-diag.val
	if(is.positive.definite(mat,,method="chol")==FALSE){
	  new.diag<-signif(1+((diag.val-trunc(diag.val))*10))
		mat<-attempt.pos.def(mat,new.diag)
	}else{
		return(mat)
	}
}

mvs.sigma.r2<-function(r2){
	diag(r2)<-1
	if(!is.positive.definite(r2,,method="chol")){
		#this recurses through various values of diag if we exceed 1 then
		#we compute the closest matrix that is positive definite.
		r2<-attempt.pos.def(r2)
	}
	r2
}

convert.snp.df<-function(fname){
	snp.df<-get(load(fname))
	ld.block.filename<-unique(snp.df$ld.block.filename)
	snp.gr<-with(snp.df,GRanges(seqnames=Rle(seqnames),ranges=IRanges(start=start,end=end)))
	mcols(snp.gr)<-snp.df[,6:ncol(snp.df)]
	list(snp.gr=snp.gr,ld.block.filename=ld.block.filename)
}


##EXPECTS

#sigma_cache_dir='/sigma/'                            
#r2_threshold=0.95                                                              
#support_files='chr10_30687816-30783325.RData'
#out_dir='/pruned/'                                  

##TODO WE MAY WANT TO IMPLEMENT A FORCED OVERWRITING OF SIGMA CACHE ?

overwrite_sigma_cache=FALSE

args<-commandArgs(TRUE)
if(length(args) < 4){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{                                      
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

print(args)
                                                                    

## prevent non complete bins messing things up
support_files<-sub("[,]+$","",support_files)

res<-lapply(unlist(strsplit(support_files,',')),function(support_file){
	print(paste("Processing",support_file))
	olist<-convert.snp.df(support_file)
	r2<-get(load(olist$ld.block.filename))
	
	## consider removing dependency on GRanges here ?
	snp.gr<-prune.snps(r2,olist$snp.gr,r2_threshold)                   
	## remove the dependenecy on GRanges libraries
	snps<-as.data.frame(mcols(snp.gr))
	rownames(snps)<-paste(seqnames(snp.gr),start(snp.gr),sep="::")
	save(snps,file=paste(out_dir,basename(olist$ld.block.filename),sep=""))

	##COMPUTE SIGMA

	sigma.file<-paste(sigma_cache_dir,basename(support_file),sep="")
	if(!file.exists(sigma.file) | overwrite_sigma_cache){
		sigma<-as.matrix(mvs.sigma.r2(r2))
		save(sigma,file=sigma.file)
	}else{
		##print("Sigma found skipping")
	}
})
	
print("Success")


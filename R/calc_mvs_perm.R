library(Matrix)
library(corpcor)
library(mvtnorm)

mvs.perm<-function(sigma,n=1000){
	if(!is.matrix(sigma))
		stop("sigma parameter is not a matrix")		
	if(!is.positive.definite(sigma,,method="chol"))
		stop("sigma is not positive definite")
	## in original paper method="chol" was not defined so I assume used eigen default
	## this is slower than the choleski decomp ! Perhaps we should contact the author ?
	rd<-rmvnorm(n,mean=rep(0,ncol(sigma)),sigma=sigma,method="chol")
	t(rd)
}

## EXPECTS 

#n_perms=1000                                                         
#sample_sigma_dir='/ds_sigma/'
#test=TRUE                                                             
#snp_dir='/pruned/'          
#sigma_cache_dir='/sigma/'                     
#run_count='6'                                                          
#chr_name='chr9'                                                       
#out_dir='/perms/'       


args<-commandArgs(TRUE)
print(args)
if(length(args) < 5){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}


#load(support.file)
pattern<-paste("^",chr_name,'_.*\\.RData',sep="")
snp_files<-list.files(pattern=pattern,path=snp_dir,full.names=TRUE)
## do all files on the same chromosome together
## otherwise too many jobs and not so efficient
dev.null<-lapply(seq_along(snp_files),function(i){
	if(!file.exists(snp_files[i])) return()
	snp_file<-snp_files[i]
	message(snp_file)
	load(snp_file)
	if(nrow(snps)==0) return ()
	#print(snps)
	snps$start<-gsub(".*::","",rownames(snps))
	## here we check to see if we have to use a downsampled version
	## of sigma
	#rownames(sigma)<-sub("^.*::","",rownames(sigma))
	sigma.file<-paste(sigma_cache_dir,basename(snp_file),sep="")
	assign('sigma',get(load(sigma.file)))
	rownames(sigma)<-sub("^.*::","",rownames(sigma))
	#colnames(sigma)<-sub("^.*::","",colnames(sigma))
	rindex<-which(rownames(sigma) %in% snps$start)
	if(length(rindex)==1){
		sname<-rownames(sigma)[rindex]
		sigma<-as.matrix(sigma[rindex,rindex])
		rownames(sigma)<-sname
	}else{
		sigma<-sigma[rindex,rindex]
	}
	#sample.sigma.file<-paste(sample_sigma_dir,basename(snp_file),sep="")
	#if(file.exists(sample.sigma.file)){
	#	assign('sigma',get(load(sample.sigma.file)))
	#}else{
	#	##otherwise we use the cannonical one
	#	sigma.file<-paste(sigma_cache_dir,basename(snp_file),sep="")
	#	assign('sigma',get(load(sigma.file)))
	#}
	### hopefully this is deprecated after testing
	if(!is.matrix(sigma))
		sigma<-as.matrix(sigma)
	perms<-mvs.perm(sigma,n_perms)
	rownames(perms)<-gsub("^.*::","",rownames(sigma))
	match.index<-which(rownames(perms) %in% snps$start)
	perms<-perms[match.index,]
	if(!is.matrix(perms)){
		 perms<-matrix(perms,ncol=n_perms)
		 rownames(perms)<-gsub("^.*::","",rownames(sigma))[match.index]
	}
	snps<-snps[match(rownames(perms), snps$start),]
	rownames(perms)<-rownames(snps)
	out.file<-gsub("RData",paste(run_count,n_perms,"RData",sep="."),basename(snp_file))
	out.file<-paste(out_dir,out.file,sep="")
	save(perms,file=out.file)
})
print(paste("Processed",chr_name,"with",length(snp_files),"files"))
print("Success")


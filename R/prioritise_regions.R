library(rtracklayer)
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


args<-commandArgs(TRUE)
if(length(args) < 6){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{                                      
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

##EXPECTS

#n_perms=20                                                               
#test='SET.EZH2.test'                                                     
#support_file='support.RData'
#snp_manifest='barrett.meta.t1d.b37.bed'       
#index_file='eur.ld.index.RData'      
#region_id='ENSG00000155729'                                              
#sigma_cache_dir='/sigma/'                      
#out_dir=

chunksize=1000;
if(chunksize>n_perms)
	chunksize=n_perms

print(args)

snp.gr<-import.bed(snp_manifest,asRangedData=FALSE)
sigma.files<-list.files(path=sigma_cache_dir,pattern="*.RData")
##get the index
load(index_file)
##get region file
load(support_file)
##subset by a test set
spec.gr<-region.gr[mcols(region.gr)[[test]],]

ztominus.p<-function(z){
	-log(2*pnorm(abs(z),lower.tail=FALSE))
}

myVegas<-function(region.id,index.gr,spec.gr,sigma_cache_dir,sample_sigma_dir,n_perms){
	gene.gr<-spec.gr[spec.gr$names==region.id,]
	##get the list of SNPs that we need
	spec.snps<-subsetByOverlaps(snp.gr,gene.gr)
	#create gene summary stats
	spec.snps$m.log.p<- -log(spec.snps$score)
	mean.log.p<-sum(spec.snps$m.log.p)/length(spec.snps)
	##next we need permutations
	##get the list of files we need
	sigma.files<-paste(sigma_cache_dir,basename(subsetByOverlaps(index.gr,spec.snps)$file),sep="")
	if(length(sigma.files)==0)
		return("NO_COVERAGE")
		
		
	ld.totals<-lapply(sigma.files,function(sfile){
			print(paste("Processing",sfile))
			#sample.sigma.file<-paste(sample_sigma_dir,basename(sfile),sep="")
			assign('sigma',get(load(sfile)))
		
			##FILTERING LIKE THIS MAKES LITTLE DIFFERENCE
			##AND SPEEDS THINGS UP A LOT
			rownames(sigma)<-sub("^.*::","",rownames(sigma))
			#colnames(sigma)<-sub("^.*::","",colnames(sigma))
			rindex<-which(rownames(sigma) %in% start(spec.snps))
			if(length(rindex)==1){
				sname<-rownames(sigma)[rindex]
				sigma<-as.matrix(sigma[rindex,rindex])
				rownames(sigma)<-sname
			}else{
				sigma<-sigma[rindex,rindex]
			}
			if(!is.matrix(sigma))
				sigma<-as.matrix(sigma)
				## do 10,000 at a time
				p.it<-lapply(1:(n_perms/chunksize),function(p.count){
				#message(p.count)
				perms<-mvs.perm(sigma,chunksize)
				rownames(perms)<-gsub("^.*::","",rownames(sigma))
				match.index<-which(rownames(perms) %in% start(spec.snps))
				perms<-perms[match.index,]
				if(!is.matrix(perms)){
					perms<-matrix(perms,ncol=chunksize)
					rownames(perms)<-gsub("^.*::","",rownames(sigma))[match.index]
				}
				p.snps<-spec.snps[match(rownames(perms), start(spec.snps)),]
				list(sum=colSums(ztominus.p(perms)),snp.count=length(p.snps))
			})
			ret<-list()
			ret$sum<-unlist(lapply(p.it,"[[",'sum'))
			ret$snp.count<-unique(unlist(lapply(p.it,"[[",'snp.count')))
			ret
		})
		total.sum<-rowSums(do.call("cbind",lapply(ld.totals,"[[","sum")))
		total.snps<-rowSums(do.call("cbind",lapply(ld.totals,"[[","snp.count")))
		if(total.snps != length(spec.snps))
			print(paste("Warning total.snps is",total.snps,"expecting",length(spec.snps)))
		null.dist<-total.sum/total.snps
		sum(ifelse(mean.log.p > null.dist, 0, 1))/length(null.dist)
}

stat<-myVegas(region_id,index.gr,spec.gr,sigma_cache_dir,sample_sigma_dir,n_perms)
write(stat,file=paste0(out_dir,region_id))
print("Success")

library(Matrix)

## some of our regions are very large and/or contain a lot of snps
## some of our regions are modest but contain a lot of snps
## this causes a slow down in computation of MVS. Need a sampling 
## method that for those regions which are large and or contain a lot 
## of SNPs downsamples so that we can more easily compute downstream MVS
## This method downsamples r2 input matrix and then spikes in the target snps

## I sampled 10,000 blocks from the genome and got this distribution
## > summary(as.numeric(counts))
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    1.0    30.0   258.0   495.8   674.0 13810.0 
## we try this (means that approx < 5% of blocks undergo sampling)
## MAX_SNPS_PER_SIGMA_BLOCK=2000

## EXPECTS 

#sigma_cache_dir<-'/sigma/'
#snp_file
#snp_dir<-'/pruned/'
#out_dir<-'sampled_sigma/'
#max_snps_per_sigma_block=2000


args<-commandArgs(TRUE)
if(length(args) < 5){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

print(args)
(load(snp_file))
#print(snps)
snps$start<-gsub(".*::","",rownames(snps))
sigma.file<-paste(sigma_cache_dir,basename(snp_file),sep="")
print(sigma.file)
assign('sigma',as.matrix(get(load(sigma.file))))

## here we implement sampling for those sigma's that are large.
## to remove mvs calculation overhead
## we spike in our target snps, so we don't break
## downstream code.
if(dim(sigma)[1] > max_snps_per_sigma_block & nrow(snps)>0){
	sigma.pos<-gsub("^.*::","",rownames(sigma))
	sample.index<-which(sigma.pos %in% snps$start)
	print(sample.index)
	sample.index<-unique(c(sample.index,sample((1:length(sigma.pos))[-sample.index],max_snps_per_sigma_block-length(sample.index))))
	sigma<-sigma[sample.index,sample.index]
	sample.sigma.file<-paste(out_dir,basename(snp_file),sep="")
	save(sigma,file=sample.sigma.file)
}
                         
print("Success")

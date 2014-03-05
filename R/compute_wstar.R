## loads in a set of blocks of perms and using SNP manifest file
## works out those in test and control groups and calculates 
## wilcoxon

library(wgsea)

## EXPECTS 

#chunk_size=1000                                                                            
#perm_dir='/perms/'                                                
#control_set='SET.control1'                                                                 
#run_count='2'                                                                              
#snp_manifest='SET.test1.SET.control1.RData'
#out_dir='/wstar/'                                         
#test_set='SET.test1'                                                          

args<-commandArgs(TRUE)
if(length(args) < 7){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
print(args)
pattern<-paste("*",run_count,chunk_size,"RData",sep="\\.")
files<-list.files(path=perm_dir,pattern=pattern,full.names=TRUE)
perms<-do.call('rbind',lapply(files,function(x) get(load(x))))

snps<-get(load(snp_manifest))

## next filter perms based on universe of SNPs

snps<-snps[snps[[test_set]] | snps[[control_set]],]
perms<-perms[which(rownames(perms) %in% rownames(snps)),]
snp.index<-which(rownames(snps) %in% rownames(perms))
if(length(snp.index) != nrow(snps)){
	## this means that we have not permuted every SNP and an error should be returned
	## indicates a failure in the sigma step
	print (snps[-snp.index,])
	print("ERROR snp manifest does not match perms matrix") 	
}else{
	snps.in<-which(rownames(perms) %in% rownames(snps[snps[[test_set]],]))
	Wstar<-wilcoxon(perms,snps.in=snps.in)
	out.file<-paste(out_dir,paste(run_count,test_set,control_set,'WSTAR','RData',sep="."),sep="")
	save(Wstar,file=out.file)
	print("Success")
}

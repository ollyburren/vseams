
library(data.table)
## loads in individual snp files, concatenates them together
## if there are snps that are in both test and controls sets
## randomly assigns them to a set

## EXPECTS 

#control_set='SET.control1'                                   
#snp_dir=/pruned/'              
#out_dir='/mytest1/snp_manifest/'
#test_set='SET.test1'                                         



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

files<-list.files(path=snp_dir,pattern="*.RData",full.names=TRUE)
snps<-lapply(files,function(x){
		tmp<-get(load(x))
		tmp<-tmp[,c('name','score',test_set,control_set)]
		tmp$name<-as.character(tmp$name)
		tmp$rnames<-rownames(tmp)
		tmp
})
snps<-rbindlist(snps)
rownames(snps)<-snps$rnames
snps$rnames<-NULL
## for cases where TEST and CONTROL sets overlap randomly assign.
d.i<-which(snps[[test_set]] & snps[[control_set]])
snps[d.i,test_set]<-sample(c(TRUE,FALSE),length(d.i),replace=TRUE)
out.file<-paste(out_dir,paste(test_set,control_set,sep="."),'.RData',sep="")
save(snps,file=out.file)
print("Success")

library(wgsea)


## EXPECTS

## wstar_dir<-'/wstar/'
## snp_dir<-'/pruned/'
## out_dir<-'/wilcoxon_params/'
## test_set<-'SET.protective'       
## control_set<-'SET.all'
## dataset<-'ds'

args<-commandArgs(TRUE)
if(length(args) < 6){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

                        
## load in and concatenate wstar
perm.pattern<-paste(test_set,control_set,'WSTAR','RData',sep=".")
wstar.files<-list.files(path=wstar_dir,pattern=paste('*',perm.pattern,sep=""),full.names=TRUE)
wstar<-unlist(lapply(wstar.files,function(x) get(load(x))))

##compute w
                                                     

#snps<-do.call('rbind',lapply(snp.files,function(x) get(load(x))))

snp<-get(load(snp_manifest))

##
snps<-snps[snps[[test_set]] | snps[[control_set]],]
snps.in<-which(snps[[test_set]])
w<-wilcoxon(snps$score,snps.in=snps.in)

n.in<-sum(snps[[test_set]])
n.out<-nrow(snps)-n.in

param<-list(w=w,wstar=wstar,n.in=n.in,n.out=n.out)

out.file<-paste(out_dir,paste(dataset,test_set,control_set,'param','RData',sep="."),sep="")
save(param,file=out.file)

print("Success")

#!/usr/bin/RScript

args<-commandArgs(TRUE)
if(length(args) < 2){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

print(args)

library(GenomicRanges)
#in.file<-''
#out.dir<-''

assign('t',get(load(in.file)))
fnames<-lapply(seq_along(t),function(y){
		gr<-t[y,]
		r2<-gr$r2[[1]]
		chr<-as.character(seqnames(gr))
		start<-start(gr)
		end<-end(gr)
		fname=paste(out.dir,chr,'_',start,'-',end,'.RData',sep="")
		print(fname)
		save(r2,file=fname)
		fname
})
print(unlist(fnames))
print("Success")


## create genome hapmap recomb regions for genome - this can then be fed into cal_sigma_1kg.pl
## for use in sandman
#source("generate_hapmap_recomb_regions.R",echo=T)

MAX_TABIX_BIN_SIZE<-(2^29)-1 ## this is the max bin size that tabix can deal with

library(GenomicRanges)

## Get HapMap Recombination regions
## curl -s http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz | tar -xvz
## set location above to be your recomb.dir location

seqnames<-paste("chr",c(1:22,'X'),sep="")
recomb.thresh=0.1
## see above
recomb.dir<-'FILL ME IN';
## valid output dir
out.dir='FILL ME IN'

calc.intervals.by.recomb<-function(r.df,thresh){
	#must be in positional order
	r.df<-r.df[order(r.df$position),]
	telo.df<-data.frame(start=c(1,max(r.df$position)+1),end=c(min(r.df$position)-1,MAX_TABIX_BIN_SIZE))
	print(telo.df)
	gmap<-r.df$genetic_map_position
	r.df$interval<-trunc(gmap/thresh)+1
	r.df.split<-split(r.df,r.df$interval)
	#here we adjust intervals slightly to use max recomb rate as within an interval this makes sense?
	interval.position<-sapply(r.df.split,function(x) x[which.max(x$recombination_rate),]$position)
	idx<-which(r.df$position %in% interval.position)
	t.df<-data.frame(start=r.df[head(idx,n=length(idx)-1),]$position+1,end=r.df[idx[-1],]$position)	
	## we will miss start and end of chromosome using this so we need to add into our data frame 
	## two extra entries
	t.df<-rbind(t.df,telo.df)
	t.df[order(t.df$start),]
}


gr.recomb<-function(seqname){
	recomb.file <- paste(recomb.dir,'genetic_map_GRCh37_',seqname,'.txt',sep="")
	print(paste("Processing",seqname))
	recomb.df<-read.table(file=recomb.file,header=TRUE,sep="\t")
	names(recomb.df)<-c('chrom','position','recombination_rate','genetic_map_position')
	regions.df<-calc.intervals.by.recomb(recomb.df,recomb.thresh)
	interval.gr<-with(regions.df,
		GRanges(seqnames=Rle(seqname),
						IRanges(start=start,end=end,name=1:nrow(regions.df)))
						)
	split.on<-trunc(as.numeric(names(ranges(interval.gr)))/10)
	sapply(split(interval.gr,split.on),function(x){
		fs<-min(start(x))
		fe<-max(end(x))
		ofile<-paste(out.dir,seqname,'_',fs,'-',fe,'.RData',sep="")
		region.gr<-x
		save(region.gr,file=ofile)
		region.gr
		}
	)
}

for(i in seqnames){
	gr.recomb(i)
}



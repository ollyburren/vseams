library(GenomicRanges)

## location of output from split_r2_1kg_vcf.pl
ind.dir<-'FILL ME IN'
## where to put index file created
out.file<-'FILL ME IN'

create_dir_index<-function(in.dir){
	files<-list.files(path=ind.dir,pattern="*.RData$",full.name=TRUE)
	## parse filenames to make a GRanges object
	tmp<-gsub("\\.RData","",basename(files))

	chr=gsub("^([^_]+).*","\\1",tmp)
	start=as.numeric(gsub(".*_([^\\-]+)\\-.*","\\1",tmp))
	end=as.numeric(gsub(".*_[^\\-]+\\-(.*)","\\1",tmp))

	GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),file=files)
}

index.gr<-create_dir_index(in.dir)
save(index.gr,file=out.file);

#!/usr/bin/RScript


library(rtracklayer)
library(GenomicRanges)

#############
###FUNCTIONS#
#############

snps2region<-function(snp.gr,region.gr,add.region.id=FALSE){
	if(!is(snp.gr,"GRanges") || !is(region.gr,"GRanges"))
		stop("Both parameters are required to be GRanges objects")
	#first step is to assign genenames to each snp
	names.gr<-names(mcols(region.gr))
	tmp.mcols<-mcols(snp.gr)
	for(r in names.gr[grep("^SET\\.",names.gr)]){
		print(paste("Proceesing ",r))
		index<-which(mcols(region.gr)[[r]]==TRUE)
		t.gr<-subsetByOverlaps(snp.gr,region.gr[index,])
		tmp.mcols[[r]]<-tmp.mcols$name %in% t.gr$name
	}
	mcols(snp.gr)<-tmp.mcols
	if(!add.region.id){	
		return(snp.gr)
	}
	
	## add gene id's 	
	
	## if region.gr$name is not defined we make it up
	
	if(length(region.gr$names)==0){
		region.gr$names<-paste("region",1:length(region.gr),sep=".")
	}
	ol<-as.data.frame(findOverlaps(snp.gr,region.gr))
	ol$region_id<-region.gr[ol$subjectHits,]$names
	#what to do if a SNP overlaps two sets of genes
	#here we just ignore it and take a first come first served.
	#for some applications this could be problematic. Take the 
	#example of multivariate sampling if a snp overlaps a test and control
	#region then we need to assign to the test region preferentially TODO !!
	nr<-ol[!duplicated(ol$queryHits),]
	##only care about snps that overlap our test regions
	snp.gr<-snp.gr[sort(nr$queryHits),]
	tmp.mcols<-mcols(snp.gr)
	tmp.mcols$region_id<-character(length=length(snp.gr))
	tmp.mcols$region_id<-as.character(nr$region_id)
	mcols(snp.gr)<-tmp.mcols
	snp.gr
}

##EXPECTS


#ld_index_file='eur.ld.index.RData'      
#region_file='support.RData'
#gwas_bed='gwas.bed'                 
#exclude_bed='exclude.regions.bed'        
#out_dir='/support/'
#max_snps_per_sigma_block=2000


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

gwas.gr<-import.bed(gwas_bed,asRangedData=FALSE)
assign('region.gr',get(load(region_file)))
assign('index.gr',get(load(ld_index_file)))

## first we remove any genes/regions that overlap
## our list of excluded regions

exclude.gr<-import.bed(exclude_bed,asRangedData=FALSE)
e<-as.character((subsetByOverlaps(region.gr,exclude.gr))$names)
if(length(e)>0)
	region.gr<-region.gr[-which(region.gr$names %in% (subsetByOverlaps(region.gr,exclude.gr))$names),]

## which GWAS snps

gwas.gr<-subsetByOverlaps(gwas.gr,region.gr)

## we may wish to assess which regions we cannot assay (no coverage)
## dump these to a file in the output directory
## note that the .tab is important.
regions.missing<-region.gr[-which(region.gr$names %in% subsetByOverlaps(region.gr,gwas.gr)$names ),]

if(length(regions.missing)>0){
	miss<-as.data.frame(regions.missing)
	write.table(miss,file=paste(out_dir,'missing.tab',sep=""),row.names=FALSE,quote=FALSE,sep="\t")
}
## which LD blocks
gwas.gr<-snps2region(gwas.gr,region.gr,TRUE)

## some of the LD blocks are excluded here we compute if
## any of our target regions are affected and return this 
## as a report. These are excluded from further analysis.

excluded.gr<-index.gr[index.gr$exclude,]
ol.excluded<-as.matrix(findOverlaps(excluded.gr,gwas.gr))
ex.df<-as.data.frame(gwas.gr[unique(ol.excluded[,2]),])
write.table(ex.df,file=paste(out_dir,'snps_overlapping_excluded_regions.tab',sep=""))

## from now on we exclude those regions that are marked as such
index.gr<-index.gr[!index.gr$exclude,]

ld.fnames<-index.gr$file
gwas.t.df<-as.data.frame(gwas.gr)
##convert factors to characters to speed things up
##considerably
for(n in names(gwas.t.df)){
	if(is.factor(gwas.t.df[[n]])){
		print(paste("Converting",n))
		gwas.t.df[[n]]<-as.character(gwas.t.df[[n]])
	}
}
ol.snp<-as.matrix(findOverlaps(index.gr,gwas.gr))

## it is expedient at this stage to create a list of all
## those regions that will need to be downsampled

ds.gr<-index.gr[unique(ol.snp[,1]),]
ds.gr<-ds.gr[ds.gr$snp.count >= max_snps_per_sigma_block,]
write(ds.gr$file,file=paste(out_dir,'downsample_files.tab',sep=""))
out.files<-lapply(unique(ol.snp[,1]),function(x){
	snp.df<-gwas.t.df[ol.snp[ol.snp[,1]==x,2],]
	ld.block.filename<-ld.fnames[x]
	snp.df$ld.block.filename<-ld.fnames[x]
	file.out<-paste(out_dir,basename(ld.block.filename),sep="")
	save(snp.df,file=file.out)
	file.out
})

## we could attempt to send back the files that have been created ??

#print(out.files)

print("Success")


#!/usr/bin/RScript


library(biomaRt)
library(GenomicRanges)

#############
###FUNCTIONS#
#############



get.bed<-function(fname){
	gr<-GRanges()
	not.bed=FALSE
	tryCatch(gr<-import.bed(fname,asRangedData=FALSE),error=function(e){
		print("Not a bed file")
		not.bed<<-TRUE
	})
	if(not.bed) return(gr)
	log.m<-lapply(unique(gr$name),function(x){
		gr$name %in% x
	})
	names(log.m)<-paste("SET",unique(gr$name),sep=".")
	mcols(gr)<-DataFrame(names=as.character(1:length(gr)),external_id=as.character(1:length(gr)))
	mcols(gr)<-cbind(mcols(gr),DataFrame(log.m))
	gr
}


convertEnsIDToGranges<-function(geneset.list,ensembl,tss.extension=NULL){
	if(length(geneset.list)<=1)
		print("Warning you have 1 or less gene sets consider using convertEnsIDtoGRanges instead")
	#first step is to unlist to get the full set
	glist<-unique(unlist(geneset.list))
  ## check paramaters
  ## attempt to guess species
  identifier<-table(gsub('[0-9]+$','',glist))
  if(length(identifier)>1)
  	stop("Trying to mix ensembl identifiers")
  external_id<-ifelse(names(identifier)=="ENSG",'hgnc_symbol','external_gene_id')
	gene.details<-getBM(
		filters= c("ensembl_gene_id"),
		attributes= c('ensembl_gene_id','chromosome_name','start_position','end_position','strand',external_id),
		values= list(ensembl_gene_id=glist),
		mart= ensembl)
	## Mitochondrial 'chromosome' is mostly named M
	
	mito.ndx<-which(gene.details$chromosome_name=="MT")
	if(length(mito.ndx)>1)
		gene.details[mito.ndx,]$chromosome_name="M"
		
	names(gene.details)<-gsub(external_id,'external_id',names(gene.details))
	## create GRange object
	gr<-with(gene.details,GRanges(
		seqnames=Rle(paste("chr",chromosome_name,sep ="")),
		ranges=IRanges(start=start_position,end=end_position),
		strand=strand,
		names=ensembl_gene_id,
		external_id=external_id))
	if(!is.null(tss.extension)){
		tmp.gr<-gr
		tss.neg<-which(strand(tmp.gr)=="-")
		tss.pos<-which(strand(tmp.gr)=="+")
		end(tmp.gr[tss.pos])<-ifelse(start(tmp.gr[tss.pos])+tss.extension>end(tmp.gr[tss.pos]),
			start(tmp.gr[tss.pos])+tss.extension,
			end(tmp.gr[tss.pos])
		)
		
		start(tmp.gr[tss.pos])<-start(tmp.gr[tss.pos])-tss.extension
		start(tmp.gr[tss.neg])<-ifelse(end(tmp.gr[tss.neg])-tss.extension < start(tmp.gr[tss.neg]),
			end(tmp.gr[tss.neg])-tss.extension,
			start(tmp.gr[tss.neg])
		)
		end(tmp.gr[tss.neg])<-end(tmp.gr[tss.neg])+tss.extension
		gr<-tmp.gr
	}
	if(length(gr) != length(unique(glist)))
		print("WARNING: Input and output gene lists have different lengths")
		
	tmp.df<-as.data.frame(mcols(gr))
	for(i in names(geneset.list)){
		tmp.df[[paste("SET",i,sep=".")]]<-gr$names %in% geneset.list[[i]]
	}
	mcols(gr)<-tmp.df
	gr
}

## THIS SCRIPT TAKES EITHER A BED FILE OR A LIST OF ENSEMBL ID's AND CONVERTS THEM TO GENOMIC RANGES

##EXPECTS

#input_file<-'file.tab'
## OR input_file<-'file.bed'
#tss_extension=200000
#mart_host<-'feb2012.archive.ensembl.org'
#output_file<-'file.RData'

args<-commandArgs(TRUE)
if(length(args) < 4){
  cat("Error incorrect number of args","\n",sep="")
  q()
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}

## check to see if we have a gene list or a region bed file
region.gr<-get.bed(input_file)
## if not we use biomart to get coords etc.
if(length(region.gr)==0){
	genes<-read.table(file=input_file,header=F,sep="\t")
	names(genes)<-c('ensid','type')
	gene.list<-split(genes$ensid,genes$type)
	mart <- useMart('ENSEMBL_MART_ENSEMBL',host=mart_host)
	ensembl.gene <- useDataset("hsapiens_gene_ensembl",mart=mart)
	region.gr<-convertEnsIDToGranges(gene.list,ensembl.gene,tss.extension=tss_extension)
}
save(region.gr,file=output_file)

print("Success")



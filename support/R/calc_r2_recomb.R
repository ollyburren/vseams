#!/usr/bin/RScript

## This is a slave script that should be called by cal_r2_1kg_vcf.pl

library(GenomicRanges)
library(snpStats)
library(VariantAnnotation)

args<-commandArgs(TRUE)
if(length(args)<1){
	cat("Error incorrect number of args","\n",sep="")
	q()
}else{
	for(i in 1:length(args)){                                    
		eval(parse(text=args[[i]]))
	}
}

## need arg region.file  

#region.file='chrX_154182953-536870911.RData'

## note that outfile is hardcoded to be a dir SIGMA
#vcfext='.phase1_release_v3.20101123.snps_indels_svs.genotypes.refpanel.EUR.vcf.bgz'

#thou_gen_data_dir<-''

##SNP QC filters
call.rate<-0.99
z.HWE.co<-25
MAF.co<-0.01

#region.file<-''

out.file<-gsub("RECOMB_0.1_REGIONS","LD/EUR",region.file)

assign('regions.gr',get(load(region.file)))

seqlevels(regions.gr)<-gsub("chr","",seqlevels(regions.gr))

chr<-gsub("^(chr[^_]+).*","\\1",basename(region.file))

vcfname<-paste(thou_gen_data_dir,chr,vcfext,sep="")
vcfp<-ScanVcfParam(which=regions.gr)
thou<-readVcf(vcfname,"hg19",vcfp)
#this code to remove largescale structural deletions which break the snpMatrix conversion code
structural <- logical(length(alt(thou)))
structural[grep("<", unlist(alt(thou)), fixed=TRUE)] <- TRUE
#remove large features as these are indicative of largescale structural issues
structural[width(thou)>100]<-TRUE
thou <- thou[!structural, ]
if(class(alt(thou)) != "DNAStringSetList"){
  alt(thou)<-VariantAnnotation:::.toDNAStringSetList(unlist(alt(thou),use.names=F))
}
lup<-data.frame(rs=names(rowData(thou)),lu=paste(names(rowData(thou)),start(rowData(thou)),sep="::"))
rd<-rowData(thou)
### do this next bit on an interval basis
regions.gr$r2<-lapply(names(ranges(regions.gr)),function(x){
  t.vcf<-thou[rd$paramRangeID == x]
	calls<-geno(t.vcf)$GT
	a0<-ref(t.vcf)
	a1<-alt(t.vcf)
	
	## supress warnings 
	region.snpm<-suppressWarnings(genotypeToSnpMatrix(calls,a0,a1))
	colnames(region.snpm$genotypes)<-lup[match(colnames(region.snpm$genotypes),lup$rs),]$lu
	## do some QC 
	## this means that there are no snps in this region
	if(nrow(region.snpm$map)==0 || sum(!region.snpm$map$ignore)==0)
		return(NA)                          
	sum<-col.summary(region.snpm$genotypes)
	ok.index<-which(with(sum, Call.rate >= call.rate & z.HWE^2 < z.HWE.co & MAF > MAF.co ))
	## now work out sigma
	print(paste(x,length(ok.index)))
	if(length(ok.index)==0)
		return(NA)
	if(length(ok.index)==1){
		snp.name<-rownames(sum[ok.index,])
		return(Matrix(1,dimnames = list(snp.name,snp.name)))
	}
	## filter out bad QC snps
	region.snpm$genotypes<-region.snpm$genotypes[,ok.index]
	r2<-forceSymmetric(ld(x=region.snpm$genotypes,depth=ncol(region.snpm$genotypes)-1,stats="R.squared",symmetric=TRUE))
	r2[which(is.na(r2))]<-0;
	r2
})

seqlevels(regions.gr)<-paste("chr",seqlevels(regions.gr),sep="")
save(regions.gr,file=out.file)
print("Success")


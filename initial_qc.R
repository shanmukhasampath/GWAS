library('ggplot2')
library('stringr')
options(stringsAsFactors=FALSE)
source('~/Scripts/gwasfunctions.r')
plinkdir <- "~/Software/plink-1.9/plink"
smartpca <- '~/Software/EIG-6.1.4/bin/smartpca'
pcasoftware <- 'EIGENSOFT'
onekgdir <- "~/1000genomes_European"

cdir <- "INPUT DATA DEPOSIT DIRECTORY"
ifile <- 'INPUTFILENAME'
ofile <- 'OUTPUTFILENAME'
dataset <- 'DATASET'
platformchip <- 'GENOTYPE PLATFORM'
system(sprintf('echo \'The STATS START HERE \n\'> %s/%s_stats.txt',cdir,ofile))
######### ------------------------------------------------------ B. Report, Annotation, Quality control ------------------------------------------------------ #########

# Steps involved:
# B1. Convert the tped,tfam files to bed,bim and fam files
# B2. Report on the initial data with number of SNPs, samples, populations, Sex/Gender, SNP stats
# B3. Update the SNP annotation
# B4. Remove SNPs with chromosomes 0, Y, XY, MT. Remove INDELS and DELETIONS. Remove SNPs with Duplicated Names and Positions
# B5. Separate Autosomal and SexChromosome SNPs
# B6. Plot MAF and Missingness before geno < 0.05 carried out in step B7
#  B6.1 Minor Allele Frequency before QC
#  B6.2 Missingness per individual before QC
#  B6.3 SNP MAF vs SNP Missingness
# B7. Keep SNPs that have missingness < 5% i.e geno < 0.05 
# B8. Plot MAF and Missingness after geno < 0.05 carried out in step B7. Remove Samples that have high individual missingness with > 2%
# B9. Remove samples that are contaminated with Fhet > 0.05 and Fhet < -0.05 --- ONLY for datasets with ONE POPULATION of either CASES or CONTROLS
# B10. Remove SexAmbiguous Samples using SexCheck calculated on Sexchromsome using PLINK
# B11. Filter for SNPs that have geno < 0.02 which is < 2% missingess
# B12. Filter for SNPs that show genotype-rate differenc between cases and controls ---- ONLY for datasets with CASES and CONTROLS present
# B13. Filter for SNPs with Hardy-Weinberg Eequilibirum p-values < 1e-06 for controls or < 1e-10 for cases  

# B1. Convert the tped,tfam files to bed,bim and fam files
if(platformchip %in% c('oe24v1_2','oe24v1_1')){
	system(sprintf('%s --noweb --tped %s/%s.tped --tfam %s/%s.tfam --keep-allele-order --allow-no-sex --make-bed --out %s/%s',plinkdir,cdir,ifile,cdir,ifile,cdir,ofile))
} else {
	system(sprintf('%s --noweb --bfile %s/%s --keep-allele-order --allow-no-sex --make-bed --out %s/%s',plinkdir,cdir,ifile,cdir,ofile))
}

# B2. Report on the initial data with number of SNPs, samples, populations, Sex/Gender, SNP stats
printstats(cdir,ofile,'Initial Dataset Stats',ofile)

# B3. Update the SNP annotation
updatesnps(cdir,ofile,platformchip)
if( platformchip %in% c('affy500','affysnp6','hhap550') ){
	system(sprintf('%s --noweb --bfile %s/%s --exclude %s/%s_update_remove.txt --update-chr %s/%s_update_chr.txt --update-map %s/%s_update_map.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%stemp',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile,cdir,ofile,cdir,ofile))
} else {
	system(sprintf('%s --noweb --bfile %s/%s --update-chr %s/%s_update_chr.txt --update-map %s/%s_update_map.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%stemp',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile,cdir,ofile))
} 

if( platformchip %in% c('oe12v1','oe12v1_1','oe24v1_2','o8','affy500','affysnp6') ){
	system(sprintf('%s --noweb --bfile %s/%stemp --update-name %s/%s_update_name.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%stemp1',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile))
	printstats(cdir,ofile,'After Updating SNP annotation',ofile)

} else if( platformchip %in% c('o25m_8v1_2','oe24v1_1','hhap550') ) {
	system(sprintf('%s --noweb --bfile %s/%stemp --keep-allele-order --allow-no-sex --make-bed --out %s/%stemp1',plinkdir,cdir,ofile,cdir,ofile))
	printstats(cdir,ofile,'After annoating SNP map and chromosome build but not name',ofile)

} else {
	system(sprintf('%s --noweb --bfile %s/%stemp --keep-allele-order --allow-no-sex --make-bed --out %s/%stemp1',plinkdir,cdir,ofile,cdir,ofile))
	printstats(cdir,ofile,'SNP annotation step could not be performed due to lack chip information',ofile)
}

# B4. Remove SNPs with chromosomes 0, Y, XY, MT. Remove INDELS and DELETIONS. Remove SNPs with Duplicated Names and Positions
bimcheck(cdir,sprintf('%stemp1',ofile))
system(sprintf('%s --noweb --bfile %s/%stemp1 --exclude %s/%stemp1_snpsremove.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%s',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile))
system(sprintf('rm %s/%stemp*.* %s/*_update_*.txt',cdir,ofile,cdir))
printstats(cdir,ofile,'Remove SNPs with chromosomes 0, Y, XY, MT. Remove INDELS and DELETIONS. Remove SNPs with Duplicated Names and Positions',ofile)

# B5. Separate Autosomal and SexChromosome SNPs
extractAutosomal(cdir,ofile)
system(sprintf('%s --noweb --bfile %s/%s --extract %s/%s_AutosomalSnps.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%sA',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile))
data <- read.table(sprintf('%s/%s_SexChromosomeSnps.txt',cdir,ofile),h=F,stringsAsFactors=F)
if( length(which(data$V1 %in% c('NoXChr'))) > 0 ){
	SexVar <- 'NoSex'
} else {
	SexVar <- 'Sex'
	system(sprintf('%s --noweb --bfile %s/%s --extract %s/%s_SexChromosomeSnps.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%sX',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile))
}
rm(data)
system(sprintf('rm %s/%s_AutosomalSnps.txt %s/%s_SexChromosomeSnps.txt',cdir,ofile,cdir,ofile))
printstats(cdir,sprintf('%sA',ofile),'Separate Autosomal and SexChromosome SNPs',ofile)

# B6. Plot MAF and Missingness before geno < 0.05 carried out in step B7
system(sprintf('%s --noweb --bfile %s/%sA --freq --allow-no-sex --out %s/%sABQCFreq',plinkdir,cdir,ofile,cdir,ofile))
system(sprintf('%s --noweb --bfile %s/%sA --missing -allow-no-sex --out %s/%sABQCmissing',plinkdir,cdir,ofile,cdir,ofile))

# B6.1 Minor Allele Frequency before QC
ggmafplot(cdir,sprintf('%sABQCFreq',ofile),'~/Plots/qualitycontrol')

# B6.2 Missingness per individual before QC
ggmissplot(cdir,sprintf('%sABQCmissing',ofile),'~/Plots/qualitycontrol')

# B6.3 SNP MAF vs SNP Missingness
ggmafvsmissplot(cdir,sprintf('%sABQCFreq',ofile),sprintf('%sABQCmissing',ofile),'~/Plots/qualitycontrol')
system(sprintf('rm %s/%sABQCFreq.* %s/%sABQCmissing.*',cdir,ofile,cdir,ofile))

# B7. Remove SNPs that have missingness of > 5% i.e geno > 0.05
system(sprintf('%s --noweb --bfile %s/%sA --geno 0.05 --keep-allele-order --allow-no-sex --make-bed --out %s/%sAG',plinkdir,cdir,ofile,cdir,ofile))
printstats(cdir,sprintf('%sAG',ofile),'Remove SNPs that have missingness of > 5% i.e geno > 0.05',ofile)

# B8. Plot MinorAlleleFrequency and Missingness after step B5 --- 4 samples failed the Missingness Test
# B8.1 Minor Allele Frequency after QC
system(sprintf('%s --noweb --bfile %s/%sAG --freq --out %s/%sAGAQCFreq',plinkdir,cdir,ofile,cdir,ofile))
ggmafplot(cdir,sprintf('%sAGAQCFreq',ofile),'~/Plots/qualitycontrol')
system(sprintf('rm %s/%sAGAQCFreq.*',cdir,ofile))

# B8.2 Missingness per individual after QC
system(sprintf('%s --noweb --bfile %s/%sAG --missing --allow-no-sex --out %s/%sAGAQCmissing',plinkdir,cdir,ofile,cdir,ofile))
ggmissplot(cdir,sprintf('%sAGAQCmissing',ofile),'~/Plots/qualitycontrol')

data <- read.table(sprintf('%s/%sAGAQCmissing.imiss',cdir,ofile),h=T,stringsAsFactors=F)
write.table(data[which(data$F_MISS > 0.02),c(1,2,1)],sprintf('%s/%s_missing_samples.txt',cdir,ofile),quote=F,row.names=F,col.names=F)
system(sprintf('%s --noweb --bfile %s/%sAG --remove %s/%s_missing_samples.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%sAGM',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile))
system(sprintf('rm %s/%sAGAQCmissing.*',cdir,ofile))
rm(data)
printstats(cdir,sprintf('%sAGM',ofile),'Remove Samples that have missingness of > 2% i.e mind > 0.02',ofile)

# B9. Remove samples that are contaminated with Fhet > 0.2 and Fhet < -0.2 --- ONLY for datasets with ONE POPULATION of either CASES or CONTROLS
famfile <- read.table(sprintf('%s/%sAGM.fam',cdir,ofile),h=F,stringsAsFactors=F)
if( length(unique(famfile$V1)) == 1 & length(unique(famfile$V6)) == 1 & length(unique(famfile$V6) %in% c(1,2)) == 1 ){
	system(sprintf('%s --noweb --bfile %s/%sAGM --het --keep-allele-order --allow-no-sex --out %s/%sAGMhet',plinkdir,cdir,ofile,cdir,ofile))
	data <- read.table(sprintf('%s/%sAGMhet.het',cdir,ofile),h=T,stringsAsFactors=F)
	write.table(data[which(data$F < -0.2 | data$F > 0.2),c(1,2,1)],sprintf('%s/%s_contaminatedsamples.txt',cdir,ofile),quote=F,row.names=F,col.names=F)
	system(sprintf('%s --noweb --bfile %s/%sAGM --remove %s/%s_contaminatedsamples.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%sAGMH',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile))
	system(sprintf('rm %s/%sAGMhet.*',cdir,ofile))
	rm(data)
	printstats(cdir,sprintf('%sAGMH',ofile),'Remove samples that are contaminated with Fhet > 0.2 and Fhet < -0.2',ofile)
} else {
	write.table(c(''),sprintf('%s/%s_contaminatedsamples.txt',cdir,ofile),quote=F,row.names=F,col.names=F)
	system(sprintf('%s --noweb --bfile %s/%sAGM --remove %s/%s_contaminatedsamples.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%sAGMH',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile))
	printstats(cdir,sprintf('%sAGMH',ofile),'Sample Contamination check is not performed due multiple populations',ofile)
}
system(sprintf('rm %s/%s_contaminatedsamples.txt',cdir,ofile))

# B10. Check for Sex Ambiguity using sexcheck from PLINK and remove samples wtih F value > 0.2 and < 0.8 --- 9 Samples failed in SexCheck
if( SexVar == 'NoSex' ){
	system(sprintf('%s --noweb --bfile %s/%sAGMH --keep-allele-order --allow-no-sex --make-bed --out %s/%sAGMHRS',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile))
	printstats(cdir,sprintf('%sAGMHRS',ofile),'Check for Sex Ambiguity Couldnot be done due to no X chromsome SNPs in the dataset',ofile)

} else {
	system(sprintf('%s --noweb --bfile %s/%sX --check-sex --allow-no-sex --out %s/%sSexCheck',plinkdir,cdir,ofile,cdir,ofile))
	data <- read.table(sprintf('%s/%sSexCheck.sexcheck',cdir,ofile),h=T,stringsAsFactors=F)
	pdf(sprintf("~/Plots/qualitycontrol/%sSexCheckAfterQC.pdf",ofile), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	print(ggplot(data,aes(x = F)) + xlab("X homozygosity Rate") + geom_histogram(colour="black",bins=100) + scale_x_continuous(breaks=seq(-0.1,1,by=0.1)) )
	dev.off()
	removesamples <- data[which( (data$F > 0.2 & data$F < 0.8) ),c(1,2,1)]
	checksamples <- data[which(data$STATUS == 'PROBLEM' & data$PEDSEX != 0 & data$SNPSEX != 0),]
	write.table(checksamples,sprintf('%s/%s_SexSamplesToCheck.txt',cdir,ofile),quote=F,sep='\t',row.names=F,col.names=T)
	write.table(removesamples,sprintf('%s/%s_sexproblem_samples.txt',cdir,ofile),quote=F,row.names=F,col.names=F)
	system(sprintf('%s --noweb --bfile %s/%sAGMH --remove %s/%s_sexproblem_samples.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%sAGMHRS',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile))
	system(sprintf('rm %s/%s_sexproblem_samples.txt %s/%sSexCheck.*',cdir,ofile,cdir,ofile))
	rm(data,checksamples,removesamples)
	printstats(cdir,sprintf('%sAGMHRS',ofile),'Check for Sex Ambiguity using sexcheck from PLINK and remove samples wtih F value > 0.2 and < 0.8',ofile)
}

# B11. Filter for SNPs that have geno < 0.02 which is < 2% missingess
system(sprintf('%s --noweb --bfile %s/%sAGMHRS --geno 0.02 --keep-allele-order --allow-no-sex --make-bed --out %s/%sAGMHRSG',plinkdir,cdir,ofile,cdir,ofile))
printstats(cdir,sprintf('%sAGMHRSG',ofile),'Filter for SNPs that have geno < 0.02 which is < 2% missingess',ofile)
# plot the missingess
system(sprintf('%s --noweb --bfile %s/%sAGMHRSG --missing --allow-no-sex --out %s/%sAGMHRSGmissing',plinkdir,cdir,ofile,cdir,ofile)) # Gentoype Call rate for every samples is 99%
ggmissplot(cdir,sprintf('%sAGMHRSGmissing',ofile),'~/Plots/qualitycontrol')
system(sprintf('rm %s/%sAGMHRSGmissing.*',cdir,ofile))

# B12. Filter for SNPs that show genotype-rate differenc between cases and controls ---- ONLY for datasets with CASES and CONTROLS present
famfile <- read.table(sprintf('%s/%sAGMHRSG.fam',cdir,ofile),h=F,stringsAsFactors=F)
if( length(unique(famfile$V6) %in% c(1,2)) == 2 ){
	system(sprintf('%s --noweb --bfile %s/%sAGMHRSG --filter-cases --missing --allow-no-sex --out %s/%sAGMHRSGcasesmissing',plinkdir,cdir,ofile,cdir,ofile))
	system(sprintf('%s --noweb --bfile %s/%sAGMHRSG --filter-controls --missing --allow-no-sex --out %s/%sAGMHRSGcontrolsmissing',plinkdir,cdir,ofile,cdir,ofile))
	casemiss <- read.table(sprintf('%s/%sAGMHRSGcasesmissing.lmiss',cdir,ofile),h=T,stringsAsFactors=F)
	controlmiss <- read.table(sprintf('%s/%sAGMHRSGcontrolsmissing.lmiss',cdir,ofile),h=T,stringsAsFactors=F) 
	write.table(casemiss[which(abs(casemiss$F_MISS - controlmiss$F_MISS) > 0.02),2],sprintf('%s/%s_snpstoremove.txt',cdir,ofile),quote=F,col.names=F,row.names=F)
	system(sprintf('%s --noweb --bfile %s/%sAGMHRSG --exclude %s/%s_snpstoremove.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%sAGMHRSGDCC',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile))
	system(sprintf('rm %s/%sAGMHRSGcasesmissing.* %s/%sAGMHRSGcontrolsmissing.*',cdir,ofile,cdir,ofile))
	printstats(cdir,sprintf('%sAGMHRSGDCC',ofile),'Filter for SNPs that show genotype-rate differenc between cases and controls',ofile)
	rm(casemiss,controlmiss)
} else {
	write.table(c(''),sprintf('%s/%s_snpstoremove.txt',cdir,ofile),quote=F,col.names=F,row.names=F)
	system(sprintf('%s --noweb --bfile %s/%sAGMHRSG --exclude %s/%s_snpstoremove.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%sAGMHRSGDCC',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile))
	printstats(cdir,sprintf('%sAGMHRSGDCC',ofile),'Dataset only contains either cases/controls so no check for SNPs with genotype-rate difference between cases and controls',ofile)
}
rm(famfile)
system(sprintf('rm %s/%s_snpstoremove.txt',cdir,ofile))

# B13. Filter for SNPs with Hardy-Weinberg Eequilibirum p-values < 1e-06 for controls or < 1e-10 for cases
famfile <- read.table(sprintf('%s/%sAGMHRSGDCC.fam',cdir,ofile),h=F,stringsAsFactors=F)
if( length(unique(famfile$V1)) == 1 ){
	system(sprintf('%s --noweb --bfile %s/%sAGMHRSGDCC --hwe 0.000001 --keep-allele-order --allow-no-sex --make-bed --out %s/%s_qc',plinkdir,cdir,ofile,cdir,ofile))
	printstats(cdir,sprintf('%s_qc',ofile),'Filter for SNPs with Hardy-Weinberg Eequilibirum p-values < 1e-06 for controls or < 1e-10 for cases',ofile)
} else {
	system(sprintf('%s --noweb --bfile %s/%sAGMHRSGDCC --keep-allele-order --allow-no-sex --make-bed --out %s/%s_qc',plinkdir,cdir,ofile,cdir,ofile))
	printstats(cdir,sprintf('%s_qc',ofile),'HWE is not performed due to different populations present in the dataset',ofile)
}
rm(famfile)
system(sprintf('rm %s/%sA.* %s/%sAG.* %s/%sAGM.* %s/%sAGMH.* %s/%sAGMHRS.* %s/%sAGMHRSG.* %s/%sAGMHRSGDCC.*',cdir,ofile,cdir,ofile,cdir,ofile,cdir,ofile,cdir,ofile,cdir,ofile,cdir,ofile))
system(sprintf('rm  %s/%sX.* %s/%s_missing_samples.txt',cdir,ofile,cdir,ofile))


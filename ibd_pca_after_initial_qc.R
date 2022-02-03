library('ggplot2')
library('stringr')
options(stringsAsFactors=FALSE)
source('~/Scripts/gwasfunctions.r')
plinkdir <- "~/Software/plink-1.9/plink"
smartpca <- '~/Software/EIG-6.1.4/bin/smartpca'
pcasoftware <- 'EIGENSOFT'
onekgdir <- "~/1000genomes_European"
onekgfile <- '1kgrsSelected' # 1000genomes with selected european populations cohort

######### ------------------------------------------------------ C. PCA - European Outliers, Population Outliers ------------------------------------------------------ #########

# Steps Involved in Ancesty Matching with 1000genomes cohort for European and Global population level
# C1. QC steps before Merging 1000genomes selected and european cohorts with Input cohort
#  C1.1 Extracting the common snps from 1000genomes cohort selected
#  C1.2 Extracting the common snps from Input
# C2. Remove SNPs that show platform diff in maf difference > 0.15, missingness diff > 0.02 and flip snps
#  C2.1 Remove SNPs that show platform diff in maf difference > 0.15, missingness diff > 0.02 and flip snps in 1000genomes selected and Input. Merge the cohorts.
# C3. Remove SNPs with geno > 0.02 , AT/CG and HEW fail
#  C3.1 Remove SNPs with geno > 0.02, AT/CG and HWE failed SNPs in merged 1000genomes selected and Input
# C4. Remove related samples
#  C4.1 Remove related samples in merged 1000genomes selected and Input
# C5. PCA/MDS 
#  C5.1 PCA/MDS on 1000genomes selected and Input
# C6. Plot MDS
#  C6.1 Plot MDS of 1000genomes selected and Input

# C1.QC steps before Merging 1000genomes selected and european cohorts with Input cohort
# C1.1 Extract commons snps between 1000genomes selcted, 1000genomes european and Input cohort using SNP annotation
system(sprintf('perl ~/Scripts/common_snps_1kg_othercohort.pl %s/%s_qc.bim %s/%s.bim',cdir,ofile,onekgdir,onekgfile))
system(sprintf('%s --noweb --bfile %s/%s_qc --update-name %s/%s_update_snps.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%sQCU',plinkdir,cdir,ofile,cdir,ofile,cdir,ofile))
commonfile <- read.table(sprintf('%s/%s_common_snps.txt',cdir,ofile),h=F,stringsAsFactors=F,blank.lines.skip = TRUE)
if(class(try(read.table(sprintf('%s/%s_update_snps.txt',cdir,ofile),h=F,stringsAsFactors=F,blank.lines.skip = TRUE)))=='try-error') { 
	updatefile <- as.data.frame(c(''))
	updatefile$V2 <- ''
} else {
	updatefile <- read.table(sprintf('%s/%s_update_snps.txt',cdir,ofile),h=F,stringsAsFactors=F,blank.lines.skip = TRUE)	
}
commonsnps <- c(as.vector(commonfile$V1),as.vector(updatefile$V2))
write.table(commonsnps,sprintf('%s/%s_common_snps.txt',cdir,ofile),quote=F,row.names=F,col.names=F)
printstats(cdir,sprintf('%sQCU',ofile),'Update the SNP information using 1000genomes reference',ofile)

# C1.2 Extracting the common snps from 1000genomes cohort selected
system(sprintf('echo \'CEU\' > %s/%s_filterpop.txt',onekgdir,ofile))
system(sprintf('%s --noweb --bfile %s/%s --extract %s/%s_common_snps.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%s_common',plinkdir,onekgdir,onekgfile,cdir,ofile,onekgdir,onekgfile))
system(sprintf('%s --noweb --bfile %s/%s_common --keep-fam %s/%s_filterpop.txt --freq --out %s/%sfreq',plinkdir,onekgdir,onekgfile,onekgdir,ofile,onekgdir,onekgfile))
system(sprintf('%s --noweb --bfile %s/%s_common --missing --out %s/%smiss',plinkdir,onekgdir,onekgfile,onekgdir,onekgfile))
system(sprintf('rm %s/%s_filterpop.txt',onekgdir,ofile))

# C1.3 Extracting the common snps from Input
system(sprintf('%s --noweb --bfile %s/%sQCU --extract %s/%s_common_snps.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%sQCU_common',plinkdir,cdir,ofile,cdir,ofile,onekgdir,ofile))
famfile <- read.table(sprintf('%s/%sQCU.fam',cdir,ofile),h=F,stringsAsFactors=F)
if( length(which(famfile$V6 %in% c(1) == TRUE)) > 0 ){
	system(sprintf('%s --noweb --bfile %s/%sQCU_common --filter-controls --freq --allow-no-sex --out %s/%sQCUfreq',plinkdir,onekgdir,ofile,onekgdir,ofile))
} else {
	system(sprintf('%s --noweb --bfile %s/%sQCU_common --freq --allow-no-sex --out %s/%sQCUfreq',plinkdir,onekgdir,ofile,onekgdir,ofile))
}
system(sprintf('%s --noweb --bfile %s/%sQCU_common --missing --allow-no-sex --out %s/%sQCUmiss',plinkdir,onekgdir,ofile,onekgdir,ofile))

# C2 Remove SNPs that show platform diff in maf difference > 0.15, missingness diff > 0.02 and flip snps
# C2.1 Remove SNPs that show platform diff in maf difference > 0.15, missingness diff > 0.02 and flip snps in 1000genomes selected and Input
famfile <- read.table(sprintf('%s/%sQCU.fam',cdir,ofile),h=F,stringsAsFactors=F)
if( length(which(famfile$V6 %in% c(1) == TRUE)) > 0 ){
	snpstoremove(onekgdir,sprintf("%sQCU",ofile),onekgdir,onekgfile,'domafdiff',0.2)	
} else {
	snpstoremove(onekgdir,sprintf("%sQCU",ofile),onekgdir,onekgfile,'dontdomafdiff',0.2)
}
system(sprintf('rm %s/%sfreq.* %s/%sQCUfreq.* %s/%smiss.* %s/%sQCUmiss.*',onekgdir,onekgfile,onekgdir,ofile,onekgdir,onekgfile,onekgdir,ofile))
system(sprintf('%s --noweb --bfile %s/%s_common --exclude %s/%sQCU_snpsremove.txt --make-bed --out %s/%s_temp',plinkdir,onekgdir,onekgfile,onekgdir,ofile,onekgdir,onekgfile))
system(sprintf('%s --noweb --bfile %s/%sQCU_common --exclude %s/%sQCU_snpsremove.txt --flip %s/%sQCU_flipsnps.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%sQCUF',plinkdir,onekgdir,ofile,onekgdir,ofile,onekgdir,ofile,onekgdir,ofile))
system(sprintf('%s --noweb --bfile %s/%sQCUF --bmerge %s/%s_temp --keep-allele-order --allow-no-sex --make-bed --out %s/%s1kgS',plinkdir,onekgdir,ofile,onekgdir,onekgfile,onekgdir,ofile))
system(sprintf('rm %s/%s_temp.* %s/%s_common.* %s/%sQCU_common.* %s/%sQCUF.* %s/%sQCU_flipsnps.txt %s/%sQCU_snpsremove.txt',onekgdir,onekgfile,onekgdir,onekgfile,onekgdir,ofile,onekgdir,ofile,onekgdir,ofile,onekgdir,ofile))
system(sprintf('rm %s/%s_common_snps.txt %s/%s_update_snps.txt',cdir,ofile,cdir,ofile))

# C3. Remove SNPs with geno > 0.02 and AT/CG
# C3.1 Remove SNPs with geno > 0.02 and AT/CG SNPs in merged 1000genomes selected and Input
system(sprintf('%s --noweb --bfile %s/%s1kgS --geno 0.02 --keep-allele-order --allow-no-sex --make-bed --out %s/%s1kgSG',plinkdir,onekgdir,ofile,onekgdir,ofile))
extractATCGsnps(onekgdir,sprintf('%s1kgSG',ofile))
system(sprintf('%s --noweb --bfile %s/%s1kgSG --exclude %s/%s1kgSG_ATorCGsnplist.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%s1kgSGS',plinkdir,onekgdir,ofile,onekgdir,ofile,onekgdir,ofile))
system(sprintf('rm %s/%s1kgSG.* %s/%s1kgSG_ATorCGsnplist.txt',onekgdir,ofile,onekgdir,ofile))

# C4 Remove related samples by pruning over maf 0.01 and removing SNPs from high LD-regions.txt
# C4.1 Remove related samples in merged 1000genomes selected and Input
system(sprintf('%s --noweb --bfile %s/%s1kgSGS --exclude ~/Scripts/high-LD-regions.txt -range --maf 0.05 --keep-allele-order --allow-no-sex --indep-pairwise 100 50 0.2 --out %s/%s1kgSLD',plinkdir,onekgdir,ofile,onekgdir,ofile))
system(sprintf('%s --noweb --bfile %s/%s1kgSGS --keep-allele-order --allow-no-sex --make-bed --extract %s/%s1kgSLD.prune.in --out %s/%s1kgSLD',plinkdir,onekgdir,ofile,onekgdir,ofile,onekgdir,ofile))
system(sprintf('%s --noweb --bfile %s/%s1kgSLD --Z-genome --min 0.1 --allow-no-sex --out %s/%s1kgSLDIBD',plinkdir,onekgdir,ofile,onekgdir,ofile))
system(sprintf('gunzip %s/%s1kgSLDIBD.genome.gz',onekgdir,ofile))
relationCheck <- read.table(sprintf("%s/%s1kgSLDIBD.genome",onekgdir,ofile),h=T,stringsAsFactors=F)
pdf(sprintf("~/Desktop/SNP/qualitycontrol/%s1kgIBDAfterQC.pdf",ofile), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
print(ggplot(relationCheck,aes(x = PI_HAT)) + xlab("RELATEDNESS") + geom_histogram(colour="black") + coord_cartesian(xlim = c(0, 1)))
dev.off()
temp <- as.data.frame(relationCheck[which(relationCheck$PI_HAT > 0.18),c(1:4,7:10)])
write.table(temp,sprintf('~/Desktop/%s_related_samples.txt',ofile),quote=F,sep='\t',row.names=F,col.names=T)
if( cohortname == 'COHORTNAME'){
	temp$SampleComb <- as.vector(paste(temp$IID1,temp$IID2,sep='-'))
	test <- as.data.frame(table(c(temp$IID1,temp$IID2)),stringsAsFactors=F)
	samplestoremove <- ''
	for(i in 1:(nrow(test)-2)){
		if( test[i,2] == 2 & test[i+1,2] == 2 & test[i+2,2] == 2 ){
			print(i)
			if( length(which( temp[,9] %in% c( paste(test[i,1],test[i+1,1],sep='-') , paste(test[i+1,1],test[i,1],sep='-'), paste(test[i,1],test[i+2,1],sep='-'), paste(test[i+2,1],test[i,1],sep='-'), paste(test[i+1,1],test[i+2,1],sep='-'), paste(test[i+2,1],test[i+1,1],sep='-')  )  )) > 1 ){
				samplestoremove <- c(samplestoremove,test[i,1],test[i+1,1])
			} else {
				i <- i + 1	
			}
		
		} else {
			i <- i + 1
		}
	}
	samplestoremove <- samplestoremove[-1]
	temp1 <- temp[-which(temp$IID1 %in% samplestoremove | temp$IID2 %in% samplestoremove),]
	test <- as.data.frame(table(c(temp1$IID1,temp1$IID2)),stringsAsFactors=F)
	temp2 <- temp1[ -which(temp1$IID1 %in% test[which(test[,2] == 2),1] | temp1$IID2 %in% test[which(test[,2] == 2),1]), ]
	samplestoremove <- c(samplestoremove,test[which(test[,2] == 2),1])
	samplestoremove <- c(samplestoremove,temp2$IID1)
	famfile <- read.table(sprintf('%s/%s1kgSLD.fam',onekgdir,ofile),h=F,stringsAsFactors=F)
	samplestoremove <- famfile[which(famfile$V2 %in% samplestoremove),c(1,2,1)]
	write.table(samplestoremove,sprintf("%s/%s_related_samples.txt",onekgdir,ofile),quote=F,row.names=F,col.names=F)
	rm(temp,relationCheck,temp1,temp2,samplestoremove,famfile)	
}
if( cohortname %in% c('3C','POPGEN','Controls','POPRES','GoNL')){
	write.table(temp[,c(3,4,3)],sprintf("%s/%s_related_samples.txt",onekgdir,ofile),quote=F,row.names=F,col.names=F)
	rm(temp,relationCheck)
}

system(sprintf('%s --noweb --bfile %s/%s1kgSGS --remove %s/%s_related_samples.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%s1kgSGSRRS',plinkdir,onekgdir,ofile,onekgdir,ofile,onekgdir,ofile))
system(sprintf('rm %s/%s1kgSLD.* %s/%s1kgSLDIBD.* %s/%s1kgSGS.*',onekgdir,ofile,onekgdir,ofile,onekgdir,ofile))

# C5. PCA/MDS by pruning over maf 0.01 and removing SNPs from high LD-regions.txt
# C5.1 PCA and MDS on 1000genomes selected and Input
system(sprintf('%s --noweb --bfile %s/%s1kgSGSRRS --exclude ~/Scripts/high-LD-regions.txt -range --maf 0.05 --indep-pairwise 100 50 0.2 --keep-allele-order --allow-no-sex --out %s/%s1kgSLD',plinkdir,onekgdir,ofile,onekgdir,ofile))
system(sprintf('%s --noweb --bfile %s/%s1kgSGSRRS --extract %s/%s1kgSLD.prune.in --keep-allele-order --allow-no-sex --make-bed --out %s/%s1kgSLD',plinkdir,onekgdir,ofile,onekgdir,ofile,onekgdir,ofile))
pcaeigensoft(onekgdir,sprintf('%s1kgSLD',ofile),smartpca,'bypop')

# C6 Plot MDS in 2D and 3D
# C6.1 Plot only 2D and the components in R
ggpcapop2d(onekgdir,ofile,cohortname,'selected',pcasoftware,sprintf('%s',ofile),'~/Desktop/SNP/europeanoutliers')



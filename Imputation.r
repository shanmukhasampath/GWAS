############################################################# Imputation #############################################
##################### IMPORTANT INFORMATION ################
Steps from 1 to 5 need to be run separately and only after the initial step is finished its job in the server.
Step 2 cannot be run before Step 1 is finished. Same rule applies till Step 5. 
############################################################
# Change the directories

plinkdir <- '~/Software/plink-1.9/plink' # Plink Directory
shapeit <- '~/Software/shapeit/bin/shapeit' # SHAPEIT directory for both strand alignment and phasing
impute <- '~/Software/impute/impute2' # IMPUTE software directory for imputation
ceudir <- '~/1000genomes/CEU' # Genetic map
onekgdir <- '~/1000genomes_European' # Directory of 1000 genomes reference panel
idir <- '~/INPUTDIRECTORY' # Input file directory
ifile <- 'INPUTFILENAME' # Input file name
odir <- '~/OUTPUTDIRECTORY' # Output file directory
ofile <- 'OUTPUTFILENAME' # Output file name
seedgen <- 13579 # Random Seed
khap <- 700     # Number of Cases in your dataset. 
infoscore <- 0.8 # Post Imputation score for selection of SNPs with good quality
system('echo "CEU" > ~/Scripts/group.list')

###################### 1. split the input file by Chromsome
coebim <- read.table(sprintf('%s/%s.bim',idir,ifile),h=F,stringsAsFactors=F)
for(i in 1:22){
	if(i == 1){
		system(sprintf('%s --bfile %s/%s --from %s --to %s --allow-no-sex --keep-allele-order --make-bed --out %s/%s_chr%s',plinkdir,idir,ifile,coebim[1,2],coebim[which(coebim$V1 == i)[length(which(coebim$V1 == i))],2],odir,ifile,i))
	} else {
		system(sprintf('%s --bfile %s/%s --from %s --to %s --allow-no-sex --keep-allele-order --make-bed --out %s/%s_chr%s',plinkdir,idir,ifile,coebim[which(coebim$V1 == i)[1],2],coebim[which(coebim$V1 == i)[length(which(coebim$V1 == i))],2],odir,ifile,i))
	}
}
rm(coebim)

######################## 2. Strand Alignment done per Chromsome using SHAPEIT software
for(i in 1:22){
	sink(sprintf('~/Scripts/%s_strandaligncheck%s.in',ifile,i))
	cat(sprintf('# FILENAME: %s_strandaligncheck%s.in\n',ifile,i))
	cat('# Strand Alignment R script\n')
	cat('options(stringsAsFactors=F)\n')
	cat(sprintf('system(\'%s -check -B %s/%s_chr%s -M %s/CEU_chr%s.txt --input-ref %s/1000GP_Phase3_chr%s.hap.gz %s/1000GP_Phase3_chr%s.legend.gz %s/1000GP_Phase3.sample --include-grp ~/Scripts/group.list --thread 20 --output-log %s/chr%s.alignments\')\n',shapeit,odir,ifile,i,ceudir,i,onekgdir,i,onekgdir,i,onekgdir,odir,i))
	cat(sprintf('system(\'awk \\\'NR > 1 {print $4}\\\' %s/chr%s.alignments.snp.strand | sort | uniq > %s/flipsnps%s.txt\')\n',odir,i,odir,i))
	cat(sprintf('system(\'%s --bfile %s/%s_chr%s --flip %s/flipsnps%s.txt --allow-no-sex --keep-allele-order --make-bed --out %s/%s_chr%s_aligned\')\n',plinkdir,odir,ifile,i,odir,i,odir,ifile,i))
	cat(sprintf('system(\'%s -check -B %s/%s_chr%s_aligned -M %s/CEU_chr%s.txt --input-ref %s/1000GP_Phase3_chr%s.hap.gz %s/1000GP_Phase3_chr%s.legend.gz %s/1000GP_Phase3.sample --include-grp ~/Scripts/group.list --thread 20 --output-log %s/chr%s.alignments\')\n',shapeit,odir,ifile,i,ceudir,i,onekgdir,i,onekgdir,i,onekgdir,odir,i))
	cat(sprintf('system(\'awk \\\'NR > 1 {print $4}\\\' %s/chr%s.alignments.snp.strand | sort | uniq > %s/flipsnps%s.txt\')\n',odir,i,odir,i))
	cat(sprintf('system(\'%s --bfile %s/%s_chr%s_aligned --exclude %s/flipsnps%s.txt --allow-no-sex --keep-allele-order --make-bed --out %s/%s_chr%s\')\n',plinkdir,odir,ifile,i,odir,i,odir,ifile,i))
	cat(sprintf('system(\'%s -check -B %s/%s_chr%s -M %s/CEU_chr%s.txt --input-ref %s/1000GP_Phase3_chr%s.hap.gz %s/1000GP_Phase3_chr%s.legend.gz %s/1000GP_Phase3.sample --include-grp ~/Scripts/group.list --thread 20 --output-log %s/chr%s.alignments\')\n',shapeit,odir,ifile,i,ceudir,i,onekgdir,i,onekgdir,i,onekgdir,odir,i))
	sink()
	sink(sprintf('~/Scripts/%s_strandaligncheck%s.sub',ifile,i))
	cat('#!/bin/sh -l\n')
	cat(sprintf('# FILENAME:  %s_strandaligncheck%s.sub\n',ifile,i))
	cat('\n')
	cat('module load r\n')
	cat('# --vanilla:\n')
	cat('# --no-save: do not save datasets at the end of an R session\n')
	cat(sprintf('R --vanilla --no-save < ~/Scripts/%s_strandaligncheck%s.in\n',ifile,i))
	sink()
	system(sprintf('qsub -q SGEqueueName -l nodes=1:ppn=1,naccesspolicy=singleuser -l walltime=4:00:00 -N StrandAlign ~/Scripts/%s_strandaligncheck%s.sub',ifile,i))
}
system(sprintf('rm %s/*aligned.* %s/chr* ~/Scripts/%s_strandaligncheck* %s/flipsnps*.txt',odir,odir,ifile,odir))

######################## 3. Phasing each chromosome file using SHAPEIT 
for(i in 1:22){
	sink(sprintf('~/Scripts/%s_chr_%s.sh',ifile,i))
	cat(sprintf('#!/bin/bash\n'))
	cat(sprintf('\n'))
	cat(sprintf('%s --seed %s -B %s/%s_chr%s -M %s/CEU_chr%s.txt -O %s/%s_chr%s_phased -T 4\n',shapeit,seedgen,odir,ifile,i,ceudir,i,odir,ifile,i))
	sink()
	system(sprintf('qsub -q SGEqueueName -l nodes=1:ppn=1,naccesspolicy=singleuser -l walltime=50:00:00 -N Phasing ~/Scripts/%s_chr_%s.sh',ifile,i))
}
system(sprintf('rm %s/%s_chr*.{bed,bim,fam,nosex,log}',odir,ifile))

######################## 4. Impute each chromosome
for(i in 1:22){
	system(sprintf('awk \'NR==1 { print $3 }\' %s/%s_chr%s_phased.haps > %s/firstline.txt',odir,ifile,i,odir))
	system(sprintf('awk \'END { print $3 }\' %s/%s_chr%s_phased.haps > %s/lastline.txt',odir,ifile,i,odir))
	fline <- read.table(sprintf('%s/firstline.txt',odir),h=F,stringsAsFactors=F)
	lline <- read.table(sprintf('%s/lastline.txt',odir),h=F,stringsAsFactors=F)
	begining <- floor(fline[,1]/1000000)
	ending <- ceiling(lline[,1]/1000000)
	k <- ceiling(abs(begining-ending)/5)
	if(begining > 0){
		m <- begining + 5
	} else {
		m <- 5
	}
	sink(sprintf('~/Scripts/%s_chr%s_impute.sh',ifile,i))
	cat('#!/bin/bash\n')
	for(j in 1:k){
		if (j == 1){
			cat(sprintf('%s -o_gz -filt_rules_l \'EUR==0\' \'TYPE!=Biallelic_SNP\' -use_prephased_g -known_haps_g %s/%s_chr%s_phased.haps -h %s/1000GP_Phase3_chr%s.hap.gz -l %s/1000GP_Phase3_chr%s.legend.gz -m %s/genetic_map_chr%s_combined_b37.txt -int %s %se6 -k_hap %s -Ne 20000 -buffer 1000 -o %s/%s_chr%s_phased_imputed_chunk%s\n',impute,odir,ifile,i,onekgdir,i,onekgdir,i,onekgdir,i,fline[1,1],m,khap,odir,ifile,i,j))
			l <- begining+5
			m <- m + 5
		} else {
			cat(sprintf('%s -o_gz -filt_rules_l \'EUR==0\' \'TYPE!=Biallelic_SNP\' -use_prephased_g -known_haps_g %s/%s_chr%s_phased.haps -h %s/1000GP_Phase3_chr%s.hap.gz -l %s/1000GP_Phase3_chr%s.legend.gz -m %s/genetic_map_chr%s_combined_b37.txt -int %s.000001e6 %se6 -k_hap %s -Ne 20000 -buffer 1000 -o %s/%s_chr%s_phased_imputed_chunk%s\n',impute,odir,ifile,i,onekgdir,i,onekgdir,i,onekgdir,i,l,m,khap,odir,ifile,i,j))
			l <- l+5
			m <- m + 5
		}
	}
	sink()
	system(sprintf('qsub -q SGEqueueName -l nodes=1:ppn=1,naccesspolicy=singleuser -l walltime=100:00:00 -N Impute ~/Scripts/%s_chr%s_impute.sh',ifile,i))
}
system(sprintf('rm ~/Scripts/%s_*chunk*.sh',ifile))

######################### 5. Post Imputation to get SNPs imputed quality information in terms of Infoscore. Visualization of imputation quality.
system(sprintf('mkdir %s/post_imputed',odir))
sink(sprintf('~/Scripts/%s_postimpute.in',ifile))
cat(sprintf('# FILENAME: %s_postimpute.in\n',ifile))
cat('# Post Imputation R script\n')
cat('options(stringsAsFactors=F)\n')
cat(sprintf('plinkdir <- \'%s\'\n',plinkdir))
cat(sprintf('shapeit <- \'%s\'\n',shapeit))
cat(sprintf('impute <- \'%s\'\n',impute))
cat(sprintf('idir <- \'%s\'\n',idir))
cat(sprintf('ifile <- \'%s\'\n',ifile))
cat(sprintf('odir <- \'%s\'\n',odir))
cat(sprintf('ofile <- \'%s\'\n',ofile))
cat(sprintf('seedgen <- %s\n',seedgen))
cat(sprintf('khap <- %s     # Number of Cases in your dataset.\n',khap))
cat(sprintf('infoscore <- %s # Post Imputation score for selection of SNPs with good quality\n',infoscore))

cat('for(i in 1:22){\n')

	cat('system(sprintf(\'cat %s/%s_chr%s_phased_imputed_chunk*.gz > %s/post_imputed/%s_chr%s_phased_imputed_chunk.gz\',odir,ifile,i,odir,ifile,i))\n')
	cat('system(sprintf(\'cp %s/%s_chr%s_phased.sample %s/post_imputed/%s_chr%s_phased.sample\',odir,ifile,i,odir,ifile,i))\n')
	cat('imputefiles <- list.files(sprintf(\'%s\',odir),pattern=sprintf(\'^(.*)_chr%s_(.*)\\\\.gz$\',i),full.names=T)\n')
	cat('temp <- as.numeric(gsub(\'^(.*)_chunk(.*?)\\\\.gz$\',\'\\\\2\',imputefiles,perl=T))\n')
	cat('imputefiles <- imputefiles[order(temp)]\n')
	
	cat('infofiles <- list.files(sprintf(\'%s\',odir),pattern=sprintf(\'^(.*)_chr%s_(.*)_info$\',i),full.names=T)\n')
	cat('temp <- as.numeric(gsub(\'^(.*)_chunk(.*?)_info$\',\'\\\\2\',infofiles,perl=T))\n')
	cat('infofiles <- infofiles[order(temp)]\n')
	cat('rm(temp)\n')
	
	cat('# Get information on number of imputed and genotyped SNPs for All, InfoScore > 0.8, InfoScore > 0.7 and InfoScore > 0.3\n')
	cat('for(j in 1:length(imputefiles)){\n')
		cat('if(j == 1){\n')
			cat('infodat <- read.table(infofiles[j],h=T)\n')
		cat('}\n')
		cat('else {\n')
			cat('temp <- read.table(infofiles[j],h=T)\n')
			cat('infodat <- rbind(infodat,temp)\n')
		cat('}\n')		
	cat('}\n')
	cat('write.table(infodat,sprintf(\'%s/post_imputed/%s_chr%s_phased_imputed_chunk_info\',odir,ifile,i),quote=F,row.names=F,col.names=T)\n')
	
	cat('if (i == 1){\n')
		cat('tot <- nrow(infodat)\n')
		cat('geno <- length(which(infodat$type > 0))\n')
		cat('imp <- length(which(infodat$type == 0))\n')
		cat('write.table(paste(tot,paste(geno,imp,sep="/"),sep=" ("),sprintf(\'%s/post_imputed/stats.txt\',odir),quote=F,row.names=F,col.names=F)\n')
		
		cat('tot <- length(which(infodat$info > 0.8))\n')
		cat('geno <- length(which(infodat$info > 0.8 & infodat$type > 0))\n')
		cat('imp <- length(which(infodat$info > 0.8 & infodat$type == 0))\n')
		cat('write.table(paste(tot,paste(geno,imp,sep="/"),sep=" ("),sprintf(\'%s/post_imputed/stats_info8.txt\',odir),quote=F,row.names=F,col.names=F)\n')
			
		cat('tot <- length(which(infodat$info > 0.7))\n')
		cat('geno <- length(which(infodat$info > 0.7 & infodat$type > 0))\n')
		cat('imp <- length(which(infodat$info > 0.7 & infodat$type == 0))\n')
		cat('write.table(paste(tot,paste(geno,imp,sep="/"),sep=" ("),sprintf(\'%s/post_imputed/stats_info7.txt\',odir),quote=F,row.names=F,col.names=F)\n')
			
		cat('tot <- length(which(infodat$info > 0.3))\n')
		cat('geno <- length(which(infodat$info > 0.3 & infodat$type > 0))\n')
		cat('imp <- length(which(infodat$info > 0.3 & infodat$type == 0))\n')
		cat('write.table(paste(tot,paste(geno,imp,sep="/"),sep=" ("),sprintf(\'%s/post_imputed/stats_info3.txt\',odir),quote=F,row.names=F,col.names=F)\n')
		cat('rm(tot,geno,imp)\n')

	cat('} else {\n')
		
		cat('tot <- nrow(infodat)\n')
		cat('geno <- length(which(infodat$type > 0))\n')
		cat('imp <- length(which(infodat$type == 0))\n')
		cat('write.table(paste(tot,paste(geno,imp,sep="/"),sep=" ("),sprintf(\'%s/post_imputed/stats.txt\',odir),quote=F,row.names=F,col.names=F,append=T)\n')
		
		cat('tot <- length(which(infodat$info > 0.8))\n')
		cat('geno <- length(which(infodat$info > 0.8 & infodat$type > 0))\n')
		cat('imp <- length(which(infodat$info > 0.8 & infodat$type == 0))\n')
		cat('write.table(paste(tot,paste(geno,imp,sep="/"),sep=" ("),sprintf(\'%s/post_imputed/stats_info8.txt\',odir),quote=F,row.names=F,col.names=F,append=T)\n')
		
		cat('tot <- length(which(infodat$info > 0.7))\n')
		cat('geno <- length(which(infodat$info > 0.7 & infodat$type > 0))\n')
		cat('imp <- length(which(infodat$info > 0.7 & infodat$type == 0))\n')
		cat('write.table(paste(tot,paste(geno,imp,sep="/"),sep=" ("),sprintf(\'%s/post_imputed/stats_info7.txt\',odir),quote=F,row.names=F,col.names=F,append=T)\n')
			
		cat('tot <- length(which(infodat$info > 0.3))\n')
		cat('geno <- length(which(infodat$info > 0.3 & infodat$type > 0))\n')
		cat('imp <- length(which(infodat$info > 0.3 & infodat$type == 0))\n')
		cat('write.table(paste(tot,paste(geno,imp,sep="/"),sep=" ("),sprintf(\'%s/post_imputed/stats_info3.txt\',odir),quote=F,row.names=F,col.names=F,append=T)\n')
		cat('rm(tot,geno,imp)\n')
	cat('}\n')
	cat('rm(temp,infodat)\n')
	cat('gc()\n')
cat('}\n')

cat('## Visualization of MAF and InfoScore\n')
cat('library(\'ggplot2\')\n')
cat('library(\'plotrix\')\n')
cat('library(\'grid\')\n')
cat('for(i in 1:22){\n')
	
	cat('if(i == 1){\n')
		cat('infodat <- read.table(sprintf(\'%s/post_imputed/%s_chr%s_phased_imputed_chunk_info\',odir,ifile,i),h=T,stringsAsFactors=F)\n')
	cat('} else {\n')
		cat('temp <- read.table(sprintf(\'%s/post_imputed/%s_chr%s_phased_imputed_chunk_info\',odir,ifile,i),h=T,stringsAsFactors=F)\n')
		cat('infodat <- rbind(infodat,temp)\n')
	cat('}\n')
cat('}\n')
cat('rm(temp)\n')
cat('infodat$MAF <- infodat$exp_freq_a1\n')
cat('infodat[which(infodat$MAF > 0.5),\'MAF\'] <- 1- infodat[which(infodat$MAF > 0.5),\'MAF\']\n')
cat(sprintf('png(\'~/Plots/%s_imputed_infoscore_histogram.png\',units="in", width=8, height=8, res=150)\n',ifile))
cat('ggplot(infodat[which(infodat$type == 0),], aes(x=info)) + geom_histogram(binwidth=.01, colour="black", fill="darkgreen") + geom_vline(aes(xintercept=0.3),color="black", linetype="dashed", size=1) + geom_vline(aes(xintercept=0.7),color="black", linetype="dashed", size=1) + geom_vline(aes(xintercept=0.8),color="black", linetype="dashed", size=1) + xlab("INFO SCORE") + ylab(sprintf("Imputed SNPs (%s)",length(which(infodat$type == 0)))) + ggtitle("INFO score distribution in Imputed SNPs in Danish Cohort")\n')
cat('dev.off()\n')

cat('mafseq <- as.data.frame(seq(0,0.5,by=0.01),stringsAsFactors=F)\n')
cat('mafseq$meaninfo <- 0\n')
cat('mafseq$imputesnps <- 0\n')
cat('colnames(mafseq)[1] <- \'MAFbin\'\n')
cat('for(i in 1:nrow(mafseq)){\n')
	cat('if(i == 1){\n')
		cat('mafseq[i,3] <- length(which(infodat$MAF < mafseq[i+1,1] & infodat$type == 0))\n')
		cat('mafseq[i,2] <- mean(infodat[which(infodat$MAF < mafseq[i+1,1] & infodat$type == 0),\'info\'])\n')
	cat('} else if(i != 1 & i < nrow(mafseq)){\n')
		cat('mafseq[i,3] <- length(which( (infodat$MAF >= mafseq[i,1] & infodat$MAF < mafseq[i+1,1]) & infodat$type == 0 ))\n')
		cat('mafseq[i,2] <- mean(infodat[which( (infodat$MAF >= mafseq[i,1] & infodat$MAF < mafseq[i+1,1]) & infodat$type == 0 ),\'info\'])\n')
	cat('} else {\n')
		cat('mafseq[i,3] <- length(which( infodat$MAF == mafseq[i,1] & infodat$type == 0 ))\n')
		cat('mafseq[i,2] <- mean(infodat[which(infodat$MAF == mafseq[i,1] & infodat$type == 0),\'info\'])\n')
	cat('}\n')
cat('}\n')

cat('p1 <- ggplot(mafseq, aes(x=MAFbin,y=meaninfo)) + geom_point( colour = "darkgreen" ) + scale_shape_identity(2) + scale_shape(solid = TRUE) + geom_line(linetype="dashed",colour=\'darkgreen\') +  xlab("MAF Bins") + ylab("Mean INFO score") + ggtitle("INFO score per MAF bin")\n')
cat('p2 <- ggplot(mafseq, aes(x=MAFbin,y=imputesnps)) + geom_point( colour = "darkblue" ) + scale_shape_identity(1) + scale_shape(solid = TRUE) + geom_line(linetype="dashed",colour=\'darkblue\') +  xlab("MAF Bins") + ylab("Number of Imputed SNPs") + ggtitle("Imputed SNPs per MAF bin")\n')
cat(sprintf('png(\'~/Plots/%s_imputed_infoscore_MAF_ggplot.png\',units="in", width=10, height=8, res=150)\n',ifile))
cat('pushViewport(viewport(layout = grid.layout(1, 2)))\n')
cat('print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))\n')
cat('print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))\n')
cat('dev.off()\n')

cat(sprintf('png(\'~/Plots/%s_imputed_infoscore_MAF_plotrix.png\',units="in", width=10, height=10, res=150)\n',ifile))
cat('twoord.plot(mafseq[,1],mafseq[,2],mafseq[,1],mafseq[,3],xlab="MAF bins of interval 0.01",ylab="mean INFO score per MAF bin",rylab="Imputed SNPs per MAF bin",lcol=4,rcol=2,lpch=16,rpch=17,lty=2,main="INFO score and Imputed SNPs per MAF bin in Danish Controls",do.first="plot_bg();grid(col=\\"white\\",lty=1)")\n')
cat('dev.off()\n')

cat('write.table(infodat[which(infodat$MAF > 0.01 & infodat$info > sprintf(\'%s\',infoscore)),2],sprintf(\'%s/snps_toinclude.txt\',odir),quote=F,row.names=F,col.names=F)\n')
sink()
sink(sprintf('~/Scripts/%s_postimpute.sub',ifile))
cat('#!/bin/sh -l\n')
cat(sprintf('# FILENAME:  %s_postimpute.sub\n',ifile))
cat('\n')
cat('module load r\n')
cat('# --vanilla:\n')
cat('# --no-save: do not save datasets at the end of an R session\n')
cat(sprintf('R --vanilla --no-save < ~/Scripts/%s_postimpute.in\n',ifile))
sink()
system(sprintf('qsub -q SGEqueueName -l nodes=1:ppn=1,naccesspolicy=singleuser -l walltime=50:00:00 -N PostImpute ~/Scripts/%s_postimpute.sub',ifile))
system(sprintf('rm ~/Scripts/%s_postimpute.sub ~/Scripts/%s_postimpute.in',ifile,ifile))

######################### 6. Converting post-imputation output GEN to PLINK PED and MAP files
for (i in 1:22){
	system(sprintf('%s --gen %s/post_imputed/%s_chr%s_phased_imputed_chunk.gz --sample %s/post_imputed/%s_chr%s_phased.sample --oxford-pheno-name plink_pheno --oxford-single-chr %s --hard-call-threshold 0.1 --keep-allele-order --allow-no-sex --make-bed --out %s/%s_chr%s',plinkdir,odir,ifile,i,odir,ifile,i,i,odir,ifile,i))
}
system(sprintf('rm %s/*chunk*',odir))

######################### 7. Merge Imputed Chromosomal PLINK files of a cohort into one PLINK file
beddat <- list.files(sprintf('%s',odir),pattern='*_chr[0-9]*.bed',full.names=T)
bimdat <- list.files(sprintf('%s',odir),pattern='*_chr[0-9]*.bim',full.names=T)
famdat <- list.files(sprintf('%s',odir),pattern='*_chr[0-9]*.fam',full.names=T) 
flist <- cbind(beddat,cbind(bimdat,famdat))
temp <- as.numeric(gsub('^(.*)_chr(.*?)\\.bed$','\\2',flist[,1],perl=T))
flist <- flist[order(temp),]
write.table(flist[2:nrow(flist),],sprintf('%s/allfiles.txt',odir),quote=F,sep=' ',row.names=F,col.names=F)

sink(sprintf('~/Scripts/%s_combine_chrs.sh',ifile))
cat('#!/bin/bash\n')
cat(sprintf('%s --bfile %s/%s_chr1 --merge-list %s/allfiles.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%stemp\n',plinkdir,odir,ifile,odir,odir,ofile))
sink()
system(sprintf('qsub -q SGEqueueName -l nodes=1:ppn=1,naccesspolicy=singleuser -l walltime=50:00:00 -N PostImpute ~/Scripts/%s_combine_chrs.sh',ifile))

# Change the FAM file Phenotype to 1 (Control) and 2 (Case) 
fam1 <- read.table(sprintf('%s/%stemp.fam',odir,ofile),h=F,stringsAsFactors=F)
fam1[which(fam1$V6 == 2),6] <- 1
fam1[which(fam1$V6 == -9),6] <- 2
write.table(fam1,sprintf('%s/%stemp.fam',odir,ofile),quote=F,row.names=F,col.names=F)

######################### 8. Filter out the SNPs that failed MAF < 0.01 & INFO score < 0.8 and Cleanup temporary imputation files
system(sprintf('grep \'^Warning\' %s/%stemp.log  | awk \'{print $5}\' > %s/snpstoremove.txt',odir,ofile,odir))
data <- read.table(sprintf('%s/snpstoremove.txt',odir),h=F,stringsAsFactors=F)
write.table(data[,1],sprintf('%s/snpstoremove.txt',odir),quote=F,row.names=F,col.names=F)
system(sprintf('%s --bfile %s/%stemp --extract %s/snps_toinclude.txt --exclude %s/snpstoremove.txt --keep-allele-order --allow-no-sex --make-bed --out %s/%s',plinkdir,odir,ofile,odir,odir,odir,ofile))
system(sprintf('rm %s/%stemp.* %s/%s_chr*.{bed,bim,fam,nosex,log}',odir,ofile,odir,ifile))

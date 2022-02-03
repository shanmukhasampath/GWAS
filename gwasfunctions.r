library('ggplot2')
options(stringsAsFactors=FALSE)

extractAutosomal <- function(currentdir,bimfile){
	data <- read.table(sprintf('%s/%s.bim',currentdir,bimfile),h=F,stringsAsFactors=F)
	write.table(data[which(data$V1 < 23),2],sprintf('%s/%s_AutosomalSnps.txt',currentdir,bimfile),quote=F,row.names=F,col.names = F)
	if( length(which(data$V1 ==23)) > 0 ) {
		write.table(data[which(data$V1 == 23),2],sprintf('%s/%s_SexChromosomeSnps.txt',currentdir,bimfile),quote=F,row.names=F,col.names = F)	
	} else {
		temp <- c('NoXChr')
		write.table(temp,sprintf('%s/%s_SexChromosomeSnps.txt',currentdir,bimfile),quote=F,row.names=F,col.names = F)
	}
}

extractATCGsnps <- function(currentdir,bimname){
	bimfile <- read.table(sprintf("%s/%s.bim",currentdir,bimname),h=F,stringsAsFactors=F)
	bimfile$Alleles <- as.vector(paste(bimfile$V5,bimfile$V6,sep=""))
	ATorCGsnps <- bimfile[which(bimfile$Alleles %in% c('AT','TA','GC','CG')),2]
	write.table(ATorCGsnps,sprintf("%s/%s_ATorCGsnplist.txt",currentdir,bimname),quote=F,row.names=F,col.names=F)
	rm(bimfile,ATorCGsnps)
}

printstats <- function(currentdir,bimname,stepname,cohortname){

	bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
	bimfile$V7 <- as.vector(paste(bimfile$V1,bimfile$V4,sep='-'))
	temp <- bimfile[which(bimfile$V1 != 0),c(2,7)]
	famfile <- read.table(sprintf('%s/%s.fam',currentdir,bimname),h=F,stringsAsFactors=F)

	sink(sprintf('%s/%s_stats.txt',currentdir,cohortname),append=T)
	cat("=============================\n")
	cat(sprintf("%s",stepname),"\n")
	cat("=============================\n")

	cat('No:of SNPs ',dim(bimfile)[1],"\n") # No:of SNPs
	print('Different SNP names:')
	print(table(gsub('^([A-za-z]*)(.*?)$','\\1',bimfile$V2,perl=T))) # Type of SNP names
	print('Type of SNPs with their Number')
	print(table(paste(bimfile$V5,bimfile$V6,sep=""))) # Type of SNPs with their number.
	cat('No:of SNPs with duplicated names ',length(which(duplicated(bimfile$V2,fromLast=T))),"\n") # No:of SNPs with duplicated names
	cat('No:of SNPs with duplicated positions ',length(which(duplicated(temp$V7,fromLast=T))),"\n") # No:of SNPs with duplicated positions
	print('No:of SNPs per Chromosome')
	print(table(bimfile$V1)) # No:of SNPs per Chromosome

	cat('No:of Samples ',dim(famfile)[1],"\n") # No:of Samples
	print('Table of Populations:')
	print(table(famfile$V1)) # No:of Populations
	print('No:of Cases/Controls : 1 - Control, 2 - Case, -9/0 - Unknown')
	print(table(famfile$V6)) # No:of Cases/Controls
	print('Table of Sex status:  1 - Male, 2 - Female, -9/0 - Unknown')
	print(table(famfile$V5)) # Table of Sex Status
	print('Table of Sex status per Population:')
	print(table(famfile$V1,famfile$V5)) # Table of Sex Status per Population
	print('Table of Case/Control status per Population:')
	print(table(famfile$V1,famfile$V6)) # Table of Sex Status per Population

	sink()

}

bimcheck <- function(currentdir,bimname){

	bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
	bimfile$V7 <- as.vector(paste(bimfile$V1,bimfile$V4,sep='-'))
	temp <- bimfile[which(bimfile$V1 != 0),c(2,7)]

	write.table(bimfile[which(bimfile$V1 > 23 | bimfile$V1 == 0),2],sprintf('%s/%s_snpsremove.txt',currentdir,bimname),quote=F,row.names=F,col.names=F) # Remove SNPs with Chromosome 0, Y, XY, MT
	write.table(bimfile[which(duplicated(bimfile$V2,fromLast=T)),2],sprintf('%s/%s_snpsremove.txt',currentdir,bimname),quote=F,row.names=F,col.names=F,append=T) # Remove SNPs with duplicated names
	write.table(bimfile[which(bimfile$V5 %in% c("I","D") | bimfile$V6 %in% c("I","D")),2],sprintf('%s/%s_snpsremove.txt',currentdir,bimname),quote=F,row.names=F,col.names=F,append=T) # Remove Indels and Deletions
	write.table(temp[which(duplicated(temp$V7,fromLast=T)),1],sprintf('%s/%s_snpsremove.txt',currentdir,bimname),quote=F,col.names=F,row.names=F,append=T) # Remove SNPs with duplicated chr-positions

}

checksnps <- function(currentdir1,bimname1,currentdir2,bimname2){

	bimfile1 <- read.table(sprintf('%s/%s.bim',currentdir1,bimname1),h=F,stringsAsFactors=F)
	bimfile2 <- read.table(sprintf('%s/%s.bim',currentdir2,bimname2),h=F,stringsAsFactors=F)

	bimfile1$Alleles <- paste(paste(paste(bimfile1$V5,bimfile1$V6,sep=''),paste(bimfile2$V5,bimfile2$V6,sep=''),sep=''))
	bimfile1$Count <- sapply(bimfile1$Alleles,function(x) length(unique(strsplit(x,"")[[1]])))
	triallelicsnpstoremove <- bimfile1[which(bimfile1$Count == 3 & !grepl('0',bimfile1$Alleles)),2] ### SNPs that need to be excluded which are tri-allelic
	write.table(triallelicsnpstoremove,sprintf("%s/%s_snpsremove.txt",currentdir1,bimname1),quote=F,row.names=F,col.names=F)	

	snpstoflip <- bimfile1[which(bimfile1$Count == 4 & grepl('0',bimfile1$Alleles)),2] ### SNPs that need to be flipped
	snpstoflip <- c(snpstoflip,bimfile1[which(bimfile1$Count == 4 & !grepl('0',bimfile1$Alleles)),2]) ### SNPs that need to be flipped
	write.table(snpstoflip,sprintf("%s/%s_flipsnps.txt",currentdir1,bimname1),quote=F,row.names=F,col.names=F)

}

checkmafdiff <- function(currentdir1,bimname1,currentdir2,bimname2,mafdifference,appendoperator){

	bimfile1 <- read.table(sprintf('%s/%s.frq',currentdir1,bimname1),h=T,stringsAsFactors=F)
	rownames(bimfile1) <- bimfile1[,2]
	colnames(bimfile1) <- c('hc_chr','hc_snp','hc_A1','hc_A2','hc_MAF','hc_NCHROBS')
	bimfile2 <- read.table(sprintf('%s/%s.frq',currentdir2,bimname2),h=T,stringsAsFactors=F)
	rownames(bimfile2) <- bimfile2[,2]
	colnames(bimfile2) <- c('onekg_chr','onekg_snp','onekg_A1','onekg_A2','onekg_MAF','onekg_NCHROBS')
	bimfile2 <- bimfile2[rownames(bimfile1),]
	bimfile <- merge(bimfile1,bimfile2,by="row.names",sort=FALSE)
	bimfile$Row.names <- NULL
	bimfile$onekg_chr <- NULL
	bimfile$onekg_snp <- NULL
	bimfile <- bimfile[,c(1:4,7:8,5,9,6,10)]
	bimfile$MAFdiff <- abs(bimfile$hc_MAF - bimfile$onekg_MAF)
	snpsfrqfail <- bimfile[which(bimfile$MAFdiff > mafdifference),2] # SNPs that need to be excluded whihc failed the MAF difference
	if(appendoperator == "Append"){
		write.table(snpsfrqfail,sprintf('%s/%s_snpsmaffail.txt',currentdir1,bimname1),quote=F,row.names=F,col.names=F,append=T)	
	} else {
		write.table(snpsfrqfail,sprintf('%s/%s_snpsmaffail.txt',currentdir1,bimname1),quote=F,row.names=F,col.names=F)
	}	

}

snpstoremove <- function(cohortdir,cohortname,onekgdir,onekg,mafdiffstatus,mafdifference){
	hcfrq <- read.table(sprintf("%s/%sfreq.frq",cohortdir,cohortname),h=T,stringsAsFactors=F)
	rownames(hcfrq) <- hcfrq[,2]
	colnames(hcfrq) <- c('hc_chr','hc_snp','hc_A1','hc_A2','hc_MAF','hc_NCHROBS')
	onekfrq <- read.table(sprintf("%s/%sfreq.frq",onekgdir,onekg),h=T,stringsAsFactors=F)
	rownames(onekfrq) <- onekfrq[,2]
	colnames(onekfrq) <- c('onekg_chr','onekg_snp','onekg_A1','onekg_A2','onekg_MAF','onekg_NCHROBS')
	onekfrq <- onekfrq[rownames(hcfrq),]
	hconekfrq <- merge(hcfrq,onekfrq,by="row.names",sort=FALSE)
	hconekfrq$Row.names <- NULL
	hconekfrq$onekg_chr <- NULL
	hconekfrq$onekg_snp <- NULL
	hconekfrq <- hconekfrq[,c(1:4,7:8,5,9,6,10)]
	hconekfrq$MAFdiff <- abs(hconekfrq$hc_MAF - hconekfrq$onekg_MAF)
	snpsfrqfail <- hconekfrq[which(hconekfrq$MAFdiff > mafdifference),2] # SNPs that need to be excluded whihc failed the MAF difference

	hconekfrq$Alleles <- paste(paste(paste(hconekfrq$hc_A1,hconekfrq$hc_A2,sep=''),paste(hconekfrq$onekg_A1,hconekfrq$onekg_A2,sep=''),sep=''))
	hconekfrq$Count <- sapply(hconekfrq$Alleles,function(x) length(unique(strsplit(x,"")[[1]])))
	triallelicsnpstoremove <- hconekfrq[which(hconekfrq$Count == 3 & !grepl('0',hconekfrq$Alleles)),2] ### SNPs that need to be excluded which are tri-allelic

	snpstoflip <- hconekfrq[which(hconekfrq$Count == 4 & grepl('0',hconekfrq$Alleles)),2] ### SNPs that need to be flipped
	snpstoflip <- c(snpstoflip,hconekfrq[which(hconekfrq$Count == 4 & !grepl('0',hconekfrq$Alleles)),2]) ### SNPs that need to be flipped
	write.table(snpstoflip,sprintf("%s/%s_flipsnps.txt",cohortdir,cohortname),quote=F,row.names=F,col.names=F)
		
	hcmiss <- read.table(sprintf("%s/%smiss.lmiss",cohortdir,cohortname),h=T,stringsAsFactors=F)
	rownames(hcmiss) <- hcmiss[,2]
	colnames(hcmiss) <- c('hc_snp','hc_chr','hc_NMISS','hc_NGENO','hc_FMISS')
	onekmiss <- read.table(sprintf("%s/%smiss.lmiss",onekgdir,onekg),h=T,stringsAsFactors=F)
	rownames(onekmiss) <- onekmiss[,2]
	colnames(onekmiss) <- c('onekg_snp','onekg_chr','onekg_NMISS','onekg_NGENO','onekg_FMISS')
	onekmiss <- onekmiss[rownames(hcmiss),]
	hconekmiss <- merge(hcmiss,onekmiss,by="row.names",sort=FALSE)
	hconekmiss$Row.names <- NULL
	hconekmiss$onekg_chr <- NULL
	hconekmiss$onekg_snp <- NULL
	hconekmiss <- hconekmiss[,c(1:4,6,7,5,8)]
	hconekmiss$hc_FMISS <- (1 - hconekmiss$hc_FMISS) * 100
	hconekmiss$onekg_FMISS <- (1 - hconekmiss$onekg_FMISS) * 100
	hconekmiss$FMISSdiff <- abs(hconekmiss$hc_FMISS - hconekmiss$onekg_FMISS)
	snpsmissfail <- hconekmiss[which(hconekmiss$FMISSdiff > 2),2]
	if(mafdiffstatus == 'domafdiff'){
		write.table(unique(c(snpsfrqfail,snpsmissfail,triallelicsnpstoremove)),sprintf("%s/%s_snpsremove.txt",cohortdir,cohortname),quote=F,row.names=F,col.names=F)	
	} else {
		write.table(unique(c(snpsmissfail,triallelicsnpstoremove)),sprintf("%s/%s_snpsremove.txt",cohortdir,cohortname),quote=F,row.names=F,col.names=F)	
	}
	
	rm(list=ls())
}

pcaeigensoft <- function(cohortdir,filename,softwaredir,popoperator){

	if(popoperator == 'bypop'){
		system(sprintf('awk \'{print $1,$2,$3,$4,$5,$1}\' %s/%s.fam > %s/%s.eig.fam',cohortdir,filename,cohortdir,filename))
		sink(sprintf('%s/%s_pca.parfile',cohortdir,filename))
		cat(sprintf('genotypename: %s/%s.bed\n',cohortdir,filename))
		cat(sprintf('snpname: %s/%s.bim\n',cohortdir,filename))
		cat(sprintf('indivname: %s/%s.eig.fam\n',cohortdir,filename))
	} else {
		sink(sprintf('%s/%s_pca.parfile',cohortdir,filename))
		cat(sprintf('genotypename: %s/%s.bed\n',cohortdir,filename))
		cat(sprintf('snpname: %s/%s.bim\n',cohortdir,filename))
		cat(sprintf('indivname: %s/%s.fam\n',cohortdir,filename))
	}
	cat(sprintf('evecoutname: %s/%sMDS.evec\n',cohortdir,filename))
	cat(sprintf('evaloutname: %s/%sMDS.eval\n',cohortdir,filename))
	cat('altnormstyle: NO\n')
	cat('numoutevec: 20\n')
	cat('numoutlieriter: 0\n')
	cat('familynames: NO\n')
	cat('numthreads: 8\n')
	sink()
	system(sprintf('%s -p %s/%s_pca.parfile',softwaredir,cohortdir,filename))
	system(sprintf('sed -i -r \'s/^\\s+//g\' %s/%sMDS.evec',cohortdir,filename))
	system(sprintf('sed -i -r \'s/\\s+/ /g\' %s/%sMDS.evec',cohortdir,filename))
}

ggmafplot <- function(cohortdir,freqfilename,plotdir){

	mafdata <- read.table(sprintf('%s/%s.frq',cohortdir,freqfilename),h=T,stringsAsFactors=F)

	pdf(sprintf("%s/%s_AllSnps.pdf",plotdir,freqfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	print(ggplot(mafdata,aes(x = MAF)) + xlab("MinorAlleleFrequency") + geom_histogram(colour="black",bins=50))
	dev.off()
	pdf(sprintf("%s/%s_RareSnps.pdf",plotdir,freqfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	print(ggplot(mafdata[which(mafdata$MAF < 0.01),],aes(x = MAF)) + xlab("MinorAlleleFrequency") + geom_histogram(colour="black",bins=20))
	dev.off()
	pdf(sprintf("%s/%s_CommonSnps.pdf",plotdir,freqfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	print(ggplot(mafdata[which(mafdata$MAF >= 0.01),],aes(x = MAF)) + xlab("MinorAlleleFrequency") + geom_histogram(colour="black",bins=49))
	dev.off()
}

ggmissplot <- function(cohortdir,missfilename,plotdir){

	imissdata <- read.table(sprintf('%s/%s.imiss',cohortdir,missfilename),h=T,stringsAsFactors=F)
	lmissdata <- read.table(sprintf('%s/%s.lmiss',cohortdir,missfilename),h=T,stringsAsFactors=F)

	pdf(sprintf("%s/%s_IndGenotypeRate.pdf",plotdir,missfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	print(ggplot(imissdata,aes(x = F_MISS)) + xlab("IndividualMissingness") + geom_histogram(colour="black",bins=60))
	dev.off()
	pdf(sprintf("%s/%s_LocusGenotypeRate.pdf",plotdir,missfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	print(ggplot(lmissdata,aes(x = F_MISS)) + xlab("SNPMissingness") + geom_histogram(colour="black",bins=100))
	dev.off()
	pdf(sprintf("%s/%s_LocusGenotypeRate_Zoomed.pdf",plotdir,missfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	print(ggplot(lmissdata,aes(x = F_MISS)) + xlab("SNPMissingness") + geom_histogram(colour="black",bins=200) + coord_cartesian(xlim = c(0,0.05)))
	dev.off()
}

ggmafvsmissplot <- function(cohortdir,freqfilename,missfilename,plotdir){

	mafdata <- read.table(sprintf('%s/%s.frq',cohortdir,freqfilename),h=T,stringsAsFactors=F)
	lmissdata <- read.table(sprintf('%s/%s.lmiss',cohortdir,missfilename),h=T,stringsAsFactors=F)

	mafdata$F_MISS <- lmissdata$F_MISS
	mafdata$CallRate <- 0
	mafdata[which(mafdata$F_MISS > 0.05),'CallRate'] <- " < 95%"
	mafdata[which(mafdata$F_MISS <= 0.05 & mafdata$F_MISS > 0.02),'CallRate'] <- "95% - 98%"
	mafdata[which(mafdata$F_MISS <= 0.02),'CallRate'] <- " > 98%"

	pdf(sprintf("%s/%s_MAFvsLocusMissingness.pdf",plotdir,freqfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	print(ggplot(mafdata, aes(x=MAF)) +  xlab("MinorAlleleFrequency") + geom_histogram(aes(fill=CallRate)))
	dev.off()

}

# ggplot of PCA for Population status Selected
ggpcapop2d <- function(cohortdir,filename,cohortname,onkgcohorttype,pcasoftware,plotfilename,plotdir){

	mycolors <- c("red","blue","green","purple","orange","black","cyan","brown","magenta","maroon","orchid","sienna",
				"orchid4","maroon4","magenta4","brown4","cyan4","orange4","purple4","green4","blue4","red4",
				"skyblue","darkturquoise", "dodgerblue2", "#FB9A99", # lt pink
                "palegreen2", "deeppink1","blue1","steelblue4","green4","green1",
				"darkorange4","skyblue2","gray70","gold1","khaki","pink","azure") # More colours "yellow4","yellow3","#CAB2D6" --# lt purple ,"#FDBF6F" --# lt orange
	mymarkers <- c(0:2,5:14,0:2,5:14)
	mycolors2 <- c("red","blue","green","purple","orange","black","cyan","brown","magenta","maroon","orchid","sienna",
				"red4","blue4","green4","purple4","orange4","cyan4","brown4","magenta4","maroon4","orchid4",
				"skyblue","darkturquoise", "dodgerblue2", "#FB9A99", # lt pink
                "palegreen2", "deeppink1","blue1","steelblue4","green4","green1",
				"darkorange4","skyblue2","gray70","gold1","khaki","pink","azure") # More colours "yellow4","yellow3","#CAB2D6" --# lt purple ,"#FDBF6F" --# lt orange

	switch(onkgcohorttype,

		selected = {
						if(pcasoftware %in% c('EIGENSOFT','eigensoft','EIGENSTRAT','eigenstrat','smartpca')){
							data <- read.table(sprintf("%s/%s1kgSLDMDS.evec",cohortdir,filename),h=F,stringsAsFactors=F)
							data <- data[,c(dim(data)[2],1:(dim(data)[2]-1))]
						}
						if(pcasoftware %in% c('plinkpca')){
							data <- read.table(sprintf("%s/%s1kgSLDMDS.eigenvec",cohortdir,filename),h=F,stringsAsFactors=F)
						}
						colnames(data)[1] <- 'FID'
						colnames(data)[2] <- 'IID'
						colnames(data)[3:22] <- c('C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')
						data$SOL <- 0
						data <- data[,c(1:2,dim(data)[2],3:22)]
						data[!(data$FID %in% c('CEU','CHB','FIN','GBR','GIH','IBS','TSI','YRI')),1] <- paste(cohortname,data[!(data$FID %in% c('CEU','CHB','FIN','GBR','GIH','IBS','TSI','YRI')),1],sep='-')
						data[which(data$FID == "CEU"),1] <- "1000g-CentralEurope"
						data[which(data$FID == "CHB"),1] <- "1000g-ChinaBeijing"
						data[which(data$FID == "FIN"),1] <- "1000g-Finland"
						data[which(data$FID == "GBR"),1] <- "1000g-GreatBritan"
						data[which(data$FID == "GIH"),1] <- "1000g-IndiaGujarat"
						data[which(data$FID == "IBS"),1] <- "1000g-SpainIberia"
						data[which(data$FID == "TSI"),1] <- "1000g-ItalyTuscan"
						data[which(data$FID == "YRI"),1] <- "1000g-NigeriaYoruba"
						mymarkersub <- rep(0,length(levels(factor(data$FID))))
						mycoloursub <- rep(0,length(levels(factor(data$FID))))
						for(i in 1:length(levels(factor(data$FID)))){
							mymarkersub[i] <- mymarkers[i %% (length(levels(factor(data$FID)))+1)]
							names(mymarkersub)[i] <- levels(factor(data$FID))[i]
							mycoloursub[i] <- mycolors[i %% (length(levels(factor(data$FID)))+1)]
							names(mycoloursub)[i] <- levels(factor(data$FID))[i]
						}
						data$FID <- as.factor(data$FID)
						colnames(data)[1] <- 'Population'
						require(RColorBrewer)
						pdf(sprintf("%s/%s_1kgS_PCAPop2D.pdf",plotdir,plotfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
						for(i in 1:10){
							for(j in i:10){
								if(i != j){
									print(ggplot(data,aes(data[,(i+3)],data[,(j+3)])) + geom_point(aes(shape = Population, colour=Population)) + scale_shape_manual(name="Cohort-Population",values=mymarkersub ) +
                             		scale_color_manual(name="Cohort-Population",values=mycoloursub)+ theme(legend.text=element_text(size=5)) + xlab(sprintf("PC%s",i)) + ylab(sprintf("PC%s",j))   )
								}
							}
						}
						dev.off()
						pdf(sprintf("%s/%s_1kgS_PCAPopIndPcs.pdf",plotdir,plotfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
						for(i in 1:20){
							print(ggplot(data,aes(c(1:nrow(data)),data[,i+3])) + geom_point(aes(shape = Population, colour=Population)) + scale_shape_manual(name="Cohort-Population",values=mymarkersub ) +
                     		scale_color_manual(name="Cohort-Population",values=mycoloursub) + theme(legend.text=element_text(size=5)) + xlab("Samples") + ylab(sprintf("PC%s",i)))
						}
						dev.off()

				   }

		)
}

# ggplot of PCA for Case-Control status
ggpcacc2d <- function(cohortdir,filename,pcasoftware,plotfilename,plotdir){

	mycolors <- c("red","blue","green","purple","orange","black","cyan","brown","magenta","gray70",
				"red4","green4","blue4","purple4","skyblue","orange4",
				"maroon","darkturquoise", "dodgerblue2", "#FB9A99", # lt pink
                "palegreen2", "deeppink1","blue1","steelblue4","green4","green1",
				"darkorange4","skyblue2","gray70","gold1","khaki","pink","azure") # More colours "yellow4","yellow3","#CAB2D6" --# lt purple ,"#FDBF6F" --# lt orange
	mymarkers <- c(0,1)
	mycolors2 <- c("red","blue","green","purple","orange","black","cyan","brown","magenta","skyblue","orchid","sienna",
				"red4","blue4","green4","purple4","orange4","cyan4","brown4","magenta4","maroon4","orchid4",
				"maroon","darkturquoise", "dodgerblue2", "#FB9A99", # lt pink
                "palegreen2", "deeppink1","blue1","steelblue4","green4","green1",
				"darkorange4","skyblue2","gray70","gold1","khaki","pink","azure") # More colours "yellow4","yellow3","#CAB2D6" --# lt purple ,"#FDBF6F" --# lt orange

	
	if(pcasoftware %in% c('EIGENSOFT','eigensoft','EIGENSTRAT','eigenstrat','smartpca')){
		data <- read.table(sprintf("%s/%sMDS.evec",cohortdir,filename),h=F,stringsAsFactors=F)
		data <- data[,c(dim(data)[2],1:(dim(data)[2]-1))]
	}
	if(pcasoftware %in% c('plinkpca')){
		data <- read.table(sprintf("%s/%sMDS.eigenvec",cohortdir,filename),h=F,stringsAsFactors=F)
	}
	#data <- read.table(sprintf('%s/%sLDMDS.eigenvec',cohortdir,filename),h=F,stringsAsFactors=F)
	colnames(data)[1] <- 'FID'
	colnames(data)[2] <- 'IID'
	colnames(data)[3:22] <- c('C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')
	data$SOL <- 0
	data <- data[,c(1:2,13,3:22)]
	famdat <- read.table(sprintf('%s/%s.fam',cohortdir,filename),h=F,stringsAsFactors=F)
	rownames(famdat) <- as.vector(famdat[,2])
	data$CCStaus <- famdat[data[,2],6]
	data[which(data$CCStaus == 1),(dim(data)[2])] <- 'Control'
	data[which(data$CCStaus == 2),(dim(data)[2])] <- 'Case'
	mymarkersub <- rep(0,length(levels(factor(data$CCStaus))))
	mycoloursub <- rep(0,length(levels(factor(data$FID))))
	for(i in 1:length(levels(factor(data$CCStaus)))){
		if(levels(factor(data$CCStaus))[i] == 'Case'){
			mymarkersub[i] <- mymarkers[2]
			names(mymarkersub)[i] <- levels(factor(data$CCStaus))[i]
		}
		else{
			mymarkersub[i] <- mymarkers[1]
			names(mymarkersub)[i] <- levels(factor(data$CCStaus))[i]
		}
	}
	for(i in 1:length(levels(factor(data$FID)))){
		mycoloursub[i] <- mycolors[i %% (length(levels(factor(data$FID)))+1)]
		names(mycoloursub)[i] <- levels(factor(data$FID))[i]
	}
	data$FID <- as.factor(data$FID)
	colnames(data)[1] <- "Population"
	data$CCStaus <- as.factor(data$CCStaus)
	require(RColorBrewer)
	pdf(sprintf("%s/%s_PCACaseControl2D.pdf",plotdir,plotfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	for(i in 1:10){
		for(j in i:10){
			if(i != j){
				print(ggplot(data,aes(data[,(i+3)],data[,(j+3)])) + geom_point(aes(shape = CCStaus, colour=Population)) + scale_shape_manual(name="AffectionStatus",values=mymarkersub ) +
                                scale_color_manual(name="Population",values=mycoloursub)+ xlab(sprintf("PC%s",i)) + ylab(sprintf("PC%s",j))   )
			}
		}
	}
	dev.off()
	pdf(sprintf("%s/%s_PCACaseControlIndPcs.pdf",plotdir,plotfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	for(i in 1:20){
		print(ggplot(data,aes(c(1:nrow(data)),data[,i+3])) + geom_point(aes(shape = CCStaus, colour=Population)) + scale_shape_manual(name="AffectionStatus",values=mymarkersub ) +
                                scale_color_manual(name="Population",values=mycoloursub) + xlab("Samples") + ylab(sprintf("PC%s",i)))
	}
	dev.off()
}

ggpcaccpop2d <- function(cohortdir,filename,pcasoftware,plotfilename,plotdir){

	mycolors <- c("red","blue","green","purple","orange","black","cyan","brown","magenta","gray70",
				"red3","green4","blue4","purple4","skyblue","orange4",
				"maroon","darkturquoise", "dodgerblue2", "#FB9A99", # lt pink
                "palegreen2", "deeppink1","blue1","steelblue4","green4","green1",
				"darkorange4","skyblue2","gray70","gold1","khaki","pink","azure") # More colours "yellow4","yellow3","#CAB2D6" --# lt purple ,"#FDBF6F" --# lt orange
	mymarkers <- c(0,1)
	mycolors2 <- c("red","blue","green","purple","orange","black","cyan","brown","magenta","skyblue","orchid","sienna",
				"red4","blue4","green4","purple4","orange4","cyan4","brown4","magenta4","maroon4","orchid4",
				"maroon","darkturquoise", "dodgerblue2", "#FB9A99", # lt pink
                "palegreen2", "deeppink1","blue1","steelblue4","green4","green1",
				"darkorange4","skyblue2","gray70","gold1","khaki","pink","azure") # More colours "yellow4","yellow3","#CAB2D6" --# lt purple ,"#FDBF6F" --# lt orange

	if(pcasoftware %in% c('EIGENSOFT','eigensoft','EIGENSTRAT','eigenstrat','smartpca')){
		data <- read.table(sprintf("%s/%sMDS.evec",cohortdir,filename),h=F,stringsAsFactors=F)
		data <- data[,c(dim(data)[2],1:(dim(data)[2]-1))]
	}
	if(pcasoftware %in% c('plinkpca')){
		data <- read.table(sprintf("%s/%sMDS.eigenvec",cohortdir,filename),h=F,stringsAsFactors=F)
	}
	#data <- read.table(sprintf('%s/%sLDMDS.eigenvec',cohortdir,filename),h=F,stringsAsFactors=F)
	colnames(data)[1] <- 'FID'
	colnames(data)[2] <- 'IID'
	colnames(data)[3:22] <- c('C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20')
	data$SOL <- 0
	data <- data[,c(1:2,13,3:22)]
	famdat <- read.table(sprintf('%s/%s.fam',cohortdir,filename),h=F,stringsAsFactors=F)
	rownames(famdat) <- as.vector(famdat[,2])
	data$AffectionStatus <- famdat[data[,2],6]
	data[which(data$AffectionStatus == 1),(dim(data)[2])] <- 'Control'
	data[which(data$AffectionStatus == 2),(dim(data)[2])] <- 'Case'
	colnames(data)[1] <- "Population"
	data$AffectionStatus <- as.factor(data$AffectionStatus)
	require(RColorBrewer)
	pdf(sprintf("%s/%s_PCACCPop2D.pdf",plotdir,plotfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	for(i in 1:10){
		for(j in i:10){
			if(i != j){
				print(ggplot(data,aes(data[,(i+3)],data[,(j+3)])) + geom_point(aes(shape = AffectionStatus, colour=AffectionStatus)) + xlab(sprintf("PC%s",i)) + ylab(sprintf("PC%s",j))   )
			}
		}
	}
	dev.off()
	
	pdf(sprintf("%s/%s_PCACCPopIndPcs.pdf",plotdir,plotfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	for(i in 1:20){
		print(ggplot(data,aes(c(1:nrow(data)),data[,i+3])) + geom_point(aes(shape = AffectionStatus, colour=AffectionStatus)) + xlab("Samples") + ylab(sprintf("PC%s",i)))
	}
	dev.off()
	pdf(sprintf("%s/%s_PCACCIndPcsPerPop.pdf",plotdir,plotfilename), onefile = TRUE, bg = "white",family = "Times", width = 6, height = 6)
	for(i in 1:2){
		for( j in 1:length(levels(as.factor(data$Population))) ){
			temp <- data[which(data$Population == levels(as.factor(data$Population))[j] | data$AffectionStatus == 'Control'),]
			temp[!(temp$Population == levels(as.factor(data$Population))[j]),'Population'] <- levels(as.factor(data$Population))[j]
			print(ggplot(temp,aes(c(1:nrow(temp)),temp[,i+3])) + geom_point(aes(shape = AffectionStatus, colour=AffectionStatus)) + xlab("Samples") + 
				ylab(sprintf("PC%s",i)) + labs( title = sprintf("%s",levels(as.factor(data$Population))[j])) )
		}
	}
	dev.off()
}

gg_qqplot <- function(xs, ci=0.95, cohortname) {
    N = length(xs)
    df = data.frame(observed=-log10(sort(xs)),
                    expected=-log10(1:N / (N+1)),
                    cupper=-log10(qbeta(ci, 1:N, N - 1:N + 1)),
                    clower=-log10(qbeta(1 - ci, 1:N, N - 1:N + 1)))
    log10Pe = expression(paste("Expected -log"[10], plain(P)))
    log10Po = expression(paste("Observed -log"[10], plain(P)))
    png(sprintf("~/Plots/pcaplots/%sqqplot.png",cohortname), width = 6, height = 6, units = "in", res = 300, colortype = "true", family = "Times")
    print( ggplot(df) + geom_point(aes(expected, observed), shape=1, size=1) + geom_abline(intercept=0, slope=1, alpha=0.5) +  geom_line(aes(expected, cupper), linetype=2) + geom_line(aes(expected, clower), linetype=2) + 
    	   xlab(log10Pe) + ylab(log10Po) )
	dev.off()
}

qqmanhatonplot <- function(currentdir,filename,plotdir,plotfilename){

	library("qqman")
	library("GWASTools")
	logisticresult <- read.table(sprintf('%s/%s.assoc.logistic',currentdir,filename),h=T,stringsAsFactors=F)
	logisticresult <- logisticresult[which(logisticresult$TEST == 'ADD'),]
	logisticresult <- logisticresult[which(!is.na(logisticresult$P)), ]
	logisticresult <- logisticresult[order(logisticresult$P), ]
	z <- qnorm(logisticresult[,"P"]/2)
	png(sprintf('%s/%s_QQ.png',plotdir,plotfilename),units="in", width=8, height=8, res=150)
	qqPlot(logisticresult$P,main=sprintf("QQ plot with %s = %s",expression(lambda),round(median(z^2)/0.456,3)),cex.axis =1.2,cex.main=1.5,cex.lab=1)
	dev.off()
	res <- logisticresult[,c(2,1,3,9)]
	snpsOfInterest <- res[which(res$P < 0.00001),1]
	png(sprintf('%s/%s_Manhattan.png',plotdir,plotfilename), units="in", width=20, height=10, res=150)
	manhattan(res,suggestiveline = -log10(0.00001),main="TS and Controls",genomewideline = -log10(0.00000001),col = c("black","#666666","#CC6600"),cex=1,cex.axis=1.5,cex.lab=1,cex.main=2,highlight = snpsOfInterest,ylim=c(0,ceiling(max(-log10(logisticresult$P))))) #
	dev.off()
}

updatesnps <- function(currentdir,bimname,beadarray){
	switch(beadarray,
		oe12v1 = { 
			data <- read.table('~/Illumina_Files/humanomniexpress-12v1/HumanOmniExpress-12-v1-0-K_Updated.csv',h=T,sep=',',stringsAsFactors=F)
			rownames(data) <- as.vector(data[,'Name'])
			bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
			temp <- data[bimfile[which(bimfile$V1 == 0),2],c('Name','Chr','MapInfo')]
			temp <- na.omit(temp)
			temp <- temp[-which(temp$Chr == 0),]
			vgprobe <- read.table('~/Illumina_Files/humanomniexpress-12v1/HumanOmniExpress-12-v1-0-K-auxilliary-file.txt',h=T,stringsAsFactors=F)
			vgprobe <- vgprobe[grep('^VG',vgprobe$Name,perl=T),]
			vgprobe[which(vgprobe$RsID == '.'),2] <- vgprobe[which(vgprobe$RsID == '.'),1]
			data[vgprobe$Name,2] <- vgprobe$RsID
			write.table(vgprobe[,c(1,2)],sprintf('%s/%s_update_name.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(temp[,c(1,3)],sprintf('%s/%s_update_map.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(temp[,c(1,2)],sprintf('%s/%s_update_chr.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
		},
		oe12v1_1 = {
			data <- read.table('~/Illumina_Files/humanomniexpress-12v1-1/HumanOmniExpress-12-v1-1-C_Updated.csv',h=T,sep=',',stringsAsFactors=F)
			rownames(data) <- as.vector(data[,'Name'])
			bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
			temp <- data[bimfile[which(bimfile$V1 == 0),2],c('Name','Chr','MapInfo')]
			temp <- na.omit(temp)
			temp <- temp[-which(temp$Chr == 0),]
			vgprobe <- read.table('~/Illumina_Files/humanomniexpress-12v1-1/HumanOmniExpress-12-v1-1-C-auxilliary-file.txt',h=T,stringsAsFactors=F)
			vgprobe <- vgprobe[grep('^VG',vgprobe$Name,perl=T),]
			vgprobe[which(vgprobe$RsID == '.'),2] <- vgprobe[which(vgprobe$RsID == '.'),1]
			data[vgprobe$Name,2] <- vgprobe$RsID
			write.table(vgprobe[,c(1,2)],sprintf('%s/%s_update_name.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(temp[,c(1,3)],sprintf('%s/%s_update_map.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(temp[,c(1,2)],sprintf('%s/%s_update_chr.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
		},
		oe24v1_1 = {
			data <- read.table('~/Illumina_Files/humanomniexpress-24v1-1/HumanOmniExpress-24v1-1_A_annotated_Updated.txt',h=T,sep='\t',stringsAsFactors=F)
			rownames(data) <- as.vector(data[,'Name'])
			bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
			temp <- data[bimfile[which(bimfile$V1 == 0),2],c('Name','Chr','MapInfo')]
			temp <- na.omit(temp)
			temp <- temp[-which(temp$Chr == 0),]
			write.table(temp[,c(1,3)],sprintf('%s/%s_update_map.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(temp[,c(1,2)],sprintf('%s/%s_update_chr.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
		},
		oe24v1_2 = {
			data <- read.table('~/Illumina_Files/infiniumomniexpress-24v1-2/InfiniumOmniExpress-24v1-2_A1_Updated.csv',h=T,sep=',',stringsAsFactors=F)
			rownames(data) <- as.vector(data[,'Name'])
			bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
			temp <- data[bimfile[which(bimfile$V1 == 0),2],c('Name','Chr','MapInfo')]
			temp <- na.omit(temp)
			temp <- temp[-which(temp$Chr == 0),]
			write.table(temp[,c(1,3)],sprintf('%s/%s_update_map.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(temp[,c(1,2)],sprintf('%s/%s_update_chr.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			test <- read.table('~/Illumina_Files/infiniumomniexpress-24v1-2/InfiniumOmniExpress-24v1-2_A1_b144_rsids.txt',h=T,stringsAsFactors=F)
			test$RsID <- gsub(",(.*?)$","",test$RsID,perl=T)
			write.table(test,sprintf('%s/%s_update_name.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
		},
		o25m_8v1_2 = {
			data <- read.table('~/Illumina_Files/omniexpressexome-25m-8v1-2/HumanOmni2-5-8-v1-2-A-Gene-Annotation-File_Updated.txt',h=T,sep='\t',stringsAsFactors=F)
			rownames(data) <- as.vector(data[,'Name'])
			bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
			temp <- data[bimfile[which(bimfile$V1 == 0),2],c('Name','Chr','MapInfo')]
			temp <- na.omit(temp)
			temp <- temp[-which(temp$Chr == 0),]
			write.table(temp[,c(1,3)],sprintf('%s/%s_update_map.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(temp[,c(1,2)],sprintf('%s/%s_update_chr.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
		},
		o8 = {
			data <- read.table('~/Illumina_Files/omniexpressexome-8/InfiniumOmniExpressExome-8v1-4_A1_annotated_Updated.txt',h=T,sep='\t',stringsAsFactors=F)
			rownames(data) <- as.vector(data[,'Name'])
			bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
			temp <- data[bimfile[which(bimfile$V1 == 0),2],c('Name','Chr','MapInfo')]
			temp <- na.omit(temp)
			temp <- temp[-which(temp$Chr == 0),]
			write.table(temp[,c(1,3)],sprintf('%s/%s_update_map.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(temp[,c(1,2)],sprintf('%s/%s_update_chr.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			test <- read.table('~/Illumina_Files/omniexpressexome-8/InfiniumOmniExpressExome-8v1-4_A1_b144_rsids.txt',h=T,stringsAsFactors=F)
			test$RsID <- gsub(",(.*?)$","",test$RsID,perl=T)
			test[which(test$RsID == '.'),2] <- test[which(test$RsID == '.'),1]
			write.table(test,sprintf('%s/%s_update_name.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
		},
		hhap550 = {
			data <- read.table('~/Illumina_Files/humanhap550/HumanHap550v3_A-b37_annotation_updated.txt',h=T,sep='\t',stringsAsFactors=F)
			bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
			write.table(data[,c(1,3)],sprintf('%s/%s_update_map.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(data[,c(1,2)],sprintf('%s/%s_update_chr.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(bimfile[!(bimfile$V2 %in% data$Name),2],sprintf('%s/%s_update_remove.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
		},
		h1mduo = {
			data <- read.table('~/Illumina_Files/human1mduo/Human1M-Duov3_annotated_Updated.txt',h=T,sep='\t',stringsAsFactors=F)
			bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
			write.table(data[,c(1,3)],sprintf('%s/%s_update_map.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(data[,c(1,2)],sprintf('%s/%s_update_chr.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(bimfile[!(bimfile$V2 %in% data$Name),2],sprintf('%s/%s_update_remove.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
		},
		h670q = {
			data <- read.table('~/Illumina_Files/human670quad/Human670-QuadCustom_v1_A-b37_Updated.txt',h=T,sep='\t',stringsAsFactors=F)
			bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
			write.table(data[,c(1,3)],sprintf('%s/%s_update_map.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(data[,c(1,2)],sprintf('%s/%s_update_chr.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(bimfile[!(bimfile$V2 %in% data$Name),2],sprintf('%s/%s_update_remove.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
		},
		affysnp6 = {
			data <- read.table('~/Affy_Files/humansnp6/AffySNP6_annotation.txt',h=T,stringsAsFactors=F)
			bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
			write.table(data[,c(1,4)],sprintf('%s/%s_update_map.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(data[,c(1,3)],sprintf('%s/%s_update_chr.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(data[,c(1,2)],sprintf('%s/%s_update_name.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(bimfile[!(bimfile$V2 %in% data$ProbeSetID),2],sprintf('%s/%s_update_remove.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
		},
		affy500 = {
			data <- read.table('~/Affy_Files/human500k/Affy500K_annotation_updated.txt',h=T,sep='\t',stringsAsFactors=F)
			bimfile <- read.table(sprintf('%s/%s.bim',currentdir,bimname),h=F,stringsAsFactors=F)
			write.table(data[,c(1,4)],sprintf('%s/%s_update_map.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(data[,c(1,3)],sprintf('%s/%s_update_chr.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(data[,c(1,2)],sprintf('%s/%s_update_name.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(bimfile[!(bimfile$V2 %in% data$ProbeSetID),2],sprintf('%s/%s_update_remove.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
		},
		nochip = {
			temp <- c('','')
			write.table(temp,sprintf('%s/%s_update_map.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
			write.table(temp,sprintf('%s/%s_update_chr.txt',currentdir,bimname),quote=F,row.names=F,col.names=F)
		}
	)
}


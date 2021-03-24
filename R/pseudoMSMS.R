#####################################################################################################################################################
### Function pseudoMSMS for raw LC-MSE chromatogram using EIC profile correlation                                                                
### Extracts and save the pseudo-MS/MS spectrum of a list of LC-MS features from LC-MSE data                                                      
### Goncalo Graca, g.gomes-da-graca@imperial.ac.uk                                                                                                                                                                                                                         
### Usage: pseudoMSMS(targetTable,filetype='CDF',method='NPC_Lipids',ESImode='POS',save.hce=FALSE,cthres=0.7)                                     
### targetTable is a .csv file containing 3 columns: fmz, frt, Sample.name (chromatograms file name)                                
###                                                                                                                                               
### Returns: pseudo-MS/MS and parent-frament EICs as well as pseudo-MS/MS peaks in .cvs format                                                    
###                                                                                                                                               
###	Notes: - correlation threshold can be changed using the parameter 'chtresh'                                                                   
### 	   - parameters 'method' and 'ESImode' are for results folders naming purposes only                                                       
###		   - filetypes accepted: 'CDF' and 'mzML'                                                                                                 
###        - MSE mzML files must be separate into two files ending in 01.mzML and 02.mzML, for no- and high colision-energy, respectively        
#####################################################################################################################################################

pseudoMSMS <- function(targetTable,filetype='mzML',method='NPC_Lipids',ESImode='POS',save.hce=FALSE,cthres=0.8){
# read XCMS options
xcmsOptions <- read.csv('XCMS_options.csv')
# Read targets table
targets <- read.csv(targetTable,header=TRUE)

# evaluate if one or more samples are to be read-----------
if(length(unique(targets[,3]))==1){
	message("Reading data...")
	xcmsF1 <- xcmsRaw(paste('./',targets[1,3],'01','.',filetype,sep=''),profstep=0.1)
	xcmsF2 <- xcmsRaw(paste('./',targets[1,3],'02','.',filetype,sep=''),profstep=0.1)
	peaksF1 <- findPeaks.centWave(xcmsF1,ppm=xcmsOptions[1,2],noise=xcmsOptions[2,2],snthresh=xcmsOptions[3,2],peakwidth=as.numeric(xcmsOptions[4,2:3]),prefilter=as.numeric(xcmsOptions[5,2:3]),mzCenterFun="wMean",integrate=2,verbose.columns=TRUE)
	peaksF2 <- findPeaks.centWave(xcmsF2,ppm=xcmsOptions[1,2],noise=xcmsOptions[2,2],snthresh=xcmsOptions[3,2],peakwidth=as.numeric(xcmsOptions[4,2:3]),prefilter=as.numeric(xcmsOptions[5,2:3]),mzCenterFun="wMean",integrate=2,verbose.columns=TRUE)
} else NULL

# Create directory to store the results
mainDir <- './pseudoMSMSspectra'
Date <- Sys.Date()
Time <- format(Sys.time(), "%X")
Time <- gsub(':','_',Time)
subDir <- paste(method,'_',ESImode,'_','Batch','_',Date,'_',Time,sep='')
dir.create(file.path(mainDir),showWarnings = FALSE)
dir.create(file.path(mainDir,subDir),showWarnings = FALSE)


for (i in 1:dim(targets)[1]){
	message(paste('####### Processing target',i,'of',dim(targets)[1]),' ########')	
	if(length(unique(targets[,3]))>1){
		message("Reading data...")
		xcmsF1 <- xcmsRaw(paste('./',targets[i,3],'01','.',filetype,sep=''),profstep=0.1)
		xcmsF2 <- xcmsRaw(paste('./',targets[i,3],'02','.',filetype,sep=''),profstep=0.1)
		peaksF1 <- findPeaks.centWave(xcmsF1,ppm=xcmsOptions[1,2],noise=xcmsOptions[2,2],snthresh=xcmsOptions[3,2],peakwidth=as.numeric(xcmsOptions[4,2:3]),prefilter=as.numeric(xcmsOptions[5,2:3]),mzCenterFun="wMean",integrate=2,verbose.columns=TRUE)
		peaksF2 <- findPeaks.centWave(xcmsF2,ppm=xcmsOptions[1,2],noise=xcmsOptions[2,2],snthresh=xcmsOptions[3,2],peakwidth=as.numeric(xcmsOptions[4,2:3]),prefilter=as.numeric(xcmsOptions[5,2:3]),mzCenterFun="wMean",integrate=2,verbose.columns=TRUE)
	} else NULL
	
	fmz <- targets[i,1]
	frt <- targets[i,2]
	SpName <- targets[i,3]
	
	# get MS spectra at feature RT --------------------
	#message("Reading high colision energy and pseudo-MS/MS spectra...")
	# get pseudo MS/MS spectrum and high colision energy spectra
	specs <- corrSpec(fmz,frt,xcmsF1,xcmsF2,peaksF1,peaksF2,cthres=cthres,inSource=FALSE) 
	highCESpec <- specs[[1]]
	pseudoSpec <- specs[[2]]

	# save pseudo-MS/MS spectrum per feature
	if(length(pseudoSpec)>1){ 
		if (length(pseudoSpec)<44) {
			pseudoSpec <- t(as.data.frame(pseudoSpec[c('mz','maxo')])) # if pseudoSpec is a vector transform it to data frame
		} else {
			pseudoSpec <- pseudoSpec[,c('mz','maxo')]
			pseudoSpec <- pseudoSpec[order(pseudoSpec[,'mz']),]
		}
		parent <- specs[[3]]
		pdf(file = paste(mainDir,"/",subDir,"/","pseudoMSMS_",SpName,"_",round(fmz,3),"mz_",round(frt,3),"s",".pdf",sep=""),height = 8,width = 5)
		par(mfrow=c(2,1))
		plot(pseudoSpec,type='h',xlim=c(50,max(pseudoSpec[,1])+100),ylim=c(0,max(pseudoSpec[,2])+max(pseudoSpec[,2])/1.5),xlab="m/z",ylab="intensity (a.u.)",col='black',lwd=1,main=paste('Pseudo-MS/MS Feature:',round(fmz,3),'m/z,',round(frt),'s'),cex.main=0.95,bty="L",xaxs="i",yaxs="i")
		text(pseudoSpec[,1]-10,pseudoSpec[,2],as.character(round(pseudoSpec[,1],3)),pos=4,cex=0.8,srt=45)
		eic.parent <- rawEIC(xcmsF1,mzrange=c(fmz-0.01,fmz+0.01),scanrange=c(parent['scmin']-6,parent['scmax']+6))
		eic.vals<-eic.parent$intensity
		eic.vals<-100*(eic.vals/eic.vals[which.max(eic.vals)])
		plot(eic.parent$scan,eic.vals,type='l',col=1,ylab='relative intensity (%)',xlab='scan number',lwd=2,lty=2,main='pseudo-MS/MS EICs',cex.main=0.95)
		eic2<-NULL
		cols<-rainbow(dim(pseudoSpec)[1])
		for(j in 1:dim(pseudoSpec)[1]){
		eic2<-rawEIC(xcmsF2,mzrange=c(pseudoSpec[j,1]-0.01,pseudoSpec[j,1]+0.01),scanrange=c(parent['scmin']-6,parent['scmax']+6))
		eic2.vals<-eic2$intensity
		eic2.vals<-100*(eic2.vals/eic2.vals[which.max(eic2.vals)])
		eic2.scan<-eic2$scan
		lines(eic2.scan,eic2.vals,col=cols[j],lwd=1,lty=1)
		legend('topright',c(paste(round(fmz,3),'m/z'),paste(round(pseudoSpec[,1],3), 'm/z')),cex = 0.8,pch=c('_'),col=c('black',cols),bty='n')
		}
		dev.off()
		write.csv(pseudoSpec,file = paste(mainDir,"/",subDir,"/","pseudoMSMS_",SpName,"_",round(fmz,3),"mz_",round(frt,3),"s",".csv",sep=""),row.names=FALSE)
	} else if (length(pseudoSpec)==1) message('no pseudo-MS/MS spectrum obtained') # if no pseudoSpec
	if(save.hce & length(specs[[1]])>1){
	highCESpec <- specs[[1]][,c('mz','into')]
	write.csv(highCESpec,file = paste(mainDir,"/",subDir,"/","high-collision_MS_",SpName,"_",round(fmz,3),"mz_",round(frt,3),"s",".csv",sep=""),row.names=FALSE)
	} else next
}
write.csv(targets,file = paste0(mainDir,"/",subDir,"/",'feature_list.csv'),row.names = FALSE)
message('Job done!')
}
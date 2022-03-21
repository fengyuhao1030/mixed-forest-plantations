rm(list = ls())
library(rstudioapi)
library(metafor)
currentPath <- getSourceEditorContext()$path
charLocations <- gregexpr('/',currentPath)[[1]]
currentPath <- substring(currentPath,1,charLocations[length(charLocations)]-1)
setwd(currentPath)

load('Data_Bio_PlotSize_Add.R')
delIDs <- which(is.na(dataMatrix$PlotSize))
if(length(delIDs) > 0){
  dataMatrix <- dataMatrix[-delIDs,]
}

##==== Random effects model ====##
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))
res <- rma.mv(Sel_yi,Sel_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
##==== Random effects model ====##

##==== Plot size ====##
res_1 <- rma.mv(Sel_yi,Sel_vi,mods = ~PlotSize,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
modelCompare_1 <- anova.rma(res,res_1)
##==== Plot size ====##

##==== Output ====##
ouputFileName <- paste0(currentPath,'/Results/Community_Bio_Sel.txt')
fileHandle <- file(ouputFileName,'w')
# Plot size
writeInfo_1 <- paste0(row.names(res_1$b)[1],': ',as.character(res_1$b[1]),', P: ',as.character(res_1$pval[1]))
writeInfo_2 <- paste0(row.names(res_1$b)[2],': ',as.character(res_1$b[2]),', P: ',as.character(res_1$pval[2]))
writeInfo_3 <- paste0('df: ',as.character(res_1$k.all - nrow(res_1$b)),'\n',
                      'Full AIC (1): ',as.character(modelCompare_1$fit.stats.f['AIC']),', Reduced AIC (null): ',as.character(modelCompare_1$fit.stats.r['AIC']),'\n',
                      'Likelihood-rato test: ',as.character(modelCompare_1$LRT),', P: ',as.character(modelCompare_1$pval))
writeLines('Plot size model',con = fileHandle)
writeLines(writeInfo_1,con = fileHandle)
writeLines(writeInfo_2,con = fileHandle)
writeLines(writeInfo_3,con = fileHandle)
writeLines('\n',con = fileHandle,sep = '')
#End
close(fileHandle)
##==== Output ====##
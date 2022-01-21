rm(list = ls())
library(rstudioapi)
library(metafor)
library(ggplot2)
currentPath <- getSourceEditorContext()$path
charLocations <- gregexpr('/',currentPath)[[1]]
currentPath <- substring(currentPath,1,charLocations[length(charLocations)]-1)
setwd(currentPath)

##==== Function ====##
theme_custom <- function(){
  myTheme <- theme(panel.background = element_rect(fill = 'white',color = 'black',size = 0.5),
                   panel.grid = element_blank(),
                   legend.position = 'none',
                   plot.margin = margin(6,6,6,6),
                   plot.background = element_blank(),
                   axis.ticks = element_line(size = 0.3),
                   axis.ticks.length = unit(-0.15,'lines'),
                   axis.title.y = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.y = element_text(size = 7,margin = margin(0,3,0,0),color = '#000000'),
                   axis.text.x = element_text(size = 7,margin = margin(5,0,0,0),color = '#000000',face = 'bold'))
  return(myTheme)
}
CMYKtoRGB <- function(C,M,Y,K){
  Rc <- (1-C)*(1-K)
  Gc <- (1-M)*(1-K)
  Bc <- (1-Y)*(1-K)
  R <- as.character(as.hexmode((ceiling(Rc*255))))
  G <- as.character(as.hexmode((ceiling(Gc*255))))
  B <- as.character(as.hexmode((ceiling(Bc*255))))
  if(nchar(R) == 1){
    R <- paste0('0',R)
  }
  if(nchar(G) == 1){
    G <- paste0('0',G)
  }
  if(nchar(B) == 1){
    B <- paste0('0',B)
  }
  RGBColor <- paste0('#',R,G,B)
  return(RGBColor)
}
##==== Function ====##

##==== Prepare data for bubble plots: Height ====##
load('HeightData.R')
dataMatrix_1 <- HeightData
# Net
dataMatrix_1_Net <- dataMatrix_1[,c('Net_yi','Net_vi')]
dataMatrix_1_Net$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix_1_Net$Net_vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix_1_Net$BubbleSize[selectIDs] <- 0.3
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix_1_Net$BubbleSize[selectIDs] <- 0.3 + (i-1)*0.3
}
selectIDs <- which((dataMatrix_1_Net$Net_yi <= quantile(dataMatrix_1_Net$Net_yi,0.98))&(dataMatrix_1_Net$Net_yi >= quantile(dataMatrix_1_Net$Net_yi,0.02)))
drawDataMatrix_1_Net <- dataMatrix_1_Net[selectIDs,c('Net_yi','BubbleSize')]
# Com
dataMatrix_1_Com <- dataMatrix_1[,c('Com_yi','Com_vi')]
dataMatrix_1_Com$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix_1_Com$Com_vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix_1_Com$BubbleSize[selectIDs] <- 0.3
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix_1_Com$BubbleSize[selectIDs] <- 0.3 + (i-1)*0.3
}
selectIDs <- which((dataMatrix_1_Com$Com_yi <= quantile(dataMatrix_1_Com$Com_yi,0.98))&(dataMatrix_1_Com$Com_yi >= quantile(dataMatrix_1_Com$Com_yi,0.02)))
drawDataMatrix_1_Com <- dataMatrix_1_Com[selectIDs,c('Com_yi','BubbleSize')]
# Sel
dataMatrix_1_Sel <- dataMatrix_1[,c('Sel_yi','Sel_vi')]
dataMatrix_1_Sel$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix_1_Sel$Sel_vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix_1_Sel$BubbleSize[selectIDs] <- 0.3
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix_1_Sel$BubbleSize[selectIDs] <- 0.3 + (i-1)*0.3
}
selectIDs <- which((dataMatrix_1_Sel$Sel_yi <= quantile(dataMatrix_1_Sel$Sel_yi,0.98))&(dataMatrix_1_Sel$Sel_yi >= quantile(dataMatrix_1_Sel$Sel_yi,0.02)))
drawDataMatrix_1_Sel <- dataMatrix_1_Sel[selectIDs,c('Sel_yi','BubbleSize')]
##==== Prepare data for bubble plots: Height ====##

##==== Prepare data for bubble plots: DBH ====##
load('DBHData.R')
dataMatrix_2 <- DBHData
# Net
dataMatrix_2_Net <- dataMatrix_2[,c('Net_yi','Net_vi')]
dataMatrix_2_Net$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix_2_Net$Net_vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix_2_Net$BubbleSize[selectIDs] <- 0.3
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix_2_Net$BubbleSize[selectIDs] <- 0.3 + (i-1)*0.3
}
selectIDs <- which((dataMatrix_2_Net$Net_yi <= quantile(dataMatrix_2_Net$Net_yi,0.98))&(dataMatrix_2_Net$Net_yi >= quantile(dataMatrix_2_Net$Net_yi,0.02)))
drawDataMatrix_2_Net <- dataMatrix_2_Net[selectIDs,c('Net_yi','BubbleSize')]
# Com
dataMatrix_2_Com <- dataMatrix_2[,c('Com_yi','Com_vi')]
dataMatrix_2_Com$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix_2_Com$Com_vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix_2_Com$BubbleSize[selectIDs] <- 0.3
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix_2_Com$BubbleSize[selectIDs] <- 0.3 + (i-1)*0.3
}
selectIDs <- which((dataMatrix_2_Com$Com_yi <= quantile(dataMatrix_2_Com$Com_yi,0.98))&(dataMatrix_2_Com$Com_yi >= quantile(dataMatrix_2_Com$Com_yi,0.02)))
drawDataMatrix_2_Com <- dataMatrix_2_Com[selectIDs,c('Com_yi','BubbleSize')]
# Sel
dataMatrix_2_Sel <- dataMatrix_2[,c('Sel_yi','Sel_vi')]
dataMatrix_2_Sel$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix_2_Sel$Sel_vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix_2_Sel$BubbleSize[selectIDs] <- 0.3
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix_2_Sel$BubbleSize[selectIDs] <- 0.3 + (i-1)*0.3
}
selectIDs <- which((dataMatrix_2_Sel$Sel_yi <= quantile(dataMatrix_2_Sel$Sel_yi,0.98))&(dataMatrix_2_Sel$Sel_yi >= quantile(dataMatrix_2_Sel$Sel_yi,0.02)))
drawDataMatrix_2_Sel <- dataMatrix_2_Sel[selectIDs,c('Sel_yi','BubbleSize')]
##==== Prepare data for bubble plots: DBH ====##

##==== Prepare data for bubble plots: Biomass ====##
load('BioData.R')
dataMatrix_3 <- BioData
# Net
dataMatrix_3_Net <- dataMatrix_3[,c('Net_yi','Net_vi')]
dataMatrix_3_Net$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix_3_Net$Net_vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix_3_Net$BubbleSize[selectIDs] <- 0.3
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix_3_Net$BubbleSize[selectIDs] <- 0.3 + (i-1)*0.3
}
selectIDs <- which((dataMatrix_3_Net$Net_yi <= quantile(dataMatrix_3_Net$Net_yi,0.98))&(dataMatrix_3_Net$Net_yi >= quantile(dataMatrix_3_Net$Net_yi,0.02)))
drawDataMatrix_3_Net <- dataMatrix_3_Net[selectIDs,c('Net_yi','BubbleSize')]
# Com
dataMatrix_3_Com <- dataMatrix_3[,c('Com_yi','Com_vi')]
dataMatrix_3_Com$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix_3_Com$Com_vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix_3_Com$BubbleSize[selectIDs] <- 0.3
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix_3_Com$BubbleSize[selectIDs] <- 0.3 + (i-1)*0.3
}
selectIDs <- which((dataMatrix_3_Com$Com_yi <= quantile(dataMatrix_3_Com$Com_yi,0.98))&(dataMatrix_3_Com$Com_yi >= quantile(dataMatrix_3_Com$Com_yi,0.02)))
drawDataMatrix_3_Com <- dataMatrix_3_Com[selectIDs,c('Com_yi','BubbleSize')]
# Sel
dataMatrix_3_Sel <- dataMatrix_3[,c('Sel_yi','Sel_vi')]
dataMatrix_3_Sel$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix_3_Sel$Sel_vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix_3_Sel$BubbleSize[selectIDs] <- 0.3
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix_3_Sel$BubbleSize[selectIDs] <- 0.3 + (i-1)*0.3
}
selectIDs <- which((dataMatrix_3_Sel$Sel_yi <= quantile(dataMatrix_3_Sel$Sel_yi,0.98))&(dataMatrix_3_Sel$Sel_yi >= quantile(dataMatrix_3_Sel$Sel_yi,0.02)))
drawDataMatrix_3_Sel <- dataMatrix_3_Sel[selectIDs,c('Sel_yi','BubbleSize')]
##==== Prepare data for bubble plots: Biomass ====##

##==== Prepare data for bubble plots: Combine ====##
drawDataMatrix_1_Net$XLocation <- 1
colnames(drawDataMatrix_1_Net)[1] <- 'yi'
drawDataMatrix_1_Com$XLocation <- 2
colnames(drawDataMatrix_1_Com)[1] <- 'yi'
drawDataMatrix_1_Sel$XLocation <- 3
colnames(drawDataMatrix_1_Sel)[1] <- 'yi'
drawDataMatrix <- rbind(drawDataMatrix_1_Net,drawDataMatrix_1_Com,drawDataMatrix_1_Sel)
##==== Prepare data for bubble plots: Combine ====##

##==== Prepare data for error bars ====##
dataMatrix_1$Level4 <- as.factor(paste(dataMatrix_1$Lat,dataMatrix_1$Lon))
dataMatrix_1$Level3 <- as.factor(paste(dataMatrix_1$Composition,dataMatrix_1$Age,dataMatrix_1$Density))
dataMatrix_1$Level2 <- as.factor(seq(1,nrow(dataMatrix_1)))
res_1 <- rma.mv(Net_yi,Net_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix_1)
wi <- weights(res_1,type = 'rowsum')
res_1 <- data.frame(Mean = as.numeric(res_1$b),LowLimit = as.numeric(res_1$ci.lb),UpLimit = as.numeric(res_1$ci.ub),
                    XLocation = 1,Class = 1)
res_2 <- rma.mv(Com_yi,Com_vi,W = wi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix_1)
res_2 <- data.frame(Mean = as.numeric(res_2$b),LowLimit = as.numeric(res_2$ci.lb),UpLimit = as.numeric(res_2$ci.ub),
                    XLocation = 2,Class = 2)
res_3 <- rma.mv(Sel_yi,Sel_vi,W = wi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix_1)
res_3 <- data.frame(Mean = as.numeric(res_3$b),LowLimit = as.numeric(res_3$ci.lb),UpLimit = as.numeric(res_3$ci.ub),
                    XLocation = 3,Class = 3)
meanDataMatrix <- rbind(res_1,res_2,res_3)
meanDataMatrix$Class <- as.factor(meanDataMatrix$Class)
##==== Prepare data for error bars ====##

##==== Draw ====##
Fig_1b_1 <- ggplot()+
  geom_jitter(data = drawDataMatrix,mapping = aes(x = XLocation,y = yi),
              width = 0.35,shape = 20,stroke = 0,color = '#666666',size = drawDataMatrix$BubbleSize,alpha = 0.15)+
  geom_hline(yintercept = 0,linetype = 'dashed',color = '#000000',size = 0.3)+
  geom_violin(data = drawDataMatrix_1_Net,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_violin(data = drawDataMatrix_1_Com,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,color = CMYKtoRGB(0.82,0.2,0.93,0.24))+
  geom_violin(data = drawDataMatrix_1_Sel,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,color = CMYKtoRGB(0,0.38,0.96,0))+
  geom_errorbar(data = meanDataMatrix,mapping = aes(x = XLocation,ymin = LowLimit,ymax = UpLimit,color = Class),size = 0.4,width = 0.1)+
  scale_color_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.2,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_x_continuous(limits = c(0.4,3.6),breaks = 2,labels = 'Height',expand = c(0,0))+
  scale_y_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),labels = c('-0.6','-0.3','0','0.3','0.6'),expand = c(0,0))+
  theme_custom()
outputFileName <- paste0(currentPath,'/Fig_1b_1.pdf')
pdf(file = outputFileName,width = 2.2,height = 1.8,colormodel = 'cmyk')
print(Fig_1b_1)
dev.off()
##==== Draw ====##
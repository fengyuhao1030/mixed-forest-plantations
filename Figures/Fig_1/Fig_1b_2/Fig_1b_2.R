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
drawDataMatrix_2_Net$XLocation <- 1
colnames(drawDataMatrix_2_Net)[1] <- 'yi'
drawDataMatrix_2_Com$XLocation <- 2
colnames(drawDataMatrix_2_Com)[1] <- 'yi'
drawDataMatrix_2_Sel$XLocation <- 3
colnames(drawDataMatrix_2_Sel)[1] <- 'yi'
drawDataMatrix_3_Net$XLocation <- 4.8
colnames(drawDataMatrix_3_Net)[1] <- 'yi'
drawDataMatrix_3_Com$XLocation <- 5.8
colnames(drawDataMatrix_3_Com)[1] <- 'yi'
drawDataMatrix_3_Sel$XLocation <- 6.8
colnames(drawDataMatrix_3_Sel)[1] <- 'yi'
drawDataMatrix <- rbind(drawDataMatrix_2_Net,drawDataMatrix_2_Com,drawDataMatrix_2_Sel,
                        drawDataMatrix_3_Net,drawDataMatrix_3_Com,drawDataMatrix_3_Sel)
##==== Prepare data for bubble plots: Combine ====##

##==== Prepare data for error bars ====##
dataMatrix_2$Level4 <- as.factor(paste(dataMatrix_2$Lat,dataMatrix_2$Lon))
dataMatrix_2$Level3 <- as.factor(paste(dataMatrix_2$Composition,dataMatrix_2$Age,dataMatrix_2$Density))
dataMatrix_2$Level2 <- as.factor(seq(1,nrow(dataMatrix_2)))
res_4 <- rma.mv(Net_yi,Net_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix_2)
wi <- weights(res_4,type = 'rowsum')
res_4 <- data.frame(Mean = as.numeric(res_4$b),LowLimit = as.numeric(res_4$ci.lb),UpLimit = as.numeric(res_4$ci.ub),
                    XLocation = 1,Class = 1)
res_5 <- rma.mv(Com_yi,Com_vi,W = wi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix_2)
res_5 <- data.frame(Mean = as.numeric(res_5$b),LowLimit = as.numeric(res_5$ci.lb),UpLimit = as.numeric(res_5$ci.ub),
                    XLocation = 2,Class = 2)
res_6 <- rma.mv(Sel_yi,Sel_vi,W = wi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix_2)
res_6 <- data.frame(Mean = as.numeric(res_6$b),LowLimit = as.numeric(res_6$ci.lb),UpLimit = as.numeric(res_6$ci.ub),
                    XLocation = 3,Class = 3)
dataMatrix_3$Level4 <- as.factor(paste(dataMatrix_3$Lat,dataMatrix_3$Lon))
dataMatrix_3$Level3 <- as.factor(paste(dataMatrix_3$Composition,dataMatrix_3$Age,dataMatrix_3$Density))
dataMatrix_3$Level2 <- as.factor(seq(1,nrow(dataMatrix_3)))
res_7 <- rma.mv(Net_yi,Net_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix_3)
wi <- weights(res_7,type = 'rowsum')
res_7 <- data.frame(Mean = as.numeric(res_7$b),LowLimit = as.numeric(res_7$ci.lb),UpLimit = as.numeric(res_7$ci.ub),
                    XLocation = 4.8,Class = 1)
res_8 <- rma.mv(Com_yi,Com_vi,W = wi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix_3)
res_8 <- data.frame(Mean = as.numeric(res_8$b),LowLimit = as.numeric(res_8$ci.lb),UpLimit = as.numeric(res_8$ci.ub),
                    XLocation = 5.8,Class = 2)
res_9 <- rma.mv(Sel_yi,Sel_vi,W = wi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix_3)
res_9 <- data.frame(Mean = as.numeric(res_9$b),LowLimit = as.numeric(res_9$ci.lb),UpLimit = as.numeric(res_9$ci.ub),
                    XLocation = 6.8,Class = 3)
meanDataMatrix <- rbind(res_4,res_5,res_6,res_7,res_8,res_9)
meanDataMatrix$Class <- as.factor(meanDataMatrix$Class)
##==== Prepare data for error bars ====##

##==== Draw ====##
Fig_1b_2 <- ggplot()+
  geom_jitter(data = drawDataMatrix,mapping = aes(x = XLocation,y = yi),
              width = 0.35,shape = 20,stroke = 0,color = '#666666',size = drawDataMatrix$BubbleSize,alpha = 0.15)+
  geom_hline(yintercept = 0,linetype = 'dashed',color = '#000000',size = 0.3)+
  geom_violin(data = drawDataMatrix_2_Net,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,fill = CMYKtoRGB(0.84,0.67,0,0),color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_violin(data = drawDataMatrix_2_Com,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),color = CMYKtoRGB(0.82,0.2,0.93,0.24))+
  geom_violin(data = drawDataMatrix_2_Sel,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,fill = CMYKtoRGB(0,0.38,0.96,0),color = CMYKtoRGB(0,0.38,0.96,0))+
  geom_violin(data = drawDataMatrix_3_Net,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,fill = CMYKtoRGB(0.84,0.67,0,0),color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_violin(data = drawDataMatrix_3_Com,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),color = CMYKtoRGB(0.82,0.2,0.93,0.24))+
  geom_violin(data = drawDataMatrix_3_Sel,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,fill = CMYKtoRGB(0,0.38,0.96,0),color = CMYKtoRGB(0,0.38,0.96,0))+
  geom_errorbar(data = meanDataMatrix,mapping = aes(x = XLocation,ymin = LowLimit,ymax = UpLimit,color = Class),size = 0.4,width = 0.1)+
  scale_color_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.2,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_x_continuous(limits = c(0.4,7.4),breaks = c(2,5.8),labels = c('DBH','Biomass'),expand = c(0,0))+
  scale_y_continuous(limits = c(-1,1.2),breaks = seq(-1,1,0.5),labels = c('-1.0','-0.5','0','0.5','1.0'),expand = c(0,0))+
  theme_custom()
outputFileName <- paste0(currentPath,'/Fig_1b_2.pdf')
pdf(file = outputFileName,width = 4.35,height = 1.8,colormodel = 'cmyk')
print(Fig_1b_2)
dev.off()
##==== Draw ====##
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

##==== Prepare data for bubble plots ====##
# Height
load('Data_Height.R')
dataMatrix_1 <- dataMatrix
dataMatrix_1$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix_1$vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix_1$BubbleSize[selectIDs] <- 0.3
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix_1$BubbleSize[selectIDs] <- 0.3 + (i-1)*0.3
}
selectIDs <- which((dataMatrix_1$yi <= quantile(dataMatrix_1$yi,0.98))&(dataMatrix_1$yi >= quantile(dataMatrix_1$yi,0.02)))
drawDataMatrix_1 <- dataMatrix_1[selectIDs,c('yi','BubbleSize')]
# DBH
load('Data_DBH.R')
dataMatrix_2 <- dataMatrix
dataMatrix_2$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix_2$vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix_2$BubbleSize[selectIDs] <- 0.3
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix_2$BubbleSize[selectIDs] <- 0.3 + (i-1)*0.3
}
selectIDs <- which((dataMatrix_2$yi <= quantile(dataMatrix_2$yi,0.98))&(dataMatrix_2$yi >= quantile(dataMatrix_2$yi,0.02)))
drawDataMatrix_2 <- dataMatrix_2[selectIDs,c('yi','BubbleSize')]
# Biomass
load('Data_Bio.R')
dataMatrix_3 <- dataMatrix
dataMatrix_3$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix_3$vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix_3$BubbleSize[selectIDs] <- 0.3
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix_3$BubbleSize[selectIDs] <- 0.3 + (i-1)*0.3
}
selectIDs <- which((dataMatrix_3$yi <= quantile(dataMatrix_3$yi,0.98))&(dataMatrix_3$yi >= quantile(dataMatrix_3$yi,0.02)))
drawDataMatrix_3 <- dataMatrix_3[selectIDs,c('yi','BubbleSize')]
# Combine
drawDataMatrix_1$XLocation <- 1
drawDataMatrix_2$XLocation <- 2
drawDataMatrix_3$XLocation <- 3
drawDataMatrix <- rbind(drawDataMatrix_1,drawDataMatrix_2,drawDataMatrix_3)
##==== Prepare data for bubble plots ====##

##==== Prepare data for error bars ====##
dataMatrix_1$Level4 <- as.factor(paste(dataMatrix_1$Lat,dataMatrix_1$Lon))
dataMatrix_1$Level3 <- as.factor(paste(dataMatrix_1$Composition,dataMatrix_1$Age,dataMatrix_1$Density))
dataMatrix_1$Level2 <- as.factor(dataMatrix_1$Species)
dataMatrix_1$Level1 <- as.factor(seq(1,nrow(dataMatrix_1)))
res_1 <- rma.mv(yi,vi,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix_1)
res_1 <- data.frame(Mean = as.numeric(res_1$b),LowLimit = as.numeric(res_1$ci.lb),UpLimit = as.numeric(res_1$ci.ub),XLocation = 1)
dataMatrix_2$Level4 <- as.factor(paste(dataMatrix_2$Lat,dataMatrix_2$Lon))
dataMatrix_2$Level3 <- as.factor(paste(dataMatrix_2$Composition,dataMatrix_2$Age,dataMatrix_2$Density))
dataMatrix_2$Level2 <- as.factor(dataMatrix_2$Species)
dataMatrix_2$Level1 <- as.factor(seq(1,nrow(dataMatrix_2)))
res_2 <- rma.mv(yi,vi,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix_2)
res_2 <- data.frame(Mean = as.numeric(res_2$b),LowLimit = as.numeric(res_2$ci.lb),UpLimit = as.numeric(res_2$ci.ub),XLocation = 2)
dataMatrix_3$Level4 <- as.factor(paste(dataMatrix_3$Lat,dataMatrix_3$Lon))
dataMatrix_3$Level3 <- as.factor(paste(dataMatrix_3$Composition,dataMatrix_3$Age,dataMatrix_3$Density))
dataMatrix_3$Level2 <- as.factor(dataMatrix_3$Species)
dataMatrix_3$Level1 <- as.factor(seq(1,nrow(dataMatrix_3)))
res_3 <- rma.mv(yi,vi,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix_3)
res_3 <- data.frame(Mean = as.numeric(res_3$b),LowLimit = as.numeric(res_3$ci.lb),UpLimit = as.numeric(res_3$ci.ub),XLocation = 3)
meanDataMatrix <- rbind(res_1,res_2,res_3)
##==== Prepare data for error bars ====##

##==== Draw ====##
Fig_1a <- ggplot()+
  geom_jitter(data = drawDataMatrix,mapping = aes(x = XLocation,y = yi),
              width = 0.35,shape = 20,stroke = 0,color = '#666666',size = drawDataMatrix$BubbleSize,alpha = 0.15)+
  geom_hline(yintercept = 0,linetype = 'dashed',color = '#000000',size = 0.3)+
  geom_violin(data = drawDataMatrix_1,mapping = aes(x = XLocation,y = yi),width = 0.85,alpha = 0,size = 0.4,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_violin(data = drawDataMatrix_2,mapping = aes(x = XLocation,y = yi),width = 0.85,alpha = 0,size = 0.4,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_violin(data = drawDataMatrix_3,mapping = aes(x = XLocation,y = yi),width = 0.85,alpha = 0,size = 0.4,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_errorbar(data = meanDataMatrix,mapping = aes(x = XLocation,ymin = LowLimit,ymax = UpLimit),width = 0.1,size = 0.4,color = CMYKtoRGB(0.84,0.67,0,0))+
  scale_x_continuous(limits = c(0.4,3.6),breaks = c(1,2,3),labels = c('Height','DBH','Biomass'),expand = c(0,0))+
  scale_y_continuous(limits = c(-1.2,1.2),breaks = seq(-1.2,1.2,0.6),labels = c('-1.2','-0.6','0','0.6','1.2'),expand = c(0,0))+
  theme_custom()
outputFileName <- paste0(currentPath,'/Fig_1a.pdf')
pdf(file = outputFileName,width = 2.2,height = 1.8,colormodel = 'cmyk')
print(Fig_1a)
dev.off()
##==== Draw ====##
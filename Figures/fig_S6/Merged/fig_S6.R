rm(list = ls())
library(rstudioapi)
library(stringr)
library(metafor)
library(ggplot2)
library(cowplot)
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

##==== Species level ====##
# Part 1: Prepare data for bubble plots
# Height
load('Data_Height_Species.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
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
load('Data_DBH_Species.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
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
load('Data_Bio_Species.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
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
# Part 2: Prepare data for error bars
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
# Part 3: Draw
fig_S6a <- ggplot()+
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
##==== Species level ====##

##==== Community level ====##
# Part 1: Prepare data for bubble plots
# Part 1.1: Height
load('Data_Height_Community.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
dataMatrix_1 <- dataMatrix
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
# Part 1.2: DBH
load('Data_DBH_Community.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
dataMatrix_2 <- dataMatrix
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
# Part 1.3: Biomass
load('Data_Bio_Community.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
dataMatrix_3 <- dataMatrix
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
# Part 1.4: Combine
drawDataMatrix_1_Net$XLocation <- 1
colnames(drawDataMatrix_1_Net)[1] <- 'yi'
drawDataMatrix_1_Com$XLocation <- 2
colnames(drawDataMatrix_1_Com)[1] <- 'yi'
drawDataMatrix_1_Sel$XLocation <- 3
colnames(drawDataMatrix_1_Sel)[1] <- 'yi'
drawDataMatrix_1 <- rbind(drawDataMatrix_1_Net,drawDataMatrix_1_Com,drawDataMatrix_1_Sel)
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
drawDataMatrix_2 <- rbind(drawDataMatrix_2_Net,drawDataMatrix_2_Com,drawDataMatrix_2_Sel,
                          drawDataMatrix_3_Net,drawDataMatrix_3_Com,drawDataMatrix_3_Sel)
# Part 2: Prepare data for error bars
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
meanDataMatrix_1 <- rbind(res_1,res_2,res_3)
meanDataMatrix_1$Class <- as.factor(meanDataMatrix_1$Class)
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
meanDataMatrix_2 <- rbind(res_4,res_5,res_6,res_7,res_8,res_9)
meanDataMatrix_2$Class <- as.factor(meanDataMatrix_2$Class)
# Part 3: Draw
fig_S6b_1 <- ggplot()+
  geom_jitter(data = drawDataMatrix_1,mapping = aes(x = XLocation,y = yi),
              width = 0.35,shape = 20,stroke = 0,color = '#666666',size = drawDataMatrix_1$BubbleSize,alpha = 0.15)+
  geom_hline(yintercept = 0,linetype = 'dashed',color = '#000000',size = 0.3)+
  geom_violin(data = drawDataMatrix_1_Net,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_violin(data = drawDataMatrix_1_Com,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,color = CMYKtoRGB(0.82,0.2,0.93,0.24))+
  geom_violin(data = drawDataMatrix_1_Sel,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,color = CMYKtoRGB(0,0.38,0.96,0))+
  geom_errorbar(data = meanDataMatrix_1,mapping = aes(x = XLocation,ymin = LowLimit,ymax = UpLimit,color = Class),size = 0.4,width = 0.1)+
  scale_color_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.2,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_x_continuous(limits = c(0.4,3.6),breaks = 2,labels = 'Height',expand = c(0,0))+
  scale_y_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),labels = c('-0.6','-0.3','0','0.3','0.6'),expand = c(0,0))+
  theme_custom()
fig_S6b_2 <- ggplot()+
  geom_jitter(data = drawDataMatrix_2,mapping = aes(x = XLocation,y = yi),
              width = 0.35,shape = 20,stroke = 0,color = '#666666',size = drawDataMatrix_2$BubbleSize,alpha = 0.15)+
  geom_hline(yintercept = 0,linetype = 'dashed',color = '#000000',size = 0.3)+
  geom_violin(data = drawDataMatrix_2_Net,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,fill = CMYKtoRGB(0.84,0.67,0,0),color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_violin(data = drawDataMatrix_2_Com,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),color = CMYKtoRGB(0.82,0.2,0.93,0.24))+
  geom_violin(data = drawDataMatrix_2_Sel,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,fill = CMYKtoRGB(0,0.38,0.96,0),color = CMYKtoRGB(0,0.38,0.96,0))+
  geom_violin(data = drawDataMatrix_3_Net,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,fill = CMYKtoRGB(0.84,0.67,0,0),color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_violin(data = drawDataMatrix_3_Com,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),color = CMYKtoRGB(0.82,0.2,0.93,0.24))+
  geom_violin(data = drawDataMatrix_3_Sel,mapping = aes(x = XLocation,y = yi),width = 0.85,size = 0.4,alpha = 0,fill = CMYKtoRGB(0,0.38,0.96,0),color = CMYKtoRGB(0,0.38,0.96,0))+
  geom_errorbar(data = meanDataMatrix_2,mapping = aes(x = XLocation,ymin = LowLimit,ymax = UpLimit,color = Class),size = 0.4,width = 0.1)+
  scale_color_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.2,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_x_continuous(limits = c(0.4,7.4),breaks = c(2,5.8),labels = c('DBH','Biomass'),expand = c(0,0))+
  scale_y_continuous(limits = c(-1,1.2),breaks = seq(-1,1,0.5),labels = c('-1.0','-0.5','0','0.5','1.0'),expand = c(0,0))+
  theme_custom()
##==== Community level ====##

##==== Combine figures ====##
fig_S6 <- ggdraw()+
  draw_plot(fig_S6a,x = 0,y = 0.5,width = 0.5,height = 0.5)+
  draw_plot(fig_S6b_1,x = 0.5,y = 0.5,width = 0.5,height = 0.5)+
  draw_plot(fig_S6b_2, x = 0,y = 0,width = 1,height = 0.5)
outputFileName <- paste0(currentPath,'/fig_S6.pdf')
pdf(file = outputFileName,width = 4.35,height = 3.6,colormodel = 'cmyk')
print(fig_S6)
dev.off()
##==== Combine figures ====##
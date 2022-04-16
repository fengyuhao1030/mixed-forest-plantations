rm(list = ls())
library(rstudioapi)
library(stringr)
library(metafor)
library(multcomp)
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
                   axis.text.x = element_text(size = 7,margin = margin(4,0,0,0),color = '#000000'))
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

##==== Part 1: Height (Species) ====##
load('Data_Height_Species.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
selectIDs <- which(dataMatrix$Density <= 4)
dataMatrix <- dataMatrix[selectIDs,]
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(dataMatrix$Species)
dataMatrix$Level1 <- as.factor(seq(1,nrow(dataMatrix)))
# Part 1.1: Collate data
dataMatrix$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix$vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix$BubbleSize <- 0.4
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix$BubbleSize[selectIDs] <- 0.4 + (i-1)*0.4
}
selectIDs <- which((dataMatrix$yi <= quantile(dataMatrix$yi,0.98))&(dataMatrix$yi >= quantile(dataMatrix$yi,0.02)))
newDataMatrix <- dataMatrix[selectIDs,c('Age','Density','yi','BubbleSize')]
# Part 1.2: Fit
res <- rma.mv(yi,vi,mods = ~Density,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix)
# Fit line
densityLine <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(densityLine) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Density),max(dataMatrix$Density),length.out = 100)
for(i in seq(1,length(tempX))){
  densityLine[i,1] <- tempX[i]
  tempResult <- predict.rma(res,newmods = tempX[i])
  densityLine[i,2] <- as.numeric(tempResult$pred)
  densityLine[i,3] <- as.numeric(tempResult$ci.ub)
  densityLine[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon
densityCIPolygon <- as.data.frame(matrix(data = 0,nrow = (2*nrow(densityLine)+1),ncol = 2))
colnames(densityCIPolygon) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(densityLine))){
  densityCIPolygon$X[count] <- densityLine$X[i]
  densityCIPolygon$Y[count] <- densityLine$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(densityLine),1,-1)){
  densityCIPolygon$X[count] <- densityLine$X[i]
  densityCIPolygon$Y[count] <- densityLine$Lower[i]
  count <- count + 1
}
densityCIPolygon$X[count] <- densityLine$X[1]
densityCIPolygon$Y[count] <- densityLine$Upper[1]
# Part 1.3: Draw
fig_S9a <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Density,y = yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = densityLine,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_polygon(data = densityCIPolygon,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = c(3,3.477,4),labels = c('1000','3000','10000'))+
  scale_y_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),labels = c('-0.6','-0.3','0','0.3','0.6'),expand = c(0,0))+
  theme_custom()
##==== Part 1: Height (Species) ====##

##==== Part 2: DBH (Species) ====##
load('Data_DBH_Species.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
selectIDs <- which(dataMatrix$Density <= 4)
dataMatrix <- dataMatrix[selectIDs,]
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(dataMatrix$Species)
dataMatrix$Level1 <- as.factor(seq(1,nrow(dataMatrix)))
# Part 2.1: Collate data
dataMatrix$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix$vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix$BubbleSize <- 0.4
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix$BubbleSize[selectIDs] <- 0.4 + (i-1)*0.4
}
selectIDs <- which((dataMatrix$yi <= quantile(dataMatrix$yi,0.98))&(dataMatrix$yi >= quantile(dataMatrix$yi,0.02)))
newDataMatrix <- dataMatrix[selectIDs,c('Age','Density','yi','BubbleSize')]
# Part 2.2: Fit
res <- rma.mv(yi,vi,mods = ~Density,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix)
# Fit line
densityLine <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(densityLine) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Density),max(dataMatrix$Density),length.out = 100)
for(i in seq(1,length(tempX))){
  densityLine[i,1] <- tempX[i]
  tempResult <- predict.rma(res,newmods = tempX[i])
  densityLine[i,2] <- as.numeric(tempResult$pred)
  densityLine[i,3] <- as.numeric(tempResult$ci.ub)
  densityLine[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon
densityCIPolygon <- as.data.frame(matrix(data = 0,nrow = (2*nrow(densityLine)+1),ncol = 2))
colnames(densityCIPolygon) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(densityLine))){
  densityCIPolygon$X[count] <- densityLine$X[i]
  densityCIPolygon$Y[count] <- densityLine$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(densityLine),1,-1)){
  densityCIPolygon$X[count] <- densityLine$X[i]
  densityCIPolygon$Y[count] <- densityLine$Lower[i]
  count <- count + 1
}
densityCIPolygon$X[count] <- densityLine$X[1]
densityCIPolygon$Y[count] <- densityLine$Upper[1]
# Part 2.3: Draw
fig_S9b <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Density,y = yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = densityLine,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0),linetype = 'dashed')+
  geom_polygon(data = densityCIPolygon,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = c(3,3.477,4),labels = c('1000','3000','10000'))+
  scale_y_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),labels = c('-0.6','-0.3','0','0.3','0.6'),expand = c(0,0))+
  theme_custom()
##==== Part 2: DBH (Species) ====##

##==== Part 3: Biomass (Species) ====##
load('Data_Bio_Species.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
selectIDs <- which(dataMatrix$Density <= 4)
dataMatrix <- dataMatrix[selectIDs,]
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(dataMatrix$Species)
dataMatrix$Level1 <- as.factor(seq(1,nrow(dataMatrix)))
# Part 3.1: Collate data
dataMatrix$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix$vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix$BubbleSize <- 0.4
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix$BubbleSize[selectIDs] <- 0.4 + (i-1)*0.4
}
selectIDs <- which((dataMatrix$yi <= quantile(dataMatrix$yi,0.98))&(dataMatrix$yi >= quantile(dataMatrix$yi,0.02)))
newDataMatrix <- dataMatrix[selectIDs,c('Age','Density','yi','BubbleSize')]
# Part 3.2: Fit
res <- rma.mv(yi,vi,mods = ~Density+Density2,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix)
# Fit line
densityLine <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(densityLine) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Density),max(dataMatrix$Density),length.out = 100)
for(i in seq(1,length(tempX))){
  densityLine[i,1] <- tempX[i]
  tempResult <- predict.rma(res,newmods = c(tempX[i],tempX[i]^2))
  densityLine[i,2] <- as.numeric(tempResult$pred)
  densityLine[i,3] <- as.numeric(tempResult$ci.ub)
  densityLine[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon
densityCIPolygon <- as.data.frame(matrix(data = 0,nrow = (2*nrow(densityLine)+1),ncol = 2))
colnames(densityCIPolygon) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(densityLine))){
  densityCIPolygon$X[count] <- densityLine$X[i]
  densityCIPolygon$Y[count] <- densityLine$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(densityLine),1,-1)){
  densityCIPolygon$X[count] <- densityLine$X[i]
  densityCIPolygon$Y[count] <- densityLine$Lower[i]
  count <- count + 1
}
densityCIPolygon$X[count] <- densityLine$X[1]
densityCIPolygon$Y[count] <- densityLine$Upper[1]
# Part 3.3: Draw
fig_S9c <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Density,y = yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = densityLine,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_polygon(data = densityCIPolygon,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = c(3,3.477,4),labels = c('1000','3000','10000'))+
  scale_y_continuous(limits = c(-1.2,1.2),breaks = seq(-1.2,1.2,0.6),labels = c('-1.2','-0.6','0','0.6','1.2'),expand = c(0,0))+
  xlab('Density')+
  theme_custom()
##==== Part 3: Biomass (Species) ====##

##==== Part 4: Height (Community) ====##
load('Data_Height_Community.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
selectIDs <- which(dataMatrix$Density <= 4)
dataMatrix <- dataMatrix[selectIDs,]
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))
# Part 4.1: Collate data
dataMatrix$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix$Net_vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix$BubbleSize <- 0.4
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix$BubbleSize[selectIDs] <- 0.4 + (i-1)*0.4
}
selectIDs <- which((dataMatrix$Net_yi <= quantile(dataMatrix$Net_yi,0.98))&(dataMatrix$Net_yi >= quantile(dataMatrix$Net_yi,0.02)))
newDataMatrix <- dataMatrix[selectIDs,c('Age','Density','Net_yi','BubbleSize')]
# Part 4.2: Fit
res_Net <- rma.mv(Net_yi,Net_vi,mods = ~Density+Density2,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Com <- rma.mv(Com_yi,Com_vi,mods = ~Density+Density2,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~Density,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
# Fit line Net
densityLineNet <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(densityLineNet) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Density),max(dataMatrix$Density),length.out = 100)
for(i in seq(1,length(tempX))){
  densityLineNet[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Net,newmods = c(tempX[i],tempX[i]^2))
  densityLineNet[i,2] <- as.numeric(tempResult$pred)
  densityLineNet[i,3] <- as.numeric(tempResult$ci.ub)
  densityLineNet[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Com
densityLineCom <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(densityLineCom) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Density),max(dataMatrix$Density),length.out = 100)
for(i in seq(1,length(tempX))){
  densityLineCom[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Com,newmods = c(tempX[i],tempX[i]^2))
  densityLineCom[i,2] <- as.numeric(tempResult$pred)
  densityLineCom[i,3] <- as.numeric(tempResult$ci.ub)
  densityLineCom[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Sel
densityLineSel <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(densityLineSel) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Density),max(dataMatrix$Density),length.out = 100)
for(i in seq(1,length(tempX))){
  densityLineSel[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Sel,newmods = tempX[i])
  densityLineSel[i,2] <- as.numeric(tempResult$pred)
  densityLineSel[i,3] <- as.numeric(tempResult$ci.ub)
  densityLineSel[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon Net
densityCIPolygonNet <- as.data.frame(matrix(data = 0,nrow = (2*nrow(densityLineNet)+1),ncol = 2))
colnames(densityCIPolygonNet) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(densityLineNet))){
  densityCIPolygonNet$X[count] <- densityLineNet$X[i]
  densityCIPolygonNet$Y[count] <- densityLineNet$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(densityLineNet),1,-1)){
  densityCIPolygonNet$X[count] <- densityLineNet$X[i]
  densityCIPolygonNet$Y[count] <- densityLineNet$Lower[i]
  count <- count + 1
}
densityCIPolygonNet$X[count] <- densityLineNet$X[1]
densityCIPolygonNet$Y[count] <- densityLineNet$Upper[1]
# CI polygon Com
densityCIPolygonCom <- as.data.frame(matrix(data = 0,nrow = (2*nrow(densityLineCom)+1),ncol = 2))
colnames(densityCIPolygonCom) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(densityLineCom))){
  densityCIPolygonCom$X[count] <- densityLineCom$X[i]
  densityCIPolygonCom$Y[count] <- densityLineCom$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(densityLineCom),1,-1)){
  densityCIPolygonCom$X[count] <- densityLineCom$X[i]
  densityCIPolygonCom$Y[count] <- densityLineCom$Lower[i]
  count <- count + 1
}
densityCIPolygonCom$X[count] <- densityLineCom$X[1]
densityCIPolygonCom$Y[count] <- densityLineCom$Upper[1]
# CI polygon Sel
densityCIPolygonSel <- as.data.frame(matrix(data = 0,nrow = (2*nrow(densityLineSel)+1),ncol = 2))
colnames(densityCIPolygonSel) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(densityLineSel))){
  densityCIPolygonSel$X[count] <- densityLineSel$X[i]
  densityCIPolygonSel$Y[count] <- densityLineSel$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(densityLineSel),1,-1)){
  densityCIPolygonSel$X[count] <- densityLineSel$X[i]
  densityCIPolygonSel$Y[count] <- densityLineSel$Lower[i]
  count <- count + 1
}
densityCIPolygonSel$X[count] <- densityLineSel$X[1]
densityCIPolygonSel$Y[count] <- densityLineSel$Upper[1]
# Part 4.3: Draw
fig_S9d <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Density,y = Net_yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = densityLineSel,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0,0.38,0.96,0),linetype = 'dashed')+
  geom_polygon(data = densityCIPolygonSel,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0,0.38,0.96,0),alpha = 0.3)+
  geom_line(data = densityLineCom,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.82,0.2,0.93,0.24))+
  geom_polygon(data = densityCIPolygonCom,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),alpha = 0.3)+
  geom_line(data = densityLineNet,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_polygon(data = densityCIPolygonNet,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = c(3,3.477,4),labels = c('1000','3000','10000'))+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = seq(-0.4,0.4,0.2),labels = c('-0.4','-0.2','0','0.2','0.4'),expand = c(0,0))+
  theme_custom()
##==== Part 4: Height (Community) ====##

##==== Part 5: DBH (Community) ====##
load('Data_DBH_Community.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
selectIDs <- which(dataMatrix$Density <= 4)
dataMatrix <- dataMatrix[selectIDs,]
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))
# Part 5.1: Collate data
dataMatrix$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix$Net_vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix$BubbleSize <- 0.4
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix$BubbleSize[selectIDs] <- 0.4 + (i-1)*0.4
}
selectIDs <- which((dataMatrix$Net_yi <= quantile(dataMatrix$Net_yi,0.98))&(dataMatrix$Net_yi >= quantile(dataMatrix$Net_yi,0.02)))
newDataMatrix <- dataMatrix[selectIDs,c('Age','Density','Net_yi','BubbleSize')]
# Part 5.2: Fit
res_Net <- rma.mv(Net_yi,Net_vi,mods = ~Density,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Com <- rma.mv(Com_yi,Com_vi,mods = ~Density+Density2,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~Density,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
# Fit line Net
densityLineNet <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(densityLineNet) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Density),max(dataMatrix$Density),length.out = 100)
for(i in seq(1,length(tempX))){
  densityLineNet[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Net,newmods = tempX[i])
  densityLineNet[i,2] <- as.numeric(tempResult$pred)
  densityLineNet[i,3] <- as.numeric(tempResult$ci.ub)
  densityLineNet[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Com
densityLineCom <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(densityLineCom) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Density),max(dataMatrix$Density),length.out = 100)
for(i in seq(1,length(tempX))){
  densityLineCom[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Com,newmods = c(tempX[i],tempX[i]^2))
  densityLineCom[i,2] <- as.numeric(tempResult$pred)
  densityLineCom[i,3] <- as.numeric(tempResult$ci.ub)
  densityLineCom[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Sel
densityLineSel <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(densityLineSel) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Density),max(dataMatrix$Density),length.out = 100)
for(i in seq(1,length(tempX))){
  densityLineSel[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Sel,newmods = tempX[i])
  densityLineSel[i,2] <- as.numeric(tempResult$pred)
  densityLineSel[i,3] <- as.numeric(tempResult$ci.ub)
  densityLineSel[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon Net
densityCIPolygonNet <- as.data.frame(matrix(data = 0,nrow = (2*nrow(densityLineNet)+1),ncol = 2))
colnames(densityCIPolygonNet) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(densityLineNet))){
  densityCIPolygonNet$X[count] <- densityLineNet$X[i]
  densityCIPolygonNet$Y[count] <- densityLineNet$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(densityLineNet),1,-1)){
  densityCIPolygonNet$X[count] <- densityLineNet$X[i]
  densityCIPolygonNet$Y[count] <- densityLineNet$Lower[i]
  count <- count + 1
}
densityCIPolygonNet$X[count] <- densityLineNet$X[1]
densityCIPolygonNet$Y[count] <- densityLineNet$Upper[1]
# CI polygon Com
densityCIPolygonCom <- as.data.frame(matrix(data = 0,nrow = (2*nrow(densityLineCom)+1),ncol = 2))
colnames(densityCIPolygonCom) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(densityLineCom))){
  densityCIPolygonCom$X[count] <- densityLineCom$X[i]
  densityCIPolygonCom$Y[count] <- densityLineCom$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(densityLineCom),1,-1)){
  densityCIPolygonCom$X[count] <- densityLineCom$X[i]
  densityCIPolygonCom$Y[count] <- densityLineCom$Lower[i]
  count <- count + 1
}
densityCIPolygonCom$X[count] <- densityLineCom$X[1]
densityCIPolygonCom$Y[count] <- densityLineCom$Upper[1]
# CI polygon Sel
densityCIPolygonSel <- as.data.frame(matrix(data = 0,nrow = (2*nrow(densityLineSel)+1),ncol = 2))
colnames(densityCIPolygonSel) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(densityLineSel))){
  densityCIPolygonSel$X[count] <- densityLineSel$X[i]
  densityCIPolygonSel$Y[count] <- densityLineSel$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(densityLineSel),1,-1)){
  densityCIPolygonSel$X[count] <- densityLineSel$X[i]
  densityCIPolygonSel$Y[count] <- densityLineSel$Lower[i]
  count <- count + 1
}
densityCIPolygonSel$X[count] <- densityLineSel$X[1]
densityCIPolygonSel$Y[count] <- densityLineSel$Upper[1]
# Part 5.3: Draw
fig_S9e <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Density,y = Net_yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = densityLineSel,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0,0.38,0.96,0))+
  geom_polygon(data = densityCIPolygonSel,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0,0.38,0.96,0),alpha = 0.3)+
  geom_line(data = densityLineCom,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.82,0.2,0.93,0.24),linetype = 'dashed')+
  geom_polygon(data = densityCIPolygonCom,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),alpha = 0.3)+
  geom_line(data = densityLineNet,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0),linetype = 'dashed')+
  geom_polygon(data = densityCIPolygonNet,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = c(3,3.477,4),labels = c('1000','3000','10000'))+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = seq(-0.4,0.4,0.2),labels = c('-0.4','-0.2','0','0.2','0.4'),expand = c(0,0))+
  theme_custom()
##==== Part 5: DBH (Community) ====##

##==== Part 6: Biomass (Community) ====##
load('Data_Bio_Community.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
selectIDs <- which(dataMatrix$Density <= 4)
dataMatrix <- dataMatrix[selectIDs,]
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))
# Part 6.1: Collate data
dataMatrix$BubbleSize <- 0
Weight <- 1/sqrt(dataMatrix$Net_vi)
Breaks <- quantile(Weight,c(0.2,0.4,0.6,0.8,1))
selectIDs <- which(Weight < Breaks[1])
dataMatrix$BubbleSize <- 0.4
for(i in seq(2,length(Breaks))){
  selectIDs <- which((Weight < Breaks[i])&(Weight >= Breaks[i-1]))
  dataMatrix$BubbleSize[selectIDs] <- 0.4 + (i-1)*0.4
}
selectIDs <- which((dataMatrix$Net_yi <= quantile(dataMatrix$Net_yi,0.98))&(dataMatrix$Net_yi >= quantile(dataMatrix$Net_yi,0.02)))
newDataMatrix <- dataMatrix[selectIDs,c('Age','Density','Net_yi','BubbleSize')]
# Part 6.2: Fit
res_Net <- rma.mv(Net_yi,Net_vi,mods = ~Density+Density2,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Com <- rma.mv(Com_yi,Com_vi,mods = ~Density+Density2,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~Density,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
# Fit line Net
densityLineNet <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(densityLineNet) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Density),max(dataMatrix$Density),length.out = 100)
for(i in seq(1,length(tempX))){
  densityLineNet[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Net,newmods = c(tempX[i],tempX[i]^2))
  densityLineNet[i,2] <- as.numeric(tempResult$pred)
  densityLineNet[i,3] <- as.numeric(tempResult$ci.ub)
  densityLineNet[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Com
densityLineCom <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(densityLineCom) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Density),max(dataMatrix$Density),length.out = 100)
for(i in seq(1,length(tempX))){
  densityLineCom[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Com,newmods = c(tempX[i],tempX[i]^2))
  densityLineCom[i,2] <- as.numeric(tempResult$pred)
  densityLineCom[i,3] <- as.numeric(tempResult$ci.ub)
  densityLineCom[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Sel
densityLineSel <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(densityLineSel) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Density),max(dataMatrix$Density),length.out = 100)
for(i in seq(1,length(tempX))){
  densityLineSel[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Sel,newmods = tempX[i])
  densityLineSel[i,2] <- as.numeric(tempResult$pred)
  densityLineSel[i,3] <- as.numeric(tempResult$ci.ub)
  densityLineSel[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon Net
densityCIPolygonNet <- as.data.frame(matrix(data = 0,nrow = (2*nrow(densityLineNet)+1),ncol = 2))
colnames(densityCIPolygonNet) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(densityLineNet))){
  densityCIPolygonNet$X[count] <- densityLineNet$X[i]
  densityCIPolygonNet$Y[count] <- densityLineNet$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(densityLineNet),1,-1)){
  densityCIPolygonNet$X[count] <- densityLineNet$X[i]
  densityCIPolygonNet$Y[count] <- densityLineNet$Lower[i]
  count <- count + 1
}
densityCIPolygonNet$X[count] <- densityLineNet$X[1]
densityCIPolygonNet$Y[count] <- densityLineNet$Upper[1]
# CI polygon Com
densityCIPolygonCom <- as.data.frame(matrix(data = 0,nrow = (2*nrow(densityLineCom)+1),ncol = 2))
colnames(densityCIPolygonCom) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(densityLineCom))){
  densityCIPolygonCom$X[count] <- densityLineCom$X[i]
  densityCIPolygonCom$Y[count] <- densityLineCom$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(densityLineCom),1,-1)){
  densityCIPolygonCom$X[count] <- densityLineCom$X[i]
  densityCIPolygonCom$Y[count] <- densityLineCom$Lower[i]
  count <- count + 1
}
densityCIPolygonCom$X[count] <- densityLineCom$X[1]
densityCIPolygonCom$Y[count] <- densityLineCom$Upper[1]
# CI polygon Sel
densityCIPolygonSel <- as.data.frame(matrix(data = 0,nrow = (2*nrow(densityLineSel)+1),ncol = 2))
colnames(densityCIPolygonSel) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(densityLineSel))){
  densityCIPolygonSel$X[count] <- densityLineSel$X[i]
  densityCIPolygonSel$Y[count] <- densityLineSel$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(densityLineSel),1,-1)){
  densityCIPolygonSel$X[count] <- densityLineSel$X[i]
  densityCIPolygonSel$Y[count] <- densityLineSel$Lower[i]
  count <- count + 1
}
densityCIPolygonSel$X[count] <- densityLineSel$X[1]
densityCIPolygonSel$Y[count] <- densityLineSel$Upper[1]
# Part 6.3: Draw
fig_S9f <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Density,y = Net_yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = densityLineSel,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0,0.38,0.96,0),linetype = 'dashed')+
  geom_polygon(data = densityCIPolygonSel,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0,0.38,0.96,0),alpha = 0.3)+
  geom_line(data = densityLineCom,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.82,0.2,0.93,0.24))+
  geom_polygon(data = densityCIPolygonCom,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),alpha = 0.3)+
  geom_line(data = densityLineNet,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_polygon(data = densityCIPolygonNet,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = c(3,3.477,4),labels = c('1000','3000','10000'))+
  scale_y_continuous(limits = c(-0.9,1.2),breaks = seq(-0.6,1.2,0.6),labels = c('-0.6','0','0.6','1.2'),expand = c(0,0))+
  theme_custom()
##==== Part 6: Biomass (Community) ====##

##==== Part 7: Combine figures ====##
fig_S9 <- ggdraw()+
  draw_plot(fig_S9a,x = 0,y = 0.5,width = 0.333,height = 0.5)+
  draw_plot(fig_S9b,x = 0.334,y = 0.5,width = 0.333,height = 0.5)+
  draw_plot(fig_S9c,x = 0.667,y = 0.5,width = 0.333,height = 0.5)+
  draw_plot(fig_S9d,x = 0,y = 0,width = 0.333,height = 0.5)+
  draw_plot(fig_S9e,x = 0.334,y = 0,width = 0.333,height = 0.5)+
  draw_plot(fig_S9f,x = 0.667,y = 0,width = 0.333,height = 0.5)
outputFileName <- paste0(currentPath,'/fig_S9.pdf')
pdf(file = outputFileName,width = 4.65,height = 2.7,colormodel = 'cmyk')
print(fig_S9)
dev.off()
##==== Part 7: Combine figures ====##
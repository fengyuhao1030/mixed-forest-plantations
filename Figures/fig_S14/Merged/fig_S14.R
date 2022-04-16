rm(list = ls())
library(rstudioapi)
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
newDataMatrix <- dataMatrix[selectIDs,c('Temperature','Precipitation','yi','BubbleSize')]
# Part 1.2: Fit
res <- rma.mv(yi,vi,mods = ~Temperature,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix)
# Fit line
temLine <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLine) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Temperature),max(dataMatrix$Temperature),length.out = 100)
for(i in seq(1,length(tempX))){
  temLine[i,1] <- tempX[i]
  tempResult <- predict.rma(res,newmods = tempX[i])
  temLine[i,2] <- as.numeric(tempResult$pred)
  temLine[i,3] <- as.numeric(tempResult$ci.ub)
  temLine[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon
temCIPolygon <- as.data.frame(matrix(data = 0,nrow = (2*nrow(temLine)+1),ncol = 2))
colnames(temCIPolygon) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(temLine))){
  temCIPolygon$X[count] <- temLine$X[i]
  temCIPolygon$Y[count] <- temLine$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(temLine),1,-1)){
  temCIPolygon$X[count] <- temLine$X[i]
  temCIPolygon$Y[count] <- temLine$Lower[i]
  count <- count + 1
}
temCIPolygon$X[count] <- temLine$X[1]
temCIPolygon$Y[count] <- temLine$Upper[1]
# Part 1.3: Draw
fig_S14a <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Temperature,y = yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = temLine,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_polygon(data = temCIPolygon,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = seq(5,25,5),labels = c('5','10','15','20','25'))+
  scale_y_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),labels = c('-0.6','-0.3','0','0.3','0.6'),expand = c(0,0))+
  theme_custom()
##==== Part 1: Height (Species) ====##

##==== Part 2: DBH (Species) ====##
load('Data_DBH_Species.R')
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
newDataMatrix <- dataMatrix[selectIDs,c('Temperature','Precipitation','yi','BubbleSize')]
# Part 2.2: Fit
res <- rma.mv(yi,vi,mods = ~Temperature,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix)
# Fit line
temLine <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLine) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Temperature),max(dataMatrix$Temperature),length.out = 100)
for(i in seq(1,length(tempX))){
  temLine[i,1] <- tempX[i]
  tempResult <- predict.rma(res,newmods = tempX[i])
  temLine[i,2] <- as.numeric(tempResult$pred)
  temLine[i,3] <- as.numeric(tempResult$ci.ub)
  temLine[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon
temCIPolygon <- as.data.frame(matrix(data = 0,nrow = (2*nrow(temLine)+1),ncol = 2))
colnames(temCIPolygon) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(temLine))){
  temCIPolygon$X[count] <- temLine$X[i]
  temCIPolygon$Y[count] <- temLine$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(temLine),1,-1)){
  temCIPolygon$X[count] <- temLine$X[i]
  temCIPolygon$Y[count] <- temLine$Lower[i]
  count <- count + 1
}
temCIPolygon$X[count] <- temLine$X[1]
temCIPolygon$Y[count] <- temLine$Upper[1]
# Part 2.3: Draw
fig_S14b <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Temperature,y = yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = temLine,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0),linetype = 'dashed')+
  geom_polygon(data = temCIPolygon,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = seq(5,25,5),labels = c('5','10','15','20','25'))+
  scale_y_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),labels = c('-0.6','-0.3','0','0.3','0.6'),expand = c(0,0))+
  theme_custom()
##==== Part 2: DBH (Species) ====##

##==== Part 3: Biomass (Species) ====##
load('Data_Bio_Species.R')
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
newDataMatrix <- dataMatrix[selectIDs,c('Temperature','Precipitation','yi','BubbleSize')]
# Part 3.2: Fit
res <- rma.mv(yi,vi,mods = ~Temperature,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix)
# Fit line
temLine <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLine) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Temperature),max(dataMatrix$Temperature),length.out = 100)
for(i in seq(1,length(tempX))){
  temLine[i,1] <- tempX[i]
  tempResult <- predict.rma(res,newmods = tempX[i])
  temLine[i,2] <- as.numeric(tempResult$pred)
  temLine[i,3] <- as.numeric(tempResult$ci.ub)
  temLine[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon
temCIPolygon <- as.data.frame(matrix(data = 0,nrow = (2*nrow(temLine)+1),ncol = 2))
colnames(temCIPolygon) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(temLine))){
  temCIPolygon$X[count] <- temLine$X[i]
  temCIPolygon$Y[count] <- temLine$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(temLine),1,-1)){
  temCIPolygon$X[count] <- temLine$X[i]
  temCIPolygon$Y[count] <- temLine$Lower[i]
  count <- count + 1
}
temCIPolygon$X[count] <- temLine$X[1]
temCIPolygon$Y[count] <- temLine$Upper[1]
# Part 3.3: Draw
fig_S14c <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Temperature,y = yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = temLine,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0),linetype = 'dashed')+
  geom_polygon(data = temCIPolygon,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = seq(5,25,5),labels = c('5','10','15','20','25'))+
  scale_y_continuous(limits = c(-1.2,1.2),breaks = seq(-1.2,1.2,0.6),labels = c('-1.2','-0.6','0','0.6','1.2'),expand = c(0,0))+
  theme_custom()
##==== Part 3: Biomass (Species) ====##

##==== Part 4: Height (Community) ====##
load('Data_Height_Community.R')
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
newDataMatrix <- dataMatrix[selectIDs,c('Temperature','Precipitation','Net_yi','BubbleSize')]
# Part 4.2: Fit
res_Net <- rma.mv(Net_yi,Net_vi,mods = ~Temperature,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Com <- rma.mv(Com_yi,Com_vi,mods = ~Temperature,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~Temperature,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
# Fit line Net
temLineNet <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLineNet) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Temperature),max(dataMatrix$Temperature),length.out = 100)
for(i in seq(1,length(tempX))){
  temLineNet[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Net,newmods = tempX[i])
  temLineNet[i,2] <- as.numeric(tempResult$pred)
  temLineNet[i,3] <- as.numeric(tempResult$ci.ub)
  temLineNet[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Com
temLineCom <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLineCom) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Temperature),max(dataMatrix$Temperature),length.out = 100)
for(i in seq(1,length(tempX))){
  temLineCom[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Com,newmods = tempX[i])
  temLineCom[i,2] <- as.numeric(tempResult$pred)
  temLineCom[i,3] <- as.numeric(tempResult$ci.ub)
  temLineCom[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Sel
temLineSel <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLineSel) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Temperature),max(dataMatrix$Temperature),length.out = 100)
for(i in seq(1,length(tempX))){
  temLineSel[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Sel,newmods = tempX[i])
  temLineSel[i,2] <- as.numeric(tempResult$pred)
  temLineSel[i,3] <- as.numeric(tempResult$ci.ub)
  temLineSel[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon Net
temCIPolygonNet <- as.data.frame(matrix(data = 0,nrow = (2*nrow(temLineNet)+1),ncol = 2))
colnames(temCIPolygonNet) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(temLineNet))){
  temCIPolygonNet$X[count] <- temLineNet$X[i]
  temCIPolygonNet$Y[count] <- temLineNet$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(temLineNet),1,-1)){
  temCIPolygonNet$X[count] <- temLineNet$X[i]
  temCIPolygonNet$Y[count] <- temLineNet$Lower[i]
  count <- count + 1
}
temCIPolygonNet$X[count] <- temLineNet$X[1]
temCIPolygonNet$Y[count] <- temLineNet$Upper[1]
# CI polygon Com
temCIPolygonCom <- as.data.frame(matrix(data = 0,nrow = (2*nrow(temLineCom)+1),ncol = 2))
colnames(temCIPolygonCom) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(temLineCom))){
  temCIPolygonCom$X[count] <- temLineCom$X[i]
  temCIPolygonCom$Y[count] <- temLineCom$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(temLineCom),1,-1)){
  temCIPolygonCom$X[count] <- temLineCom$X[i]
  temCIPolygonCom$Y[count] <- temLineCom$Lower[i]
  count <- count + 1
}
temCIPolygonCom$X[count] <- temLineCom$X[1]
temCIPolygonCom$Y[count] <- temLineCom$Upper[1]
# CI polygon Sel
temCIPolygonSel <- as.data.frame(matrix(data = 0,nrow = (2*nrow(temLineSel)+1),ncol = 2))
colnames(temCIPolygonSel) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(temLineSel))){
  temCIPolygonSel$X[count] <- temLineSel$X[i]
  temCIPolygonSel$Y[count] <- temLineSel$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(temLineSel),1,-1)){
  temCIPolygonSel$X[count] <- temLineSel$X[i]
  temCIPolygonSel$Y[count] <- temLineSel$Lower[i]
  count <- count + 1
}
temCIPolygonSel$X[count] <- temLineSel$X[1]
temCIPolygonSel$Y[count] <- temLineSel$Upper[1]
# Part 4.3: Draw
fig_S14d <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Temperature,y = Net_yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = temLineSel,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0,0.38,0.96,0),linetype = 'dashed')+
  geom_polygon(data = temCIPolygonSel,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0,0.38,0.96,0),alpha = 0.3)+
  geom_line(data = temLineCom,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.82,0.2,0.93,0.24))+
  geom_polygon(data = temCIPolygonCom,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),alpha = 0.3)+
  geom_line(data = temLineNet,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_polygon(data = temCIPolygonNet,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = seq(5,25,5),labels = c('5','10','15','20','25'))+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = seq(-0.4,0.4,0.2),labels = c('-0.4','-0.2','0','0.2','0.4'),expand = c(0,0))+
  theme_custom()
##==== Part 4: Height (Community) ====##

##==== Part 5: DBH (Community) ====##
load('Data_DBH_Community.R')
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
newDataMatrix <- dataMatrix[selectIDs,c('Temperature','Precipitation','Net_yi','BubbleSize')]
# Part 5.2: Fit
res_Net <- rma.mv(Net_yi,Net_vi,mods = ~Temperature,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Com <- rma.mv(Com_yi,Com_vi,mods = ~Temperature,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~Temperature,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
# Fit line Net
temLineNet <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLineNet) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Temperature),max(dataMatrix$Temperature),length.out = 100)
for(i in seq(1,length(tempX))){
  temLineNet[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Net,newmods = tempX[i])
  temLineNet[i,2] <- as.numeric(tempResult$pred)
  temLineNet[i,3] <- as.numeric(tempResult$ci.ub)
  temLineNet[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Com
temLineCom <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLineCom) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Temperature),max(dataMatrix$Temperature),length.out = 100)
for(i in seq(1,length(tempX))){
  temLineCom[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Com,newmods = tempX[i])
  temLineCom[i,2] <- as.numeric(tempResult$pred)
  temLineCom[i,3] <- as.numeric(tempResult$ci.ub)
  temLineCom[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Sel
temLineSel <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLineSel) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Temperature),max(dataMatrix$Temperature),length.out = 100)
for(i in seq(1,length(tempX))){
  temLineSel[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Sel,newmods = tempX[i])
  temLineSel[i,2] <- as.numeric(tempResult$pred)
  temLineSel[i,3] <- as.numeric(tempResult$ci.ub)
  temLineSel[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon Net
temCIPolygonNet <- as.data.frame(matrix(data = 0,nrow = (2*nrow(temLineNet)+1),ncol = 2))
colnames(temCIPolygonNet) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(temLineNet))){
  temCIPolygonNet$X[count] <- temLineNet$X[i]
  temCIPolygonNet$Y[count] <- temLineNet$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(temLineNet),1,-1)){
  temCIPolygonNet$X[count] <- temLineNet$X[i]
  temCIPolygonNet$Y[count] <- temLineNet$Lower[i]
  count <- count + 1
}
temCIPolygonNet$X[count] <- temLineNet$X[1]
temCIPolygonNet$Y[count] <- temLineNet$Upper[1]
# CI polygon Com
temCIPolygonCom <- as.data.frame(matrix(data = 0,nrow = (2*nrow(temLineCom)+1),ncol = 2))
colnames(temCIPolygonCom) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(temLineCom))){
  temCIPolygonCom$X[count] <- temLineCom$X[i]
  temCIPolygonCom$Y[count] <- temLineCom$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(temLineCom),1,-1)){
  temCIPolygonCom$X[count] <- temLineCom$X[i]
  temCIPolygonCom$Y[count] <- temLineCom$Lower[i]
  count <- count + 1
}
temCIPolygonCom$X[count] <- temLineCom$X[1]
temCIPolygonCom$Y[count] <- temLineCom$Upper[1]
# CI polygon Sel
temCIPolygonSel <- as.data.frame(matrix(data = 0,nrow = (2*nrow(temLineSel)+1),ncol = 2))
colnames(temCIPolygonSel) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(temLineSel))){
  temCIPolygonSel$X[count] <- temLineSel$X[i]
  temCIPolygonSel$Y[count] <- temLineSel$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(temLineSel),1,-1)){
  temCIPolygonSel$X[count] <- temLineSel$X[i]
  temCIPolygonSel$Y[count] <- temLineSel$Lower[i]
  count <- count + 1
}
temCIPolygonSel$X[count] <- temLineSel$X[1]
temCIPolygonSel$Y[count] <- temLineSel$Upper[1]
# Part 5.3: Draw
fig_S14e <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Temperature,y = Net_yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = temLineSel,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0,0.38,0.96,0),linetype = 'dashed')+
  geom_polygon(data = temCIPolygonSel,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0,0.38,0.96,0),alpha = 0.3)+
  geom_line(data = temLineCom,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.82,0.2,0.93,0.24),linetype = 'dashed')+
  geom_polygon(data = temCIPolygonCom,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),alpha = 0.3)+
  geom_line(data = temLineNet,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0),linetype = 'dashed')+
  geom_polygon(data = temCIPolygonNet,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = seq(5,25,5),labels = c('5','10','15','20','25'))+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = seq(-0.4,0.4,0.2),labels = c('-0.4','-0.2','0','0.2','0.4'),expand = c(0,0))+
  theme_custom()
##==== Part 5: DBH (Community) ====##

##==== Part 6: Biomass (Community) ====##
load('Data_Bio_Community.R')
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
newDataMatrix <- dataMatrix[selectIDs,c('Temperature','Precipitation','Net_yi','BubbleSize')]
# Part 6.2: Fit
res_Net <- rma.mv(Net_yi,Net_vi,mods = ~Temperature,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Com <- rma.mv(Com_yi,Com_vi,mods = ~Temperature,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~Temperature,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
# Fit line Net
temLineNet <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLineNet) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Temperature),max(dataMatrix$Temperature),length.out = 100)
for(i in seq(1,length(tempX))){
  temLineNet[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Net,newmods = tempX[i])
  temLineNet[i,2] <- as.numeric(tempResult$pred)
  temLineNet[i,3] <- as.numeric(tempResult$ci.ub)
  temLineNet[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Com
temLineCom <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLineCom) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Temperature),max(dataMatrix$Temperature),length.out = 100)
for(i in seq(1,length(tempX))){
  temLineCom[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Com,newmods = tempX[i])
  temLineCom[i,2] <- as.numeric(tempResult$pred)
  temLineCom[i,3] <- as.numeric(tempResult$ci.ub)
  temLineCom[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Sel
temLineSel <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLineSel) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Temperature),max(dataMatrix$Temperature),length.out = 100)
for(i in seq(1,length(tempX))){
  temLineSel[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Sel,newmods = tempX[i])
  temLineSel[i,2] <- as.numeric(tempResult$pred)
  temLineSel[i,3] <- as.numeric(tempResult$ci.ub)
  temLineSel[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon Net
temCIPolygonNet <- as.data.frame(matrix(data = 0,nrow = (2*nrow(temLineNet)+1),ncol = 2))
colnames(temCIPolygonNet) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(temLineNet))){
  temCIPolygonNet$X[count] <- temLineNet$X[i]
  temCIPolygonNet$Y[count] <- temLineNet$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(temLineNet),1,-1)){
  temCIPolygonNet$X[count] <- temLineNet$X[i]
  temCIPolygonNet$Y[count] <- temLineNet$Lower[i]
  count <- count + 1
}
temCIPolygonNet$X[count] <- temLineNet$X[1]
temCIPolygonNet$Y[count] <- temLineNet$Upper[1]
# CI polygon Com
temCIPolygonCom <- as.data.frame(matrix(data = 0,nrow = (2*nrow(temLineCom)+1),ncol = 2))
colnames(temCIPolygonCom) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(temLineCom))){
  temCIPolygonCom$X[count] <- temLineCom$X[i]
  temCIPolygonCom$Y[count] <- temLineCom$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(temLineCom),1,-1)){
  temCIPolygonCom$X[count] <- temLineCom$X[i]
  temCIPolygonCom$Y[count] <- temLineCom$Lower[i]
  count <- count + 1
}
temCIPolygonCom$X[count] <- temLineCom$X[1]
temCIPolygonCom$Y[count] <- temLineCom$Upper[1]
# CI polygon Sel
temCIPolygonSel <- as.data.frame(matrix(data = 0,nrow = (2*nrow(temLineSel)+1),ncol = 2))
colnames(temCIPolygonSel) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(temLineSel))){
  temCIPolygonSel$X[count] <- temLineSel$X[i]
  temCIPolygonSel$Y[count] <- temLineSel$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(temLineSel),1,-1)){
  temCIPolygonSel$X[count] <- temLineSel$X[i]
  temCIPolygonSel$Y[count] <- temLineSel$Lower[i]
  count <- count + 1
}
temCIPolygonSel$X[count] <- temLineSel$X[1]
temCIPolygonSel$Y[count] <- temLineSel$Upper[1]
# Part 6.3: Draw
fig_S14f <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Temperature,y = Net_yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = temLineSel,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0,0.38,0.96,0),linetype = 'dashed')+
  geom_polygon(data = temCIPolygonSel,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0,0.38,0.96,0),alpha = 0.3)+
  geom_line(data = temLineCom,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.82,0.2,0.93,0.24),linetype = 'dashed')+
  geom_polygon(data = temCIPolygonCom,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),alpha = 0.3)+
  geom_line(data = temLineNet,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_polygon(data = temCIPolygonNet,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = seq(5,25,5),labels = c('5','10','15','20','25'))+
  scale_y_continuous(limits = c(-0.9,1.2),breaks = seq(-0.6,1.2,0.6),labels = c('-0.6','0','0.6','1.2'),expand = c(0,0))+
  theme_custom()
##==== Part 6: Biomass (Community) ====##

##==== Part 7: Combine figures ====##
fig_S14 <- ggdraw()+
  draw_plot(fig_S14a,x = 0,y = 0.5,width = 0.333,height = 0.5)+
  draw_plot(fig_S14b,x = 0.334,y = 0.5,width = 0.333,height = 0.5)+
  draw_plot(fig_S14c,x = 0.667,y = 0.5,width = 0.333,height = 0.5)+
  draw_plot(fig_S14d,x = 0,y = 0,width = 0.333,height = 0.5)+
  draw_plot(fig_S14e,x = 0.334,y = 0,width = 0.333,height = 0.5)+
  draw_plot(fig_S14f,x = 0.667,y = 0,width = 0.333,height = 0.5)
outputFileName <- paste0(currentPath,'/fig_S14.pdf')
pdf(file = outputFileName,width = 4.65,height = 2.7,colormodel = 'cmyk')
print(fig_S14)
dev.off()
##==== Part 7: Combine figures ====##
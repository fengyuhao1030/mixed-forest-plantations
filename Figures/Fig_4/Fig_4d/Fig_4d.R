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
                   axis.title.x = element_text(size = 7,margin = margin(2,0,0,0),face = 'bold'),
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

load('HeightData.R')
dataMatrix <- get('HeightData')
selectIDs <- which(dataMatrix$Density <= 4)
dataMatrix <- dataMatrix[selectIDs,]

##==== Part 1: Collate data ====##
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

dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))
##==== Part 1: Collate data ====##

##==== Part 2: Fit Density ====##
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
##==== Part 2: Fit Density ====##

##==== Part 3: Draw ====##
Fig_4d <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Density,y = Net_yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = densityLineSel,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0,0.38,0.96,0),linetype = 'dashed')+
  geom_polygon(data = densityCIPolygonSel,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0,0.38,0.96,0),alpha = 0.3)+
  geom_line(data = densityLineCom,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.82,0.2,0.93,0.24),linetype = 'dashed')+
  geom_polygon(data = densityCIPolygonCom,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),alpha = 0.3)+
  geom_line(data = densityLineNet,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_polygon(data = densityCIPolygonNet,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = c(3,3.477,4),labels = c('1000','3000','10000'))+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = seq(-0.4,0.4,0.2),labels = c('-0.4','-0.2','0','0.2','0.4'),expand = c(0,0))+
  xlab('Density')+
  theme_custom()
outputFileName <- paste0(currentPath,'/Fig_4d.pdf')
pdf(file = outputFileName,width = 1.55,height = 1.46,colormodel = 'cmyk')
print(Fig_4d)
dev.off()
##==== Part 3: Draw ====##
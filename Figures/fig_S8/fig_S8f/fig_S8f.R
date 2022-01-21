rm(list = ls())
library(rstudioapi)
library(multcomp)
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

load('BioData.R')
dataMatrix <- get('BioData')
selectIDs <- which(!is.na(dataMatrix$PlotSize))
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
newDataMatrix <- dataMatrix[selectIDs,c('Age','Density','AridityIndex','TEM','PRE','PlotSize','Net_yi','BubbleSize')]

dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))
##==== Part 1: Collate data ====##

##==== Part 2: Fit plot size ====##
res_Net <- rma.mv(Net_yi,Net_vi,mods = ~PlotSize,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Com <- rma.mv(Com_yi,Com_vi,mods = ~PlotSize,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~PlotSize,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
# Fit line Net
plotsizeLineNet <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(plotsizeLineNet) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$PlotSize),max(dataMatrix$PlotSize),length.out = 100)
for(i in seq(1,length(tempX))){
  plotsizeLineNet[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Net,newmods = tempX[i])
  plotsizeLineNet[i,2] <- as.numeric(tempResult$pred)
  plotsizeLineNet[i,3] <- as.numeric(tempResult$ci.ub)
  plotsizeLineNet[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Com
plotsizeLineCom <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(plotsizeLineCom) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$PlotSize),max(dataMatrix$PlotSize),length.out = 100)
for(i in seq(1,length(tempX))){
  plotsizeLineCom[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Com,newmods = tempX[i])
  plotsizeLineCom[i,2] <- as.numeric(tempResult$pred)
  plotsizeLineCom[i,3] <- as.numeric(tempResult$ci.ub)
  plotsizeLineCom[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Sel
plotsizeLineSel <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(plotsizeLineSel) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$PlotSize),max(dataMatrix$PlotSize),length.out = 100)
for(i in seq(1,length(tempX))){
  plotsizeLineSel[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Sel,newmods = tempX[i])
  plotsizeLineSel[i,2] <- as.numeric(tempResult$pred)
  plotsizeLineSel[i,3] <- as.numeric(tempResult$ci.ub)
  plotsizeLineSel[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon Net
plotsizeCIPolygonNet <- as.data.frame(matrix(data = 0,nrow = (2*nrow(plotsizeLineNet)+1),ncol = 2))
colnames(plotsizeCIPolygonNet) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(plotsizeLineNet))){
  plotsizeCIPolygonNet$X[count] <- plotsizeLineNet$X[i]
  plotsizeCIPolygonNet$Y[count] <- plotsizeLineNet$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(plotsizeLineNet),1,-1)){
  plotsizeCIPolygonNet$X[count] <- plotsizeLineNet$X[i]
  plotsizeCIPolygonNet$Y[count] <- plotsizeLineNet$Lower[i]
  count <- count + 1
}
plotsizeCIPolygonNet$X[count] <- plotsizeLineNet$X[1]
plotsizeCIPolygonNet$Y[count] <- plotsizeLineNet$Upper[1]
# CI polygon Com
plotsizeCIPolygonCom <- as.data.frame(matrix(data = 0,nrow = (2*nrow(plotsizeLineCom)+1),ncol = 2))
colnames(plotsizeCIPolygonCom) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(plotsizeLineCom))){
  plotsizeCIPolygonCom$X[count] <- plotsizeLineCom$X[i]
  plotsizeCIPolygonCom$Y[count] <- plotsizeLineCom$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(plotsizeLineCom),1,-1)){
  plotsizeCIPolygonCom$X[count] <- plotsizeLineCom$X[i]
  plotsizeCIPolygonCom$Y[count] <- plotsizeLineCom$Lower[i]
  count <- count + 1
}
plotsizeCIPolygonCom$X[count] <- plotsizeLineCom$X[1]
plotsizeCIPolygonCom$Y[count] <- plotsizeLineCom$Upper[1]
# CI polygon Sel
plotsizeCIPolygonSel <- as.data.frame(matrix(data = 0,nrow = (2*nrow(plotsizeLineSel)+1),ncol = 2))
colnames(plotsizeCIPolygonSel) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(plotsizeLineSel))){
  plotsizeCIPolygonSel$X[count] <- plotsizeLineSel$X[i]
  plotsizeCIPolygonSel$Y[count] <- plotsizeLineSel$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(plotsizeLineSel),1,-1)){
  plotsizeCIPolygonSel$X[count] <- plotsizeLineSel$X[i]
  plotsizeCIPolygonSel$Y[count] <- plotsizeLineSel$Lower[i]
  count <- count + 1
}
plotsizeCIPolygonSel$X[count] <- plotsizeLineSel$X[1]
plotsizeCIPolygonSel$Y[count] <- plotsizeLineSel$Upper[1]
##==== Part 2: Fit plot size ====##

##==== Part 3: Draw ====##
fig_S8f <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = PlotSize,y = Net_yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = plotsizeLineSel,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0,0.38,0.96,0))+
  geom_polygon(data = plotsizeCIPolygonSel,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0,0.38,0.96,0),alpha = 0.3)+
  geom_line(data = plotsizeLineCom,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.82,0.2,0.93,0.24),linetype = 'dashed')+
  geom_polygon(data = plotsizeCIPolygonCom,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),alpha = 0.3)+
  geom_line(data = plotsizeLineNet,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0),linetype = 'dashed')+
  geom_polygon(data = plotsizeCIPolygonNet,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = seq(0,1500,500),labels = c('0','500','1000','1500'))+
  scale_y_continuous(limits = c(-0.9,1.2),breaks = seq(-0.6,1.2,0.6),labels = c('-0.6','0','0.6','1.2'),expand = c(0,0))+
  xlab('Plot size')+
  theme_custom()
outputFileName <- paste0(currentPath,'/fig_S8f.pdf')
pdf(file = outputFileName,width = 1.55,height = 1.46,colormodel = 'cmyk')
print(fig_S8f)
dev.off()
##==== Part 3: Draw ====##
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

##==== Part 2: Fit TEM ====##
res_Net <- rma.mv(Net_yi,Net_vi,mods = ~TEM,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Com <- rma.mv(Com_yi,Com_vi,mods = ~TEM,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~TEM,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
# Fit line Net
temLineNet <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(temLineNet) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$TEM),max(dataMatrix$TEM),length.out = 100)
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
tempX <- seq(min(dataMatrix$TEM),max(dataMatrix$TEM),length.out = 100)
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
tempX <- seq(min(dataMatrix$TEM),max(dataMatrix$TEM),length.out = 100)
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
##==== Part 2: Fit TEM ====##

##==== Part 3: Draw ====##
fig_S14f <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = TEM,y = Net_yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = temLineSel,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0,0.38,0.96,0),linetype = 'dashed')+
  geom_polygon(data = temCIPolygonSel,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0,0.38,0.96,0),alpha = 0.3)+
  geom_line(data = temLineCom,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.82,0.2,0.93,0.24),linetype = 'dashed')+
  geom_polygon(data = temCIPolygonCom,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),alpha = 0.3)+
  geom_line(data = temLineNet,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0),linetype = 'dashed')+
  geom_polygon(data = temCIPolygonNet,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = seq(5,25,5),labels = c('5','10','15','20','25'))+
  scale_y_continuous(limits = c(-0.9,1.2),breaks = seq(-0.6,1.2,0.6),labels = c('-0.6','0','0.6','1.2'),expand = c(0,0))+
  xlab(expression('TEM ('*degree*'C)'))+
  theme_custom()
outputFileName <- paste0(currentPath,'/fig_S14f.pdf')
pdf(file = outputFileName,width = 1.55,height = 1.46,colormodel = 'cmyk')
print(fig_S14f)
dev.off()
##==== Part 3: Draw ====##
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

load('HeightData.R')
dataMatrix <- get('HeightData')

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

##==== Part 2: Fit age ====##
res_Net <- rma.mv(Net_yi,Net_vi,mods = ~Age+Age2-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Com <- rma.mv(Com_yi,Com_vi,mods = ~Age+Age2-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~Age-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
# Fit line Net
ageLineNet <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(ageLineNet) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Age),max(dataMatrix$Age),length.out = 100)
for(i in seq(1,length(tempX))){
  ageLineNet[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Net,newmods = c(tempX[i],tempX[i]^2))
  ageLineNet[i,2] <- as.numeric(tempResult$pred)
  ageLineNet[i,3] <- as.numeric(tempResult$ci.ub)
  ageLineNet[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Com
ageLineCom <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(ageLineCom) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Age),max(dataMatrix$Age),length.out = 100)
for(i in seq(1,length(tempX))){
  ageLineCom[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Com,newmods = c(tempX[i],tempX[i]^2))
  ageLineCom[i,2] <- as.numeric(tempResult$pred)
  ageLineCom[i,3] <- as.numeric(tempResult$ci.ub)
  ageLineCom[i,4] <- as.numeric(tempResult$ci.lb)
}
# Fit line Sel
ageLineSel <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(ageLineSel) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Age),max(dataMatrix$Age),length.out = 100)
for(i in seq(1,length(tempX))){
  ageLineSel[i,1] <- tempX[i]
  tempResult <- predict.rma(res_Sel,newmods = tempX[i])
  ageLineSel[i,2] <- as.numeric(tempResult$pred)
  ageLineSel[i,3] <- as.numeric(tempResult$ci.ub)
  ageLineSel[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon Net
ageCIPolygonNet <- as.data.frame(matrix(data = 0,nrow = (2*nrow(ageLineNet)+1),ncol = 2))
colnames(ageCIPolygonNet) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(ageLineNet))){
  ageCIPolygonNet$X[count] <- ageLineNet$X[i]
  ageCIPolygonNet$Y[count] <- ageLineNet$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(ageLineNet),1,-1)){
  ageCIPolygonNet$X[count] <- ageLineNet$X[i]
  ageCIPolygonNet$Y[count] <- ageLineNet$Lower[i]
  count <- count + 1
}
ageCIPolygonNet$X[count] <- ageLineNet$X[1]
ageCIPolygonNet$Y[count] <- ageLineNet$Upper[1]
# CI polygon Com
ageCIPolygonCom <- as.data.frame(matrix(data = 0,nrow = (2*nrow(ageLineCom)+1),ncol = 2))
colnames(ageCIPolygonCom) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(ageLineCom))){
  ageCIPolygonCom$X[count] <- ageLineCom$X[i]
  ageCIPolygonCom$Y[count] <- ageLineCom$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(ageLineCom),1,-1)){
  ageCIPolygonCom$X[count] <- ageLineCom$X[i]
  ageCIPolygonCom$Y[count] <- ageLineCom$Lower[i]
  count <- count + 1
}
ageCIPolygonCom$X[count] <- ageLineCom$X[1]
ageCIPolygonCom$Y[count] <- ageLineCom$Upper[1]
# CI polygon Sel
ageCIPolygonSel <- as.data.frame(matrix(data = 0,nrow = (2*nrow(ageLineSel)+1),ncol = 2))
colnames(ageCIPolygonSel) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(ageLineSel))){
  ageCIPolygonSel$X[count] <- ageLineSel$X[i]
  ageCIPolygonSel$Y[count] <- ageLineSel$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(ageLineSel),1,-1)){
  ageCIPolygonSel$X[count] <- ageLineSel$X[i]
  ageCIPolygonSel$Y[count] <- ageLineSel$Lower[i]
  count <- count + 1
}
ageCIPolygonSel$X[count] <- ageLineSel$X[1]
ageCIPolygonSel$Y[count] <- ageLineSel$Upper[1]
##==== Part 2: Fit age ====##

##==== Part 3: Draw ====##
fig_S11d <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Age,y = Net_yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = ageLineSel,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0,0.38,0.96,0),linetype = 'dashed')+
  geom_polygon(data = ageCIPolygonSel,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0,0.38,0.96,0),alpha = 0.3)+
  geom_line(data = ageLineCom,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.82,0.2,0.93,0.24))+
  geom_polygon(data = ageCIPolygonCom,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.82,0.2,0.93,0.24),alpha = 0.3)+
  geom_line(data = ageLineNet,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_polygon(data = ageCIPolygonNet,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = seq(0,60,20),labels = c('0','20','40','60'))+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = seq(-0.4,0.4,0.2),labels = c('-0.4','-0.2','0','0.2','0.4'),expand = c(0,0))+
  xlab('Age')+
  theme_custom()
outputFileName <- paste0(currentPath,'/fig_S11d.pdf')
pdf(file = outputFileName,width = 1.55,height = 1.46,colormodel = 'cmyk')
print(fig_S11d)
dev.off()
##==== Part 3: Draw ====##
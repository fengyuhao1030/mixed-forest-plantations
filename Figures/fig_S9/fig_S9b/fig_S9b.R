rm(list = ls())
library(rstudioapi)
library(stringr)
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

load('Data_DBH.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
dataMatrix <- dataMatrix[selectIDs,]
selectIDs <- which(dataMatrix$Density <= 4)
dataMatrix <- dataMatrix[selectIDs,]

##==== Part 1: Collate data ====##
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

dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(dataMatrix$Species)
dataMatrix$Level1 <- as.factor(seq(1,nrow(dataMatrix)))
##==== Part 1: Collate data ====##

##==== Part 2: Fit density ====##
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
##==== Part 2: Fit density ====##

##==== Part 3: Draw ====##
fig_S9b <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Density,y = yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = densityLine,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0),linetype = 'dashed')+
  geom_polygon(data = densityCIPolygon,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = c(3,3.477,4),labels = c('1000','3000','10000'))+
  scale_y_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),labels = c('-0.6','-0.3','0','0.3','0.6'),expand = c(0,0))+
  xlab('Density')+
  theme_custom()
outputFileName <- paste0(currentPath,'/fig_S9b.pdf')
pdf(file = outputFileName,width = 1.55,height = 1.46,colormodel = 'cmyk')
print(fig_S9b)
dev.off()
##==== Part 3: Draw ====##
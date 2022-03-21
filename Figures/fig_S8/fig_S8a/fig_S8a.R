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

load('Data_Height.R')
dataMatrix$SpeciesRichness <- str_count(dataMatrix$Composition,'\\+') + 1
selectIDs <- which(dataMatrix$SpeciesRichness == 2)
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

##==== Part 2: Fit age ====##
res <- rma.mv(yi,vi,mods = ~Age+Age2-1,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix)
# Fit line
ageLine <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
colnames(ageLine) <- c('X','Y','Upper','Lower')
tempX <- seq(min(dataMatrix$Age),max(dataMatrix$Age),length.out = 100)
for(i in seq(1,length(tempX))){
  ageLine[i,1] <- tempX[i]
  tempResult <- predict.rma(res,newmods = c(tempX[i],tempX[i]^2))
  ageLine[i,2] <- as.numeric(tempResult$pred)
  ageLine[i,3] <- as.numeric(tempResult$ci.ub)
  ageLine[i,4] <- as.numeric(tempResult$ci.lb)
}
# CI polygon
ageCIPolygon <- as.data.frame(matrix(data = 0,nrow = (2*nrow(ageLine)+1),ncol = 2))
colnames(ageCIPolygon) <- c('X','Y')
count <- 1
for(i in seq(1,nrow(ageLine))){
  ageCIPolygon$X[count] <- ageLine$X[i]
  ageCIPolygon$Y[count] <- ageLine$Upper[i]
  count <- count + 1
}
for(i in seq(nrow(ageLine),1,-1)){
  ageCIPolygon$X[count] <- ageLine$X[i]
  ageCIPolygon$Y[count] <- ageLine$Lower[i]
  count <- count + 1
}
ageCIPolygon$X[count] <- ageLine$X[1]
ageCIPolygon$Y[count] <- ageLine$Upper[1]
##==== Part 2: Fit age ====##

##==== Part 3: Draw ====##
fig_S8a <- ggplot()+
  geom_point(newDataMatrix,mapping = aes(x = Age,y = yi),
             shape = 20,stroke = 0,color = '#666666',size = newDataMatrix$BubbleSize,alpha = 0.4)+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_line(data = ageLine,mapping = aes(x = X,y = Y),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  geom_polygon(data = ageCIPolygon,mapping = aes(x = X,y = Y),size = 0.5,fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.3)+
  scale_x_continuous(breaks = seq(0,60,20),labels = c('0','20','40','60'))+
  scale_y_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),labels = c('-0.6','-0.3','0','0.3','0.6'),expand = c(0,0))+
  xlab('Age')+
  theme_custom()
outputFileName <- paste0(currentPath,'/fig_S8a.pdf')
pdf(file = outputFileName,width = 1.55,height = 1.46,colormodel = 'cmyk')
print(fig_S8a)
dev.off()
##==== Part 3: Draw ====##
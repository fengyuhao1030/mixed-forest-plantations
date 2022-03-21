rm(list = ls())
library(rstudioapi)
library(stringr)
library(metafor)
library(multcomp)
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
                   axis.title.y = element_text(size = 7,margin = margin(0,4,0,0),face = 'bold'),
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

load('Data_DBH.R')

##==== Part 1: Analysis ====##
dataMatrix$SR <- str_count(dataMatrix$Composition,'\\+') + 1
dataMatrix$SRG <- 0
selectIDs_1 <- which(dataMatrix$SR == 2)
selectIDs_2 <- which((dataMatrix$SR == 3) | (dataMatrix$SR == 4))
selectIDs_3 <- which(dataMatrix$SR > 4)
dataMatrix$SRG[selectIDs_1] <- 1
dataMatrix$SRG[selectIDs_2] <- 2
dataMatrix$SRG[selectIDs_3] <- 3
dataMatrix$SRG <- as.factor(dataMatrix$SRG)
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(dataMatrix$Species)
dataMatrix$Level1 <- as.factor(seq(1,nrow(dataMatrix)))
res_null <- rma.mv(yi,vi,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix)
res <- rma.mv(yi,vi,mods = ~SRG-1,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix)
modelCompare <- anova.rma(res_null,res)
multiCompare <- summary(glht(res,linfct = rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))),test = adjusted('none'))
##==== Part 1: Analysis ====##

##==== Part 2: Prepare draw data ====##
figMatrix <- data.frame()
for(i in seq(1,nrow(res$b))){
  tempRowName <- row.names(res$b)[i]
  selectID <- gregexpr('SRG',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 3),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res$b[i]),
                           UpLimit = as.numeric(res$ci.ub[i]),LowLimit = as.numeric(res$ci.lb[i]))
  figMatrix <- rbind(figMatrix,tempMatrix)
}
figMatrix$Class <- factor(figMatrix$Class,levels = c('1','2','3'),ordered = TRUE)
# Points
figMatrix$XLocation <- 0
selectID <- which(figMatrix$Class == '1')
figMatrix$XLocation[selectID] <- 3
selectID <- which(figMatrix$Class == '2')
figMatrix$XLocation[selectID] <- 2
selectID <- which(figMatrix$Class == '3')
figMatrix$XLocation[selectID] <- 1
# Polygons
width <- 0.33
polygons <- data.frame()
for(i in seq(1,nrow(figMatrix))){
  point_1 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$LowLimit[i])
  point_2 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$UpLimit[i])
  point_3 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$UpLimit[i])
  point_4 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$LowLimit[i])
  tempPolygon <- rbind(point_1,point_2,point_3,point_4,point_1)
  tempPolygon$ID <- i
  polygons <- rbind(polygons,tempPolygon)
}
##==== Part 2: Prepare draw data ====##

##==== Part 3: Draw ====##
fig_S11b <- ggplot()+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID),fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean),size = 0.7,color = CMYKtoRGB(0.84,0.67,0,0))+
  scale_x_continuous(limits = c(0.4,3.6),breaks = c(1,2,3),labels = c('> 4','3 or 4','2'),expand = c(0,0))+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = seq(-0.4,0.4,0.2),labels = c('-0.4','-0.2','0','0.2','0.4'),expand = c(0,0))+
  xlab('Species richness')+
  coord_flip()+
  theme_custom()
outputFileName <- paste0(currentPath,'/fig_S11b.pdf')
pdf(file = outputFileName,width = 2,height = 1.35,colormodel = 'cmyk')
print(fig_S11b)
dev.off()
##==== Part 3: Draw ====##
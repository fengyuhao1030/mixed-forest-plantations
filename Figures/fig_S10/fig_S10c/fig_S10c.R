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
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(dataMatrix$Species)
dataMatrix$Level1 <- as.factor(seq(1,nrow(dataMatrix)))

##==== Part 1: Needleleaf ====##
dataMatrix$NeedleleafClass <- as.factor(dataMatrix$NeedleleafClass)
res_1 <- rma.mv(yi,vi,mods = ~NeedleleafClass-1,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix)
figMatrix_1 <- data.frame()
for(i in seq(1,nrow(res_1$b))){
  tempRowName <- row.names(res_1$b)[i]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_1$b[i]),
                           UpLimit = as.numeric(res_1$ci.ub[i]),LowLimit = as.numeric(res_1$ci.lb[i]))
  figMatrix_1 <- rbind(figMatrix_1,tempMatrix)
}
figMatrix_1$Class <- factor(figMatrix_1$Class,levels = c('Broad (Broad)','Broad (Broad*Needle)','Needle (Broad*Needle)','Needle (Needle)'),ordered = TRUE)
##==== Part 1: Needleleaf ====##

##==== Part 2: Evergreen ====##
dataMatrix$EvergreenClass <- as.factor(dataMatrix$EvergreenClass)
res_2 <- rma.mv(yi,vi,mods = ~EvergreenClass-1,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix)
figMatrix_2 <- data.frame()
for(i in seq(1,nrow(res_2$b))){
  tempRowName <- row.names(res_2$b)[i]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_2$b[i]),
                           UpLimit = as.numeric(res_2$ci.ub[i]),LowLimit = as.numeric(res_2$ci.lb[i]))
  figMatrix_2 <- rbind(figMatrix_2,tempMatrix)
}
figMatrix_2$Class <- factor(figMatrix_2$Class,levels = c('Deciduous (Deciduous)','Deciduous (Deciduous*Evergreen)','Evergreen (Deciduous*Evergreen)','Evergreen (Evergreen)'),ordered = TRUE)
##==== Part 2: Evergreen ====##

##==== Part 3: Nitrogen ====##
dataMatrix$NitrogenClass <- as.factor(dataMatrix$NitrogenClass)
res_3 <- rma.mv(yi,vi,mods = ~NitrogenClass-1,random = ~1|Level4/Level3/Level2/Level1,method = 'ML',data = dataMatrix)
figMatrix_3 <- data.frame()
for(i in seq(1,nrow(res_3$b))){
  tempRowName <- row.names(res_3$b)[i]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_3$b[i]),
                           UpLimit = as.numeric(res_3$ci.ub[i]),LowLimit = as.numeric(res_3$ci.lb[i]))
  figMatrix_3 <- rbind(figMatrix_3,tempMatrix)
}
figMatrix_3$Class <- factor(figMatrix_3$Class,levels = c('N (N)','N (Non-N*N)','Non-N (Non-N*N)','Non-N (Non-N)'),ordered = TRUE)
##==== Part 3: Nitrogen ====##

##==== Part 4: Multiple comparison ====##
multiCompare_1 <- summary(glht(res_1, linfct = rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))),test = adjusted('none'))
multiCompare_2 <- summary(glht(res_2, linfct = rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))),test = adjusted('none'))
multiCompare_3 <- summary(glht(res_3, linfct = rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))),test = adjusted('none'))
##==== Part 4: Multiple comparison ====##

##==== Part 5: Draw ====##
# BL vs. NL
figMatrix_1$XLocation <- 0
selectID <- which(figMatrix_1$Class == 'Broad (Broad)')
figMatrix_1$XLocation[selectID] <- 13
selectID <- which(figMatrix_1$Class == 'Broad (Broad*Needle)')
figMatrix_1$XLocation[selectID] <- 12
selectID <- which(figMatrix_1$Class == 'Needle (Broad*Needle)')
figMatrix_1$XLocation[selectID] <- 11
selectID <- which(figMatrix_1$Class == 'Needle (Needle)')
figMatrix_1$XLocation[selectID] <- 10
# DE vs. EV
figMatrix_2$XLocation <- 0
selectID <- which(figMatrix_2$Class == 'Deciduous (Deciduous)')
figMatrix_2$XLocation[selectID] <- 8.5
selectID <- which(figMatrix_2$Class == 'Deciduous (Deciduous*Evergreen)')
figMatrix_2$XLocation[selectID] <- 7.5
selectID <- which(figMatrix_2$Class == 'Evergreen (Deciduous*Evergreen)')
figMatrix_2$XLocation[selectID] <- 6.5
selectID <- which(figMatrix_2$Class == 'Evergreen (Evergreen)')
figMatrix_2$XLocation[selectID] <- 5.5
# NF vs. non-NF
figMatrix_3$XLocation <- 0
selectID <- which(figMatrix_3$Class == 'N (N)')
figMatrix_3$XLocation[selectID] <- 4
selectID <- which(figMatrix_3$Class == 'N (Non-N*N)')
figMatrix_3$XLocation[selectID] <- 3
selectID <- which(figMatrix_3$Class == 'Non-N (Non-N*N)')
figMatrix_3$XLocation[selectID] <- 2
selectID <- which(figMatrix_3$Class == 'Non-N (Non-N)')
figMatrix_3$XLocation[selectID] <- 1
figMatrix <- rbind(figMatrix_1,figMatrix_2,figMatrix_3)
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
# Draw
fig_S10c <- ggplot()+
  geom_vline(xintercept = 4.75,size = 0.3,color = '#000000')+
  geom_vline(xintercept = 9.25,size = 0.3,color = '#000000')+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID),fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  scale_x_continuous(breaks = c(1,2,3,4,5.5,6.5,7.5,8.5,10,11,12,13),
                     labels = c('non-NF (non-NF x non-NF)','non-NF (NF x non-NF)','NF (NF x non-NF)','NF (NF x NF)',
                                'EV (EV x EV)','EV (DE x EV)','DE (DE x EV)','DE (DE x DE)',
                                'NL (NL x NL)','NL (BL x NL)','BL (BL x NL)','BL (BL x BL)'))+
  scale_y_continuous(limits = c(-0.7,0.7),breaks = seq(-0.6,0.6,0.3),labels = c('-0.6','-0.3','0','0.3','0.6'),expand = c(0,0))+
  ylab('Effect size')+
  coord_flip()+
  theme_custom()
outputFileName <- paste0(currentPath,'/fig_S10c.pdf')
pdf(file = outputFileName,width = 2.45,height = 3,colormodel = 'cmyk')
print(fig_S10c)
dev.off()
##==== Part 5: Draw ====##
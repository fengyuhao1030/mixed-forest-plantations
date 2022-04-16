rm(list = ls())
library(rstudioapi)
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
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(dataMatrix$Species)
dataMatrix$Level1 <- as.factor(seq(1,nrow(dataMatrix)))
# Part 1.1: Needleleaf
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
# Part 1.2: Evergreen
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
# Part 1.3: Nitrogen
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
# Part 1.4: Multiple comparison
multiCompare_1 <- summary(glht(res_1, linfct = rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))),test = adjusted('none'))
multiCompare_2 <- summary(glht(res_2, linfct = rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))),test = adjusted('none'))
multiCompare_3 <- summary(glht(res_3, linfct = rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))),test = adjusted('none'))
# Part 1.5: Draw
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
Fig_2a <- ggplot()+
  geom_vline(xintercept = 4.75,size = 0.3,color = '#000000')+
  geom_vline(xintercept = 9.25,size = 0.3,color = '#000000')+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID),fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  scale_x_continuous(breaks = c(1,2,3,4,5.5,6.5,7.5,8.5,10,11,12,13),
                     labels = c('non-NF (non-NF x non-NF)','non-NF (NF x non-NF)','NF (NF x non-NF)','NF (NF x NF)',
                                'EV (EV x EV)','EV (DE x EV)','DE (DE x EV)','DE (DE x DE)',
                                'NL (NL x NL)','NL (BL x NL)','BL (BL x NL)','BL (BL x BL)'))+
  scale_y_continuous(limits = c(-0.22,0.22),breaks = seq(-0.2,0.2,0.1),labels = c('-0.2','-0.1','0','0.1','0.2'),expand = c(0,0))+
  coord_flip()+
  theme_custom()
##==== Part 1: Height (Species) ====##

##==== Part 2: DBH (Species) ====##
load('Data_DBH_Species.R')
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(dataMatrix$Species)
dataMatrix$Level1 <- as.factor(seq(1,nrow(dataMatrix)))
# Part 2.1: Needleleaf
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
# Part 2.2: Evergreen
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
# Part 2.3: Nitrogen
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
# Part 2.4: Multiple comparison
multiCompare_1 <- summary(glht(res_1, linfct = rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))),test = adjusted('none'))
multiCompare_2 <- summary(glht(res_2, linfct = rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))),test = adjusted('none'))
multiCompare_3 <- summary(glht(res_3, linfct = rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))),test = adjusted('none'))
# Part 2.5: Draw
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
Fig_2b <- ggplot()+
  geom_vline(xintercept = 4.75,size = 0.3,color = '#000000')+
  geom_vline(xintercept = 9.25,size = 0.3,color = '#000000')+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID),fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean),size = 0.5,color = CMYKtoRGB(0.84,0.67,0,0))+
  scale_x_continuous(breaks = c(1,2,3,4,5.5,6.5,7.5,8.5,10,11,12,13),
                     labels = c('non-NF (non-NF x non-NF)','non-NF (NF x non-NF)','NF (NF x non-NF)','NF (NF x NF)',
                                'EV (EV x EV)','EV (DE x EV)','DE (DE x EV)','DE (DE x DE)',
                                'NL (NL x NL)','NL (BL x NL)','BL (BL x NL)','BL (BL x BL)'))+
  scale_y_continuous(limits = c(-0.22,0.22),breaks = seq(-0.2,0.2,0.1),labels = c('-0.2','-0.1','0','0.1','0.2'),expand = c(0,0))+
  coord_flip()+
  theme_custom()+
  theme(axis.text.y = element_blank())
##==== Part 2: DBH (Species) ====##

##==== Part 3: Biomass (Species) ====##
load('Data_Bio_Species.R')
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(dataMatrix$Species)
dataMatrix$Level1 <- as.factor(seq(1,nrow(dataMatrix)))
# Part 3.1: Needleleaf
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
# Part 3.2: Evergreen
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
# Part 3.3: Nitrogen
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
# Part 3.4: Multiple comparison
multiCompare_1 <- summary(glht(res_1, linfct = rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))),test = adjusted('none'))
multiCompare_2 <- summary(glht(res_2, linfct = rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))),test = adjusted('none'))
multiCompare_3 <- summary(glht(res_3, linfct = rbind(c(-1,1,0,0),c(-1,0,1,0),c(-1,0,0,1),c(0,-1,1,0),c(0,-1,0,1),c(0,0,-1,1))),test = adjusted('none'))
# Part 3.5: Draw
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
Fig_2c <- ggplot()+
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
  coord_flip()+
  theme_custom()+
  theme(axis.text.y = element_blank())
##==== Part 3: Biomass (Species) ====##

##==== Part 4: Height (Community) ====##
load('Data_Height_Community.R')
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))
# Part 4.1: Needleleaf
dataMatrix$NeedleleafClass <- as.factor(dataMatrix$NeedleleafClass)
res_1_Net <- rma.mv(Net_yi,Net_vi,mods = ~NeedleleafClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_1_Com <- rma.mv(Com_yi,Com_vi,mods = ~NeedleleafClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_1_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~NeedleleafClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
figMatrix_1 <- data.frame()
for(i in seq(1,nrow(res_1_Net$b))){
  # Net
  tempRowName <- row.names(res_1_Net$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_1_Net$b[i]),
                           UpLimit = as.numeric(res_1_Net$ci.ub[i]),LowLimit = as.numeric(res_1_Net$ci.lb[i]),
                           Group = 'Net')
  figMatrix_1 <- rbind(figMatrix_1,tempMatrix)
  # Com
  tempRowName <- row.names(res_1_Com$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_1_Com$b[i]),
                           UpLimit = as.numeric(res_1_Com$ci.ub[i]),LowLimit = as.numeric(res_1_Com$ci.lb[i]),
                           Group = 'Com')
  figMatrix_1 <- rbind(figMatrix_1,tempMatrix)
  # Sel
  tempRowName <- row.names(res_1_Sel$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID+5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_1_Sel$b[i]),
                           UpLimit = as.numeric(res_1_Sel$ci.ub[i]),LowLimit = as.numeric(res_1_Sel$ci.lb[i]),
                           Group = 'Sel')
  figMatrix_1 <- rbind(figMatrix_1,tempMatrix)
}
# Part 4.2: Evergreen
dataMatrix$EvergreenClass <- as.factor(dataMatrix$EvergreenClass)
res_2_Net <- rma.mv(Net_yi,Net_vi,mods = ~EvergreenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_2_Com <- rma.mv(Com_yi,Com_vi,mods = ~EvergreenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_2_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~EvergreenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
figMatrix_2 <- data.frame()
for(i in seq(1,nrow(res_2_Net$b))){
  # Net
  tempRowName <- row.names(res_2_Net$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_2_Net$b[i]),
                           UpLimit = as.numeric(res_2_Net$ci.ub[i]),LowLimit = as.numeric(res_2_Net$ci.lb[i]),
                           Group = 'Net')
  figMatrix_2 <- rbind(figMatrix_2,tempMatrix)
  # Com
  tempRowName <- row.names(res_2_Com$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_2_Com$b[i]),
                           UpLimit = as.numeric(res_2_Com$ci.ub[i]),LowLimit = as.numeric(res_2_Com$ci.lb[i]),
                           Group = 'Com')
  figMatrix_2 <- rbind(figMatrix_2,tempMatrix)
  # Sel
  tempRowName <- row.names(res_2_Sel$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_2_Sel$b[i]),
                           UpLimit = as.numeric(res_2_Sel$ci.ub[i]),LowLimit = as.numeric(res_2_Sel$ci.lb[i]),
                           Group = 'Sel')
  figMatrix_2 <- rbind(figMatrix_2,tempMatrix)
}
# Part 4.3: Nitrogen
dataMatrix$NitrogenClass <- as.factor(dataMatrix$NitrogenClass)
res_3_Net <- rma.mv(Net_yi,Net_vi,mods = ~NitrogenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_3_Com <- rma.mv(Com_yi,Com_vi,mods = ~NitrogenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_3_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~NitrogenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
figMatrix_3 <- data.frame()
for(i in seq(1,nrow(res_3_Net$b))){
  # Net
  tempRowName <- row.names(res_3_Net$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_3_Net$b[i]),
                           UpLimit = as.numeric(res_3_Net$ci.ub[i]),LowLimit = as.numeric(res_3_Net$ci.lb[i]),
                           Group = 'Net')
  figMatrix_3 <- rbind(figMatrix_3,tempMatrix)
  # Com
  tempRowName <- row.names(res_3_Com$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_3_Com$b[i]),
                           UpLimit = as.numeric(res_3_Com$ci.ub[i]),LowLimit = as.numeric(res_3_Com$ci.lb[i]),
                           Group = 'Com')
  figMatrix_3 <- rbind(figMatrix_3,tempMatrix)
  # Sel
  tempRowName <- row.names(res_3_Sel$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_3_Sel$b[i]),
                           UpLimit = as.numeric(res_3_Sel$ci.ub[i]),LowLimit = as.numeric(res_3_Sel$ci.lb[i]),
                           Group = 'Sel')
  figMatrix_3 <- rbind(figMatrix_3,tempMatrix)
}
# Part 4.4: Draw
# BL vs. NL
figMatrix_1$XLocation <- 0
selectID <- which(figMatrix_1$Class == 'Needle/Broad')
figMatrix_1$XLocation[selectID] <- 7
selectID <- which(figMatrix_1$Class == 'Needle*Broad')
figMatrix_1$XLocation[selectID] <- 6
# DE vs. EV
figMatrix_2$XLocation <- 0
selectID <- which(figMatrix_2$Class == 'Deciduous/Evergreen')
figMatrix_2$XLocation[selectID] <- 4.5
selectID <- which(figMatrix_2$Class == 'Deciduous*Evergreen')
figMatrix_2$XLocation[selectID] <- 3.5
# NF vs. non-NF
figMatrix_3$XLocation <- 0
selectID <- which(figMatrix_3$Class == 'Non-N/N')
figMatrix_3$XLocation[selectID] <- 2
selectID <- which(figMatrix_3$Class == 'Non-N*N')
figMatrix_3$XLocation[selectID] <- 1
figMatrix <- rbind(figMatrix_1,figMatrix_2,figMatrix_3)
# Polygons
width <- 0.2
polygons <- data.frame()
for(i in seq(1,nrow(figMatrix))){
  point_1 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$LowLimit[i])
  point_2 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$UpLimit[i])
  point_3 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$UpLimit[i])
  point_4 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$LowLimit[i])
  tempPolygon <- rbind(point_1,point_2,point_3,point_4,point_1)
  tempPolygon$ID <- i
  tempPolygon$Group <- figMatrix$Group[i]
  polygons <- rbind(polygons,tempPolygon)
}
dodge <- 0.3
selectID_1 <- which(figMatrix$Group == 'Net')
selectID_2 <- which(figMatrix$Group == 'Sel')
figMatrix$XLocation[selectID_1] <- figMatrix$XLocation[selectID_1] + dodge
figMatrix$XLocation[selectID_2] <- figMatrix$XLocation[selectID_2] - dodge
figMatrix$Group <- factor(figMatrix$Group,levels = c('Net','Com','Sel'),ordered = TRUE)
selectID_1 <- which(polygons$Group == 'Net')
selectID_2 <- which(polygons$Group == 'Sel')
polygons$X[selectID_1] <- polygons$X[selectID_1] + dodge
polygons$X[selectID_2] <- polygons$X[selectID_2] - dodge
polygons$Group <- factor(polygons$Group,levels = c('Net','Com','Sel'),ordered = TRUE)
# Draw
Fig_2d <- ggplot()+
  geom_vline(xintercept = 2.75,size = 0.3,color = '#000000')+
  geom_vline(xintercept = 5.25,size = 0.3,color = '#000000')+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID,fill = Group),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean,group = Group,color = Group),size = 0.5)+
  scale_color_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.20,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_fill_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.20,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_x_continuous(breaks = c(1,2,3.5,4.5,6,7),labels = c('NF x non-NF','NF / non-NF','DE x EV','DE / EV','BL x NL','BL / NL'))+
  scale_y_continuous(limits = c(-0.22,0.22),breaks = seq(-0.2,0.2,0.1),labels = c('-0.2','-0.1','0','0.1','0.2'),expand = c(0,0))+
  coord_flip()+
  theme_custom()
##==== Part 4: Height (Community) ====##

##==== Part 5: DBH (Community) ====##
load('Data_DBH_Community.R')
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))
# Part 5.1: Needleleaf
dataMatrix$NeedleleafClass <- as.factor(dataMatrix$NeedleleafClass)
res_1_Net <- rma.mv(Net_yi,Net_vi,mods = ~NeedleleafClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_1_Com <- rma.mv(Com_yi,Com_vi,mods = ~NeedleleafClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_1_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~NeedleleafClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
figMatrix_1 <- data.frame()
for(i in seq(1,nrow(res_1_Net$b))){
  # Net
  tempRowName <- row.names(res_1_Net$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_1_Net$b[i]),
                           UpLimit = as.numeric(res_1_Net$ci.ub[i]),LowLimit = as.numeric(res_1_Net$ci.lb[i]),
                           Group = 'Net')
  figMatrix_1 <- rbind(figMatrix_1,tempMatrix)
  # Com
  tempRowName <- row.names(res_1_Com$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_1_Com$b[i]),
                           UpLimit = as.numeric(res_1_Com$ci.ub[i]),LowLimit = as.numeric(res_1_Com$ci.lb[i]),
                           Group = 'Com')
  figMatrix_1 <- rbind(figMatrix_1,tempMatrix)
  # Sel
  tempRowName <- row.names(res_1_Sel$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID+5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_1_Sel$b[i]),
                           UpLimit = as.numeric(res_1_Sel$ci.ub[i]),LowLimit = as.numeric(res_1_Sel$ci.lb[i]),
                           Group = 'Sel')
  figMatrix_1 <- rbind(figMatrix_1,tempMatrix)
}
# Part 5.2: Evergreen
dataMatrix$EvergreenClass <- as.factor(dataMatrix$EvergreenClass)
res_2_Net <- rma.mv(Net_yi,Net_vi,mods = ~EvergreenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_2_Com <- rma.mv(Com_yi,Com_vi,mods = ~EvergreenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_2_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~EvergreenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
figMatrix_2 <- data.frame()
for(i in seq(1,nrow(res_2_Net$b))){
  # Net
  tempRowName <- row.names(res_2_Net$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_2_Net$b[i]),
                           UpLimit = as.numeric(res_2_Net$ci.ub[i]),LowLimit = as.numeric(res_2_Net$ci.lb[i]),
                           Group = 'Net')
  figMatrix_2 <- rbind(figMatrix_2,tempMatrix)
  # Com
  tempRowName <- row.names(res_2_Com$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_2_Com$b[i]),
                           UpLimit = as.numeric(res_2_Com$ci.ub[i]),LowLimit = as.numeric(res_2_Com$ci.lb[i]),
                           Group = 'Com')
  figMatrix_2 <- rbind(figMatrix_2,tempMatrix)
  # Sel
  tempRowName <- row.names(res_2_Sel$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_2_Sel$b[i]),
                           UpLimit = as.numeric(res_2_Sel$ci.ub[i]),LowLimit = as.numeric(res_2_Sel$ci.lb[i]),
                           Group = 'Sel')
  figMatrix_2 <- rbind(figMatrix_2,tempMatrix)
}
# Part 5.3: Nitrogen
dataMatrix$NitrogenClass <- as.factor(dataMatrix$NitrogenClass)
res_3_Net <- rma.mv(Net_yi,Net_vi,mods = ~NitrogenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_3_Com <- rma.mv(Com_yi,Com_vi,mods = ~NitrogenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_3_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~NitrogenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
figMatrix_3 <- data.frame()
for(i in seq(1,nrow(res_3_Net$b))){
  # Net
  tempRowName <- row.names(res_3_Net$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_3_Net$b[i]),
                           UpLimit = as.numeric(res_3_Net$ci.ub[i]),LowLimit = as.numeric(res_3_Net$ci.lb[i]),
                           Group = 'Net')
  figMatrix_3 <- rbind(figMatrix_3,tempMatrix)
  # Com
  tempRowName <- row.names(res_3_Com$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_3_Com$b[i]),
                           UpLimit = as.numeric(res_3_Com$ci.ub[i]),LowLimit = as.numeric(res_3_Com$ci.lb[i]),
                           Group = 'Com')
  figMatrix_3 <- rbind(figMatrix_3,tempMatrix)
  # Sel
  tempRowName <- row.names(res_3_Sel$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_3_Sel$b[i]),
                           UpLimit = as.numeric(res_3_Sel$ci.ub[i]),LowLimit = as.numeric(res_3_Sel$ci.lb[i]),
                           Group = 'Sel')
  figMatrix_3 <- rbind(figMatrix_3,tempMatrix)
}
# Part 5.4: Draw
# BL vs. NL
figMatrix_1$XLocation <- 0
selectID <- which(figMatrix_1$Class == 'Needle/Broad')
figMatrix_1$XLocation[selectID] <- 7
selectID <- which(figMatrix_1$Class == 'Needle*Broad')
figMatrix_1$XLocation[selectID] <- 6
# DE vs. EV
figMatrix_2$XLocation <- 0
selectID <- which(figMatrix_2$Class == 'Deciduous/Evergreen')
figMatrix_2$XLocation[selectID] <- 4.5
selectID <- which(figMatrix_2$Class == 'Deciduous*Evergreen')
figMatrix_2$XLocation[selectID] <- 3.5
# NF vs. non-NF
figMatrix_3$XLocation <- 0
selectID <- which(figMatrix_3$Class == 'Non-N/N')
figMatrix_3$XLocation[selectID] <- 2
selectID <- which(figMatrix_3$Class == 'Non-N*N')
figMatrix_3$XLocation[selectID] <- 1
figMatrix <- rbind(figMatrix_1,figMatrix_2,figMatrix_3)
# Polygons
width <- 0.2
polygons <- data.frame()
for(i in seq(1,nrow(figMatrix))){
  point_1 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$LowLimit[i])
  point_2 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$UpLimit[i])
  point_3 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$UpLimit[i])
  point_4 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$LowLimit[i])
  tempPolygon <- rbind(point_1,point_2,point_3,point_4,point_1)
  tempPolygon$ID <- i
  tempPolygon$Group <- figMatrix$Group[i]
  polygons <- rbind(polygons,tempPolygon)
}
dodge <- 0.3
selectID_1 <- which(figMatrix$Group == 'Net')
selectID_2 <- which(figMatrix$Group == 'Sel')
figMatrix$XLocation[selectID_1] <- figMatrix$XLocation[selectID_1] + dodge
figMatrix$XLocation[selectID_2] <- figMatrix$XLocation[selectID_2] - dodge
figMatrix$Group <- factor(figMatrix$Group,levels = c('Net','Com','Sel'),ordered = TRUE)
selectID_1 <- which(polygons$Group == 'Net')
selectID_2 <- which(polygons$Group == 'Sel')
polygons$X[selectID_1] <- polygons$X[selectID_1] + dodge
polygons$X[selectID_2] <- polygons$X[selectID_2] - dodge
polygons$Group <- factor(polygons$Group,levels = c('Net','Com','Sel'),ordered = TRUE)
# Draw
Fig_2e <- ggplot()+
  geom_vline(xintercept = 2.75,size = 0.3,color = '#000000')+
  geom_vline(xintercept = 5.25,size = 0.3,color = '#000000')+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID,fill = Group),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean,group = Group,color = Group),size = 0.5)+
  scale_color_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.20,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_fill_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.20,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_x_continuous(breaks = c(1,2,3.5,4.5,6,7),labels = c('NF x non-NF','NF / non-NF','DE x EV','DE / EV','BL x NL','BL / NL'))+
  scale_y_continuous(limits = c(-0.22,0.22),breaks = seq(-0.2,0.2,0.1),labels = c('-0.2','-0.1','0','0.1','0.2'),expand = c(0,0))+
  coord_flip()+
  theme_custom()+
  theme(axis.text.y = element_blank())
##==== Part 5: DBH (Community) ====##

##==== Part 6: Biomass (Community) ====##
load('Data_Bio_Community.R')
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))
# Part 6.1: Needleleaf
dataMatrix$NeedleleafClass <- as.factor(dataMatrix$NeedleleafClass)
res_1_Net <- rma.mv(Net_yi,Net_vi,mods = ~NeedleleafClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_1_Com <- rma.mv(Com_yi,Com_vi,mods = ~NeedleleafClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_1_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~NeedleleafClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
figMatrix_1 <- data.frame()
for(i in seq(1,nrow(res_1_Net$b))){
  # Net
  tempRowName <- row.names(res_1_Net$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_1_Net$b[i]),
                           UpLimit = as.numeric(res_1_Net$ci.ub[i]),LowLimit = as.numeric(res_1_Net$ci.lb[i]),
                           Group = 'Net')
  figMatrix_1 <- rbind(figMatrix_1,tempMatrix)
  # Com
  tempRowName <- row.names(res_1_Com$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_1_Com$b[i]),
                           UpLimit = as.numeric(res_1_Com$ci.ub[i]),LowLimit = as.numeric(res_1_Com$ci.lb[i]),
                           Group = 'Com')
  figMatrix_1 <- rbind(figMatrix_1,tempMatrix)
  # Sel
  tempRowName <- row.names(res_1_Sel$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID+5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_1_Sel$b[i]),
                           UpLimit = as.numeric(res_1_Sel$ci.ub[i]),LowLimit = as.numeric(res_1_Sel$ci.lb[i]),
                           Group = 'Sel')
  figMatrix_1 <- rbind(figMatrix_1,tempMatrix)
}
# Part 6.2: Evergreen
dataMatrix$EvergreenClass <- as.factor(dataMatrix$EvergreenClass)
res_2_Net <- rma.mv(Net_yi,Net_vi,mods = ~EvergreenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_2_Com <- rma.mv(Com_yi,Com_vi,mods = ~EvergreenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_2_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~EvergreenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
figMatrix_2 <- data.frame()
for(i in seq(1,nrow(res_2_Net$b))){
  # Net
  tempRowName <- row.names(res_2_Net$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_2_Net$b[i]),
                           UpLimit = as.numeric(res_2_Net$ci.ub[i]),LowLimit = as.numeric(res_2_Net$ci.lb[i]),
                           Group = 'Net')
  figMatrix_2 <- rbind(figMatrix_2,tempMatrix)
  # Com
  tempRowName <- row.names(res_2_Com$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_2_Com$b[i]),
                           UpLimit = as.numeric(res_2_Com$ci.ub[i]),LowLimit = as.numeric(res_2_Com$ci.lb[i]),
                           Group = 'Com')
  figMatrix_2 <- rbind(figMatrix_2,tempMatrix)
  # Sel
  tempRowName <- row.names(res_2_Sel$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_2_Sel$b[i]),
                           UpLimit = as.numeric(res_2_Sel$ci.ub[i]),LowLimit = as.numeric(res_2_Sel$ci.lb[i]),
                           Group = 'Sel')
  figMatrix_2 <- rbind(figMatrix_2,tempMatrix)
}
# Part 6.3: Nitrogen
dataMatrix$NitrogenClass <- as.factor(dataMatrix$NitrogenClass)
res_3_Net <- rma.mv(Net_yi,Net_vi,mods = ~NitrogenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_3_Com <- rma.mv(Com_yi,Com_vi,mods = ~NitrogenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_3_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~NitrogenClass-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
figMatrix_3 <- data.frame()
for(i in seq(1,nrow(res_3_Net$b))){
  # Net
  tempRowName <- row.names(res_3_Net$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_3_Net$b[i]),
                           UpLimit = as.numeric(res_3_Net$ci.ub[i]),LowLimit = as.numeric(res_3_Net$ci.lb[i]),
                           Group = 'Net')
  figMatrix_3 <- rbind(figMatrix_3,tempMatrix)
  # Com
  tempRowName <- row.names(res_3_Com$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_3_Com$b[i]),
                           UpLimit = as.numeric(res_3_Com$ci.ub[i]),LowLimit = as.numeric(res_3_Com$ci.lb[i]),
                           Group = 'Com')
  figMatrix_3 <- rbind(figMatrix_3,tempMatrix)
  # Sel
  tempRowName <- row.names(res_3_Sel$b)[[i]]
  selectID <- gregexpr('Class',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 5),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Mean = as.numeric(res_3_Sel$b[i]),
                           UpLimit = as.numeric(res_3_Sel$ci.ub[i]),LowLimit = as.numeric(res_3_Sel$ci.lb[i]),
                           Group = 'Sel')
  figMatrix_3 <- rbind(figMatrix_3,tempMatrix)
}
# Part 6.4: Draw
# BL vs. NL
figMatrix_1$XLocation <- 0
selectID <- which(figMatrix_1$Class == 'Needle/Broad')
figMatrix_1$XLocation[selectID] <- 7
selectID <- which(figMatrix_1$Class == 'Needle*Broad')
figMatrix_1$XLocation[selectID] <- 6
# DE vs. EV
figMatrix_2$XLocation <- 0
selectID <- which(figMatrix_2$Class == 'Deciduous/Evergreen')
figMatrix_2$XLocation[selectID] <- 4.5
selectID <- which(figMatrix_2$Class == 'Deciduous*Evergreen')
figMatrix_2$XLocation[selectID] <- 3.5
# NF vs. non-NF
figMatrix_3$XLocation <- 0
selectID <- which(figMatrix_3$Class == 'Non-N/N')
figMatrix_3$XLocation[selectID] <- 2
selectID <- which(figMatrix_3$Class == 'Non-N*N')
figMatrix_3$XLocation[selectID] <- 1
figMatrix <- rbind(figMatrix_1,figMatrix_2,figMatrix_3)
# Polygons
width <- 0.2
polygons <- data.frame()
for(i in seq(1,nrow(figMatrix))){
  point_1 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$LowLimit[i])
  point_2 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$UpLimit[i])
  point_3 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$UpLimit[i])
  point_4 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$LowLimit[i])
  tempPolygon <- rbind(point_1,point_2,point_3,point_4,point_1)
  tempPolygon$ID <- i
  tempPolygon$Group <- figMatrix$Group[i]
  polygons <- rbind(polygons,tempPolygon)
}
dodge <- 0.3
selectID_1 <- which(figMatrix$Group == 'Net')
selectID_2 <- which(figMatrix$Group == 'Sel')
figMatrix$XLocation[selectID_1] <- figMatrix$XLocation[selectID_1] + dodge
figMatrix$XLocation[selectID_2] <- figMatrix$XLocation[selectID_2] - dodge
figMatrix$Group <- factor(figMatrix$Group,levels = c('Net','Com','Sel'),ordered = TRUE)
selectID_1 <- which(polygons$Group == 'Net')
selectID_2 <- which(polygons$Group == 'Sel')
polygons$X[selectID_1] <- polygons$X[selectID_1] + dodge
polygons$X[selectID_2] <- polygons$X[selectID_2] - dodge
polygons$Group <- factor(polygons$Group,levels = c('Net','Com','Sel'),ordered = TRUE)
# Draw
Fig_2f <- ggplot()+
  geom_vline(xintercept = 2.75,size = 0.3,color = '#000000')+
  geom_vline(xintercept = 5.25,size = 0.3,color = '#000000')+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID,fill = Group),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean,group = Group,color = Group),size = 0.5)+
  scale_color_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.20,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_fill_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.20,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_x_continuous(breaks = c(1,2,3.5,4.5,6,7),labels = c('NF x non-NF','NF / non-NF','DE x EV','DE / EV','BL x NL','BL / NL'))+
  scale_y_continuous(limits = c(-0.5,0.5),breaks = seq(-0.4,0.4,0.2),labels = c('-0.4','-0.2','0','0.2','0.4'),expand = c(0,0))+
  coord_flip()+
  theme_custom()+
  theme(axis.text.y = element_blank())
##==== Part 6: Biomass (Community) ====##

##==== Part 7: Combine figures ====##
Fig_2 <- ggdraw()+
  draw_plot(Fig_2a,x = 0,y = 0.5,width = 0.52,height = 0.5)+
  draw_plot(Fig_2b,x = 0.52,y = 0.5,width = 0.24,height = 0.5)+
  draw_plot(Fig_2c,x = 0.76,y = 0.5,width = 0.24,height = 0.5)+
  draw_plot(Fig_2d,x = 0.14,y = 0,width = 0.38,height = 0.5)+
  draw_plot(Fig_2e,x = 0.52,y = 0,width = 0.24,height = 0.5)+
  draw_plot(Fig_2f,x = 0.76,y = 0,width = 0.24,height = 0.5)
outputFileName <- paste0(currentPath,'/Fig_2.pdf')
pdf(file = outputFileName,width = 4.68,height = 5.6,colormodel = 'cmyk')
print(Fig_2)
dev.off()
##==== Part 7: Combine figures ====##
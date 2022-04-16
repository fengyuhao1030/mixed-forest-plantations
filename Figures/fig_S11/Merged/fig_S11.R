rm(list = ls())
library(rstudioapi)
library(stringr)
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
# Part 1.1: Analysis
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
# Part 1.2: Prepare draw data
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
# Part 1.3: Draw
fig_S11a <- ggplot()+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID),fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean),size = 0.7,color = CMYKtoRGB(0.84,0.67,0,0))+
  scale_x_continuous(limits = c(0.4,3.6),breaks = c(1,2,3),labels = c('> 4','3 or 4','2'),expand = c(0,0))+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = seq(-0.4,0.4,0.2),labels = c('-0.4','-0.2','0','0.2','0.4'),expand = c(0,0))+
  coord_flip()+
  theme_custom()
##==== Part 1: Height (Species) ====##

##==== Part 2: DBH (Species) ====##
load('Data_DBH_Species.R')
# Part 2.1: Analysis
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
# Part 2.2: Prepare draw data
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
# Part 2.3: Draw
fig_S11b <- ggplot()+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID),fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean),size = 0.7,color = CMYKtoRGB(0.84,0.67,0,0))+
  scale_x_continuous(limits = c(0.4,3.6),breaks = c(1,2,3),labels = c('> 4','3 or 4','2'),expand = c(0,0))+
  scale_y_continuous(limits = c(-0.4,0.4),breaks = seq(-0.4,0.4,0.2),labels = c('-0.4','-0.2','0','0.2','0.4'),expand = c(0,0))+
  coord_flip()+
  theme_custom()
##==== Part 2: DBH (Species) ====##

##==== Part 3: Biomass (Species) ====##
load('Data_Bio_Species.R')
# Part 3.1: Analysis
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
# Part 3.2: Prepare draw data
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
# Part 3.3: Draw
fig_S11c <- ggplot()+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID),fill = CMYKtoRGB(0.84,0.67,0,0),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean),size = 0.7,color = CMYKtoRGB(0.84,0.67,0,0))+
  scale_x_continuous(limits = c(0.4,3.6),breaks = c(1,2,3),labels = c('> 4','3 or 4','2'),expand = c(0,0))+
  scale_y_continuous(limits = c(-0.8,0.8),breaks = seq(-0.8,0.8,0.4),labels = c('-0.8','-0.4','0','0.4','0.8'),expand = c(0,0))+
  coord_flip()+
  theme_custom()
##==== Part 3: Biomass (Species) ====##

##==== Part 4: Height (Community) ====##
load('Data_Height_Community.R')
# Part 4.1: Analysis
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
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))
# Net
res_null_Net <- rma.mv(Net_yi,Net_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Net <- rma.mv(Net_yi,Net_vi,mods = ~SRG-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
modelCompare_Net <- anova.rma(res_null_Net,res_Net)
multiCompare_Net <- summary(glht(res_Net,linfct = rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))),test = adjusted('none'))
# Com
res_null_Com <- rma.mv(Com_yi,Com_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Com <- rma.mv(Com_yi,Com_vi,mods = ~SRG-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
modelCompare_Com <- anova.rma(res_null_Com,res_Com)
multiCompare_Com <- summary(glht(res_Com,linfct = rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))),test = adjusted('none'))
# Sel
res_null_Sel <- rma.mv(Sel_yi,Sel_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~SRG-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
modelCompare_Sel <- anova.rma(res_null_Sel,res_Sel)
multiCompare_Sel <- summary(glht(res_Sel,linfct = rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))),test = adjusted('none'))
# Part 4.2: Prepare draw data
figMatrix_Net <- data.frame()
for(i in seq(1,nrow(res_Net$b))){
  tempRowName <- row.names(res_Net$b)[i]
  selectID <- gregexpr('SRG',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 3),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Group = 'Net',Mean = as.numeric(res_Net$b[i]),
                           UpLimit = as.numeric(res_Net$ci.ub[i]),LowLimit = as.numeric(res_Net$ci.lb[i]))
  figMatrix_Net <- rbind(figMatrix_Net,tempMatrix)
}
figMatrix_Com <- data.frame()
for(i in seq(1,nrow(res_Com$b))){
  tempRowName <- row.names(res_Com$b)[i]
  selectID <- gregexpr('SRG',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 3),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Group = 'Com',Mean = as.numeric(res_Com$b[i]),
                           UpLimit = as.numeric(res_Com$ci.ub[i]),LowLimit = as.numeric(res_Com$ci.lb[i]))
  figMatrix_Com <- rbind(figMatrix_Com,tempMatrix)
}
figMatrix_Sel <- data.frame()
for(i in seq(1,nrow(res_Sel$b))){
  tempRowName <- row.names(res_Sel$b)[i]
  selectID <- gregexpr('SRG',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 3),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Group = 'Sel',Mean = as.numeric(res_Sel$b[i]),
                           UpLimit = as.numeric(res_Sel$ci.ub[i]),LowLimit = as.numeric(res_Sel$ci.lb[i]))
  figMatrix_Sel <- rbind(figMatrix_Sel,tempMatrix)
}
figMatrix <- rbind(figMatrix_Net,figMatrix_Com,figMatrix_Sel)
# Points
figMatrix$XLocation <- 0
selectIDs <- which(figMatrix$Class == '1')
figMatrix$XLocation[selectIDs] <- 3
selectIDs <- which(figMatrix$Class == '2')
figMatrix$XLocation[selectIDs] <- 2
selectIDs <- which(figMatrix$Class == '3')
figMatrix$XLocation[selectIDs] <- 1
# Polygons
width <- 0.18
polygons <- data.frame()
for(i in seq(1,nrow(figMatrix))){
  point_1 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$LowLimit[i])
  point_2 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$UpLimit[i])
  point_3 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$UpLimit[i])
  point_4 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$LowLimit[i])
  tempPolygon <- rbind(point_1,point_2,point_3,point_4,point_1)
  tempPolygon$ID <- i
  tempPolygon$Group <- figMatrix$Group[i]
  polygons = rbind(polygons,tempPolygon)
}
# Dodge
dodge <- 0.3
selectIDs_1 <- which(figMatrix$Group == 'Net')
selectIDs_2 <- which(figMatrix$Group == 'Sel')
figMatrix$XLocation[selectIDs_1] <- figMatrix$XLocation[selectIDs_1] + dodge
figMatrix$XLocation[selectIDs_2] <- figMatrix$XLocation[selectIDs_2] - dodge
figMatrix$Group <- factor(figMatrix$Group,levels = c('Net','Com','Sel'),ordered = TRUE)
selectIDs_1 <- which(polygons$Group == 'Net')
selectIDs_2 <- which(polygons$Group == 'Sel')
polygons$X[selectIDs_1] <- polygons$X[selectIDs_1] + dodge
polygons$X[selectIDs_2] <- polygons$X[selectIDs_2] - dodge
polygons$Group <- factor(polygons$Group,levels = c('Net','Com','Sel'),ordered = TRUE)
# Part 4.3: Draw
fig_S11d <- ggplot()+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID,fill = Group),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean,group = Group,color = Group),size = 0.7)+
  scale_color_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.20,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_fill_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.20,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_x_continuous(limits = c(0.4,3.6),breaks = c(1,2,3),labels = c('> 4','3 or 4','2'),expand = c(0,0))+
  scale_y_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),labels = c('-0.6','-0.3','0','0.3','0.6'),expand = c(0,0))+
  coord_flip()+
  theme_custom()
##==== Part 4: Height (Community) ====##

##==== Part 5: DBH (Community) ====##
load('Data_DBH_Community.R')
# Part 5.1: Analysis
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
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))
# Net
res_null_Net <- rma.mv(Net_yi,Net_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Net <- rma.mv(Net_yi,Net_vi,mods = ~SRG-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
modelCompare_Net <- anova.rma(res_null_Net,res_Net)
multiCompare_Net <- summary(glht(res_Net,linfct = rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))),test = adjusted('none'))
# Com
res_null_Com <- rma.mv(Com_yi,Com_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Com <- rma.mv(Com_yi,Com_vi,mods = ~SRG-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
modelCompare_Com <- anova.rma(res_null_Com,res_Com)
multiCompare_Com <- summary(glht(res_Com,linfct = rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))),test = adjusted('none'))
# Sel
res_null_Sel <- rma.mv(Sel_yi,Sel_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~SRG-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
modelCompare_Sel <- anova.rma(res_null_Sel,res_Sel)
multiCompare_Sel <- summary(glht(res_Sel,linfct = rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))),test = adjusted('none'))
# Part 5.2: Prepare draw data
figMatrix_Net <- data.frame()
for(i in seq(1,nrow(res_Net$b))){
  tempRowName <- row.names(res_Net$b)[i]
  selectID <- gregexpr('SRG',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 3),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Group = 'Net',Mean = as.numeric(res_Net$b[i]),
                           UpLimit = as.numeric(res_Net$ci.ub[i]),LowLimit = as.numeric(res_Net$ci.lb[i]))
  figMatrix_Net <- rbind(figMatrix_Net,tempMatrix)
}
figMatrix_Com <- data.frame()
for(i in seq(1,nrow(res_Com$b))){
  tempRowName <- row.names(res_Com$b)[i]
  selectID <- gregexpr('SRG',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 3),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Group = 'Com',Mean = as.numeric(res_Com$b[i]),
                           UpLimit = as.numeric(res_Com$ci.ub[i]),LowLimit = as.numeric(res_Com$ci.lb[i]))
  figMatrix_Com <- rbind(figMatrix_Com,tempMatrix)
}
figMatrix_Sel <- data.frame()
for(i in seq(1,nrow(res_Sel$b))){
  tempRowName <- row.names(res_Sel$b)[i]
  selectID <- gregexpr('SRG',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 3),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Group = 'Sel',Mean = as.numeric(res_Sel$b[i]),
                           UpLimit = as.numeric(res_Sel$ci.ub[i]),LowLimit = as.numeric(res_Sel$ci.lb[i]))
  figMatrix_Sel <- rbind(figMatrix_Sel,tempMatrix)
}
figMatrix <- rbind(figMatrix_Net,figMatrix_Com,figMatrix_Sel)
# Points
figMatrix$XLocation <- 0
selectIDs <- which(figMatrix$Class == '1')
figMatrix$XLocation[selectIDs] <- 3
selectIDs <- which(figMatrix$Class == '2')
figMatrix$XLocation[selectIDs] <- 2
selectIDs <- which(figMatrix$Class == '3')
figMatrix$XLocation[selectIDs] <- 1
# Polygons
width <- 0.18
polygons <- data.frame()
for(i in seq(1,nrow(figMatrix))){
  point_1 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$LowLimit[i])
  point_2 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$UpLimit[i])
  point_3 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$UpLimit[i])
  point_4 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$LowLimit[i])
  tempPolygon <- rbind(point_1,point_2,point_3,point_4,point_1)
  tempPolygon$ID <- i
  tempPolygon$Group <- figMatrix$Group[i]
  polygons = rbind(polygons,tempPolygon)
}
# Dodge
dodge <- 0.3
selectIDs_1 <- which(figMatrix$Group == 'Net')
selectIDs_2 <- which(figMatrix$Group == 'Sel')
figMatrix$XLocation[selectIDs_1] <- figMatrix$XLocation[selectIDs_1] + dodge
figMatrix$XLocation[selectIDs_2] <- figMatrix$XLocation[selectIDs_2] - dodge
figMatrix$Group <- factor(figMatrix$Group,levels = c('Net','Com','Sel'),ordered = TRUE)
selectIDs_1 <- which(polygons$Group == 'Net')
selectIDs_2 <- which(polygons$Group == 'Sel')
polygons$X[selectIDs_1] <- polygons$X[selectIDs_1] + dodge
polygons$X[selectIDs_2] <- polygons$X[selectIDs_2] - dodge
polygons$Group <- factor(polygons$Group,levels = c('Net','Com','Sel'),ordered = TRUE)
# Part 5.3: Draw
fig_S11e <- ggplot()+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID,fill = Group),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean,group = Group,color = Group),size = 0.7)+
  scale_color_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.20,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_fill_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.20,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_x_continuous(limits = c(0.4,3.6),breaks = c(1,2,3),labels = c('> 4','3 or 4','2'),expand = c(0,0))+
  scale_y_continuous(limits = c(-0.6,0.6),breaks = seq(-0.6,0.6,0.3),labels = c('-0.6','-0.3','0','0.3','0.6'),expand = c(0,0))+
  coord_flip()+
  theme_custom()
##==== Part 5: DBH (Community) ====##

##==== Part 6: Biomass (Community) ====##
load('Data_Bio_Community.R')
# Part 6.1: Analysis
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
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))
# Net
res_null_Net <- rma.mv(Net_yi,Net_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Net <- rma.mv(Net_yi,Net_vi,mods = ~SRG-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
modelCompare_Net <- anova.rma(res_null_Net,res_Net)
multiCompare_Net <- summary(glht(res_Net,linfct = rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))),test = adjusted('none'))
# Com
res_null_Com <- rma.mv(Com_yi,Com_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Com <- rma.mv(Com_yi,Com_vi,mods = ~SRG-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
modelCompare_Com <- anova.rma(res_null_Com,res_Com)
multiCompare_Com <- summary(glht(res_Com,linfct = rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))),test = adjusted('none'))
# Sel
res_null_Sel <- rma.mv(Sel_yi,Sel_vi,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
res_Sel <- rma.mv(Sel_yi,Sel_vi,mods = ~SRG-1,random = ~1|Level4/Level3/Level2,method = 'ML',data = dataMatrix)
modelCompare_Sel <- anova.rma(res_null_Sel,res_Sel)
multiCompare_Sel <- summary(glht(res_Sel,linfct = rbind(c(-1,1,0),c(-1,0,1),c(0,-1,1))),test = adjusted('none'))
# Part 6.2: Prepare draw data
figMatrix_Net <- data.frame()
for(i in seq(1,nrow(res_Net$b))){
  tempRowName <- row.names(res_Net$b)[i]
  selectID <- gregexpr('SRG',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 3),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Group = 'Net',Mean = as.numeric(res_Net$b[i]),
                           UpLimit = as.numeric(res_Net$ci.ub[i]),LowLimit = as.numeric(res_Net$ci.lb[i]))
  figMatrix_Net <- rbind(figMatrix_Net,tempMatrix)
}
figMatrix_Com <- data.frame()
for(i in seq(1,nrow(res_Com$b))){
  tempRowName <- row.names(res_Com$b)[i]
  selectID <- gregexpr('SRG',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 3),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Group = 'Com',Mean = as.numeric(res_Com$b[i]),
                           UpLimit = as.numeric(res_Com$ci.ub[i]),LowLimit = as.numeric(res_Com$ci.lb[i]))
  figMatrix_Com <- rbind(figMatrix_Com,tempMatrix)
}
figMatrix_Sel <- data.frame()
for(i in seq(1,nrow(res_Sel$b))){
  tempRowName <- row.names(res_Sel$b)[i]
  selectID <- gregexpr('SRG',tempRowName)[[1]]
  tempClass <- substring(tempRowName,(selectID + 3),nchar(tempRowName))
  tempMatrix <- data.frame(Class = tempClass,Group = 'Sel',Mean = as.numeric(res_Sel$b[i]),
                           UpLimit = as.numeric(res_Sel$ci.ub[i]),LowLimit = as.numeric(res_Sel$ci.lb[i]))
  figMatrix_Sel <- rbind(figMatrix_Sel,tempMatrix)
}
figMatrix <- rbind(figMatrix_Net,figMatrix_Com,figMatrix_Sel)
# Points
figMatrix$XLocation <- 0
selectIDs <- which(figMatrix$Class == '1')
figMatrix$XLocation[selectIDs] <- 3
selectIDs <- which(figMatrix$Class == '2')
figMatrix$XLocation[selectIDs] <- 2
selectIDs <- which(figMatrix$Class == '3')
figMatrix$XLocation[selectIDs] <- 1
# Polygons
width <- 0.18
polygons <- data.frame()
for(i in seq(1,nrow(figMatrix))){
  point_1 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$LowLimit[i])
  point_2 <- data.frame(X = figMatrix$XLocation[i] - width/2,Y = figMatrix$UpLimit[i])
  point_3 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$UpLimit[i])
  point_4 <- data.frame(X = figMatrix$XLocation[i] + width/2,Y = figMatrix$LowLimit[i])
  tempPolygon <- rbind(point_1,point_2,point_3,point_4,point_1)
  tempPolygon$ID <- i
  tempPolygon$Group <- figMatrix$Group[i]
  polygons = rbind(polygons,tempPolygon)
}
# Dodge
dodge <- 0.3
selectIDs_1 <- which(figMatrix$Group == 'Net')
selectIDs_2 <- which(figMatrix$Group == 'Sel')
figMatrix$XLocation[selectIDs_1] <- figMatrix$XLocation[selectIDs_1] + dodge
figMatrix$XLocation[selectIDs_2] <- figMatrix$XLocation[selectIDs_2] - dodge
figMatrix$Group <- factor(figMatrix$Group,levels = c('Net','Com','Sel'),ordered = TRUE)
selectIDs_1 <- which(polygons$Group == 'Net')
selectIDs_2 <- which(polygons$Group == 'Sel')
polygons$X[selectIDs_1] <- polygons$X[selectIDs_1] + dodge
polygons$X[selectIDs_2] <- polygons$X[selectIDs_2] - dodge
polygons$Group <- factor(polygons$Group,levels = c('Net','Com','Sel'),ordered = TRUE)
# Part 6.3: Draw
fig_S11f <- ggplot()+
  geom_hline(yintercept = 0,size = 0.3,color = '#000000',linetype = 'dashed')+
  geom_polygon(data = polygons,mapping = aes(x = X,y = Y,group = ID,fill = Group),alpha = 0.4)+
  geom_point(data = figMatrix,mapping = aes(x = XLocation,y = Mean,group = Group,color = Group),size = 0.7)+
  scale_color_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.20,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_fill_manual(values = c(CMYKtoRGB(0.84,0.67,0,0),CMYKtoRGB(0.82,0.20,0.93,0.24),CMYKtoRGB(0,0.38,0.96,0)))+
  scale_x_continuous(limits = c(0.4,3.6),breaks = c(1,2,3),labels = c('> 4','3 or 4','2'),expand = c(0,0))+
  scale_y_continuous(limits = c(-1.6,1.6),breaks = seq(-1.6,1.6,0.8),labels = c('-1.6','-0.8','0','0.8','1.6'),expand = c(0,0))+
  coord_flip()+
  theme_custom()
##==== Part 6: Biomass (Community) ====##

##==== Combine figures ====##
fig_S11 <- ggdraw()+
  draw_plot(fig_S11a,x = 0,y = 0.63,width = 0.333,height = 0.37)+
  draw_plot(fig_S11b,x = 0.334,y = 0.63,width = 0.333,height = 0.37)+
  draw_plot(fig_S11c,x = 0.667,y = 0.63,width = 0.333,height = 0.37)+
  draw_plot(fig_S11d,x = 0,y = 0,width = 0.333,height = 0.63)+
  draw_plot(fig_S11e,x = 0.334,y = 0,width = 0.333,height = 0.63)+
  draw_plot(fig_S11f,x = 0.667,y = 0,width = 0.333,height = 0.63)
outputFileName <- paste0(currentPath,'/fig_S11.pdf')
pdf(file = outputFileName,width = 5.07,height = 3.3,colormodel = 'cmyk')
print(fig_S11)
dev.off()
##==== Combine figures ====##
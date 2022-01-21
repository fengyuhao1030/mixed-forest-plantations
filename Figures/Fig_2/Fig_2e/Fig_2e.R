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

load('DBHData.R')
dataMatrix <- get('DBHData')
dataMatrix$Level4 <- as.factor(paste(dataMatrix$Lat,dataMatrix$Lon))
dataMatrix$Level3 <- as.factor(paste(dataMatrix$Composition,dataMatrix$Age,dataMatrix$Density))
dataMatrix$Level2 <- as.factor(seq(1,nrow(dataMatrix)))

##==== Part 1: Needleleaf ====##
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
##==== Part 1: Needleleaf ====##

##==== Part 2: Evergreen ====##
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
##==== Part 2: Evergreen ====##

##==== Part 3: Nitrogen ====##
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
##==== Part 3: Nitrogen ====##

##==== Part 4: Draw ====##
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
  ylab('Effect size')+
  coord_flip()+
  theme_custom()
outputFileName <- paste0(currentPath,'/Fig_2e.pdf')
pdf(file = outputFileName,width = 1.79,height = 3,colormodel = 'cmyk')
print(Fig_2e)
dev.off()
##==== Part 4: Draw ====##
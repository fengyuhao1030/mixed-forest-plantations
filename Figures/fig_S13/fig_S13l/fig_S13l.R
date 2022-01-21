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

load('BioData.R')
dataMatrix <- get('BioData')

##==== Part 1: Prepare data ====##
dataMatrix <- dataMatrix[,c('Sel_yi','Sel_vi')]
dataMatrix <- dataMatrix[order(dataMatrix$Sel_yi),]
dataMatrix$Sel_vi <- sqrt(dataMatrix$Sel_vi)
selectIDs <- which((dataMatrix$Sel_yi <= quantile(dataMatrix$Sel_yi,0.98))&(dataMatrix$Sel_yi >= quantile(dataMatrix$Sel_yi,0.02))&(dataMatrix$Sel_vi < 0.4))
figData <- data.frame(Sample = seq(1,length(selectIDs)),Mean = dataMatrix$Sel_yi[selectIDs],SD = dataMatrix$Sel_vi[selectIDs])
extraFigData <- data.frame()
for(i in seq(1,nrow(figData))){
  if((figData[i,'Mean']+qnorm(0.975)*figData[i,'SD']) < 0){
    tempData <- figData[i,]
    tempData$Color <- CMYKtoRGB(0.1,0.91,0.91,0.35)
    extraFigData <- rbind(extraFigData,tempData)
  }
  if((figData[i,'Mean']-qnorm(0.975)*figData[i,'SD']) > 0){
    tempData <- figData[i,]
    tempData$Color <- CMYKtoRGB(0.84,0.67,0,0)
    extraFigData <- rbind(extraFigData,tempData)
  }
}
##==== Part 1: Prepare data ====##

##==== Part 2: Draw ====##
# Draw
fig_S13l <- ggplot()+
  geom_errorbar(data = figData,mapping = aes(x = Sample,ymin = Mean - qnorm(0.975)*SD,ymax = Mean + qnorm(0.975)*SD),
                color = '#C8C8C8',size = 0.1,width = 0)+
  geom_errorbar(data = extraFigData,mapping = aes(x = Sample,ymin = Mean - qnorm(0.975)*SD,ymax = Mean + qnorm(0.975)*SD),
                color = extraFigData$Color,size = 0.1,width = 0)+
  geom_point(data = extraFigData,mapping = aes(x = Sample,y = Mean),color = extraFigData$Color,size = 0.1)+
  geom_hline(yintercept = 0,linetype = 'dashed',color = '#000000',size = 0.3)+
  # axis
  scale_x_continuous(limits = c(1,nrow(figData)),breaks = c(1,nrow(figData)),labels = c('1',as.character(nrow(dataMatrix))),expand = c(0.03,0.03))+
  scale_y_continuous(limits = c(-1.6,1.6),breaks = seq(-1.6,1.6,0.8),labels = c('-1.6','-0.8','0','0.8','1.6'),expand = c(0,0))+
  theme_custom()
outputFileName <- paste0(currentPath,'/fig_S13l.pdf')
pdf(file = outputFileName,width = 1.55,height = 1.46,colormodel = 'cmyk')
print(fig_S13l)
dev.off()
# Insert
fig_S13l_Insert <- ggplot()+
  geom_histogram(data = figData,mapping = aes(x = Mean),binwidth = 0.02,fill = CMYKtoRGB(0,0.38,0.96,0))+
  geom_vline(xintercept = 0,size = 0.2,color = '#000000',linetype = 'dashed')+
  scale_x_continuous(limits = c(-0.4,0.4),breaks = seq(-0.3,0.3,0.3),labels = c('-0.3','0','0.3'),expand = c(0,0))+
  scale_y_continuous(limits = c(0,70),breaks = seq(0,60,30))+
  theme_custom()+
  theme(plot.margin = margin(2,2,2,2),
        axis.ticks.length = unit(-0.1,'lines'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6,margin = margin(2,0,0,0),color = '#000000'),
        axis.text.y = element_text(size = 6,margin = margin(0,2,0,0),color = '#000000'))
outputFileName <- paste0(currentPath,'/fig_S13l_Insert.pdf')
pdf(file = outputFileName,width = 0.56,height = 0.52,colormodel = 'cmyk')
print(fig_S13l_Insert)
dev.off()
##==== Part 2: Draw ====##
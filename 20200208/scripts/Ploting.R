##Plot todtomato cells 

rm(list = ls())
library(ggplot2)



df <- read.csv("results/Final_rotated_moved_coord.csv")


p <- ggplot(df)+theme_gray()+
  geom_segment(aes(x = Rotated_Moved_Nuclei_X,y = Rotated_Moved_Nuclei_Y,xend = Rotated_Moved_GM130_X ,yend = Rotated_Moved_GM130_Y,color=Tissue),arrow=arrow(length = unit(0.02, "npc")))+
  coord_equal(ratio=1)+ #make sure the x, y axis scale the same
  #scale_x_continuous(name='',breaks = seq(-50,50,by =20))+ # modify the x axis breaks and range
  #scale_y_continuous(name='',breaks=seq(-100, 150, by =20))+
  #change the color
  #scale_color_manual(name='Tip Types',values = c('Parent tip'='green', 'Right tip'='black', 'Left tip'='blue', 'Cleft'='red'))#+
  #facet_wrap(~sample, nrow=2)
  xlab('')+
  ylab('')+
  xlim(-100,100)+
  ylim(-100,100)+
  #theme(legend.position="bottom", legend.box = "horizontal")+
  ggtitle('E11.5 Mammary Area Cells Oritation')+
  geom_segment(x=40, y=75, xend=80, yend=75,color='grey34',size=0.5, arrow=arrow(length = unit(0.02, "npc")))+
  geom_segment(x=60, y=60, xend=60, yend=90,color='grey34',size=0.5, arrow=arrow(length = unit(0.02, "npc")))+
  geom_point(aes(x=0, y=0), color='black', inherit.aes = FALSE,shape=8)+
  annotate('text', x=85, y=75, label='A')+
  annotate('text', x=35, y=75, label='P')+
  annotate('text', x=60, y=55, label='V')+
  annotate('text', x=60, y=96, label='D')+
  facet_wrap(~Sample, nrow=2)
p

vector_df <- df %>% 
  mutate(angle_X_axis=atan2((Rotated_Moved_GM130_Y-Rotated_Moved_Nuclei_Y), (Rotated_Moved_GM130_X-Rotated_Moved_Nuclei_X)))

rose.diag(circular(abs(vector_df$angle_X_axis)), units = 'degrees', prop=5, main = 'Tdtomato E11.5 Mammary Area orintation_abs')


library(RColorBrewer)
my.color <- brewer.pal(nlevels(vector_df$Tissue), name = 'Set2')

# names(my.color) <- levels(vector_df$Tissue) 
# 
# color_df <- my.color%>%as.data.frame() %>%  rownames_to_column(var='Tissue')
# names(color_df)[2] <- 'Color'
# vector_df <- vector_df %>% merge(color_df, by='Tissue')


for(i in c(1:3)){
  #i=2
  current_tissue <- levels(vector_df$Tissue)[i]
  temp <- vector_df %>% dplyr::filter(Tissue==current_tissue)
  if (i==1){
    rose.diag(circular(temp$angle_X_axis), 
              units = 'degrees', 
              prop=3, 
              main = 'Tdtomato E11.5 Mammary Area orintation_abs', 
              col=my.color[i])
  }else{
    rose.diag(circular(temp$angle_X_axis), 
              units = 'degrees', 
              prop=3, 
              #main = 'Tdtomato E11.5 Mammary Area orintation_abs', 
              col=my.color[i], add=T)
  }
}


for(i in c(1:3)){
  #i=2
  current_tissue <- levels(vector_df$Tissue)[i]
  temp <- vector_df %>% dplyr::filter(Tissue==current_tissue)
    rose.diag(circular(abs(temp$angle_X_axis)), 
              units = 'degrees', 
             prop=1.6, 
              bins=12,
              main = paste0('E11.5 Mammary Area orintation ',current_tissue), 
              col=my.color[i])
}


##########

library(spatstat)

for(i in c(1:3)){
  #i=2
  current_tissue <- levels(vector_df$Tissue)[i]
  temp <- vector_df %>% dplyr::filter(Tissue==current_tissue)
  rose(abs(temp$angle_X_axis), 
       unit = 'radian', 
            # prop=3, 
           # bins=12,
            main = paste0('E11.5 Mammary Area orintation ',current_tissue), 
            col=my.color[i])
}

rose(vector_df$angle_X_axis, unit = 'radian', labels=F)


arrows.1 <- function(x0, y0, length.ar, angle.ar, ...){
  
  ab <- cos(angle.ar) * length.ar
  bc <- sign(sin(angle.ar)) * sqrt(length.ar^2 - ab^2)
  
  x1 <- x0 + ab
  y1 <- y0 + bc
  
  arrows(x0, y0, x1, y1, ...)
}

## Circle function
circle <- function(xorig, yorig, radius, add, ...){
  
  x <- seq(-radius, radius, length.out = 1000)
  
  y <- sapply(x, function(z) sqrt(radius^2 - z^2))
  
  if(add == TRUE){
    
    lines(xorig + c(x, rev(x)), c(yorig + y, yorig + rev(-y)), type = "l", ...)
    
  } else {
    
    plot(xorig + c(x, rev(x)), c(yorig + y, yorig + rev(-y)), type = "l", ...)
    
  }
  
}

circle(1, 1, 1, add = FALSE)
for(i in 1:nrow(vector_df)){
  arrows.1(x0 = 1, y0 = 1, length.ar = 1, angle.ar = abs(vector_df$angle[i]))
}




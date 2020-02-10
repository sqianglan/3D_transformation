##data ploting
rm(list=ls())
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)

df_orignal_all <- read_csv('/Users/qianglan/Projects/3D_transformation/RawData/Cell polarity in R from Ewe/Tdtomato_moved_df_all.csv') %>% as.data.frame()

p <- ggplot(df_orignal_all)+
  theme_gray()+
  #theme(axis.title.x = element_blank(),
  #      axis.title.y = element_blank())+
  labs(x='µm', y ='µm')+
  geom_segment(aes(x = Nuclei_X,y =Nuclei_Y,xend =GM130_X ,yend =GM130_Y,color=Sample, group=Tissue),
               arrow=arrow(length = unit(0.03, "npc")), show.legend = F)+
  coord_equal(ratio=1)+
  ggtitle('Project in 2D')+
  facet_wrap(Sample ~Tissue)
p

# read thea .csv file

## chang the file location if necessary
folder <- "/Users/qianglan/Projects/3D_transformation/RawData/Cell polarity in R from Ewe/"

#fullpath <- paste(current_pwd,folder,sep = '/')

#### chang the file name pattern if necessary
pattern <- "*thea.csv" # regular expression pattern for the files containing "Green"


file_names <-  list.files(folder, pattern=pattern, recursive = T)
setwd(folder) # set the working dir into the folder place to facilitate the read_csv, otherwise the function get some problem.

##combine the data
for (i in 1:length(file_names)){
  name <- strsplit(file_names[i], '\\/|_|\\.') %>% unlist()
  tissue <- name[2]
  stage <- paste(name[4], name[5],sep = '.')
  orignal_file_name <- str_replace(file_names[i],'_thea', '')
  thea <- read_csv(file_names[i]) %>% t(.)
  colnames(thea)[1] <- 'thea_3D'
  
  orignal_file <- read_csv(orignal_file_name) # read the orignal file 
  # merged the orignal file with the thea data to matche the sample and tissue
  temp <- cbind(orignal_file, thea) %>%
    mutate(stage=stage) %>% 
    mutate(tissue= case_when(
      tissue =='Dermal' ~'de',
      tissue =='MB' ~'mb3'))
  if (i==1){
    df <- temp
  } else{df <- rbind(df,temp)}
}
df2 <-df %>%  mutate(angle=thea_3D*180/pi) %>% 
  mutate(point_to_center=ifelse(angle > 90, 1, 0))

df3 <- df2 %>%  group_by(Sample, tissue, stage, point_to_center)%>% summarise(Count_number=length(point_to_center))
total_number <- df2 %>% group_by(Sample, tissue, stage)%>% summarise(total_number=length(point_to_center))
Percent_sample <- left_join(df3, total_number, by=c('Sample', 'tissue', 'stage')) %>% mutate(perc = 100*Count_number/total_number )

Percent_sample_test <- Percent_sample %>% group_by(tissue, stage)
# t.test




q <- ggplot(df2, aes(x=point_to_center))+
  theme_gray()+
  geom_bar(aes(y=..prop.., fill = factor(..x..)), stat = 'count', position = 'dodge')+
  facet_grid(~tissue)
q

t <- ggplot(Percent_sample, aes(x=as.factor(point_to_center), y=perc, color=Sample))+
  theme_gray()+
  geom_point()+
  facet_grid(tissue~stage)
t


s <- ggboxplot(ungroup(Percent_sample), 
               x='point_to_center', 
               y='perc',
               color='point_to_center',
               add = "jitter")+
  theme_gray()+
  #ggboxplot()+
  ylab('Percentage (%)')+
  xlab('Relative cell orintation towards center ')+
  scale_x_discrete(labels=c("0" = " <90°", "1" = ">=90°"))+
  ggtitle('Cell_orintation_in_3D')+
  facet_grid(~tissue)
s <- s+stat_compare_means(method = "t.test", label.y = 90, label.x=1.2)+rremove('legend')



grid.arrange(grobs=list(p,s),
             layout_matrix=rbind(c(1,1),
                                 c(1,1),
                                 c(1,1),
                                 c(2,2)))




# try to use plotly to plot the 3D arrow

# temp <- df %>% filter(stage=='E11.25'& tissue=='mb3') 
# library(plot3D)
# 
# arrows3D (temp$Nuclei_X, temp$Nuclei_Y, temp$Nuclei_Z, x1 = temp$GM130_X, y1 = temp$GM130_Y, z1 = temp$GM130_Z,  
#           colvar = temp$Sample, phi = 0, theta = 50,
#           col = NULL, NAcol = "white", breaks = NULL,
#           colkey = NULL, panel.first = NULL,
#           clim = NULL, clab = NULL, bty = "g", type = "triangle", 
#           add = FALSE, plot = TRUE)
# scatter3D(0,0,0,add=T, col = 'black', )
# 
# arrows2D(temp$Nuclei_X, temp$Nuclei_Y, x1 = temp$GM130_X, y1 = temp$GM130_Y,
#          colvar = temp$ID, bty ="n", xlim=c(-100,100), ylim=c(-100,100))
# points2D(0, 0, add = TRUE, col="black", 
#          colkey = FALSE, pch = 19, cex = 1)

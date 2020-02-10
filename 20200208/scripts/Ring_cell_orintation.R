
rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggplot2)
library(readxl)
library(gtools)
library(spdep) # for cooridnate rotation



current_pwd <- getwd() # save the current working direct path to rest it back later

## chang the file location if necessary
folder <- "/Users/qianglan/Projects/3D_transformation/RawData/Cell polarity in R from Ewe/Ring cells orientation E12.5"

#fullpath <- paste(current_pwd,folder,sep = '/')

#### chang the file name pattern if necessary
pattern <- "*.csv" # regular expression pattern for the files containing "Green"


file_names <-  list.files(folder, pattern=pattern, recursive = T)
setwd(folder) # set the working dir into the folder place to facilitate the read_csv, otherwise the function get some problem.


for (i in 1:length(file_names)){
  #read each sheet and process seperately and then combine together.
  name <- strsplit(file_names[i], '\\/|_|\\.') %>% unlist() # split the string by '/', or '_' or '.
  sample_name <- paste(name[3], name[4], sep = '_') #extract the sample name
  data_label <- name[length(name)-1]# extact the type of data (center, axix or start point or end point)
  # better way to get type
  type <- regmatches(file_names[i], gregexpr('GM130|Nuclei|center|APDV',file_names[i], ignore.case = T)) %>% unlist()
  if (grepl("Anterior|Dorsal|Posterior|Ventral", data_label, ignore.case=T)) next #those information already in the ADPCV files, skip
  extension <- name[length(name)]
  
  if (type=='APDV') {
    temp <- read_csv(file_names[i]) %>% as.data.frame(.)
    names(temp)[1] <- "Direction"
    temp <- gather(temp, "X","Y","Z", key="axis", value="coord") %>% 
      unite(Direction, axis, col = "Direc_axis", sep = "_") %>% as.data.frame(.)
    rownames(temp) <- temp[,1]
    temp <- temp %>% select(-1) %>% t(.)
    rownames(temp) <- sample_name
    data_type =type
    
  } else if (type=='center'){
    temp <- read_csv(file_names[i], skip=1) %>% as.data.frame(.)
    center_X <- temp$Mean[which(temp$Variable=='Position X')]
    center_Y <- temp$Mean[which(temp$Variable=='Position Y')]
    center_Z <- temp$Mean[which(temp$Variable=='Position Z')]
    
    temp <- data.frame('center_X'=center_X,
                                  'center_Y'= center_Y,
                                  'center_Z'= center_Z)
    rownames(temp) <- sample_name
    data_type =type
  } else if (type %in% c('GM130','Nuclei')){
    temp <-read_csv(file_names[i], skip=1) %>% as.data.frame(.)%>%  # read the files and tranform into dataframe
      select(Variable, Min,'Surpass Object') %>% #selected interested data
      filter(Variable %in% c('Position X', 'Position Y', 'Position Z')) %>% #fileter the coordinate data
      separate(Variable, into=c('Variable','Coord')) %>% # split the name
      separate('Surpass Object', into=c('Marker','1','2','ID')) %>%  # split the object ID.
      select(Coord, ID,Min) %>% 
      mutate(Marker=type)%>% #label the data with data type (GM130 or nuclear)
      unite('Marker',Marker,Coord, sep = '_')%>% 
      spread(Marker, Min) %>%
      mutate(sample=sample_name) 
    data_type ="spots"
  }
  ## Merging the data into thw same files
  
  ##check if the variables exist or not
if (exists(data_type)){ 
  if (data_type =='APDV'){ 
      APDV <- rbind(APDV, temp)  # rbind for axis cooridinate data directly
    } 
    else if (data_type == "center"){center <- rbind(center, temp)}
    else if(type %in% c('GM130', 'Nuclei')){ # if the data is the GM130 dataset,
      if (exists(sample_name)){
        df <- merge(eval(as.name(sample_name)), temp, by =c('ID', 'sample'))
        spots <- rbind(spots,df)
      } else{assign(sample_name, temp)}
    } 
  }else if(type %in% c('GM130', 'Nuclei') & !exists(sample_name)){
    assign(sample_name, temp)
    }else if(type %in% c('GM130', 'Nuclei') & exists(sample_name)){
    spots <- merge(eval(as.name(sample_name)), temp, by =c('ID', 'sample'))
  }else{assign(data_type, temp)}
}
  
  # save the temperte files. A lazy way.
 # write.csv(temp, paste(sample_name,marker,'csv', sep = '.'))  

## get the center coordinate.
center_point <- center %>% rownames_to_column(var='sample')
# 
#write.csv(center, 'Ring_cells_center_coordinate.csv')
# #cacultet the of angel of parent tip to the x axis 
# write.csv(spots, 'Ring_cell_direction_coord.csv')
# write.csv(APDV, 'Ring_cell_APDV_axis_coord.csv')

thea <- APDV  %>% as.data.frame(.) %>% rownames_to_column(var='sample') %>% select(1:9) %>% 
  left_join(center_point, by='sample')
rotated_thea <- thea
# for (i in c(2:5)){
#   rotated_thea[,i] <- thea[,i]-thea[["center_X"]]
#   rotated_thea[,i+4] <- thea[,i+4]-thea[["center_Y"]]
# }

rotated_thea <- rotated_thea %>% select(1:9) %>% 
  mutate(Length=sqrt((Anterior_X-Posterior_X)^2+ (Anterior_Y-Posterior_Y)^2)) %>% 
  mutate(thea=acos((Anterior_X-Posterior_X)/Length)) %>% mutate('rotation_angel' =ifelse(Anterior_Y>=Posterior_Y, 2*pi-thea, thea))
##rotation angle done!!

for (i in 1:length(rownames(rotated_thea))){
  rotated_thea[i, c('Anterior_X', 'Anterior_Y')] <-as.data.frame(Rotation(rotated_thea[i, c('Anterior_X', 'Anterior_Y')], rotated_thea$rotation_angel[i]))
  rotated_thea[i, c('Posterior_X', 'Posterior_Y')] <- as.data.frame(Rotation(rotated_thea[i, c('Posterior_X', 'Posterior_Y')], rotated_thea$rotation_angel[i]))
  rotated_thea[i, c('Dorsal_X', 'Dorsal_Y')] <- as.data.frame(Rotation(rotated_thea[i, c('Dorsal_X', 'Dorsal_Y')], rotated_thea$rotation_angel[i]))
  rotated_thea[i, c('Ventral_X', 'Ventral_Y')] <- as.data.frame(Rotation(rotated_thea[i, c('Ventral_X', 'Ventral_Y')], rotated_thea$rotation_angel[i]))
  #temp <- cbind(temp_Anterior, temp_Posterior,temp_Dorsal, temp_Ventral )
  # if (i==1){
  #   rotated_AP <- temp
  # }else {rotated_AP <- rbind(rotated_AP,temp)}
}
rotated_AP <- rotated_thea



library(spdep)




  
df <- left_join(spots,center_point, by ="sample")

#math sbustrate GM130 with center-cord
n <-  c(3:5) 
new_df <- mapply(`-`, df[n], df[n+6])
colnames(new_df) <- paste0("Moved_", colnames(new_df))
df2 <- cbind(df, new_df)
# mathe substrate Nuclei with center_cord
n <-  c(6:8) 
new_df2 <- mapply(`-`, df[n], df[n+3])
colnames(new_df2) <- paste0("Moved_", colnames(new_df2))
df2 <- cbind(df2,new_df2)

coord_df <- df2 %>% select(c('ID','sample',c(12:17)))


rotation_df <- left_join(coord_df, rotated_thea[,c('sample','rotation_angel')], by='sample') %>% select(-c(5,8))


for (i in 1:length(rownames(rotation_df))){
  temp_GM130 <-as.data.frame(Rotation(rotation_df[i, c(3,4)], rotation_df$rotation_angel[i]))
  temp_Nuclei <- as.data.frame(Rotation(rotation_df[i, c(5,6)], rotation_df$rotation_angel[i]))
  temp <- cbind(temp_GM130, temp_Nuclei)
  if (i==1){
    rotated_df <- temp
  }else {rotated_df <- rbind(rotated_df,temp)}
}
names(rotated_df) <- paste('Rotated', names(rotation_df)[3:6], sep = '_')

rotated_df <- cbind(rotation_df[,c(1,2)], rotated_df)

#write.csv(rotated_df, "Final_rotated_moved_coord.csv" )
# rotated_df <- read.csv('/Users/qianglan/Projects/3D_transformation/RawData/Cell polarity in R from Ewe/Ring cells orientation E12.5/Final_rotated_moved_coord.csv') %>%
#   select(-1)

p <- ggplot(rotated_df)+theme_gray()+
  geom_segment(aes(x = Rotated_Moved_Nuclei_X,y = Rotated_Moved_Nuclei_Y,xend = Rotated_Moved_GM130_X ,yend = Rotated_Moved_GM130_Y,color=sample),arrow=arrow(length = unit(0.02, "npc")))+
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
  ggtitle('Ring Cell Oritation')+
  geom_segment(x=40, y=75, xend=80, yend=75,color='grey34',size=0.5, arrow=arrow(length = unit(0.02, "npc")))+
  geom_segment(x=60, y=60, xend=60, yend=90,color='grey34',size=0.5, arrow=arrow(length = unit(0.02, "npc")))+
  geom_point(aes(x=0, y=0), color='red', inherit.aes = FALSE,shape=8)+
  annotate('text', x=85, y=75, label='A')+
  annotate('text', x=35, y=75, label='P')+
  annotate('text', x=60, y=55, label='V')+
  annotate('text', x=60, y=96, label='D')

  p


q <- ggplot(rotated_df)+theme_gray()+
  geom_segment(aes(x = Rotated_Moved_Nuclei_X,y = Rotated_Moved_Nuclei_Y,xend = Rotated_Moved_GM130_X ,
                   yend = Rotated_Moved_GM130_Y,color=sample),arrow=arrow(length = unit(0.02, "npc")),
               show.legend = F)+
  coord_equal(ratio=1)+ #make sure the x, y axis scale the same
  #scale_x_continuous(name='',breaks = seq(-50,50,by =20))+ # modify the x axis breaks and range
  #scale_y_continuous(name='',breaks=seq(-100, 150, by =20))+
  #change the color
  #scale_color_manual(name='Tip Types',values = c('Parent tip'='green', 'Right tip'='black', 'Left tip'='blue', 'Cleft'='red'))#+
  #facet_wrap(~sample, nrow=2)
  xlab('')+
  ylab('')+
  xlim(-100,110)+
  ylim(-100,110)+
  geom_point(aes(x=0, y=0), color='red', inherit.aes = FALSE,shape=8,size=1)+
  #theme(legend.position="bottom", legend.box = "horizontal")+
  #ggtitle('Ring Cell Oritation')+
  geom_segment(x=45, y=80, xend=75, yend=80,color='grey34',size=0.5, arrow=arrow(length = unit(0.02, "npc")))+
  geom_segment(x=60, y=65, xend=60, yend=95,color='grey34',size=0.5, arrow=arrow(length = unit(0.02, "npc")))+
  annotate('text', x=85, y=80, label='A', size=3)+
  annotate('text', x=35, y=80, label='P', size=3)+
  annotate('text', x=60, y=57, label='V', size=3)+
  annotate('text', x=60, y=105, label='D', size=3)+
  facet_grid(~sample)
q


## Caculta the angel of Nuclei to GM130

moving_angle <- rotated_df %>% 
  mutate(dot= Rotated_Moved_Nuclei_X*Rotated_Moved_GM130_X+Rotated_Moved_GM130_Y*Rotated_Moved_Nuclei_Y) %>% 
  mutate(det=Rotated_Moved_Nuclei_X*Rotated_Moved_GM130_Y-Rotated_Moved_Nuclei_Y*Rotated_Moved_GM130_X) %>% 
  mutate(angle=atan2(det,dot))%>% 
  mutate(Direction = ifelse(angle<=0, "clockwise", 'counterclockwise')) %>%  
  group_by(sample, Direction) %>% summarise(Count=n())
#If the orientation of the coordinate system is mathematical with y up, you get counter-clockwise angles as is the convention in mathematics. 
#cacultet the clockwise angel based on the matrxi caculateion. (Detail not sure)

## alternative way, vecter(M,N), angel betwen M,N, is atan2(N[2],N[1]) - atan2(M[2],M[1])


t_test <- t.test(data=moving_angle, Count~Direction)



  
bar_angle <- ggplot(moving_angle)+
  theme_gray()+
  xlab('')+
  geom_bar(aes (x=Direction, y=Count, fill=Direction),stat = "summary", fun.y = "mean", show.legend = F)+
  annotate('text', x =1.5,y=9, label=paste0('p.value = ',round(t_test$p.value, 4)))#+
  #theme(axis.text.x = element_text(angle = 45,hjust=1))
  #theme(legend.position="bottom", legend.box = "horizontal")
bar_angle

bar_sample <- ggplot(moving_angle)+
  theme_gray()+
  xlab('')+
  geom_bar(aes (x=sample, y=Count, fill=Direction), stat = "identity", position="dodge")+
  #theme(axis.text.x = element_text(hjust=1))+
  theme(legend.position="bottom", legend.box = "horizontal", legend.title = element_text(size = 10),
        legend.text = element_text( size = 8), 
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        legend.box.spacing=unit(-0.5, 'cm'))

bar_sample


#
library(ggpubr)
ggarrange(ggarrange(p,bar_angle, labels = c('A','B'), ncol=2),
          ggarrange(q, bar_sample, labels = c('C','D'), nrow  = 2), 
          nrow =2,
          common.legend = TRUE, legend="bottom")




mixed_sample_plot <- ggarrange(p,bar_angle, labels = c('A','B'), ncol=2, widths = c(2,1), heights = c(1,1))
mixed_sample_plot
each_sample_plot <- ggarrange(q, bar_sample, labels = c('C','D'), nrow  = 2, heights = c(1,1.2), widths = c(2,1))
each_sample_plot
ggarrange(mixed_sample_plot,each_sample_plot, ncol = 1, nrow = 2, heights = c(1,1))

library(gridExtra)
grid.arrange(grobs=list(mixed_sample_plot, each_sample_plot))

# grid.arrange(grobs=list(p,bar_angle,q, bar_sample),
#              widths=c(2,1),
#              layout_matrix=rbind(c(1,2),
#                                  c(3,3),
#                                  c(4,4)))



# need to vector caculation to caculate the matrix
##end point need to moved to the original v3 = v1 + v2, v2 is our vector, v1 is original point. v3, is the end point
vector_df <- rotated_df %>% 
  mutate(v2_X=Rotated_Moved_GM130_X-Rotated_Moved_Nuclei_X) %>% 
  mutate(v2_Y=Rotated_Moved_GM130_Y-Rotated_Moved_Nuclei_Y) %>% 
  mutate(angle=atan2(v2_Y, v2_X)-atan2(Rotated_Moved_Nuclei_Y, Rotated_Moved_Nuclei_X)) %>% # angle negative means coutercolockwise
  mutate(simple_angle=case_when(
    angle>pi ~ angle-2*pi,
    angle< -pi ~ 2*pi+ angle,
    TRUE ~angle
    )) %>% 
  mutate(degree=180/pi*simple_angle) 

hist(abs(vector_df$degree))

#write.csv(vector_df, 'results/Ring_cell_cooridiate_Final_with_Angle_towards_Center.csv', row.names = F)

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


#  take the coordiate and shrikage to r to 1. 

arrows.2 <- function(x0, y0, x1, y1, length.ar, angle.ar, ...){
  r= sqrt(x0^2+y0^2)
  x=x0/r
  y=y0/r
  
  ab <- cos(angle.ar) * length.ar
  bc <- sign(sin(angle.ar)) * sqrt(length.ar^2 - ab^2)
  
  x2 <- x + ab
  y2 <- y + bc
  
  arrows(x, y, x2, y2, ...)
}


circle(0, 0, 1, add = FALSE)
for(i in 1:nrow(vector_df)){
  arrows.2(x0 = vector_df$Rotated_Moved_Nuclei_X[i], 
           y0 = vector_df$Rotated_Moved_Nuclei_Y[i], 
           x1 = vector_df$Rotated_Moved_GM130_X[i],
           y1 = vector_df$Rotated_Moved_Nuclei_Y[i],
           length.ar = 0.2, angle.ar = vector_df$simple_angle[i])
}


library(circular)

ray_test <- rayleigh.test(vector_df$simple_angle)


rose.diag(circular(vector_df$angle), units = 'degrees', prop=2.5, main = 'Ring cells orintation', 
          sub =paste0('Rayleigh Test  p-value=', ray_test$p.value ))

ray_test_abs <- rayleigh.test(abs(vector_df$simple_angle))

rose.diag(circular(abs(vector_df$simple_angle)), units = 'degrees', prop=2.5, main = 'Ring cells orintation_abs', 
          sub =paste0('Rayleigh Test p-value=', ray_test_abs$p.value ))


#library(dplyr)
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


#cacultet the of angel of parent tip to the x axis 
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
  #theme(legend.position="bottom", legend.box = "horizontal")+
  ggtitle('Ring Cell Oritation')+
  geom_segment(x=40, y=75, xend=80, yend=75,color='grey34',size=0.5, arrow=arrow(length = unit(0.02, "npc")))+
  geom_segment(x=60, y=60, xend=60, yend=90,color='grey34',size=0.5, arrow=arrow(length = unit(0.02, "npc")))+
  annotate('text', x=85, y=75, label='A')+
  annotate('text', x=35, y=75, label='P')+
  annotate('text', x=60, y=55, label='V')+
  annotate('text', x=60, y=95, label='D')

  p


# q <- ggplot(coord_df)+theme_gray()
# q+geom_segment(aes(x = Moved_Nuclei_X,y = Moved_Nuclei_Y,xend = Moved_GM130_X ,yend = Moved_GM130_Y,color=sample),arrow=arrow(length = unit(0.02, "npc")))+
#   coord_equal(ratio=1)+ #make sure the x, y axis scale the same
#   #scale_x_continuous(name='',breaks = seq(-50,50,by =20))+ # modify the x axis breaks and range
#   #scale_y_continuous(name='',breaks=seq(-100, 150, by =20))+
#   #change the color
#   #scale_color_manual(name='Tip Types',values = c('Parent tip'='green', 'Right tip'='black', 'Left tip'='blue', 'Cleft'='red'))#+
#   facet_wrap(~sample, nrow=2)


## Caculta the angel of Nuclei to GM130

moving_angle <- rotated_df %>% mutate(dot= Rotated_Moved_Nuclei_X*Rotated_Moved_GM130_X+Rotated_Moved_GM130_Y*Rotated_Moved_Nuclei_Y) %>% 
  mutate(det=Rotated_Moved_Nuclei_X*Rotated_Moved_GM130_Y-Rotated_Moved_Nuclei_Y*Rotated_Moved_GM130_X) %>% 
  mutate(angle=atan2(det,dot)) %>% mutate(Direction = ifelse(angle<=0, "clockwise", 'counterclockwise')) %>%  
  group_by(sample, Direction) %>% summarise(Count=n())
#If the orientation of the coordinate system is mathematical with y up, you get counter-clockwise angles as is the convention in mathematics. 

t_test <- t.test(data=moving_angle, Count~Direction)
#cacultet the clockwise angel based on the matrxi caculateion. (Detail not sure)
  
bar_angle <- ggplot(moving_angle)+
  theme_gray()+
  xlab('')+
  geom_bar(aes (x=Direction, y=Count, fill=Direction),stat = "summary", fun.y = "mean", show.legend = F)+
  annotate('text', x =1.5,y=9, label=paste0('p.value = ',round(t_test$p.value, 4)))+
  theme(axis.text.x = element_text(angle = 45,hjust=1))
  #theme(legend.position="bottom", legend.box = "horizontal")
bar_angle

bar_sample <- ggplot(moving_angle)+
  theme_gray()+
  xlab('')+
  geom_bar(aes (x=sample, y=Count, fill=Direction), stat = "identity", position="dodge")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))
  #theme(legend.position="bottom", legend.box = "horizontal")

bar_sample


#
library(ggpubr)
ggarrange(p,ggarrange(bar_angle, bar_sample, labels = c('B','C'), ncol = 2, widths = c(1,2)), labels = c('A'), heights = c(1.7,1),ncol = 1, nrow =2)











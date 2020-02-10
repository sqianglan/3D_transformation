##Rotation and change the coordination
rm(list=ls())
library(dplyr)
library(tidyverse)

library(spdep) # for rotation

# center point
center <- read.csv('/Users/qianglan/Projects/3D_transformation/RawData/Cell polarity in R from Ewe/Tdtomato_bud_center_coord_20200107.csv')%>%
  separate(X, into = c('mice','sample','mbs')) %>% 
  unite('Sample',mice, sample, sep='_') %>% 
  select(-2)

APDV <- read.csv('/Users/qianglan/Projects/3D_transformation/RawData/Cell polarity in R from Ewe/Tdtomato_APDV_axis_coord_20200107.csv') %>% 
  select(1:9) %>% 
  separate(X, into = c('mice','sample','extra')) %>% 
  unite('Sample',mice, sample, sep='_') %>% 
  select(-2)

#caculate the roataion thea and degree
rotated_thea <- APDV %>% select(1:9) %>% 
  mutate(Length=sqrt((Anterior_X-Posterior_X)^2+ (Anterior_Y-Posterior_Y)^2)) %>% 
  mutate(thea=acos((Anterior_X-Posterior_X)/Length)) %>% mutate('rotation_angel' =ifelse(Anterior_Y>=Posterior_Y, 2*pi-thea, thea))

#rotation
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

#load the coord data
coord <- read.csv('results/Tdtomato_coord_original_20200209.csv') %>% select(-1)


df <- left_join(coord,center, by ="Sample")

#math sbustrate GM130 with center-cord
n <-  c(4:6) 
new_df <- mapply(`-`, df[n], df[n+6])
colnames(new_df) <- paste0("Moved_", colnames(new_df))
df2 <- cbind(df, new_df)
# mathe substrate Nuclei with center_cord
n <-  c(7:9) 
new_df2 <- mapply(`-`, df[n], df[n+3])
colnames(new_df2) <- paste0("Moved_", colnames(new_df2))


df2 <- cbind(df2,new_df2)

coord_df <- df2 %>% select(c('ID','Sample','Tissue', c(13:18))) ## coord moved into the 0,0 center. 

rotation_df <- left_join(coord_df, rotated_thea[,c('Sample','rotation_angel')], by='Sample') #%>% select(-c(6,9)) # remove the z coordinate


for (i in 1:length(rownames(rotation_df))){
  temp_GM130 <-as.data.frame(Rotation(rotation_df[i, c(4,5)], rotation_df$rotation_angel[i]))
  temp_Nuclei <- as.data.frame(Rotation(rotation_df[i, c(7,8)], rotation_df$rotation_angel[i]))
  temp <- cbind(temp_GM130, temp_Nuclei)
  if (i==1){
    rotated_df <- temp
  }else {rotated_df <- rbind(rotated_df,temp)}
}
names(rotated_df) <- paste('Rotated', names(rotation_df)[c(4,5,7,8)], sep = '_')

rotated_df <- cbind(rotation_df[,c(1:3,6,9)], rotated_df)

write.csv(rotated_df, "results/Final_rotated_moved_coord_with_Z.csv", row.names = F )

for(i in 1:3){
  current_tissue=levels(coord_df$Tissue)[i]
  temp <- filter(coord_df, Tissue==current_tissue)
  write.csv(temp, paste0("results/",current_tissue,"_moved_coord_with_Z.csv"), row.names = F)
}











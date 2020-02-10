# read the pooled dataset and do some simple cleaning for ploting.
rm(list=ls())
library(tidyverse)
library(dplyr)

### spots center data is not complete. Try to get other data.
####   ETr919_1R ETr919_4L ETr919_5L *ETr919_6R* KEK62hom_3L KEK62hom_3R

###Placode :"ETr919_1R_mb3"   "ETr919_4L_mb3"   "ETr919_5L_mb3"   "KEK62hom_3L_mb3" "KEK62hom_3R_mb3"


spots <- read_csv('../Tdtomato_spots_20200107.csv')

placode <- read_csv('../Tdtomato_bud_center_coord_20200107.csv')

spots2 <- spots %>% separate(Sample, into = c('Sample','add','origin'))

#check the origin equal to Tissue

sum(spots2$origin!=spots2$Tissue)  # if 0, pooling date process is right.

spots2 <- spots2 %>% unite('Sample', c(Sample, add), sep = '_') %>% select(-c(X1,origin, Stage.x))


placode_center <- placode %>% separate(X1, into = c('Sample','add','origin')) %>% 
  unite(Sample, c(Sample, add), sep = '_') %>% 
  select(-origin)

df <- left_join(spots2, placode_center, by='Sample') %>% na.omit()

moved_df <- df
for (i in c(4:6)){
  moved_df[,i] <- moved_df[,i]-moved_df[,i+7]
  moved_df[,i+3] <- moved_df[,i+3]-moved_df[,i+7]
}

Dermal_df_E11_25 <- moved_df %>% filter(Tissue=='de') %>% filter(Stage.y=='E11.25') %>% select(-c(Tissue,10:13))
Dermal_df_E11_5 <- moved_df %>% filter(Tissue=='de') %>% filter(Stage.y=='E11.5')%>% select(-c(Tissue,10:13))
MB_df_E11_25 <- moved_df %>% filter(Tissue=='mb3') %>%  filter(Stage.y=='E11.25')%>% select(-c(Tissue,10:13))
MB_df_E11_5 <- moved_df %>% filter(Tissue=='mb3') %>%  filter(Stage.y=='E11.5')%>% select(-c(Tissue,10:13))

Dermal_df <- moved_df %>% filter(Tissue=='de') %>%  select(-c(Tissue,10:13))
MB_df <- moved_df %>% filter(Tissue=='mb3') %>%  select(-c(Tissue,10:13))

write_csv(Dermal_df_E11_25, '../Tdtomato_Dermal_df_E11_25.csv')
write_csv(Dermal_df_E11_5, '../Tdtomato_Dermal_df_E11_5.csv')
write_csv(MB_df_E11_25, '../Tdtomato_MB_df_E11_25.csv')
write_csv(MB_df_E11_5, '../Tdtomato_MB_df_E11_5.csv')
write_csv(Dermal_df, '../Tdtomato_Dermal_df.csv')
write_csv(MB_df, '../Tdtomato_MB_df.csv')

write_csv(moved_df, '../Tdtomato_moved_df_all.csv')

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
folder <- "RawData/Ring_cells_directionality_181019"

fullpath <- paste(current_pwd,folder,sep = '/')

#### chang the file name pattern if necessary
pattern <- "*.csv" # regular expression pattern for the files containing "Green"


file_names <-  list.files(fullpath, pattern=pattern, recursive = T)
setwd(fullpath) # set the working dir into the folder place to facilitate the read_csv


for (i in 1:length(file_names)){
  #read each sheet and process seperately and then combine together.
  name <- strsplit(file_names[i], '_') %>% unlist()
  sample_name <- paste(name[1], name[2], sep = '_')
  
  temp <- read_delim(file_names[i], skip=1,delim=',') %>%  select(Variable, Min,'Surpass Object') %>% 
    filter(Variable %in% c('Position X', 'Position Y')) %>% separate(Variable, into=c('Variable','Coord')) %>%
    separate(4, into=c('Marker','1','2','ID')) %>% 
    select(Coord, ID,Min) 
  marker <- regmatches(file_names[i], gregexpr('GM130|Nuclei',file_names[i], ignore.case = T)) %>% unlist()
  temp <- temp %>% mutate(Marker=marker)%>% unite('Marker',Marker,Coord, sep = '_')%>% 
    spread(Marker, Min) %>% mutate(sample=sample_name) 
  
  if (i %% 2 == 1){
    df1 <- temp
  }
  else{
    df1 <- merge(df1, temp, by =c('ID', 'sample'))
    if (exists('df2')){
      df2 <- rbind(df2, df1)
    }
    else{
      df2 <- df1
    }
  }
}
setwd(current_pwd)

p <- ggplot(df2)+theme_gray()
p+geom_segment(aes(x = Nuclei_X,y = Nuclei_Y,xend = GM130_X ,yend = GM130_Y,color=sample),arrow=arrow(length = unit(0.02, "npc")))+
  coord_equal(ratio=1)+ #make sure the x, y axis scale the same
  #scale_x_continuous(name='',breaks = seq(-50,50,by =20))+ # modify the x axis breaks and range
  #scale_y_continuous(name='',breaks=seq(-100, 150, by =20))+
  #change the color
  #scale_color_manual(name='Tip Types',values = c('Parent tip'='green', 'Right tip'='black', 'Left tip'='blue', 'Cleft'='red'))#+
  facet_grid(~sample)



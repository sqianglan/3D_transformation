rm(list=ls())
library(tidyverse)
library(dplyr)



current_pwd <- getwd() # save the current working direct path to rest it back later

## chang the file location if necessary
folder <- "/Users/qianglan/Projects/3D_transformation/20200208/rawData/E11.5 Placode cell polarity"

#fullpath <- paste(current_pwd,folder,sep = '/')

#### chang the file name pattern if necessary
pattern <- "*.csv" # regular expression pattern for the files containing "Green"


file_names <-  list.files(folder, pattern=pattern, recursive = T)
setwd(folder) # set the working dir into the folder place to facilitate the read_csv, otherwise the function get some problem.

df_GM130 <- data.frame()  ##create an empty data frame
df_Nuclei <- data.frame()
n=0
for (i in 1:length(file_names)){
  print(i)
  #read each sheet and process seperately and then combine together
  name <- strsplit(file_names[i], '\\/|_|\\.') %>% unlist() # split the string by '/', or '_' or '.
  sample_name <- paste(name[5], name[6], sep = '_') #extract the sample name
  data_label <- name[length(name)-2]# extact the type of data (center, axix or start point or end point)
  #tissue <- name[length(name)-3]
  #stage <- strsplit(file_names[i], '\\/|_') %>% unlist() %>% .[1] %>% strsplit(' ') %>% unlist %>% .[1]
  #sample_name <- paste(sample_name,tissue, sep = '_')
  # better way to get type
  type <- regmatches(file_names[i], gregexpr('GM130|Nuclei',file_names[i], ignore.case = T)) %>% unlist()
  tissue <- gsub(paste0(type,' '), '', data_label)
  # if (grepl("Anterior|Dorsal|Posterior|Ventral", data_label, ignore.case=T)) next #those information already in the ADPCV files, skip
  #extension <- name[length(name)]
 
  
  temp <-read_csv(file_names[i], skip=3) %>% as.data.frame(.)%>%  # read the files and tranform into dataframe
      select(Variable, Min,'Surpass Object') %>% #selected interested data
      filter(Variable %in% c('Position X', 'Position Y', 'Position Z')) %>% #fileter the coordinate data
      separate(Variable, into=c('Variable','Coord')) # split the name
    
    # if (i %in% c(16,17,26,27,46,47,56,57)){temp <- temp %>% separate('Surpass Object', into=c('Marker','ID','1'))
    # }else {temp <- temp %>% separate('Surpass Object', into=c('Marker','1','ID','2'))}# split the object ID.
    
  temp <- temp %>% separate('Surpass Object', into=c('1','ID','2'), sep='\\[|\\]|\\{|\\}') %>% ## split by [ or ] or {or }
    select(Coord, ID,Min) %>% 
      mutate(Type=type)%>%#label the data with data type (GM130 or nuclear)
      mutate(Tissue=tissue) %>% 
      unite('Marker',Type,Coord, sep = '_')%>% 
      spread(Marker, Min) %>%
      mutate(Sample=sample_name) 
    #data_type ="spots"
n=n+nrow(temp)
  ## Merging the data into thw same files
  ##Merge the data
  if(type =='Nuclei'){df_Nuclei <- rbind(df_Nuclei, temp)}else if (type =='GM130'){df_GM130 <- rbind(df_GM130, temp)}
  
  if(i==length(file_names)){df <- merge(df_GM130,df_Nuclei, by=c('ID', 'Sample', 'Tissue'))}
}
  ##check if the variables exist or not
#   if (exists(data_type)){ 
#     if (data_type =='APDV'){ 
#       APDV <- rbind(APDV, temp)  # rbind for axis cooridinate data directly
#     } 
#     else if (data_type == "Placode"){Placode <- rbind(Placode, temp)}
#     else if(type %in% c('GM130', 'Nuclei')){ # if the data is the GM130 dataset,
#       if (exists(sample_name)){
#         df <- merge(eval(as.name(sample_name)), temp, by =c('ID', 'Sample', 'Tissue'))
#         spots <- rbind(spots,df)
#       } else{assign(sample_name, temp)}
#     } 
#   }else if(type %in% c('GM130', 'Nuclei') & !exists(sample_name)){
#     assign(sample_name, temp)
#   }else if(type %in% c('GM130', 'Nuclei') & exists(sample_name)){
#     spots <- merge(eval(as.name(sample_name)), temp, by =c('ID', 'Sample', 'Tissue'))
#   }else{assign(data_type, temp)}
#}
  
#df <- unique(df)

setwd(current_pwd)

write.csv(df, 'results/Tdtomato_coord_original_20200209.csv') # 


  





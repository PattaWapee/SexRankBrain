
###############################
##### main code ###############
###############################
#import file path
file.list <- scan("file.path.txt", what = "character", sep="\n")

typeof(file.list)

file.list

#set data name list
# Separate elements by /
path.file <- strsplit(file.list, "/")
path.file
name.file <- vector()
library("stringr")
for (i in 1:length(path.file)){
	if(str_detect(path.file[[i]][11],"GPL|array|RNA")){

		name.file[i] <- paste0(path.file[[i]][10],"_",path.file[[i]][11])
    print(name.file[i])
	}

	else{
		name.file[i] <- path.file[[i]][10]
    print(path.file[[i]][10])
	}

}

name.file

names(file.list) <- name.file
#import data table from file list

destination_dir <- "/Users/pattama/Desktop/UiB/sex_differences/web_apps/shinydashboard/Test/data/CC/"


for (i in 1:length(file.list)){
  file_des <- paste0(destination_dir,name.file[i],".txt")
  print(file.list[i])
  file.copy(file.list[i],  file_des)
}

library(EBImage)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


#str(celltype_probability_maps)
#a <-imageData(celltype_probability_maps)
convert_probmap_to_csv <- function(probmap, savedir){
  if(!dir.exists(savedir)) dir.create(savedir)
  celltypes <- names(probmap[1,1,])
  for (name in celltypes) {
    write.csv(probmap[,,name],file = paste0(savedir,"/",name,".csv"), quote=FALSE)
  }
}

## These are saved as .rdata where every item becomes the variable "celltype_probability_maps"

load("~/Dropbox (HMS)/Tessa_DSM/TNBC/Barber_outputs_072623/Leap001_ROI_002_celltype_probability_maps.RData")
convert_probmap_to_csv(celltype_probability_maps, "LEAP001_ROI_002")


load("~/Dropbox (HMS)/Tessa_DSM/TNBC/Barber_outputs_072623/Leap001_ROI_003_celltype_probability_maps.RData")
convert_probmap_to_csv(celltype_probability_maps, "LEAP001_ROI_003")



load("~/Dropbox (HMS)/Tessa_DSM/TNBC/Barber_outputs_072623/Leap003_ROI_003_celltype_probability_maps.RData")
convert_probmap_to_csv(celltype_probability_maps, "LEAP003_ROI_003")

### If Paul changes the variable name at some point, you can use this to load to a variable of your choice
d <- loadRData("~/Dropbox (HMS)/Tessa_DSM/TNBC/Barber_outputs_072623/Leap001_ROI_002_celltype_probability_maps.RData")


## example output https://www.dropbox.com/scl/fo/v2xsh53k9qq9nk0agjsjv/h?rlkey=gt2f44fx7ugh2scq6wg3xn0l2&dl=0
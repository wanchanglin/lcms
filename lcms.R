###### Z.Hall, Aug 2014 #####
## wl-12-03-2018, Mon: commence 
## wl-15-03-2018, Thu: tidy R codes

library(xcms)
library(reshape2)

############ Script to process lcms data using xcms #################
### Functions for making lipid library, deisotoping and annotating were from
### Image Scope (Nick Bond) and modified for use with lcms data #####
############ Z.Hall, Oct 2014 #########################################

############################# Make lipid library #######################
makelibrary <- function(sel.class, fixed=F, fixed_FA, lookup_lipid_class,
                        lookup_FA, lookup_element){

  print("Making library of lipid masses...")
  lookup_lipid_class <- cbind(lookup_lipid_class, sel.class)

  library<- numeric()
  for(i in 1:nrow(lookup_lipid_class)){
    if(lookup_lipid_class[i,"sel.class"] == T){
      ## key variables
      rounder = 3 # number of decimals the rounded masses are rounded to.
      ## lipidclass = "TG"
      lipidclass <- row.names(lookup_lipid_class[i,])

      ##determine how many FAs places to be used for combination and generate combination of FAs
      FA_number <- as.numeric(lookup_lipid_class[lipidclass,"FA_number"])
      if(fixed==TRUE)FAnum <- FA_number - 1 else FAnum <- FA_number
      s1<-combn(FA_expt,FAnum)

      ## if one place is fixed add this FA to the matrix
      if(fixed==TRUE){
        s1 <- rbind(s1,"fixed"=fixed_FA)
        FAnum <- FAnum +1
      }

      ## if sn2 or sn3 does not have FA bind 'empty' FA channel.
      if(FAnum == 1){ s1<-rbind(s1,sn2<-vector(mode="numeric",length=ncol(s1)),sn3<-vector(mode="numeric",length=ncol(s1)))
        FAnum = FAnum +2}
      if(FAnum == 2){ s1<-rbind(s1,sn3<-vector(mode="numeric",length=ncol(s1)))
        FAnum = FAnum +1}


      ## label the matrix
      if(FAnum == 3) row.names(s1) <-c("FA1", "FA2","FA3")

      ## add rows to matrix for massofFAs and formula
      massofFAs<-vector(mode="numeric",length=ncol(s1))
      s1 <- rbind(s1,massofFAs)
      formula<-vector(mode="numeric",length=ncol(s1))
      s1 <- rbind(s1,formula)
      ##row.names(s1) <-c("FA1", "FA2","FA3", "massofFAs")
      for(i in 1:ncol(s1)){

        ## for 3 FAs
        if(FAnum == 3){
          FA_1 <- as.character((s1[1,i]))
          FA_2 <- as.character((s1[2,i]))
          FA_3 <- as.character((s1[3,i]))
          s1["massofFAs",i] <- as.numeric((lookup_FA[FA_1, "FAmass"]))+as.numeric((lookup_FA[FA_2, "FAmass"]))+as.numeric((lookup_FA[FA_3, "FAmass"]))
          ##determine the formula
          temp_carbon <- as.numeric((lookup_FA[FA_1, "FAcarbon"]))+as.numeric((lookup_FA[FA_2, "FAcarbon"]))+as.numeric((lookup_FA[FA_3, "FAcarbon"]))
          temp_doublebond <- as.numeric((lookup_FA[FA_1, "FAdoublebond"]))+as.numeric((lookup_FA[FA_2, "FAdoublebond"]))+as.numeric((lookup_FA[FA_3, "FAdoublebond"]))
          s1["formula",i] <- paste(lipidclass,"(",temp_carbon,":",temp_doublebond,")", sep = "")
        }
      }

      ## calculate total mass
      totalmass<-vector(mode="numeric",length=ncol(s1))
      s1 <- rbind(s1,totalmass)

      for(i in 1:ncol(s1)){
        s1["totalmass",i] <- as.numeric(s1["massofFAs",i]) + as.numeric(as.character(lookup_lipid_class[lipidclass,"headgroup_mass"])) - (as.numeric(lookup_lipid_class[lipidclass,"FA_number"])*as.numeric(lookup_element["H","mass"]))
      }

      ##make rows for charged lipids masses
      protonated<-vector(mode="numeric",length=ncol(s1))
      ammoniated<-vector(mode="numeric",length=ncol(s1))
      sodiated<-vector(mode="numeric",length=ncol(s1))
      potassiated<-vector(mode="numeric",length=ncol(s1))
      deprotonated<-vector(mode="numeric",length=ncol(s1))
      chlorinated<-vector(mode="numeric",length=ncol(s1))
      acetate<-vector(mode="numeric",length=ncol(s1))
      s1 <- rbind(s1,protonated, ammoniated, sodiated, potassiated, deprotonated, chlorinated, acetate)

      ##calculate charged lipids masses
      for(i in 1:ncol(s1)){
        s1["protonated",i] <- round((as.numeric(s1["totalmass",i]) + as.numeric(lookup_element["H","mass"])),digits = 4)
        s1["ammoniated",i] <- round((as.numeric(s1["totalmass",i]) + as.numeric(lookup_element["NH4","mass"])),digits = 4)
        s1["sodiated",i] <- round((as.numeric(s1["totalmass",i]) + as.numeric(lookup_element["Na","mass"])),digits = 4)
        s1["potassiated",i] <- round((as.numeric(s1["totalmass",i]) + as.numeric(lookup_element["K","mass"])),digits = 4)
        s1["deprotonated",i] <- round((as.numeric(s1["totalmass",i]) - as.numeric(lookup_element["H","mass"])),digits = 4)
        s1["chlorinated",i] <- round((as.numeric(s1["totalmass",i]) + as.numeric(lookup_element["Cl","mass"])),digits = 4)
        s1["acetate",i] <- round((as.numeric(s1["totalmass",i]) + as.numeric(lookup_element["CH3COO","mass"])),digits = 4)
      }

      ##make rows for rounded charged lipids masses
      round.protonated<-vector(mode="numeric",length=ncol(s1))
      round.ammoniated<-vector(mode="numeric",length=ncol(s1))
      round.sodiated<-vector(mode="numeric",length=ncol(s1))
      round.potassiated<-vector(mode="numeric",length=ncol(s1))
      round.deprotonated<-vector(mode="numeric",length=ncol(s1))
      round.chlorinated<-vector(mode="numeric",length=ncol(s1))
      round.acetate<-vector(mode="numeric",length=ncol(s1))
      s1 <- rbind(s1,round.protonated, round.ammoniated, round.sodiated, round.potassiated, round.deprotonated, round.chlorinated, round.acetate)

      ##calculate rounded charged lipids masses
      for(i in 1:ncol(s1)){
        s1["round.protonated",i] <- round(as.numeric(s1["protonated",i]),digits = rounder)
        s1["round.ammoniated",i] <- round(as.numeric(s1["ammoniated",i]),digits = rounder)
        s1["round.sodiated",i] <- round(as.numeric(s1["sodiated",i]),digits = rounder)
        s1["round.potassiated",i] <- round(as.numeric(s1["potassiated",i]),digits = rounder)
        s1["round.deprotonated",i] <- round(as.numeric(s1["deprotonated",i]),digits = rounder)
        s1["round.chlorinated",i] <- round(as.numeric(s1["chlorinated",i]),digits = rounder)
        s1["round.acetate",i] <- round(as.numeric(s1["acetate",i]),digits = rounder)
      }

      library<-cbind(library,s1)
    }
  }
  return(library)
}

########################  Deisotoping #################################
deisotoping <- function(ppm=5, no_isotopes=2, prop.1=0.9, prop.2=0.5, 
                        spectra=spectra){

  C13_1 <- 1.003355
  C13_2 <- C13_1*2

  k = 0
  m = 0

  ## run loop to find isotopes for each ion.
  for(i in (1:(nrow(spectra)-1))){
    ## values of search
    mass <- as.numeric(spectra[i,1])
    intensity <- as.numeric(spectra[i,2])
    ## calculated values
    offset = (ppm * mass) / 1000000

    ## find isotope with ppm filter on isotpe
    search <- round((mass+C13_1),digits = 3)
    top <- search + offset
    bottom <- search - offset
    result <- spectra[as.numeric(spectra[,"intensity"]) <= (intensity*prop.1) & spectra[,1] >= bottom & spectra[,1] <= top & spectra[,"isotope"] == "",]
    result <- rbind(result,blank1 = "", blank2 = "")


    if(no_isotopes ==2){
      ## find isotope with ppm filter on isotpe
      search <- round((mass+C13_2),digits = 3)
      top <- search + offset
      bottom <- search - offset
      result_2 <- spectra[as.numeric(spectra[,"intensity"]) <= (intensity*prop.2) & spectra[,1] >= bottom & spectra[,1] <= top & spectra[,"isotope"] == "",]
      result_2 <- rbind(result_2,blank1 = "", blank2 = "")

    }


    if(nrow(result)>2){
      k = k +1
      spectra[i, "isotope"] <- paste(spectra[i, "isotope"]," ", "[",k,"]","[M]", sep="")
      for(j in 1:(nrow(result)-2)){
        indices <- which(spectra == result[j,1], arr.ind=TRUE)
        spectra[indices[,"row"], "isotope"] <- paste(spectra[indices[,"row"],"isotope"]," ","[",k,"]","[M+1]", sep="")
      }
      if(no_isotopes ==2 && nrow(result_2)>2){
        for(j in 1:(nrow(result_2)-2)){
          indices <- which(spectra == result_2[j,1], arr.ind=TRUE)
          spectra[indices[,"row"], "isotope"] <- paste(spectra[indices[,"row"],"isotope"]," ","[",k,"]","[M+2]", sep="")
        }
      }
    }

  }

  allpeaks <- as.data.frame(spectra)
  deisotoped  <- allpeaks[(grep("\\[M\\+", allpeaks$isotope, invert = T)),]
  isotopes  <- allpeaks[(grep("\\[M\\+", allpeaks$isotope, invert = F)),]
  results<-list(allpeaks, deisotoped, isotopes)
  summary <- paste(length(as.vector(deisotoped$mz.obs)),"monoisotopic peaks retained and", length(as.vector(isotopes$mz.obs)), "C13 isotopes discarded from", length(as.vector(allpeaks$mz.obs)), "detected ions", sep=" ")
  print(summary); log<-c(log,summary)
  return(deisotoped)
}

############################# Annotating ##############################
annotating <- function(deisotoped,
                       adducts=c(H=T, NH4=F, Na=T, K=F, dH=F, Cl=F, OAc=F),
                       ppm.annotate=10, dbase){
  print("Starting annotation")
  d.finalmz   <- as.vector(deisotoped$mz.obs)
  s1          <- dbase
  spectra     <- cbind(round(as.numeric(d.finalmz), digits=3),d.finalmz)
  combined    <- vector()
  sel.adducts <- vector()
  index       <- 13 # offset to search only rounded masses in library

  for(a in 1:length(adducts)){
    if(adducts[a]==T) sel.adducts<- c(sel.adducts, index+a)
  }

  for(i in 1:nrow(spectra)){
    search<-as.numeric(spectra[i,1])
    offset = (ppm.annotate * search) / 1000000
    top <- search + offset
    bottom <- search - offset
    result<-which( s1[sel.adducts,] >= bottom & s1[sel.adducts,] <= top ,arr.ind=TRUE)
    if(nrow(result)>0){
      for(j in 1:nrow(result))
      {
        col<-result[j,"col"]
        row<-result[j,"row"]; row<-sel.adducts[row]
        ## determine the adduct that was matched, summarising match information from library for matched mass (as 'data')
        ##determine which adduct
        if(row == "14"){adduct <- "protonated"; name.adduct <- "H"}
        if(row == "15"){adduct <- "ammoniated"; name.adduct <- "NH4"}
        if(row == "16"){adduct <- "sodiated"; name.adduct <- "Na"}
        if(row == "17"){adduct <- "potassiated"; name.adduct <- "K"}
        if(row == "18"){adduct <- "deprotonated"; name.adduct <- "-H"}
        if(row == "19"){adduct <- "chlorinated"; name.adduct <- "Cl"}
        if(row == "20"){adduct <- "acetate"; name.adduct <- "OAc"}

        a.ppm = round(abs(((as.numeric(spectra[i,2]) - as.numeric(s1[adduct,col])) / as.numeric(spectra[i,2]))* 1000000), digits=1)

        ##make vector with summary of match and paired match
        data<- c(s1[row,col],  s1[adduct,col],spectra[i,2], a.ppm, s1["formula",col], name.adduct,  s1["protonated",col],  s1["FA1",col], s1["FA2",col], s1["FA3",col])

        ##make matrix of search results
        combined <- rbind(combined,unlist(data, use.names=F))
      }
    }
  }
  if(length(combined)>0){
    colnames(combined) <- c("mz.matched",  "mz.matched.lib", "mz.observed", "ppm", "formula", "adduct", "mz.lib.protonated", "FA1", "FA2", "FA3")

    ids <- unique.matrix(combined[,c(3, 5,6)])
    annotations<-cbind(d.finalmz, "")
    for(i in 1:nrow(annotations)){
      result <- which(ids[,1] == annotations[i,1], arr.ind = T)
      if(length(result)>0){
        for(j in 1:length(result)){
          annotations[i,2]<-paste(annotations[i,2],"[",ids[result[j],"formula"],"+",ids[result[j],"adduct"],"]", sep="")
        }
      }
    }

    summary <- paste(length(annotations[annotations[,2] !="",2]),"from",length(as.vector(deisotoped$mz.obs)),"monoisotopic peaks were annoated (using accuract mass) with a",ppm.annotate,"ppm tollerance", sep=" "); log<-c(log,summary)
    print(summary)

    return(annotations[,2])
  }
}

## ======================================================================== 
## parameter setting

spectra_dir     <- "C:/Users/hallz/Documents/Processing_temp/" # where your mzML files are stored
lib_dir         <- "O:/Groups/LPS/BioComputing/lcms_processing/libraries" # the folder where your library files are stored

FWHM            <- 3 # set approximate FWHM (in seconds) of chromatographic peaks
snthresh        <- 5 # set the signal to noise threshold
minfrac         <- 0.25 # minimum fraction of samples necessary for it to be a valid peak group
profmethod      <- "binlin"  # use either "bin" (better for centroid, default), "binlin" (better for profile)

adducts         <- c(H = T, NH4 = T, Na = T, K = F, dH = F, Cl = F, OAc = F) # choose which adduct you want to include in library search
ppm.annotate    <- 15   # choose ppm for annotations
ionisation_mode <- "positive"  # either "positive" or "negative" - determines which classes of lipids to search for

## ======================================================================== 
## Read in library files
setwd(lib_dir)

read                 <- read.csv('lib_FA.csv', sep=",", header=T)
lookup_FA            <- read[,2:4]
row.names(lookup_FA) <- read[,1]

read                          <- read.csv('lib_class.csv', sep=",", headeq=T)
lookup_lipid_class            <- read[,2:3]
row.names(lookup_lipid_class) <- read[,1]

read                      <- read.csv('lib_element.csv', sep=",", header=T)
lookup_element            <- read[,2:3]
row.names(lookup_element) <- read[,1]

read                   <- read.csv('lib_modification.csv', sep=",", header=T)
lookup_mod             <- read[,2:ncol(read)]
row.names(lookup_mod ) <- read[,1]

FA_expt <- list('10','12','14','15','16','16.1','17','17.1','18','18.1',
                '18.2','18.3','20.3','20.4','20.5','21','22','22.5',
                '22.6','24.1')

## fatty acids to be included in library

if (ionisation_mode == "positive"){
  sel.class <- c( ## lipid classes to include - will need to change for negative ion mode
      T, #TG
      T, #DG
      T, #PC
      F, #PA
      T, #PE
      T, #PS
      F, #PG
      F, #PI
      F, #PIP
      F, #PIP2
      F, #PIP3
      T, #LysoPC
      T, #DG-H20
      T, #CE
      F, #FFA
      T, #SM
      T  #Cer
  )
}

if (ionisation_mode == "negative"){
  sel.class <- c( ## lipid classes to include - will need to change for negative ion mode
      F, #TG
      F, #DG
      T, #PC
      T, #PA
      T, #PE
      T, #PS
      T, #PG
      T, #PI
      T, #PIP
      T, #PIP2
      T, #PIP3
      F, #LysoPC
      F, #DG-H20
      F, #CE
      T, #FFA
      F, #SM
      F  #Cer
  )
}

dbase <- makelibrary(sel.class, fixed=F, fixed_FA, lookup_lipid_class, 
                     lookup_FA, lookup_element)

## wl-15-03-2018, 周四: 'fixed-FA' will be one of 'FA_expt'.

## ======================================================================== 
## XCMS

## place converted mzML files in working directory
setwd(spectra_dir)
path  <- paste(getwd())
files <- list.files(path, recursive=T, full.names=F, pattern=".mzML")

## ------------------------------------------------------------------------ 
## xcms scripts for peak picking and integration
xset <- xcmsSet(files, method="matchedFilter", step = 0.1, sigma=FWHM/2.3548, 
                snthresh=snthresh, profmethod=profmethod) 
## feature recognition and convert group of files to xcms raw object

xset <- group(xset) # groups peaks, temporarily
## xset2 <- retcor(xset, family = "symmetric", plottype = "mdevden", 
##                 smooth = "loess", span = 0.2, missing = 5, extra = 1) 

## corrects retention times
xset2 <- retcor(xset, method="obiwarp", profStep=0.1, plottype="deviation")  
## different method for retention time correction - more robust but slower

xset2 <- group(xset2, bw = 5,  minfrac = minfrac, mzwid = 0.025) 
## groups data based on corrected retention times

## ------------------------------------------------------------------------ 
## check for retention time outliers and exclude samples if necessary
rmsd <- sapply(1:length(xset2@rt$corrected), function(x) {
                 sqrt(sum((xset2@rt$raw[[x]]-xset2@rt$corrected[[x]]) ** 2))
               })
plot(rmsd)

rmsd <- cbind(files, rmsd)
write.csv(rmsd, "RT_rmsd.csv")

## make output based on unique peak groups (mz and RT groups)
names <- groupnames(xset2, rtdec=2, mzdec=4)  
## define decimal places for mz and RT values
names <- substr(names, 2, 18)

out <- groupval(xset2, method="maxint", value="maxo", intensity="maxo") 
## if multiple peaks choose one with highest intensity based on maxo (peak
## hieght raw)
rownames(out) <- names
out[is.na(out)] <- 0 # replace NAs with 0
write.csv(out, "xcms_peak_height_raw.csv")

out <- groupval(xset2, method="maxint", value="into", intensity="maxo")
rownames(out) <- names
out[is.na(out)] <- 0
write.csv(out, "xcms_peak_area_raw.csv")

## =======================================================================
## deisotoping - working with raw peak area - can also change "out" to use
## for peak height or the fitered results (value = into/intf/maxo/maxf)

out <- groupval(xset2, method="maxint", value="into", intensity="maxo") 
## if multiple peaks choose one with highest intensity based on maxo (peak
## hieght raw)
rownames(out) <- names
out[is.na(out)] <- 0 # replace NAs with 0
names_split <- colsplit(as.vector(as.character(names)), "T", c("mz", "RT"))

spectra <- cbind(names_split$mz, out[,1], "", "")
colnames(spectra) <- c("mz.obs", "intensity", "isotope", "modification" )

## -------------------------------------------------------------------- 
deisotoped <- deisotoping(ppm=5, no_isotopes=2, prop.1=0.9, prop.2=0.5, 
                          spectra=spectra)

indices <- match(rownames(deisotoped), rownames(out))
deisotoped_out <- out[indices, ]

names_deiso <- rownames(deisotoped_out)
names_deiso <- colsplit(as.vector(as.character(names_deiso)), "T", c("mz", "RT"))
deisotoped_out <- cbind(names_deiso$mz, names_deiso$RT, 
                        deisotoped_out[,1:ncol(deisotoped_out)])

write.csv(deisotoped_out, "xcms_peak_area_raw_deisotoped.csv")

## =======================================================================
## annotating
annotated <- annotating(deisotoped, adducts, ppm.annotate, dbase)
final_out <- cbind(annotated, deisotoped_out)
write.csv(final_out, "final_annotated_desiotoped.csv")

## =======================================================================
f# normalising based on TIC
data <- read.csv("final_annotated_desiotoped.csv")
rownames(data) <- data[,1]
sums <- as.vector(colSums(data[,5:ncol(data)], na.rm=T))
factor <- sums/mean(sums)
TIC_normalised <- cbind(data[,2:4], t(t(data[,5:ncol(data)])/factor))

write.csv(TIC_normalised, "final_annotated_deisotoped_norm.csv")

##hist(factor, breaks=100)

#########################################################################
## =======================================================================
## fill in missing peaks by integrating noise
xset3 <- fillPeaks(xset2) # fills in missing peaks
names <- groupnames(xset3, rtdec=2, mzdec=4)  
## define decimal places for mz and RT values
names <- substr(names, 2, 18)

## =======================================================================
## deisotope
out <- groupval(xset3, method="maxint", value="into", intensity="maxo") # if multiple peaks choose one with highest intensity based on maxo (peak hieght raw)
rownames(out) <- names
names_split <- colsplit(as.vector(as.character(names)), "T", c("mz", "RT"))

spectra <- cbind(names_split$mz, out[,1], "", "")
colnames(spectra) <- c("mz.obs", "intensity", "isotope", "modification" )

deisotoped <- deisotoping(ppm=5, no_isotopes=2, prop.1=0.9, prop.2=0.5, 
                          spectra=spectra)

indices <- match(rownames(deisotoped), rownames(out))
deisotoped_out <- out[indices, ]

names_deiso <- rownames(deisotoped_out)
names_deiso <- colsplit(as.vector(as.character(names_deiso)), "T", c("mz", "RT"))
deisotoped_out <- cbind(names_deiso$mz, names_deiso$RT, 
                        deisotoped_out[,1:ncol(deisotoped_out)])

write.csv(deisotoped_out, "xcms_fill_peak_area_raw_deisotoped.csv")

## =======================================================================
## annotate
annotated <- annotating(deisotoped, adducts, ppm.annotate, dbase)
final_out <- cbind(annotated, deisotoped_out)
write.csv(final_out, "final_annotated_desiotoped_fill.csv")

## =======================================================================
## normalising based on TIC
data <- read.csv("final_annotated_desiotoped_fill.csv")
rownames(data) <- data[,1]
sums <- as.vector(colSums(data[,5:ncol(data)], na.rm=T))
factor <- sums/mean(sums)
TIC_normalised <- cbind(data[,2:4], t(t(data[,5:ncol(data)])/factor))

write.csv(TIC_normalised, "final_annotated_deisotoped_fill_norm.csv")

##hist(factor, breaks=100)

##########################################################################

## plot extracted ion chromatograms for each peak group (corrected RT) if
## desired, can choose specific samples and/or mz

## xic.corr <- getEIC(xset3, rt = "corrected", groupidx = 1:nrow(xset3@groups))
## plot(xic.corr, xset3, groupidx=1:nrow(xset3@groups), sampleidx=1:length(files))


################################################################################################################################


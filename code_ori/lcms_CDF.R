## wl-24-08-2017, Thu:

## ########## Script to process lcms data using xcms #####

## Functions for making lipid library, deisotoping and annotating were from
## Image Scope (Nick Bond) and modified for use with lcms data
## ########## Z.Hall, Oct 2014####

## ########### Make lipid library #######################
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


## ########################  Deisotoping ###########################
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

## ################################# Annotating #######################
annotating <- function(deisotoped,
                       adducts=c(H=T, NH4=F, Na=T, K=F, dH=F, Cl=F, OAc=F),
                       ppm.annotate=10, dbase){
  print("Starting annotation")
  d.finalmz <-as.vector(deisotoped$mz.obs)
  s1<-dbase
  spectra<-cbind(round(as.numeric(d.finalmz), digits=3),d.finalmz)
  combined <- vector()
  sel.adducts<- vector()
  index<- 13 # offset to search only rounded masses in library
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


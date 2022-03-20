library("readr")
library("RMassBank")
library("enviPat")
library("OrgMassSpecR")
data("isotopes")

#Function that provides m/z given a molecular formula and an adduct----
adductsTOprecursors<-function(screened_adducts,formula=Compound$Formula){
  library("enviPat")
  data("isotopes")
  data("adducts")
  if(!is.na(formula) & !is.na(screened_adducts)){
    i<-1
    for(i in 1:dim(adducts)[1]){
      if(adducts$Ion_mode[i]=="positive"){
        if(adducts$Charge[i]==1 || adducts$Charge[i]==-1)  adducts$Name[i]<-paste0("[",  adducts$Name[i],"]+")
        else adducts$Name[i]<-paste0("[",  adducts$Name[i],"]",adducts$Charge[i],"+")
      } else {
        if(adducts$Charge[i]==1 || adducts$Charge[i]==-1)  adducts$Name[i]<-paste0("[",  adducts$Name[i],"]-")
        else adducts$Name[i]<-paste0("[",  adducts$Name[i],"]",adducts$Charge[i],"-")
      }
    }
    
    adducts[c(dim(adducts)[1]+1),]<-c("[M-H2O+H]+","M-18.010565",1,1,-18.010565,"positive","H2O", FALSE,1)
    adducts[c(dim(adducts)[1]+1),]<-c("[M+CH3OH+H]+","M+33.033491",1,1,33.033491,"positive","C1H3O1H2",FALSE,1)
    adducts[c(dim(adducts)[1]+1),]<-c("[M]+","M-0.00054858",1,1,-0.00054858,"positive",FALSE,FALSE,1)
    
    adducts$Name[adducts$Name=="[M-]-"]<-c("[M]-")
    
    adducts$Mass<-as.numeric(adducts$Mass)
    adducts$Charge<-as.numeric(adducts$Charge)
    adducts$Mult<-as.numeric(adducts$Mult)
    adducts$Multi<-as.numeric(adducts$Multi)
    # adducts[c(dim(adducts)[1]+1),]<-c("[M+H]2+","M+0.503364",2,1,0.503364,"positive","H1",FALSE,2)
    
    i<-1; mz<-rep(NA,length(screened_adducts)); 
    
    Compound2<-check_chemform(isotopes,formula); i<-1;
    if(length(adducts[which(as.character(adducts$Name)==screened_adducts[i]),5])>0){
      for(i in 1:length(screened_adducts)) mz[i] <- print(x=c(Compound2$monoisotopic_mass + adducts[which(as.character(adducts$Name)==screened_adducts[i]),5]),
                                                          digits=6)
      names(mz)<-screened_adducts
    } else {
      mz<-c(-50)
    }
  } else {
    mz<-c(-50)
  }
  
  return(mz)
  
}

#Download MassBank data repository and unzip----
files<-list.files("MassBank-data-main",pattern=".txt", recursive = T)

#Processing MassBank records to get the spectra in a list object named "output" and their m/z in a vector named "mz_output"----
i<-1; output<-list()
for(i in 1:length(files)){
  tmp <- read_csv(paste0("MassBank-data-main/",files[i]))
  tmp <- as.data.frame(tmp)
  spectra<-tmp[c(which(tmp=="PK$PEAK: m/z int. rel.int.")+1):c(nrow(tmp)-1),]
  spectra2<-do.call(rbind.data.frame, strsplit(spectra," "))
  names(spectra2)<-c("mz","int","rel.int")
  
  spectra2<-spectra2[,c("mz","rel.int")]
  spectra2$mz<-as.numeric(spectra2$mz)
  spectra2$rel.int<-as.numeric(spectra2$rel.int)
  names(spectra2)<-c("mz","int")
  
  output[[i]]<-spectra2
  
  screened_adduct<-gsub(x=tmp[which(grepl(x=tmp[,1], pattern="PRECURSOR_TYPE")),1],pattern="MS\\$FOCUSED_ION: PRECURSOR_TYPE ", replacement="")
  
  if(length(screened_adduct)>0){
    mz_theor<-adductsTOprecursors(screened_adducts=screened_adduct,
                      formula=gsub(x=tmp[which(grepl(x=tmp[,1], pattern="FORMULA")),1],pattern="CH\\$FORMULA: ", replacement="")
    )
  } else {
    mz_theor<-c(-50)
  }
  
  
  names(output)[i]<-paste0(tmp[1,],";;",mz_theor)
  print(i)
}
mz_output <-  unlist(lapply(strsplit(names(output),";;"),function(x){ x[2]}))
mz_output <- as.numeric(mz_output)
mz_output <- data.frame(names(output),mz_output)
names(mz_output)<-c("name","mz")

save.image("environment_MassBank.RData")



#Load the directory that contains the MSMS spectra stored as csv----
queries<-list.files("HRMSMS17",pattern=".csv",full.names = TRUE)
output_spectra_similarity<-data.frame("row"=queries,"similarity"=NA, "matched record"=NA)


#Perform the spectral match between the experimental HRMS/MS and the library spectra----
i<-1
for(i in 1:length(queries)){
  spectra_query <- read_csv(queries[i])
  spectra_query <- as.data.frame(spectra_query); names(spectra_query)<-c("mz","int"); spectra_query$mz; spectra_query$int
  mz_query <- as.numeric(strsplit(strsplit(queries[i],"mz=")[[1]]," RTI")[[2]][1])
  
  selection<-mz_output[which(abs(mz_output$mz-mz_query)<0.02),]
  
  if(nrow(selection)>0 & nrow(spectra_query)>0){
    j<-1; a<-c()
    for(j in 1:nrow(selection)){
      a[j]<-OrgMassSpecR::SpectrumSimilarity(spectra_query,  
                                             output[[selection$name[j]]],
                                             print.graphic=FALSE,t=0.05)
    }
    names(a)<-selection$name

    if(!all(is.nan(a))){
    output_spectra_similarity$similarity[i]<-a[which.max(a)]
    output_spectra_similarity$matched.record[i]<-names(a)[which.max(a)]
    }
  }
  print(i)
}

write.csv(output_spectra_similarity,"y2017_annotated.csv", row.names=FALSE)

#Add number of fragments evaluated for spectral similarity----
Y2017 <- read_csv("y2017_annotated.csv")
Y2017 <- as.data.frame(Y2017)
Y2017$row <- gsub(x=Y2017$row,pattern="HRMSMS17/",replacement="")

#Load the component list to get the origin----
Dataset_Y2017 <- read_csv("Dataset_Y2017_map data_1.csv")

#Add number of fragments and source origin columns----
i<-1
Y2017$num_of_fragments<-NA
Y2017$origin<-NA
for(i in 1:nrow(Y2017)){
  tmp <- read_csv(paste0("HRMSMS17/",Y2017$row[i]))
  tmp <- as.data.frame(tmp)
  head(tmp,4)
  
  Y2017$num_of_fragments[i]<-nrow(tmp)
  Y2017$origin[i]<-  Dataset_Y2017$`Predicted label`[which(Dataset_Y2017$Feature==gsub(x=Y2017$row[i],pattern=".csv",replacement=""))]
  print(i)
}

write.csv(Y2017,"y2017_annotated_v2.csv", row.names = FALSE)
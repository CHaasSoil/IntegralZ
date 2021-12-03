#R-Script to calculate the value of the Integral Z. See C.Haas, D.Heibach, P.Pértile, D.Holthusen, and R.Horn (2021). Comparison of different approaches for the calculation of the integrated zone and the cross-over value in amplitude sweep tests. 
#Create a folder named "data" and store the examplary data files (i.e., Examplary_data_1.txt and Examplary_data_2.txt) there. The script processes all *.txt-files stored in the folder "data" and saves the results (i.e., the value of the Integral Z; the shear strain and the shear stress where G' == G'') to the file "IntZ_real.txt". Check formatting of the input files before you run the script.
#See l.38 ff and l.83-85 for script settings
#require(stats)
require(sfsmisc)
#UPDATE:
#Additional calculation of YP_tau
Int_Z_real<-function(){#This is the basic function for the calculations.
  #data$loss_modulus<-data$Verlustmodul
 
  tau<-0
  quit<-0
  G<-0
  IntZ<-0
  schnitt<-0
  #calculate the first segment (Fig. 1) first
  sub_data<-subset(data, data$Deformation >= 0.001) #will be used for all other segments...
  #find first measurement point (MP)..
  MP_1<-min(sub_data$MP)
  #data$Verlustfaktor[MP_1] #first MP with loss factor >0.001
  #data$Verlustfaktor[MP_1-1] #last MP with loss factor <0.001
  d_y<-(data$Verlustfaktor[MP_1]-data$Verlustfaktor[MP_1-1]) #calculate the slope between x and x+1
  d_x<-(data$Deformation[MP_1]-data$Deformation[MP_1-1])
  slope<-d_y/d_x
  #calculate loss factor where Deformation equals 0.001
  # f(x) = m*Deformation + b
  #here: loss factor = (slope * 0.01) + data$Verlustfaktor[MP_1-1]
  loss_factor_01<-(slope* (0.001-data$Deformation[MP_1-1]))+ data$Verlustfaktor[MP_1-1] # Eq. 12
 
  #continue with data for Deformation >0.01
 # data<-sub_data #x<-1
  for (x in 2:length(data$Deformation)-1){ #for x==1 to x_max
    if(data$Verlustfaktor[x+1] >= 1 && quit==0){
      d_y<-(data$Verlustfaktor[x+1]-data$Verlustfaktor[x]) #calculate the slope between x and x+1
      d_x<-(data$Deformation[x+1]-data$Deformation[x])
      slope<-d_y/d_x#this is to add the last segment (Fig. 2): If the loss factor of x+1 is >=1...
      schnitt<-((1-data$Verlustfaktor[x])/slope)+data$Deformation[x] #Eq. 15; Deformation where the loss factor equals 1
      quit<-1 #needed to match the next else-condition
      d_y<-(data$loss_modulus[x+1]-data$loss_modulus[x]) #calculate the slope between x and x+1
      slope<-d_y/d_x
      
      G<-slope*(schnitt-data$Deformation[x]) + data$loss_modulus[x]#G'==G'', where the loss factor equals 1
    }
  }

  
  deformations<-c(0.001,subset(sub_data$Deformation,sub_data$Deformation<=schnitt),schnitt)
  loss_factors<-c(loss_factor_01,subset(sub_data$Verlustfaktor,sub_data$Verlustfaktor<=1),1)
  IntZ<-schnitt-integrate.xy(deformations,loss_factors,  use.spline=F)
  
  ret<-paste(IntZ,G,schnitt,sep=";")
  return(ret)
}


####Config:
#getwd()
read_from_line<-28 #line number with MP 1
out_file<-"IntZ_real.txt" #name of the output file
#Write header to file
write.table("Org_Name;IntZ;G;gamma;tau", file=paste(getwd(),out_file, sep="/"),quote=F, row.names=FALSE, na="",col.names=FALSE, append=F)
input_folder<-paste(getwd(),"//data//",sep="") #folder with input data


input_file<-(list.files(path=input_folder, pattern="*.txt",recursive=F))
print(input_file) #all input files
#y<-1
for(y in 1:length(input_file)){#calculate IZ for data in all files...
  name<-read.table(file=paste(input_folder,input_file[y],sep=""),skip=1,nrow=1, head=F, sep="\t", dec=",")
  print(paste(input_folder,input_file[y],sep=""))
  data<-read.table(file=paste(input_folder,input_file[y],sep=""),skip=read_from_line, head=F, sep="\t", dec=",", fill=T, na.strings=c("ungültiger Messpunkt","1.000.000.000.000.000.000.000.000.000.000"))
  #check data; make sure that MP1 is in the uppermost row. Define names of used columns 
  data$Verlustfaktor<-as.numeric(data$V7)#variable where the data of the loss factor is stored
  data$Deformation<-as.numeric(data$V2)#variable where the data of the shear strain factor is stored
  data$loss_modulus<-as.numeric(data$V6)#variable where the data of the loss modulus is stored
  data$MP<-as.numeric(data$V1)
  if(data$MP[1]!=1){print("Fatal error occured. Check variable: read_from_line and check input files for homogeneity")}
#  str(data)
  if(is.numeric(data$V2)){
    Z_real<-Int_Z_real()
    print(Z_real)
    out<-paste(name$V4,Z_real,sep=";")
    out<-gsub(" ","_",out)
   # out<-gsub("_",";",out)
    write.table(out, file=paste(getwd(),out_file, sep="//"), quote=F,row.names=FALSE, na="",col.names=FALSE, append=T)
  }
  else{
    out<-paste(name$V4,"NA;NA;NA",sep=";")
    out<-gsub(" ","_",out)
    #out<-gsub("_",";",out)
    write.table(out, file=paste(getwd(),out_file, sep="//"),quote=F, row.names=FALSE, na="",col.names=FALSE, append=T)
  }
}


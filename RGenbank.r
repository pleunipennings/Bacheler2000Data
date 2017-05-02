library(ape)
library(seqinr)
library(pegas)

source("GetConsensusB.r")

#NOTE: remove sequence P00083-14-03 
#NOTE: remove sequence P00102-85-03
#NOTE: remove sequence P00077-761-01

#so i downloaded the whole dataset from genbank. lets see what I can do with that!

#find number of patient in file.txt
#find the first AY number in file.txt

source("RResistanceMutations.r")

text <- readLines("sequence.gb")
grep("LOCUS       AY",text)->listlinesAY
grep("patient ",text)->listlinespatient
grep("study day ",text)->listlinesday
grep("/translation=",text)->listlinesAA
grep("ORIGIN",text)->listlinesDNA
grep("in study",text)->listlinesSt
grep("was treated with",text)->listlinesTx
totalseqs=length(listlinesAY)

patientinfo<-data.frame("seqID"=numeric(totalseqs),"seqName"=numeric(totalseqs),"patient"=numeric(totalseqs),"day"=numeric(totalseqs),"Tx"=numeric(totalseqs),"studyname"=numeric(totalseqs))
patientinfo[,(length(patientinfo[1,])+1):(length(patientinfo[1,])+NumImportantMuts)]<-0
for (i in 1:length(NNRTImuts$pos)){
	names(patientinfo)[min(which(substr(names(patientinfo),1,1)=="V"))]<-
	paste("NN",NNRTImuts$pos[i],sep="")}
for (i in 1:length(NRTImuts$pos)){
	names(patientinfo)[min(which(substr(names(patientinfo),1,1)=="V"))]<-
	paste("NR",NRTImuts$pos[i],sep="")}
for (i in 1:length(Lamimuts$pos)){
	names(patientinfo)[min(which(substr(names(patientinfo),1,1)=="V"))]<-
	paste("LA",Lamimuts$pos[i],sep="")}
for (i in 1:length(Indinavirmuts$pos)){
	names(patientinfo)[min(which(substr(names(patientinfo),1,1)=="V"))]<-
	paste("PI",Indinavirmuts$pos[i],sep="")}

DNAseqs<-data.frame("seqID"=numeric(totalseqs),"seqName"=numeric(totalseqs),"DNA"=numeric(totalseqs))
AAseqs<-data.frame("seqID"=numeric(totalseqs),"seqName"=numeric(totalseqs),"AA"=numeric(totalseqs))

for (i in 1:length(listlinesAY)){
#find seqID
l = listlinesAY[i];t=text[l];patientinfo$seqID[i]<-substr(t,regexpr("AY",t)[1],regexpr("AY",t)[1]+7)
#find patientID
lpa = listlinespatient[which(listlinespatient>l)[1]];t=text[lpa]
patientinfo$patient[i]=substr(t,regexpr("patient",t)[1]+8,regexpr("patient",t)[1]+13)
#find day of collection
lday = listlinesday[which(listlinesday>(l))[1]];t=text[lday]
patientinfo$day[i]=	substr(t,regexpr("study day ",t)[1]+10,regexpr("study day ",t)[1]+12)
for (c in 2:3)if (substr(patientinfo$day[i],c,c)==";"){patientinfo$day[i]=substr(patientinfo$day[i],1,c-1)}	
if (patientinfo$day[i]=="-1"){patientinfo$day[i]=0}
patientinfo$day[i]<-gsub(pattern="\\D","",patientinfo$day[i])
#get the amino acids
	lAA = listlinesAA[which(listlinesAA>(l))[1]]
	t=paste(text[lAA:(lAA+5)],collapse="")
	t<-gsub(pattern = "/translation=", replacement = "",  x = t)	
	t<-gsub('\\s', "",t, perl = TRUE)#get rid of whitespaces
	t<-sub('\\W', "",t, perl = TRUE)#and the first characted
	t<-strsplit(t,'\\W',  perl = TRUE)[[1]][1]
	AAseqs$AA[i]<-t;AAseqs$seqID[i]<-patientinfo$seqID[i]	
#get the DNA
	lDNA = listlinesDNA[which(listlinesDNA>(l))[1]]
	lDNAend =  min(length(text),listlinesAY[i+1]-2,na.rm=TRUE)
	t=paste(text[(lDNA+1):lDNAend],collapse="")
	t<-gsub(pattern = " ", replacement = "",  x = t)	
	t<-gsub(pattern = "//", replacement = "",  x = t)	
	t<-gsub(pattern = "[[:digit:]]", replacement = "",  x = t)	
	DNAseqs$DNA[i]<-t;DNAseqs$seqID[i]<-patientinfo$seqID[i]	
#create seqName	
#for (i in 1:length(listlinesAY)){
	nameforseqID<-"000"
	lengthday<-nchar(patientinfo$day[i])
	substr(nameforseqID,4-lengthday,3)<-patientinfo$day[i]
	patientinfo$seqName[i]<-paste(substr(patientinfo$patient[i],4,10),"_",nameforseqID,"_00",sep="")
#what Tx were they getting?
	lSt = listlinesSt[which(listlinesSt>(l))[1]]
	St<-paste(as.character(text[lSt]),as.character(text[lSt+1]))
	patientinfo$studyname[i]<-gsub(" ","",strsplit(strsplit(St,"study")[[1]][2],"on")[[1]])[1]
	lTx = listlinesTx[which(listlinesTx>(l))[1]]
	Tx<-paste(as.character(text[lTx]),as.character(text[lTx+1]),as.character(text[lTx+2]))
	Tx<-gsub("[: :][: :]"," ",strsplit(strsplit(Tx,"was ")[[1]][2],"gene")[[1]])[1]
	for (k in 1:5){Tx<-gsub("[: :][: :]"," ",Tx)}#remove spaces
	patientinfo$Tx[i]<-substr(Tx,1,nchar(Tx)-2)
#determine whether they have NNRTImuts
NNRTIcolumns<-which(substr(names(patientinfo),1,2)=="NN")
for (p in 1:length(NNRTImuts$pos)){
	patAA<-substr(AAseqs$AA[i],99+NNRTImuts$pos[p],99+NNRTImuts$pos[p])
	if (nchar(patAA)>0){
		if(grepl(patAA,NNRTImuts$mut[p])){
			patientinfo[i,NNRTIcolumns[p]]<-1
		}}}
#other NNRTImuts
NRTIcolumns<-which(substr(names(patientinfo),1,2)=="NR")
for (p in 1:length(NRTImuts$pos)){
	patAA<-substr(AAseqs$AA[i],99+NRTImuts$pos[p],99+NRTImuts$pos[p])
	if (nchar(patAA)>0){
		if(grepl(patAA,NRTImuts$mut[p])){
			patientinfo[i,NRTIcolumns[p]]<-1
		}}}
#Lamimuts
Lamicolumns<-which(substr(names(patientinfo),1,2)=="LA")
for (p in 1:length(Lamimuts$pos)){
	patAA<-substr(AAseqs$AA[i],99+Lamimuts$pos[p],99+Lamimuts$pos[p])
	if (nchar(patAA)>0){
		if(grepl(patAA,Lamimuts$mut[p])){
			patientinfo[i,Lamicolumns[p]]<-1
		}}}
#PImuts
Indinavircolumns<-which(substr(names(patientinfo),1,2)=="PI")
for (p in 1:length(Indinavirmuts$pos)){
	patAA<-substr(AAseqs$AA[i],Indinavirmuts$pos[p],Indinavirmuts$pos[p])
	if (nchar(patAA)>0){
		if(grepl(patAA,Indinavirmuts$mut[p])){
			patientinfo[i,Indinavircolumns[p]]<-1
		}}}
}

#remove incomplete sequences
if(length(which(nchar(DNAseqs$DNA)!=984))>0){
patientinfo<-patientinfo[-which(nchar(DNAseqs$DNA)!=984),]
AAseqs<-AAseqs[-which(nchar(DNAseqs$DNA)!=984),]
DNAseqs<-DNAseqs[-which(nchar(DNAseqs$DNA)!=984),]}
nchar(AAseqs$AA[1])->lengthAA
if(length(which(nchar(AAseqs$AA)!=lengthAA))>0){
patientinfo<-patientinfo[-which(nchar(AAseqs$AA)!=lengthAA),]
DNAseqs<-DNAseqs[-which(nchar(AAseqs$AA)!=lengthAA),]
AAseqs<-AAseqs[-which(nchar(AAseqs$AA)!=lengthAA),]}

#make sure all seqNames names are same length and unique
for (n in unique(patientinfo$seqName)){
	listseqs<-which(patientinfo$seqName==n)
	for (j in 1:length(listseqs)){
		pat<-listseqs[j]
		substr(patientinfo$seqName[pat],11-nchar(j),10)<-paste(j)}}

#make sure that DNAseqs and AAseqs have right seqNames
for (i in 1:length(DNAseqs[,1])){
	DNAseqs$seqName[i]<-patientinfo$seqName[which(patientinfo$seqID==DNAseqs$seqID[i])]
	AAseqs$seqName[i]<-patientinfo$seqName[which(patientinfo$seqID==AAseqs$seqID[i])]}

#remove unaligned seqs
patientinfo<-patientinfo[which(patientinfo$seqName!="102_085_03"),]
AAseqs<-AAseqs[which(AAseqs$seqName!="102_085_03"),]
DNAseqs<-DNAseqs[which(DNAseqs$seqName!="102_085_03"),]
patientinfo<-patientinfo[which(patientinfo$seqName!="083_014_03"),]
AAseqs<-AAseqs[which(AAseqs$seqName!="083_014_03"),]
DNAseqs<-DNAseqs[which(DNAseqs$seqName!="083_014_03"),]
patientinfo<-patientinfo[which(patientinfo$seqName!="077_671_01"),]
AAseqs<-AAseqs[which(AAseqs$seqName!="077_671_01"),]
DNAseqs<-DNAseqs[which(DNAseqs$seqName!="077_671_01"),]

#sort the dataframes
#patientoverview<-patientoverview[order(patientoverview$seqName),]
patientinfo<-patientinfo[order(patientinfo$seqName),]
AAseqs<-AAseqs[order(AAseqs$seqName),]
DNAseqs<-DNAseqs[order(DNAseqs$seqName),]

#get overview for patients per day
#for each patient 
patientoverview<-patientinfo[which(substr(patientinfo$seqName,9,10)=="01"),]
patientoverview$day<-as.numeric(patientoverview$day)

#add numclones
for (i in 1:length(patientoverview[,1])){
	patientoverview$numclones[i]<-length(which(patientinfo$patient==patientoverview$patient[i]&patientinfo$day==patientoverview$day[i]))}

NNRTIcolumns<-which(substr(names(patientinfo),1,2)=="NN")
NRTIcolumns<-which(substr(names(patientinfo),1,2)=="NR")
Lamicolumns<-which(substr(names(patientinfo),1,2)=="LA")
Indinavircolumns<-which(substr(names(patientinfo),1,2)=="PI")

for (p in c(NNRTIcolumns,NRTIcolumns,Lamicolumns,Indinavircolumns)){
	for (i in 1:length(patientoverview[,1])){
		patientoverview[i,p]<-mean(patientinfo[which(patientinfo$patient==patientoverview$patient[i]&patientinfo$day==patientoverview$day[i]),p])}}

#take only patients with baseline (less than day 5)+ at least one time point, and at least a sample of 5 at bl and 10 after. 
#I removed the need for a BL sample
AllPatients<-unique(patientoverview$patient)
PatientList<-AllPatients
for (p in AllPatients){
	p=AllPatients[1]
	NumDays<-length(which(patientoverview$patient==p))
	BL_OK<-patientoverview$day[which(patientoverview$patient==p)][1]<5&patientoverview$numclones[which(patientoverview$patient==p)][1]>4
	ENOUGHSAMPLES<-FALSE
	if (NumDays>1){ENOUGHSAMPLES<-sum(patientoverview$numclones[which(patientoverview$patient==p)][2:NumDays])>9}
	if (!ENOUGHSAMPLES)PatientList<-PatientList[-which(PatientList==p)]}
	
#write fasta files
if(FALSE){
for (patname in PatientList){
patDNA=DNAseqs[which(substr(DNAseqs$seqName,1,3)==substr(patname,4,6)),]	
numseqs=length(patDNA$seqID)
patfasta<-patDNA$DNA
names(patfasta)<-patDNA$seqName
filename=paste(patname,".fasta",sep="")	
write.dna(patfasta, filename, format = "fasta", append = FALSE, nbcol = 6, colsep = " ", colw = nchar(DNAseqs$DNA[1]), indent = NULL, blocksep = 1)}
}

write.csv(patientoverview,"PatientOverview.csv")#includes all patients
write(PatientList, "PatientListFastaFiles.txt")	#only patients with enough samples 
write.csv(patientinfo, "PatientInfo.csv")	
write.csv(AAseqs,"AAseqs.csv")


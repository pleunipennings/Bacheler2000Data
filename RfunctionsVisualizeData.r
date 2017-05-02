source("/Users/pleunipennings/Documents/Research/HIV/SoftSweepsInHIV/Bacheler2000/RResistanceMutations.r")
options(stringsAsFactors= FALSE)

#"N"= not segregating , "P" = polymorphic, "S" = singleton
if(FALSE){
	FASTAfile<-patfasta
	DAYtipsdayzero<-DayTipsDay0
#	DAYtipsdayzero<-1:(length(patfasta[,1]))
	AASEQS<-AAseqs
}

get.ListOfSegSites <- function(FASTAfile,DAYtipsdayzero,AASEQS) {
	counts=data.frame("nucleotide"=c("a","c","g","t"),"count"=c(0,0,0,0),"countday0"=c(0,0,0,0),"countafterday0"=c(0,0,0,0))
	ListFourFold<-c("gc","cg","gg","ct","cc","tc", "ac","gt")
	len<-length(FASTAfile[1,])	
	ListOfSegSites<-data.frame("site"=numeric(len),"majoritydayzero"=character(len),"nucleotide2"=character(len),"nucleotide3"=character(len),"nucleotide4"=character(len),"morethantwo"=numeric(len),"morethanthree"=numeric(len),"status"=character(len),"codon"=numeric(len),"triplet"=numeric(len),"fourfold"=numeric(len),"codon1"=character(len),"codon2"=character(len),"codon3"=character(len),"codon4"=character(len),"resistancecodon"=numeric(len),"RTcodon"=numeric(len),"polyafterday0"=numeric(len),syn2=numeric(len),syn3=numeric(len),syn4=numeric(len))
	ListOfSegSites$status<-"N"
	ListOfSegSites$polyafterday0<-0
	ListOfSegSites$strongmultipleton<-0
	for (site in 1:len){
	ListOfSegSites$site[site]<-site
	ListOfSegSites$codon[site]<-floor((site-1)/3)+1		
	counts$count[1]=length(which(FASTAfile[,site]=="a"))
	counts$count[2]=length(which(FASTAfile[,site]=="c"))			
	counts$count[3]=length(which(FASTAfile[,site]=="g"))
	counts$count[4]=length(which(FASTAfile[,site]=="t"))
	counts$countday0[1]=length(which(FASTAfile[DAYtipsdayzero,site]=="a"))
	counts$countday0[2]=length(which(FASTAfile[DAYtipsdayzero,site]=="c"))	
	counts$countday0[3]=length(which(FASTAfile[DAYtipsdayzero,site]=="g"))
	counts$countday0[4]=length(which(FASTAfile[DAYtipsdayzero,site]=="t"))
	l<-1:length(FASTAfile[,1]);l<-l[-DAYtipsdayzero]
	counts$countafterday0[1]=length(which(FASTAfile[l,site]=="a"))
	counts$countafterday0[2]=length(which(FASTAfile[l,site]=="c"))	
	counts$countafterday0[3]=length(which(FASTAfile[l,site]=="g"))
	counts$countafterday0[4]=length(which(FASTAfile[l,site]=="t"))
	ListOfSegSites$majoritydayzero[site]<-counts$nucleotide[order(counts$countday0,decreasing=TRUE)][1]
	or<-counts$nucleotide[order(counts$count,decreasing=TRUE)]
	or<-or[-which(or==ListOfSegSites$majoritydayzero[site])]
		for (x in or){if(counts$count[which(counts$nucleotide==x)]==0)or<-or[-which(or==x)]}	
	if (max(counts$count)==sum(counts$count)-1){#its a singleton
		ListOfSegSites$status[site]<-"S"
		ListOfSegSites$nucleotide2[site]<-or[1]}
	if (max(counts$count)<sum(counts$count)-1){#its a multipleton
		ListOfSegSites$status[site]<-"P"
		ListOfSegSites$nucleotide2[site]<-or[1]
		if(max(counts$countafterday0)<sum(counts$countafterday0))ListOfSegSites$polyafterday0[site]<-1
		if(sort(counts$count)[3]==2&sort(counts$countday0)[3]==1&sort(counts$countafterday0)[3]==1){ListOfSegSites$strongmultipleton[site]<-0}else{ListOfSegSites$strongmultipleton[site]<-1}	
#i want to ignore sites that are really two singletons
	}
	if (sort(counts$count,decreasing=TRUE)[3]>0){#there are more than 2 nucls
		ListOfSegSites$morethantwo[site]<-1
		ListOfSegSites$nucleotide3[site]<-or[2]}
	if (sort(counts$count,decreasing=TRUE)[4]>0){#there are more than 3 nucls
		ListOfSegSites$morethanthree[site]<-1
		ListOfSegSites$nucleotide4[site]<-or[3]}
	}
	for (site in 1:len){
#get triplet& check whether it is a fourfold degenerate site
		ListOfSegSites$triplet[site]<-paste(ListOfSegSites$majoritydayzero[which(ListOfSegSites$codon==ListOfSegSites$codon[site])],collapse="")
		ListOfSegSites$fourfold[site]<-0
		if (site/3==floor(site/3)&length(grep(paste(ListOfSegSites$majoritydayzero[c(site-2,site-1)],collapse=""),ListFourFold))>0)ListOfSegSites$fourfold[site]<-1
	}
	ListOfSegSites$synnonsyn<-0
	for (site in which(ListOfSegSites$status!="N")){
		codon = ListOfSegSites$codon[site]	
		#which inds carry the minority allele? 
		allele1<-ListOfSegSites$majoritydayzero[site]	
		allele2<-ListOfSegSites$nucleotide2[site]	
		allele3<-ListOfSegSites$nucleotide3[site]	
		allele4<-ListOfSegSites$nucleotide4[site]	
		carriesallele1<-names(FASTAfile[which(FASTAfile[,site]==allele1),1])[1]
		carriesallele2<-names(FASTAfile[which(FASTAfile[,site]==allele2),1])[1]
		carriesallele3<-names(FASTAfile[which(FASTAfile[,site]==allele3),1])[1]
		carriesallele4<-names(FASTAfile[which(FASTAfile[,site]==allele4),1])[1]
		ListOfSegSites$codon1[site]<-substr(AASEQS$AA[which(AASEQS$seqName==carriesallele1)],codon,codon)
		ListOfSegSites$codon2[site]<-substr(AASEQS$AA[which(AASEQS$seqName==carriesallele2)],codon,codon)
		if 	(ListOfSegSites$codon2[site]==ListOfSegSites$codon1[site]) {ListOfSegSites$syn2[site]<-"syn"} else {ListOfSegSites$syn2[site]<-"non"}
		if (ListOfSegSites$morethantwo[site]==1){
			ListOfSegSites$codon3[site]<-substr(AASEQS$AA[which(AASEQS$seqName==carriesallele3)],codon,codon)
			if (ListOfSegSites$codon3[site]==ListOfSegSites$codon1[site]) {ListOfSegSites$syn3[site]<-"syn"} else {ListOfSegSites$syn3[site]<-"non"}}
		if (ListOfSegSites$morethanthree[site]==1){ListOfSegSites$codon4[site]<-substr(AASEQS$AA[which(AASEQS$seqName==carriesallele4)],codon,codon)
			if 	(ListOfSegSites$codon4[site]==ListOfSegSites$codon1[site]) {ListOfSegSites$syn4[site]<-"syn"} else {ListOfSegSites$syn4[site]<-"non"}}
		if (ListOfSegSites$syn2[site]=="non"|ListOfSegSites$syn3[site]=="non"|ListOfSegSites$syn4[site]=="non"){ListOfSegSites$synnonsyn[site]<-"non"}
	}
	ListOfSegSites$resistancecodon<-0
	ListOfSegSites$resistanceaa<-0
#see if one of the two codons at this site are in the list of resistance mutations. 
#Protein=PRO 
	for (site in which(ListOfSegSites$synnonsyn=="non")){
	for (Rsite in 1:length(Indinavirmuts$pos)){
	if (ListOfSegSites$codon[site]==Indinavirmuts$pos[Rsite]){
	if (IsThereAMatch(paste(ListOfSegSites$codon1[site],ListOfSegSites$codon2[site],ListOfSegSites$codon3[site],ListOfSegSites$codon4[site],sep=""),Indinavirmuts$mut[Rsite])){ListOfSegSites$resistancecodon[site]<-"PI"
			ListOfSegSites$resistanceaa[site]<-Indinavirmuts$mut[Rsite]}}}}		
#repeat for NNRTI sites
	for (site in which(ListOfSegSites$synnonsyn=="non")){
	for (Rsite in 1:length(NNRTImuts$pos)){
	if (ListOfSegSites$codon[site]==(NNRTImuts$pos[Rsite]+99)){
	if (IsThereAMatch(paste(ListOfSegSites$codon1[site],ListOfSegSites$codon2[site],ListOfSegSites$codon3[site],ListOfSegSites$codon4[site],sep=""),NNRTImuts$mut[Rsite]))
		{ListOfSegSites$resistancecodon[site]<-"NN"
		ListOfSegSites$resistanceaa[site]<-NNRTImuts$mut[Rsite]}}}}		
#repeat for NRTI sites
	for (site in which(ListOfSegSites$synnonsyn=="non")){
	for (Rsite in 1:length(NRTImuts$pos)){
	if (ListOfSegSites$codon[site]==(NRTImuts$pos[Rsite]+99)){
	if(IsThereAMatch(paste(ListOfSegSites$codon1[site],ListOfSegSites$codon2[site],ListOfSegSites$codon3[site],ListOfSegSites$codon4[site],sep=""),NRTImuts$mut[Rsite])){ListOfSegSites$resistancecodon[site]<-"ZDV"
			ListOfSegSites$resistanceaa[site]<-NRTImuts$mut[Rsite]}}}}	
#repeat for Lami sites
	for (site in which(ListOfSegSites$synnonsyn=="non")){
	for (Rsite in 1:length(Lamimuts$pos)){
	if (ListOfSegSites$codon[site]==(Lamimuts$pos[Rsite]+99)){
	if(IsThereAMatch(paste(ListOfSegSites$codon1[site],ListOfSegSites$codon2[site],ListOfSegSites$codon3[site],ListOfSegSites$codon4[site],sep=""),Lamimuts$mut[Rsite])){ListOfSegSites$resistancecodon[site]<-"3TC"
			ListOfSegSites$resistanceaa[site]<-Lamimuts$mut[Rsite]}}}}
	for (site in 298:984){ListOfSegSites$RTcodon[site]<-ListOfSegSites$codon[site]-99}
return<-ListOfSegSites}



IsThereAMatch <- function(string1,string2) {
	length1<-nchar(string1)
	length2<-nchar(string2)
	match=0
	for (i in 1:length1){
		for (j in 1:length2){
			if (substr(string1,i,i)==substr(string2,j,j)){match=1}}}
	return(match)}

#                                  Author: Nadia Talent, 2010 
#              The author waives copyright and places this work in the Public Domain.
#                        http://creativecommons.org/publicdomain/zero/1.0/
#
#                                     Bestchar source code
#                      version in the R language (http://www.r-project.org/)
#
# A program for calculating some representative "Best character" coefficients for
# categorial data. (The same procedure would be used for numeric ranges, after they are 
# converted to categories.) A single character is processed.
#
#------------------------------------------------------------------------------------------------------
# The input file, called charinput.txt, is a single column representing one character.
# The first line is ignored, and could be used as a label for the character. 
# Subsequent lines list character states for each taxon, with '/' to separate states.
# Spaces and tab characters are ignored.
#
# For example, a character called 'Petal colour' with the states for four taxa is coded:
# Petal colour
# red
# white/red
# white
# white
#
# Two output files are created in the same workspace called bestcharoutput.txt
# and bestchardetails.txt.
#------------------------------------------------------------------------------------------------------
# Coefficients based on information theory
# 1. Information coefficient reverse engineered from Intkey, logarithms to base 2 are used.
#    using {} to represent subscripts
#    H=sigma(H{i}), H{i}=-(k{i}/t)*log2(k{i}/t)*SIG{i}/k{i}
#    where t=number of remaining taxa,
#    j varies from 1 to t,
#    s is the total number of states allowed for this character (including all taxa),
#    i varies from 1 to s,
#    k{i} is the number of taxa that have state i, 
#    SIG{i}=sigma(1/n{j})
# 2. The "Normalized information coefficient" as described by 
#    Talent, Dickinson, and Dickinson (submitted); 
#    logarithms to base t are used.
#
# Coefficients that use pairwise comparisons of taxa, i.e. distance measures:
# 1. The Separation Coefficient counts only taxa that are completely separable, so if any state
# is allowed by two taxa they are not separable
# 2. Pairwise average Jaccard Coefficient
#
#------------------------------------------------------------------------------------------------------
#  define the files
#------------------------------------------------------------------------------------------------------

f<-file(description="charinput.txt",open="r")
logfile<-file(description="bestcharoutput.txt",open="w") #explicit open, to keep the file open
detailfile<-file(description="bestchardetails.txt",open="w")

cat("Bestchar program, R version, diagnostic output\n", file=detailfile)
cat("Output from the Bestchar program, R version\n", file=logfile)

#------------------------------------------------------------------------------------------------------
#  input the data
#------------------------------------------------------------------------------------------------------
# step 1, Read in the data for all taxa (a single taxonomic character) and store them as a vector

line<-readLines(f,1)
cat(line, "\n", file=logfile,append=TRUE)
cat(line, "\n", file=detailfile,append=TRUE)
inputDat<-readLines(f,-1)
write(inputDat, file=logfile)
write(inputDat, file=detailfile)

totaltaxa<-inputDat
numtaxa<-length(totaltaxa) # called t above
cat("totaltaxa=", totaltaxa, " length=", numtaxa, "\n\n", file=detailfile,append=TRUE)
cat("Number of taxa=",numtaxa, "\n", file=logfile,append=TRUE)

#------------------------------------------------------------------------------------------------------
# For the Information Statistic, build three parallel lists.
# statelist: a list of the possible states i, is used only during the first pass through the data
# listofstatecounts: is a list of tallies of taxa that allow each state (=k{i})
# listofsigmas: a list of the sums for each state i, of 1/n{j}, where j is the number of states
# allowed for each taxon that allows state i
#------------------------------------------------------------------------------------------------------
statelist<-list()
listofstatecounts<-list() # these are called k{i} above, in H{i}=log2(k{i}/t)*SIG{i}/t)
listofsigmas<-list() # these are called SIG{i} above 
for(taxon1 in totaltaxa) {
   taxon2<-c(strsplit(taxon1,"/"))[[1]]
   numstates<-length(taxon2) #the number of states for this taxon
   for(l in taxon2) {
      index1<-grep(l,statelist)
      cat("index1=", index1, "\n", file=detailfile,append=TRUE)
   	  if (length(index1)<1) { #state not already present in the list
    	   statelist<-append (statelist, l, after = length(statelist))
    	   listofstatecounts<-append (listofstatecounts, 1, after = length(listofstatecounts))
    	   listofsigmas<-append (listofsigmas, 1/numstates, after = length(listofsigmas))
           cat("appending new state =",l,"\n", file=detailfile,append=TRUE)
           cat("numstates=", numstates, " set sigma to ", 1/numstates,"\n",
              file=detailfile,append=TRUE)
    	   }
      else {
    	   index2<-grep(l,statelist)
           cat("incrementing at index2=",index2,"(state=",l,")\n", file=detailfile,append=TRUE)
           temp<-listofstatecounts[[index2]]
           listofstatecounts[[index2]]<-temp+1 # increment the counter of taxa that have this state
           temp<-listofsigmas[[index2]]
           listofsigmas[[index2]]<-temp+(1/numstates)
           cat("numstates=", numstates, " added= ", 1/numstates,"to sigma\n", 
              file=detailfile,append=TRUE) }
    	}
     }
 
write("\n", file=detailfile,append=TRUE)
write(paste("statelist=",statelist), file=detailfile,append=TRUE)
write(paste("listofstatecounts=",listofstatecounts), file=detailfile,append=TRUE)
write(paste("listofsigmas=",listofsigmas), file=detailfile,append=TRUE)
write("\n", file=detailfile,append=TRUE)

Totalnumstates<-length(statelist)
cat("Totalnumstates =",Totalnumstates, "\n", file=detailfile,append=TRUE)
cat("\nNumber of states=",Totalnumstates, "\n\n", file=logfile,append=TRUE)

IntKeyH<-0
PankhurstH<-0
H<-0
mindex<-1 #index for both lists, should be integer, was zero in Python
for (m in listofstatecounts) {
   cat("m=", m, "\n", file=detailfile, append=TRUE)
   mtemp<-m/numtaxa
   m1=log(mtemp, base=2) #log2(k{i}/t)
   m2=log(mtemp, base=Totalnumstates) #log base Totalnumstates
   m3=log(mtemp, base=numtaxa) #log base numtaxa
   cat("m1=",m1, "\n",file=detailfile,append=TRUE)
   cat("m2=",m2,"\n", file=detailfile,append=TRUE)
   cat("m3=",m3,"\n", file=detailfile,append=TRUE)
   n<-listofsigmas[[mindex]] # the equivalent entry in listofsigmas
   mindex<-mindex+1
   n1<-n/numtaxa
   cat("n1=",n,"/",numtaxa,"=",n1,"\n", file=detailfile,append=TRUE)
   IntKeyH<-IntKeyH+n1*m1 #base 2 logarithms, as in IntKey
   cat("IntKeyH=",IntKeyH,"\n", file=detailfile,append=TRUE)
   PankhurstH<-PankhurstH+n1*m2 #base numstates, as recommended by Pankhurst
   cat("PankhurstH=",PankhurstH,"\n", file=detailfile,append=TRUE)
   H<-H+n1*m3 #base numtaxa logarithms (to normalize the range as the number of taxa (t) decreases)
   cat("H=",H,"\n\n", file=detailfile,append=TRUE)
   }

#getcontext().prec=10
IntKeyH1<-round(IntKeyH*(-1),digits=2)
cat("Results to two decimal places:\n", file=detailfile,append=TRUE)
cat("Intkey-style information coefficient=",IntKeyH1, "\n", file=detailfile,append=TRUE)
cat("Results to two decimal places:\n", file=logfile,append=TRUE)
cat("Intkey-style information coefficient=",IntKeyH1, "\n", file=logfile,append=TRUE)
#getcontext().prec=2

PankhurstH1<-round(PankhurstH*(-1),digits=2)
cat("Pankhurst's information coefficient=",PankhurstH1, "\n", file=detailfile,append=TRUE)
cat("Pankhurst's information coefficient=",PankhurstH1, "\n", file=logfile,append=TRUE)

H1<-round(H*(-1),digits=2)
cat("Normalized information coefficient=",H1, "\n\n", file=detailfile,append=TRUE)
cat("Normalized information coefficient=",H1, "\n", file=logfile,append=TRUE)

#------------------------------------------------------------------------------------------------------
# For Jaccard and separation coefficients
# calculate number of possible pairs of taxa n!/2*(n-2)!=n*(n-1)/2
#------------------------------------------------------------------------------------------------------
cat("numtaxa=",numtaxa,"\n", file=detailfile,append=TRUE)
total_num_pairs<-numtaxa*(numtaxa-1)/2
cat("total_num_pairs=",total_num_pairs,"\n\n", file=detailfile,append=TRUE)

#------------------------------------------------------------------------------------------------------
# For Jaccard and Separation coefficients
# destructively compare taxon character-states from the top row to each of the other rows (taxa)
# taxon1 holds the current top-most taxon
#------------------------------------------------------------------------------------------------------
remainingtaxa<-list()
if (length(totaltaxa)>1) {
   taxon1<-totaltaxa[1] #pop off the top line (data for the first taxon)
   cat("popped taxon=", taxon1, "\n", file=detailfile,append=TRUE)
   remainingtaxa<-totaltaxa[2:length(totaltaxa)]
   totaltaxa<-remainingtaxa
   }

numdiffs<-0
Jaccard_similarity<-0

while (length(remainingtaxa)>0) { # no need to match anything with the last taxon in the queue
   write(paste("remainingtaxa at start of while loop=",remainingtaxa), file=detailfile,append=TRUE)
   for (taxon in remainingtaxa) {
#------------------------------------------------------------------------------------------------------
#  The separation coefficient counts only taxa that are completely separable, so if any element of
#  taxon1's list is included in taxon's list, they don't count as separable
#------------------------------------------------------------------------------------------------------
       match<-FALSE
#------------------------------------------------------------------------------------------------------
#  Average pairwise Jaccard distance uses the count of elements (character states) in the intersection 
#  (between two taxa) divided by the count of elements in the union. Real numbers are required.
#------------------------------------------------------------------------------------------------------
       unionlist<-list()
       intersection_list<-list()
       taxon2<-c(strsplit(taxon1,"/"))[[1]]
       for (k in taxon2) {
         index1<-grep(k,unionlist)
         if (length(index1)<1) unionlist<-append(unionlist,k,after=length(unionlist))
         taxon3<-c(strsplit(taxon,"/"))[[1]]
         for (j in taxon3) {
           if (k==j) {
             match<-TRUE # for the separation coefficient
             intersection_list<-append(intersection_list,k,after=length(intersection_list))
             }
           else {
             index2<-grep(j,unionlist)
             if (length(index2)<1) unionlist<-append(unionlist,j,after=length(unionlist))
             }
           }
         }
#------------------------------------------------------------------------------------------------------
       write(paste("unionlist=",unionlist), file=detailfile,append=TRUE)
       write(paste("intersection_list=",intersection_list), file=detailfile,append=TRUE)
       Jaccard_coefficient<-length(intersection_list)/length(unionlist)
       cat("Jaccard_coefficient= ",Jaccard_coefficient, "\n\n", file=detailfile,append=TRUE)
       Jaccard_similarity<-Jaccard_similarity+Jaccard_coefficient
#------------------------------------------------------------------------------------------------------
       if (match==FALSE) numdiffs<-numdiffs+1 #increment numdiffs
       }
   if (length(totaltaxa)>1) {   
      taxon1<-totaltaxa[1] #pop off the top line (data for the first taxon)
      cat("popped taxon=", taxon1, "\n", file=detailfile,append=TRUE)
      remainingtaxa<-totaltaxa[2:length(totaltaxa)]
      totaltaxa<-remainingtaxa
      }
   else {if (length(remainingtaxa)==1) {
      taxon1<-totaltaxa[1] #pop off the top line (data for the first taxon)
      cat("popped taxon=", taxon1, "\n", file=detailfile,append=TRUE)
      remainingtaxa<-list()
      totaltaxa<-remainingtaxa
      }}
}

#------------------------------------------------------------------------------------------------------
#  separation coefficient = total separable pairs/total possible pairs
#------------------------------------------------------------------------------------------------------
separation_coefficient<-numdiffs/total_num_pairs
cat("\nSeparation coefficient=",separation_coefficient, "\n", file=detailfile,append=TRUE)
cat("\nSeparation coefficient=",separation_coefficient, "\n", file=logfile,append=TRUE)

#------------------------------------------------------------------------------------------------------
#  Average pairwise Jaccard distance
#------------------------------------------------------------------------------------------------------
Jaccard_distance<- 1-(Jaccard_similarity/total_num_pairs)
cat("Jaccard similarity=",Jaccard_similarity, file=detailfile, "\n",append=TRUE)
cat("Average pairwise Jaccard distance=",Jaccard_distance, file=detailfile, "\n",append=TRUE)
cat("Average pairwise Jaccard distance=",Jaccard_distance, "\n", file=logfile,append=TRUE)

#logfile.write("\n")

close(f)
close(logfile)
close(detailfile)

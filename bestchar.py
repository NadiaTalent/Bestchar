#                                  Author: Nadia Talent, 2010 
#                 The author waives copyright and places this work in the Public Domain.
#                        http://creativecommons.org/publicdomain/zero/1.0/
#
#                                     Bestchar source code
#                     version in the Python language (http://www.python.org/) 
#
# A program for calculating some representative 'Best character' coefficients for
# categorial data (the same procedure would be used for numeric ranges
# after they are converted to categories). A single character is processed.
#
#------------------------------------------------------------------------------------------------------
# The input file called charinput.txt is a single column representing one character. 
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
# and bestchardetails.txt
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
# 2. The 'Normalized information coefficient' as described by 
#    Talent, Dickinson, and Dickinson (submitted); 
#    logarithms to base t are used.
#
# Coefficients that use pairwise comparisons of taxa, i.e. distance measures:
# 1. The Separation Coefficient counts only taxa that are completely separable, so if any state
#    is possible for two taxa they are not separable
# 2. Pairwise average Jaccard Coefficient

#------------------------------------------------------------------------------------------------------
#  set up to use real numbers with two decimal places, as well as logarithms
#------------------------------------------------------------------------------------------------------
from decimal import *
getcontext().prec=16
import math

#------------------------------------------------------------------------------------------------------
#  define the files
#------------------------------------------------------------------------------------------------------
f=open('charinput.txt','U')
logfile=open('bestcharoutput.txt','w')
detailfile=open('bestchardetails.txt','w')

detailfile.write('Bestchar program, Python version, diagnostic output\n\n')
logfile.write('Output from the Bestchar program, Python version\n\n')

#------------------------------------------------------------------------------------------------------
#  input the data
#------------------------------------------------------------------------------------------------------
line=f.readline() # the name of the character is the first line of the input file, with an EOL
detailfile.write(line)
logfile.write(line) # just echo that first line

# step 1, Read in the data for all taxa (a single character) and store them as a list (queue, fifo)

totaltaxa=[]
for line in f:
   detailfile.write('input line='+line)
   logfile.write(line) # echo the input line
   line1=line.split('\n',1) # split off the EOL
   line2=line1[0] # pick up the line minus the EOL
   line3=line2.replace('\t','') # remove tab characters
   line4=line3.replace(' ','') # remove spaces
   elements=line4.split('/') # split up the entries by the slashes, if any
   totaltaxa.append(elements) # add the line for this taxon to the whole list

numtaxa=len(totaltaxa) # called t in the comments above
detailfile.write('\ntotaltaxa='+repr(totaltaxa)+', length='+repr(numtaxa)+'\n\n')
logfile.write('\nNumber of taxa='+repr(numtaxa)+'\n')

#------------------------------------------------------------------------------------------------------
# For the Information Statistic, build three parallel lists:
# statelist: a list of the possible states i, is used only during the first pass through the data
# listofstatecounts: a list of tallies of taxa that allow each state (=k{i})
# listofsigmas: a list of the sums for each state i, of 1/n{j} where j is the number of states
# allowed for each taxon that allows state i
#------------------------------------------------------------------------------------------------------
statelist=[]
listofstatecounts=[] # these are called k{i} above, in H{i}=log2(k{i}/t)*SIG{i}/t)
listofsigmas=[] # these are called SIG{i} above
for taxon2 in totaltaxa:
   numstates=len(taxon2) #the number of states for this taxon
   for l in taxon2:
      index1=statelist.count(l) #returns 1 if already present, zero if not
      detailfile.write('index1='+str(index1)+'\n')
      if index1==0: # the state is not already present
         statelist.append(l)
         detailfile.write('appending new state='+l+'\n')
         listofstatecounts.append(1) #initialize the counter to 1 (# of taxa that have this state)
         listofsigmas.append(float(1)/float(numstates)) #initialize the SIG{i}
         detailfile.write('numstates='+str(numstates)+' set sigma to '+
            str(float(1)/float(numstates))+'\n')
      else: # this state was previously found
         index2=statelist.index(l)
         temp=listofstatecounts[index2]
         detailfile.write('incrementing at index2='+str(index2)+' (state='+l+') count was'+temp+'\n')
         listofstatecounts[index2]=temp+1 # increment the counter of taxa that have this state
         temp2=listofsigmas[index2]+float(1)/float(numstates) # add 1/(the number of states for that taxon)
         listofsigmas[index2]=temp2
         detailfile.write('numstates='+str(numstates)+' added='+str(float(1)/float(numstates))+
            ' to sigma giving'+temp2+'\n')

detailfile.write('\nstatelist='+repr(statelist)+'\n')
detailfile.write('listofstatecounts ='+str(listofstatecounts)+'\n')
detailfile.write('listofsigmas ='+str(listofsigmas)+'\n')

Totalnumstates=len(statelist)
detailfile.write('\nTotalnumstates ='+str(Totalnumstates)+'\n')
logfile.write('\nNumber of states='+repr(Totalnumstates)+'\n\n')

IntKeyH=float(0)
PankhurstH=float(0)
H=float(0)
mindex=0 #index for both lists
for m in listofstatecounts:
   detailfile.write('m='+str(m)+'\n')
   mtemp=float(m)/float(numtaxa)
   m1=math.log(mtemp,2) #log2(k{i}/t)
   m2=math.log(mtemp,Totalnumstates) #log base Totalnumstates
   m3=math.log(mtemp,numtaxa) #log base numtaxa
   detailfile.write('m1='+str(m1)+'\n')
   detailfile.write('m2='+str(m2)+'\n')
   detailfile.write('m3='+str(m3)+'\n')
   n=listofsigmas[mindex] # the equivalent entry in listofsigmas
   mindex+=1 #increment mindex
   n1=float(n)/float(numtaxa)
   detailfile.write('n1='+str(n)+'/'+str(numtaxa)+'='+str(n1)+'\n')
   IntKeyH=IntKeyH+n1*m1 #base 2 logarithms, as in IntKey
   detailfile.write('IntKeyH='+str(IntKeyH)+'\n')
   PankhurstH=PankhurstH+n1*m2 #base numstates, as recommended by Pankhurst
   detailfile.write('PankhurstH='+str(PankhurstH)+'\n')
   H=H+n1*m3 #base numtaxa logarithms (to normalize the range as the number of taxa (t) decreases)
   detailfile.write('H='+str(H)+'\n\n')

IntKeyH1=-Decimal(str(IntKeyH)).quantize(Decimal('.01'), rounding=ROUND_HALF_DOWN)
detailfile.write('Results to two decimal places:\n')
detailfile.write('Intkey-style information coefficient='+str(IntKeyH1)+'\n')
logfile.write('Results to two decimal places:\n')
logfile.write('Intkey-style information coefficient='+str(IntKeyH1)+'\n')

PankhurstH1=-Decimal(str(PankhurstH)).quantize(Decimal('.01'), rounding=ROUND_HALF_DOWN)
detailfile.write("Pankhurst's information coefficient="+str(PankhurstH1)+'\n')
logfile.write("Pankhurst's information coefficient="+str(PankhurstH1)+'\n')

H1=-Decimal(str(H)).quantize(Decimal('.01'), rounding=ROUND_HALF_DOWN)
detailfile.write('Normalized information coefficient='+str(H1)+'\n\n')
logfile.write('Normalized information coefficient='+str(H1)+'\n\n')

#------------------------------------------------------------------------------------------------------
# For Jaccard and Separation coefficients
# calculate number of possible pairs of taxa n!/2*(n-2)!=n*(n-1)/2
#------------------------------------------------------------------------------------------------------
detailfile.write('numtaxa='+str(numtaxa)+'\n')
total_num_pairs=numtaxa*(numtaxa-1)/2
detailfile.write('total_num_pairs='+str(total_num_pairs)+'\n\n')

#------------------------------------------------------------------------------------------------------
# For Jaccard and Separation coefficients
# destructively compare taxon character-states from the top row to each of the other rows (taxa)
# taxon1 holds the current top-most taxon
#------------------------------------------------------------------------------------------------------
taxon1=totaltaxa.pop(0) #pop off the data for the first taxon
detailfile.write('popped taxon='+repr(taxon1)+'\n')
numdiffs=0
Jaccard_similarity=float(0)

while len(totaltaxa)>0: # no need to match anything with the last taxon in the queue, but see 'else' suite
   detailfile.write('remaining taxa at start of while loop='+repr(totaltaxa)+'\n')
   for taxon in totaltaxa:
#------------------------------------------------------------------------------------------------------
#  The separation coefficient counts only taxa that are completely separable, so if any element of
#  taxon1's list is included in taxon's list, they don't count as separable
#------------------------------------------------------------------------------------------------------
       match=False
#------------------------------------------------------------------------------------------------------
#  Average pairwise Jaccard distance uses the count of elements (character states) in the intersection 
#  (between two taxa) divided by the count of elements in the union. Real numbers are required.
#------------------------------------------------------------------------------------------------------
       unionlist=[]
       intersection_list=[]
       for k in taxon1:
         if unionlist.count(k)==0: #not present
            unionlist.append(k)
         for j in taxon:
           if k==j:
             match=True # for the separation coefficient
             intersection_list.append(k)
           else:
             if unionlist.count(j)==0:
                unionlist.append(j)
#------------------------------------------------------------------------------------------------------
       detailfile.write('unionlist= '+str(unionlist)+'\n')
       detailfile.write('intersection_list= '+str(intersection_list)+'\n')
       Jaccard_coefficient=float(len(intersection_list))/float(len(unionlist))
       detailfile.write('Jaccard_coefficient= '+str(Jaccard_coefficient)+'\n\n')
       Jaccard_similarity=Jaccard_similarity+Jaccard_coefficient
#------------------------------------------------------------------------------------------------------
       if not(match): # separation coefficient
         numdiffs+=1 #increment numdiffs
   taxon1=totaltaxa.pop(0) #pop off the next line (the next least recently queued)
   detailfile.write('popped taxon='+repr(taxon1)+'\n')
else:
   detailfile.write('\ntotal unmatched pairs='+str(numdiffs)+'\n')

#------------------------------------------------------------------------------------------------------
#  separation coefficient = total separable pairs/total possible pairs
#------------------------------------------------------------------------------------------------------
separation_coefficient=float(numdiffs)/float(total_num_pairs)
detailfile.write('\nSeparation coefficient='+str(separation_coefficient)+'\n')
separation_coefficient1=-Decimal(str(separation_coefficient)).quantize(Decimal('.01'), rounding=ROUND_HALF_DOWN)
logfile.write('Separation coefficient='+str(separation_coefficient1)+'\n')

#------------------------------------------------------------------------------------------------------
#  Average pairwise Jaccard distance
#------------------------------------------------------------------------------------------------------
detailfile.write('Jaccard similarity='+str(Jaccard_similarity)+'\n')
Jaccard_distance=1-(Jaccard_similarity/float(total_num_pairs))
detailfile.write('Average pairwise Jaccard distance='+str(Jaccard_distance))
Jaccard_distance1=-Decimal(str(Jaccard_distance)).quantize(Decimal('.01'), rounding=ROUND_HALF_DOWN)
logfile.write('Average pairwise Jaccard distance='+str(Jaccard_distance1))

logfile.write('\n')

f.close()
logfile.close()
detailfile.close()

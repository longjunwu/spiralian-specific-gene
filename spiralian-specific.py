# Longjun Wu (longjun.wu@yale.edu)
# This script is to analyze BLAST results of genes from a query (spiralian) genome against target genomes
# and identify genes that do not have any hit in outgroup, while have hits in all three spiralian phyla (spiralian-specific genes).
# The input is a file containing BLAST results of a query spiralian genome (such as oyster, Crassostrea gigas)against target genomes
# The output file contains a list of gene ID: genes in the query genome that are identified as spiralian specific genes
# requirement: input file in the same directory
# use python 2.7


# Note 1: E eavlue threshould can be adjusted
# Note 2: An example of blast results (input file) is attached in the end  


import sys
import os
import time


# set input file name (a file containing blast results)
input_file=raw_input('please enter the input file name \n') 
# set output file name 
output_file=raw_input('please enter the final output file name \n')
# set e value thredhold (This threshold is for hit to target spiralian genomes only)
e_value=raw_input('please enter the e value threshold \n')

# open input file and store content of input file in "all"
all=open(input_file,"r")

# split blast results of each query gene in the query genome, and store results from each query gene in a list "all_split"
# the blast results of each query gene is seperated by '&>'
all_split=all.split('&>')[1:]

# set up output file 
output_file_handler = open (output_file, "a")

# loop through blast results of each query gene 
# analyze one query gene at a time

for one_query in all_split: 

# start to analyze the blast result for this query gene  
# the purpose is to figure out in what target genomes does this query gene have blast hit 
# set variables to store whether the the query gene has hit in certain target genome

# varialble k stores whether the query gene has hit to outgroups (non-spiralian genomes). 
# Initially, set variable k=0, when a hit is found in any of the outgroup, k will be changed to 1.  
 k=0

# Intially, set the following variables for each spiralian genome to 1, when a hit is found, the variable will be changed to 0.

 #lgi refer to Lottia gigantea
 lgi=1
 #cte refer to Capitella teleta
 cte=1 
 # hro refer to Helobdella robusta 
 hro=1
 # sme refer to Schmidtea mediterranea
 sme=1
 # emu refer to Echinococcus multilocularis
 emu=1
 # hmi refer to Hymenolepis microstoma
 hmi=1
 # sma refer to Schistosoma mansonia
 sma=1
 

# in the blast results, all result elements are seperated by '@'
# the first element stores the equery gene id, thus query gene id is acquired by:
 query_id=one_query.split('@')[0]

# other than the first element (query id), one element stores result of blast search to one target genome
# store these hits into a list "query_results"
 query_results=one_query.split('@')[1:]

# process result of blast search to one target genome per time
# the results are called "hit" here

 for hit in query_results:
 
# in the blast result format we use, the result only store the best hit (only one best hit per target genome)  
# in the blast result format we have, if no hit if found , '|' will not show up in the blast result
# thus if we don't find '|' in the blast result, no hit is found in the target genome  
  
  if hit.count('|')<1:
# if no hit is found, then we don't need to further analyze the blast result of that target genome (also do not need to change the value of variables)
# so we continue to the next hit   
   continue 

# blast results pass to this step has a hit to the target genome (otherwise will be passed and switched to other hit because of "continue")
# then we need to process the hit

# the target genome name of the hit is the first element seperated by '\n' 
# we store the name of the target genome for that hit "seq_species"
# make the variable empty before assigning variable, so it is empty each time before assigned, or it could cause issues
  seq_species='empty'
  seq_species=hit.split('\n')[0]
 
# next we need to know is the e value of the hit
# if e value > e value thredhold, this is not considered as "significant" hit, in other words, the hit is too weak to be considered as a real hit. So we contiune to the next hit
# get the e value and store it in "hit_e"
# in the blast result format we have, the e value is stored in the fourth from the last element seperated by "|" 
  hit_e=hit.split('|')[-4] 

# then we need to see if the query gene has hit to non-spiralian outgroup
# here are the non-spiralian outgroups we use in the study  
  
  #cel: Caenorhabditis elegans
  #ame: Apis mellifera
  #ath: Arabidopsis thaliana
  #cin: Ciona intestinalis
  #ddi: Dictyostelium discoideum
  #dme: Drosophila melanogaster
  #hma: Hydra magnipapillata
  #mbr: Monosiga brevicollis
  #mmu: Mus musculus
  #xla: Xenopus laevis
  
# note: if one (or more） significant hit is found in one (or more) genomes from a phylum, this query gene is considered existing in the phylum 
  
# if a hit (with e value smaller than the thredhold) is found in any of the outgroup target genome, change the value of k to 1. 
  if (seq_species=='cel'  or  seq_species=='ame'  or seq_species=='ath' or seq_species=='cin' or seq_species=='ddi' or seq_species=='dme'  or seq_species=='hma'  or seq_species=='mbr' or seq_species=='mmu'  or seq_species=='xla' ) and (float(hit_e)<0.00001): # note the e value threshould to outgroup can be adjusted here

# when hit is found in any outgroup, k is changed to 1. If no hit in outgroup, k is always 0
   k=1

 
# analyze the hits of query gene to spiralian genomes
# if signficant hit in found the targeting spiralian genome, the variable for that genome will be change to 0. Otherwise it always is 1.  
  if seq_species=='lgi'and (float(hit_e)<float(e_value)):
   lgi=0
  if seq_species=='cte'and (float(hit_e)<float(e_value)):
   cte=0
  if seq_species=='hro'and (float(hit_e)<float(e_value)):
   hro=0
  if seq_species=='sme'and (float(hit_e)<float(e_value)):
   sme=0
  if seq_species=='emu'and (float(hit_e)<float(e_value)):
   emu=0
  if seq_species=='hmi'and (float(hit_e)<float(e_value)):
   hmi=0
  if seq_species=='sma'and (float(hit_e)<float(e_value)):
   sma=0
  
# then use the value of these variables to find genes that fit our spirlaian specific gene criteria

# for spiralian specific genes: 

# first requirement: k==0 means the query gene has not hit at all in the outgroup genome

# next we need to know if the query gene has hit to all three spiralian phyla

# second requirement:  cte or hro ==0 means the query gene has at least one hit in annelid phylum 
# third requirement: sme==0 or emu==0 or hmi==0 or sma==0 means the query has at least one hit in flatworm phylum 
# so the query gene has hits in all three spiralian phyla examined (mollusc, annelid and flatworm)
# because our query species is oyster, a mollusc, hit to other mollusc is not need to show the gene exist in mollusc

# if the query fit this requirement, this query gene is considered as spiralian-specific gene 
 if  k==0 and ( (hro==0 or cte==0)  and (sme==0 or emu==0 or hmi==0 or sma==0) ) : 
  
# write the query id into the output file 

   output_file_handler.write(query_name + '\n')
 
  
  
'''blast result example:



#each query (gene) starts with "&>", followed by gene ID
#each element (seperated by '@')shows the result of a blast search against one target genome
#if hit is found, best hit is followed in the next line (two line for the element)
#best hit include gene ID, e value, score, sperated by "|", and ending with "%" 
#if no hit is found nothing is followed by the target genome symbol, (only one line for the element) 

&>CGI_10000780
@cte
gnl|BL_ORD_ID|29702 jgi|Capca1|191939|fgenesh1_pg.C_scaffold_805000002|1.83571e-27|278.0|1|%
@hro
gnl|BL_ORD_ID|17407 jgi|Helro1|178217|2.46349e-05|106.0|1|%
@lgi
gnl|BL_ORD_ID|10108 jgi|Lotgi1|157147|fgenesh2_pg.C_sca_13000081|6.68925e-64|552.0|1|%
@cgi
gnl|BL_ORD_ID|0 CGI_10000780|0.0|2427.0|1|%
@cel
@ame
@ath
@cin
@ddi
@dme
@hma
gnl|BL_ORD_ID|12815 Hma2.213275|PACid:15579114|5.15202e-26|277.0|2|%
@mbr
@mmu
@xla

&>CGI_10000456
@cte
gnl|BL_ORD_ID|9180 jgi|Capca1|175568|estExt_Genewise1Plus.C_3690031|1.63874e-63|493.0|1|%
@hro
gnl|BL_ORD_ID|7161 jgi|Helro1|71670|2.81618e-39|329.0|1|%
@lgi
gnl|BL_ORD_ID|17494 jgi|Lotgi1|189273|estExt_Genewise1.C_sca_280271|4.08298e-69|529.0|1|%
@cgi
gnl|BL_ORD_ID|1 CGI_10000456|7.17874e-103|753.0|1|%
@cel
gnl|BL_ORD_ID|19240 T27A3.6#CE14227#WBGene00020842#status:Partially_confirmed#UniProt:P91500#protein_id:CCD72006.1|1.71202e-36|326.0|1|%
@ame
gnl|BL_ORD_ID|3375 gi|66564230|ref|XP_397469.2| PREDICTED: molybdopterin synthase catalytic subunit [Apis mellifera]|1.9114e-60|490.0|1|%
@ath
gnl|BL_ORD_ID|7371 AT2G43760.1 | Symbols:  | molybdopterin biosynthesis MoaE family protein | chr2:18133276-18133872 FORWARD LENGTH=198|1.32042e-41|352.0|1|%
@cin
gnl|BL_ORD_ID|900 jgi|Cioin2|241841|estExt_genewise1.C_chr_07q0325|7.6029e-49|394.0|1|%
@ddi
gnl|BL_ORD_ID|6155 DDB0267170|DDB_G0271864 |Protein|gene: mocs2l on chromosome: 2 position 980070 to 980546|2.78331e-45|372.0|1|%
@dme
gnl|BL_ORD_ID|20450 FBpp0084157 pep:known chromosome:BDGP5:3R:20951354:20953081:-1 gene:FBgn0039280 transcript:FBtr0084782 gene_biotype:protein_coding transcript_biotype:protein_coding|3.15175e-49|416.0|1|%
@hma
gnl|BL_ORD_ID|8566 Hma2.208824|PACid:15574865|2.93574e-42|351.0|1|%
@mbr
gnl|BL_ORD_ID|2361 jgi|Monbr1|14411|e_gw1.3.549.1|1.0183e-44|367.0|1|%
@mmu
gnl|BL_ORD_ID|22628 ENSMUSP00000015680 pep:known chromosome:GRCm38:13:114818277:114828723:1 gene:ENSMUSG00000015536 transcript:ENSMUST00000015680 gene_biotype:protein_coding transcript_biotype:protein_coding|6.2952e-61|482.0|1|%
@xla
"



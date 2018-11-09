import pandas as pd
import numpy as np
import re
import argparse
from collections import defaultdict
from collections import OrderedDict
from itertools import chain
import datetime

################################### ADD OPTIONS ##########################################################
parser = argparse.ArgumentParser(description="Converts anvi'o ouput to VCF format")
parser.add_argument('-i','--input', help='Input file : Output from anvio', required=True)
parser.add_argument('-o','--output', help='VCF output file name', required=True)
args = parser.parse_args()
anviTabInput = args.input
#VCFoutput = args.output
VCFoutput = open(args.output, "w")
anviTab = pd.read_table(anviTabInput)
#######################   SELECTING COLUMNS OF INTEREST ##################################################
anviTab=anviTab.filter(items=['unique_pos_identifier','split_name', 'pos','sample_id','coverage','reference','competing_nts'])
##################### ADDING ALLELE COLUMNS FOR FINDING GENOTYPE #########################################
anviTab['competing_nts'] = anviTab.competing_nts.astype(str)
anviTab['allele1'] = anviTab.competing_nts.str[0]
anviTab['allele2'] = anviTab.competing_nts.str[1]
anviTab["allele2"], anviTab["allele1"] = np.where(anviTab['allele2']==anviTab['reference'], [anviTab["allele1"], anviTab["allele2"]], [anviTab["allele2"], anviTab["allele1"] ])


genotype = defaultdict(dict)
alt_alleleDict = defaultdict(list)

for index, row in anviTab.iterrows() :
    key = row['unique_pos_identifier']
    Sample = row['sample_id']
    Ref_allele = row['reference']
    a1 = row['allele1']
    a2 = row['allele2']

    if a1 not in alt_alleleDict[key] and a1 != row['reference']:
        alt_alleleDict[key].append(a1)
    if a2 not in alt_alleleDict[key] and a2 != row['reference']:
        alt_alleleDict[key].append(a2)
        
    if row['reference'] in alt_alleleDict[key]:
        alt_alleleDict[key].remove( row['reference'])
 ## SET UP A GENOTYPE dict.       
    if a2 in row['reference']:
        genotype[key][Sample]='0/0'
    elif a1 in row['reference']:
        genotype[key][Sample]='0/n'
    else:
        genotype[key][Sample]='p/q'

##########################################################################################################
## FOR KEEPING TRACK OF THE IDS WHILE PRINTING DIRECTLY, not used when dictionary printed 
current_ID=[]
finalVCF = defaultdict(dict)
sampleNames=list(OrderedDict.fromkeys(anviTab['sample_id'])) ## To maintain order of the sample 
header=["#CHROM" ,"POS", "ID", "REF", "ALT" ,"QUAL" ,"FILTER", "INFO","FORMAT"]
##########################################################################################################
for index, row in anviTab.iterrows():
  #  chromosome = row['#CHROM']
    key = row['unique_pos_identifier'] 
    Sample = row['sample_id']
    Ref_allele = row['reference']
    a1 = row['allele1']
    a2 = row['allele2']
    if genotype[key][Sample] == '0/n':
        genotype[key][Sample] = re.sub('n',  str(alt_alleleDict[key].index(a2)+1), genotype[key][Sample])    
      
    elif genotype[key][Sample] == 'p/q':
        genotype[key][Sample] = re.sub('p',  str(alt_alleleDict[key].index(a1)+1), genotype[key][Sample])
        genotype[key][Sample] = re.sub('q',  str(alt_alleleDict[key].index(a2)+1), genotype[key][Sample])

    sampleCol=str(row['coverage'])+':'+genotype[key][Sample]
    #PRINT  DIRECTLY  
    ''' 
    if key in current_ID:
        print(sampleCol,end="\t")    
    else:
        print()
        print(row['split_name'],row['pos'],key,Ref_allele,','.join(alt_alleleDict[key]),99,'PASS','.','DP:GT',sep="\t",end="\t")
        current_ID.append(key)
        print(sampleCol,end="\t")
    '''
    
    if key in finalVCF:
        TableFields.append(sampleCol)
        finalVCF[key] = TableFields
    else :
        TableFields = [row['split_name'],row['pos'],key,Ref_allele,','.join(alt_alleleDict[key]),99,'PASS','.','DP:GT',sampleCol]
        finalVCF[key] = TableFields
######################################### PRINT THE VCF ################################################
###### PRINT HEADERS AND INFO

now = datetime.datetime.now()
VCFoutput.write('##fileformat=VCFv4.0' + '\n')
VCFoutput.write('##fileDate='+now.strftime("%Y%m%d") + '\n')
VCFoutput.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + '\n')
VCFoutput.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">'+ '\n')
VCFoutput.write('\t'.join(header+sampleNames) +'\n')

###### sort?  
sortedfinalVCF = sorted(finalVCF.values())
for i in range(0, len(sortedfinalVCF)):
	VCFoutput.write('\t'.join([str(x) for x in sortedfinalVCF[i]]) + '\n')

   

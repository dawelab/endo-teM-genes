#This script filters a vcf file produced by CGmapTools snv to remove variants with extreme coverage values and non-heterozygous genotypes


##########################################################################################################
#USER INPUT 
##########################################################################################################

#directory of input files
directory = '/scratch/gent/imprinting'
 
#input file endings
ending = '.vcf'

#min depth 
minDP = 15

#max depth 
maxDP = 70

##########################################################################################################
#initialize variables and open files
##########################################################################################################

import gzip
import glob
import os
os.chdir(directory) 

for filename in glob.glob(directory + '/*' + ending):
	inFile = open(filename)
	outFile = open(filename[:-len(ending)] + '_' + 'filtered' + ending, 'w')
	
	
	
	##########################################################################################################
	#Run through entire file and filter each line
	##########################################################################################################
	

	for line in inFile:
		cols = line.strip().split('\t')
		
		try:
			if len(cols[7].split(':')) == 2: #skip lines with more than one alternative genotype.
				DP = int(cols[7].split(':')[-1].split('=')[-1])
				if DP >= minDP:
					if DP <= maxDP: #skip lines with extreme depths
						if cols[9].split(':')[0] == '0/1':  #skip lines with non-heterozygous genotypes.
							outFile.write(line)
		except IndexError:
			outFile.write(line)

	##########################################################################################################
	#close files
	##########################################################################################################
		
	inFile.close()
	outFile.close()

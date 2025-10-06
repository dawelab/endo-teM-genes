#This script converts CGmapTools Allele-Specific methylated Site files (.ass) to .CGmap.gz files


##########################################################################################################
#USER INPUT 
##########################################################################################################

#directory of input files
directory = '/scratch/gent/imprinting/ASSs'
 
#input file endings
ending = '.CG2.ass'

#cytosine context
context = 'CG'

#parental genomes
Allele1 = 'W22'
Allele2 = 'B73'

##########################################################################################################
#initialize variables and open files
##########################################################################################################

import gzip
import glob
import os
os.chdir(directory) 

for filename in glob.glob(directory + '/*' + ending):
	inFile = open(filename)
	CGmap1File = gzip.open(filename[:-len(ending)] + '_' + Allele1 + '_' + context + '.CGmap.gz', 'wb')
	CGmap2File = gzip.open(filename[:-len(ending)] + '_' + Allele2 + '_' + context + '.CGmap.gz', 'wb')
	
	
	
	##########################################################################################################
	#Run through entire input file converting each line to two CGmap lines
	##########################################################################################################

	#Check file format, and skip first header line
	print(inFile.readline())
	
	#Print info from first line
	cols = inFile.readline().strip().split('\t')
	chr = cols[0]
	strand = 'C' #This is fake. It could actually be a G, since ASS file doesn't give strand info.
	pos = int(float(cols[5])) #C position, not SNP position. Must be converted to int because ASS format sometimes uses scientific notation, which will cause trouble in CGmaps.
	meth1 = cols[8]
	methylC1 = cols[6].split('-')[0]
	totalC1 = str(int(cols[6].split('-')[0]) + int(cols[6].split('-')[1]))
	meth2 = cols[9]
	methylC2 = cols[7].split('-')[0]
	totalC2 = str(int(cols[7].split('-')[0]) + int(cols[7].split('-')[1]))
	CGmap1File.write((chr + '\t' + strand + '\t' + str(pos) + '\t' + context + '\tunknown' + '\t' + meth1 + '\t' + methylC1 + '\t' + totalC1).encode())
	CGmap2File.write((chr + '\t' + strand + '\t' + str(pos) + '\t' + context + '\tunknown' + '\t' + meth2 + '\t' + methylC2 + '\t' + totalC2).encode())

	#for all other lines in file write '\n' first
	for line in inFile:
		cols = line.strip().split('\t')
		if cols[6].split('-')[0] != 'Allele1_linked_C': #skip header line
			chr = cols[0]
			strand = 'C' #This is fake. It could actually be a G, since ASS file doesn't give strand info.
			pos = int(float(cols[5])) #C position, not SNP position. Must be converted to int because ASS format sometimes uses scientific notation, which will cause trouble in CGmaps.
			meth1 = cols[8]
			methylC1 = cols[6].split('-')[0]
			totalC1 = str(int(cols[6].split('-')[0]) + int(cols[6].split('-')[1]))
			meth2 = cols[9]
			methylC2 = cols[7].split('-')[0]
			totalC2 = str(int(cols[7].split('-')[0]) + int(cols[7].split('-')[1]))
			CGmap1File.write(('\n' + chr + '\t' + strand + '\t' + str(pos) + '\t' + context + '\tunknown' + '\t' + meth1 + '\t' + methylC1 + '\t' + totalC1).encode())
			CGmap2File.write(('\n' + chr + '\t' + strand + '\t' + str(pos) + '\t' + context + '\tunknown' + '\t' + meth2 + '\t' + methylC2 + '\t' + totalC2).encode())

	##########################################################################################################
	#close files
	##########################################################################################################
		
	inFile.close()
	CGmap1File.close()
	CGmap2File.close()

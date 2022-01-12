from pandas import DataFrame
import logging
import sys

from numpy import arange, histogram
import matplotlib.pyplot as plt

import argparse

def read_vcf(fh):
	'''
	Read the VCF file obtained from any program included into Parliment2. Adds columns to the records if they are lacking.

	Args:
		fh (file): a VCF file.

	Returns:
		DF (pandas.DataFrame)   : It contains all the records of SV
		metadata (array)        : It contains all the annotation lines begining with '##'
	'''
	metadata = []
	DF = None

	for i in fh:
		if(i[0:2]=='##'):
			metadata.append(i.strip())
		elif(i[0:6]=='#CHROM'):
			DF = DataFrame( columns=i.strip().split() )
		else:
			try:
				DF = DF.append( DataFrame( [[j for j in i.strip().split()]] , columns= list(DF.columns) ) )
			except:
				new_columns = [ 'col'+str(j+1) for j in range(len(DF.columns), len(DF.columns)+len(i.split())-len(DF.columns) ) ]
				for j in new_columns: DF[j] = None
				DF = DF.append( DataFrame( [[j for j in i.strip().split()]] , columns= list(DF.columns) ) )
				logging.warning('Extra column(s) were added as they were missing in the generated VCF file.')

	return DF, metadata


def analyze_chromosome( chr ,records, metadata , minlen=0 ,bins=100, output_figure='output.png'):
	'''
	Analyze chromosome function divides the chromosome into 'n' number of the equal-sized bin (parameter bins) and counts the different types of structural variants (SVs) to plot them as a figure. It allows the user to analyze what parts of the chromosomes are most affected by what type of variation.

	Args:
		chr (str)                       : Chromosome name user wishes to analyze.
		records (read_vcf output [0])   : Output obtained from read_vcf function.
		metadata (read_vcf output [01]) : Output obtained from read_vcf function.
		minlen (int)                    : Cutoff length of SV (Please note that SVTYPE BND does not have length defined)
		bins (int)                      : Number of parts user wishes to chop chromosomes into (read the description)
		output_figure (str)             : Name for the output figure.
	'''
	
	try:
		records = records.loc[ records['#CHROM']==chr ]
	except:
		logging.error('Given chromosome not found in the records.')

	chr_lengths = {}
	for i in metadata:
		try:
			temp = i.split(',')
			if(temp[1].split('=')[0]=='length'):
				chr_lengths[temp[0].split('=')[2]] = arange( 0 , float(temp[1].split('=')[1][:-1])+1 , float(temp[1].split('=')[1][:-1])/bins )
		except:
			None
	
	#SVTYPE: DEL, INV, DUP, UNK, BND, INS
	freq_data = { 'DEL':[], 'INV':[], 'DUP':[], 'UNK':[], 'BND':[], 'INS':[] }

	for index, row in records.iterrows():
		flag = False
		INFO = {i.split('=')[0]:i.split('=')[1] for i in row['INFO'].split(';') if len(i.split('=')) == 2 }
		try:
			if( minlen<=abs(int(INFO['SVLEN'])) ):
				flag = True
			if( flag or INFO['SVTYPE'] == 'BND' ):
				freq_data[ INFO['SVTYPE'] ].append(float(row['POS']))
				flag = False
		except:
			if( INFO['SVTYPE'] == 'BND' ):
				freq_data[ INFO['SVTYPE'] ].append(float(row['POS']))
				flag = False


	fig, ( ax0, ax1, ax2, ax3, ax4, ax5 ) = plt.subplots(6,figsize=(20,20))
	fig.suptitle('Chromosome SV locations for '+chr)
	
	ax0.hist(freq_data['DEL'],bins = chr_lengths[chr] )
	ax0.set_title('DEL', loc='left')
	ax1.hist(freq_data['INV'],bins = chr_lengths[chr] )
	ax1.set_title('INV', loc='left')
	ax2.hist(freq_data['DUP'],bins = chr_lengths[chr] )
	ax2.set_title('DUP', loc='left')
	ax3.hist(freq_data['BND'],bins = chr_lengths[chr] )
	ax3.set_title('BND', loc='left')
	ax4.hist(freq_data['INS'],bins = chr_lengths[chr] )
	ax4.set_title('INS', loc='left')
	ax5.hist(freq_data['UNK'],bins = chr_lengths[chr] )
	ax5.set_title('UNK', loc='left')

	fig.tight_layout()
	plt.savefig(output_figure)
	logging.info('Task completed successfully!. Please check the output image in the present working directory.')

def IO():
	parser=argparse.ArgumentParser(description='VCFanalysis.py for SV analysis')
	parser.add_argument('-vcf',type=str, help='Provide the name of the VCF file to parse.')
	parser.add_argument('-chr',type=str, help='Provide the name of the chromosome (eg... chr1)')
	parser.add_argument('-bins',type=int, help='Provide the number of bins chromosome need to be divided (Please red the description for more details; Default: 100 ie.. percentile)',default=100)
	parser.add_argument('-minlen',type=int, help='Provide the minimum length of the SV')
	args = parser.parse_args()
	return args

def main():
	#Example
	args = IO()

	fh = open(args.vcf)
	inchr = args.chr
	bins = args.bins
	minlen = args.minlen

	records , metadata = read_vcf(fh)

	analyze_chromosome( inchr ,records, metadata ,bins= bins, minlen= minlen )
	

if(__name__=='__main__'):
	main()
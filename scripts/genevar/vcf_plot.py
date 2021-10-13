# usage: python3 vcf_plot.py --path '/Users/priyalakra/Desktop/hack_oct/graph/final' --infile 'test.output.vcf'
# total time to process 1 vcf few seconds

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import os, sys
import argparse 
from time import strftime
pd.options.mode.chained_assignment = None



def open_vcf(filename):
    '''reads vcf files
    input: vcf file path
    returns: metadata and vcf information as separate files'''
    
    metadata = []
    data = []
    
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("##"): metadata.append(line.strip())
            elif line.startswith("#CHROM"): headers = [x for x in line.strip("##\n").split('\t')]
            else: data.append(line.split('\t'))
                
    data_final = pd.DataFrame(data, columns=headers)
    # drop extra columns 
    data_req = data_final.drop([ "ALT", "REF", "QUAL"], axis = 1) 
    
    return data_req, metadata


def extract_sv(infile, svtype):
    '''extract a type of SV specified by user
    
    inputs: VCF data, SV type
    
    returns: VCF file as dataframe containing only user specified SV type
    '''

    req_sv = infile.loc[infile["INFO"].str.contains(f"SVTYPE={svtype}")]
    req_sv["AF"] = req_sv["INFO"].str.extract(r"(AF=\d*\.\d*;|AF=0;)", expand=False) 
    
    req_sv = req_sv.dropna(subset=['AF'])
   # try:
    req_sv["allele freq"] = req_sv["AF"].str[3:-1].astype(float) #.astype(int)    
   # except ValueError: None
        
    return req_sv
    
    
def extract_chr(data_req, chrnum):
    
    req_chr = data_req.loc[ data_req['CHROM']==chrnum ]
    req_chr["sv"] = req_chr["INFO"].str.extract(r"(SVTYPE=[A-Z]{3})", expand=False)
    req_chr["POS"] = req_chr["POS"].astype(int)
    
    req_chr["AF"] = req_chr["INFO"].str.extract(r"(AF=\d*\.\d*;|AF=0;)", expand=False)  
    req_chr = req_chr.dropna(subset=['AF'])
    req_chr["allele freq"] = req_chr["AF"].str[3:-1].astype(float) #.astype(int)    
    req_chr["svtype"] = req_chr["sv"].str[7:]
   
    return req_chr


def edit_chr(req_sv):
    
    req_sv["chr"] = req_sv["CHROM"].str.extract(r"([0-9]+)", expand=False)
    req_sv["chr"] = req_sv["chr"].replace(['X', 'Y'], ['23', '24'])
    
    req_sv = req_sv.dropna(subset=['chr'])
    req_sv["chr"] = req_sv["chr"].astype(int)
    
    return req_sv


def edit_pos(file_new,position,chromosome):

    # add new position for plotting across x-axis
    positions=file_new[[chromosome, position]] 
    positions[position] = positions[position].astype(int)
    new_pos = []
    add=0
    for chro,posi in positions.groupby(chromosome):
        new_pos.append(posi[[position]]+add)
        add+=posi[position].max() # maximum position within each chromosome
        
    # append new positions to the data frame     
    file_new['new_pos'] = pd.concat(new_pos)
   
    return file_new


def bin_plot_all(svfile,col,chromosome,bins):
	'''divides each chromosome into n bins as specified by user and counts SV in that bin. Represent which
	   part of chromosome is affected by a particular SV
	'''
	
	fig = plt.figure(figsize=(15, 8))
	ax = fig.add_subplot(111)
	ax_label = []
	ax_pos = []
	leg = []
	if col == 'y': colors = itertools.cycle(['red','blue','green','black'])
	else: colors = itertools.cycle(['gray','black'])
	
	for (chrname,detail) in list(svfile.groupby(chromosome)):
	    ax.hist( x = detail["new_pos"], bins=bins, color=next(colors) )
	    ax_label.append(chrname)
	    ax_pos.append((detail['new_pos'].iloc[-1] + detail['new_pos'].iloc[0])/2)
	    
	ax.set_xticks(ax_pos)
	ax.set_xticklabels(ax_label,rotation='vertical', fontsize=14)
	ax.set_xlabel("Chromosomes", fontsize=18)
	ax.set_ylabel("SV Counts", fontsize=18)
	
	for (chrname,detail) in list(svfile.groupby(chromosome)): leg.append(f"chr {chrname}")
	
	ax.set_title("Histogram of DELETION SV across chromosomes", fontsize=20)
	ax.legend(leg, bbox_to_anchor=(1, 1), fontsize=10)
	
	plt.savefig(f'hist_plot.jpeg')
    
    
def freq_plot(svfile,col,chromosome):
    '''Plots allele frequency as given in VCF file for each SV across all chromosomes
    '''
    
    fig = plt.figure(figsize=(22,8))
    ax = fig.add_subplot(111)
    ax_label = []
    ax_pos = []
    leg = []
    if col == 'y': colors = itertools.cycle(['red','blue','green','black'])
    else: colors = itertools.cycle(['gray','black'])
    
    
    for (chrname,detail) in list(svfile.groupby(chromosome)): 
        detail.plot(kind = 'scatter', x = "new_pos", y = "allele freq", color = next(colors), ax=ax,  s=25) #, marker='s');
        ax_label.append(chrname)
        ax_pos.append((detail['new_pos'].iloc[-1] + detail['new_pos'].iloc[0])/2)
        
    ax.set_xticks(ax_pos)
    ax.set_xticklabels(ax_label,rotation='vertical', fontsize=18)
    ax.set_xlabel("Chromosomes", fontsize=22)
    ax.set_ylabel("Allele frequency", fontsize=22)
    
    for (chrname,detail) in list(svfile.groupby(chromosome)):
        leg.append(f"chr {chrname}")
    
    ax.set_title("DELETION SV", fontsize=18)
    ax.legend(leg, bbox_to_anchor=(1.1, 1), fontsize=12)
   
  #  plt.axhline(y = 0.45, color = "blue", linewidth = 0.5)
  
    labels = svfile[svfile["allele freq"] >= 0.5]
    name = []
    x = []
    y = []
    for i, j in labels.iterrows():
        name.append(j["ID"][0:-9])
        x.append(j["new_pos"])
        y.append(j["allele freq"])
    
    for i, label in enumerate(name):
        plt.annotate(label, (x[i], y[i]), fontsize=14)
        
    plt.savefig(f'freq_plot.jpeg')
        
        
def freq_plot_perchr(svfile,col,chromosome,svtype,chrnum):
	'''Plots allele frequency as given in VCF file for all SVs across chromosome specified by the user.
	'''
	
	fig = plt.figure(figsize=(10, 5))
	ax = fig.add_subplot(111)
	ax_label = []
	ax_pos = []
	if col == 'y': colors = itertools.cycle(['red','blue','green', 'black'])
	else: colors = itertools.cycle(['gray','black'])
	
	for (svtype,detail) in list(svfile.groupby(svtype)):
	    detail.plot(kind = 'scatter', x = "POS", y = "allele freq", color = next(colors), ax=ax,  s = 50) #, marker='s');
	    ax_label.append(svtype)
	    
	ax_pos.append(min(svfile["POS"]))
	ax_pos.append(int(min(svfile["POS"]) + max(svfile["POS"]))//2)
	ax_pos.append(max(svfile["POS"]))
	ax.set_xticks(ax_pos)
	ax.set_xticklabels(ax_pos,rotation='45')
	ax.set_xlabel("Chromosomal position", fontsize=14)
	ax.set_ylabel("Allele frequency", fontsize=14)
	ax.legend(ax_label, loc="upper left")
	ax.set_title(f"Structural variants on {chrnum}", fontsize=14)
	ax.set_ylim(-0.05,svfile["allele freq"].max()+0.2)
	
	labels = svfile[svfile["allele freq"] >= 0.5]
	name = []
	x = []
	y = []
	for i, j in labels.iterrows():
	    name.append(j["ID"][0:-9])
	    x.append(j["POS"])
	    y.append(j["allele freq"])
	    
	for i, label in enumerate(name):
	    plt.annotate(label, (x[i], y[i]), fontsize=10)
	    
	plt.savefig(f'freq_plot_perchr_plot.jpeg', bbox_inches = 'tight')
	
	
if __name__ == "__main__":

	# initialize your parser
    parser = argparse.ArgumentParser(description = "Script for designing a manhattan plot")   

	# parse the arguments
    parser.add_argument("--path", type=str, help="Path of the directory where your VCF files are stored")
    parser.add_argument("--infile", type=str, help="Name of the file containing SNPs/GWAS data to be analysed")
    #parser.add_argument("--chromosome",  type=str, help="How chromosomes are indicated in your input file")
    #parser.add_argument("--info",  type=str, help="How INFO fields are indicated in your input file")
    #parser.add_argument("--position", type=str, help="How chromosomal positions are indicated in your input file")
  #  parser.add_argument("--svtype", type=str, help="Which you want to analyze ")
    #parser.add_argument("--col", type=str, help="Do you want a colored graph? If yes, Write 'y' ")

    args = parser.parse_args()
    
    start_time = strftime("%m-%d-%Y %H:%M:%S")
    sys.stdout.write("Program started at " + start_time + '\n') # can also use print() here
    os.chdir(args.path)

    # call functions
    file, metadata = open_vcf(args.infile) 
    
    filename = args.infile.rsplit(".",1)[0]
    new_path = f'{args.path}/{filename}_result'
    try: os.makedirs(new_path)
    except OSError: print('directory exists')
    
    os.chdir(new_path)
        
    # graphs for all chromosomes
    dsv1 = extract_sv(file, svtype="DEL")
    dsv2 = edit_chr(dsv1)
    dsv3 = edit_pos(dsv2,position="POS",chromosome="chr")
    bin_plot_all(dsv3,col="y",chromosome="chr",bins=100)
    freq_plot(dsv3,col="y",chromosome="chr")
    
    
    # graph per chromosome
    dch1 = extract_chr(file, 'chr1')
    dch2 = edit_chr(dch1)
    freq_plot_perchr(dch2,col="y",chromosome="ch",svtype="svtype", chrnum="chr1")

    end_time = strftime("%m-%d-%Y %H:%M:%S")
    sys.stdout.write("Program ended at " + end_time + '\n')
    sys.exit(0)
    
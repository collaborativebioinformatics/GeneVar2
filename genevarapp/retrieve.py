import os
file = open("input_location.txt")
line = file.read()
file.close()
wholecommand = 'python VCFanalyse.py  -chr chr2 -vcf '+line+ ' -bins 100 -minlen 100'
print(wholecommand)
os.system(wholecommand)

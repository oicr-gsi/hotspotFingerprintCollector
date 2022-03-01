import gzip
import re
import argparse

parser = argparse.ArgumentParser(description='fincreator')
parser.add_argument('inputVcf', help='Input .bam file for analysis sample')
parser.add_argument('inputCoverage', help='Optional input .bam file for control sample')
parser.add_argument('hotspotSNPs',help='Hotspot vcf, used as a reference')
parser.add_argument('sampleID')
parser.add_argument('chroms')
args = parser.parse_args()

flags = {"zerocov": "N",  # No coverage
        "nosnp": "M"}  # no SNP, base matches the reference

# Read lines from the file with reference SNPs
refFile = gzip.open(args.hotspotSNPs, mode='r')
refLines = refFile.readlines()
refFile.close()

# Retrieve reference SNPsb"abcde".decode("utf-8")
refs = {}
for rLine in refLines:
    refLine = rLine.decode('utf-8')
    if refLine.strip().startswith('#'):
        continue
    temp = refLine.split("\t")
    if not refs.get(temp[0]):
        refs[temp[0]] = {}
    refs[temp[0]][temp[1]] = dict(flag=flags['zerocov'], dbsnp=temp[2], ref=temp[3])

# Read lines from the file with Depth data
depthFile = open(args.inputCoverage, mode='r')
depthLines = depthFile.readlines()
depthFile.close()

# Retrieve Depth values
for dLine in depthLines:
    if dLine.startswith("Target"):
        continue
    #temp = dLine.split("\t")
    temp = dLine.split(",")
    #coords = temp[0].split(":")
    coords = temp[0].replace("-",":").split(":")
    if len(coords) != 2:
        continue
    if refs[coords[0]].get(coords[1]) and int(temp[1]) > 0:
        refs[coords[0]][coords[1]]['flag'] = flags['nosnp']
    else:
        refs[coords[0]][coords[1]]['flag'] = flags['zerocov']

# Read lines from the file with Depth data
callFile = open(args.inputVcf, mode='r')
callLines = callFile.readlines()
callFile.close()

for cLine in callLines:
    if cLine.startswith('#'):
        continue
    temp = cLine.split("\t")
    if refs[temp[0]] and refs[temp[0]].get(temp[1]):
        variant = temp[4].upper()
        if re.fullmatch('[ACGT]{1}', variant) and temp[4] != temp[3]:
            refs[temp[0]][temp[1]]["flag"] = variant
        else:
            refs[temp[0]][temp[1]]["flag"] = flags['nosnp']

# Prepare lines for printing
finLines = []
for chr in [args.chroms]:
    if not refs.get(chr):
        continue
    for start in sorted(refs[chr].keys()):
        snp = [refs[chr][start]['ref']]
        if refs[chr][start]['flag'] == 'M':
            snp.append(refs[chr][start]['ref'])
        else:
            snp.append(refs[chr][start]['flag'])
        if refs[chr][start]['flag'] == 'N':
            snp = [""]
        snpflag = "".join(snp)
        finLines.append("\t".join([chr, start, refs[chr][start]['dbsnp'], snpflag, refs[chr][start]['flag']]))

# Print into a .fin file
finHeader = ["CHROM", "POS", "ID", "SNP", "FLAG"]
f = open(args.sampleID + ".fin", "w+")
f.write('\t'.join(finHeader) + '\n')
f.write('\n'.join(finLines) + '\n')
f.close()


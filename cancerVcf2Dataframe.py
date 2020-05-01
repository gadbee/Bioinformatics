#!/usr/bin/python

import os
import gzip

# Get absolute path of vcf files.
#files = os.popen('ls /home/lideyang/projects/Liaoning_cancer_hospital/breast_cancer/KY073/annotatedMutation/*.vcf')
files = os.popen('ls /home/lideyang/projects/Liaoning_cancer_hospital/breast_cancer/KY073/*.vcf')
files = [i.strip() for i in files]

# Get variant locations of a patient from a vcf file.
def getAllVarLoc(fileAbsPath,LocList):
    fileHandle = open(fileAbsPath)
    for i in fileHandle:
        if i[0] != '#':
            i = i.strip().split('\t')
            chrom = i[0][3:]
            pos = i[1]
            loc = chrom + '_' + pos
            LocList.append(loc)
    fileHandle.close()
    return LocList

LocList = []
for i in files:
    LocList = getAllVarLoc(i,LocList)
LocList = list(set(LocList))

# Process the fna file and transform the seq string to one line and write into a file. Do not write into a dictionary directly because it will be very time-consuming. Once this script had runned, the file of processed genome would have been written to disk and this part would not need to run again.
'''
f = gzip.open('/home/lideyang/projects/Liaoning_cancer_hospital/breast_cancer/KY073/GRCh37_latest_genomic.fna.gz')
fout = open('processedRefGenome.txt','w')
for i in f:
    i = i.strip()
    if i[0] == '>':
        fout.write('\n'+i+'\n')
    else:
        fout.write(i)
fout.close()
f.close()
'''

# Transform the processed genome file into a dictionary.
f = open('processedRefGenome.txt')
WGS_Dict = {}
for i in f:
    i = i.strip()
    if len(i) > 1:
        if i[0] == '>':
            chrom = i.split(',')[0].split(' ')[-1]
            WGS_Dict[chrom] = ''
        else:
            WGS_Dict[chrom] += i.upper()
f.close()

# Transform vcf file to a dataframe in which a row is a patient and a column is a location in genome.
def cancerVcf2Dataframe(fileAbsPath,patientLocDict):
    patient = fileAbsPath.split('/')[-1].split('-')[0]
    patientLocDict[patient] = {}
    fileHandle = open(fileAbsPath)
    for i in fileHandle:
        if i[0] != '#':
            i = i.strip().split('\t')
            chrom = i[0][3:]
            pos = i[1]
            ref = i[3]
            alt = i[4]
            loc = chrom + '_' + pos
            gt = i[-1].split(':')[0]
            if gt == '0/0':
                base = ref
            elif gt == '0/1' or gt == '1/0':
                base = ref + '/' + alt
            elif gt == '1/1':
                base = alt
            patientLocDict[patient][loc] = base
    fileHandle.close()
    return patientLocDict

patientLocDict = {}
for i in files:
    patientLocDict = cancerVcf2Dataframe(i,patientLocDict)

# All the locations of every patient must have a base. Otherwise, get the base from reference genome.
for k,v in patientLocDict.items():
    for i in LocList:
        if i not in v:
            iList = i.split('_')
            chrom = iList[0]
            pos = int(iList[1]) - 1
            v[i] = WGS_Dict[chrom][pos]

# Output the dataframe to stdout.
print 'participant\t',
print '\t'.join(LocList)
for k,v in patientLocDict.items():
    print k,'\t',
    for i in LocList:
        print v[i],'\t',
    print

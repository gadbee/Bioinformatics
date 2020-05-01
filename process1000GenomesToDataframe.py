#!/usr/bin/python

import os

'''def annoteInfo2Dict(annoteInfo):
    annoteInfoDict = {}
    annoteInfoList = annoteInfo.split(';')
    for i in annoteInfoList:
        if '=' in i:
            iList = i.split('=')
            annoteInfoDict[iList[0]] = iList[1]

participants = []
locs = []

def get_participantsAndLocs(fileAbsPath):
    fileHandle = open(fileAbsPath)
    for i in fileHandle:
        if i[0:2] != '##':
            i = i.strip().split('\t')
            if i[0][0] == '#':
                participants.extend(i[9:])
            else:
                loc = i[0] + '_' + i[1]
                locs.append(loc)
    fileHandle.close()
    return participants,locs

vcfFiles = os.popen('ls /home/lideyang/public_data/1000Genomes_vcf/HanWomenData_vcf/HanWomenData_annotated/425CancerPanelAnnotated/*.vcf')
for i in vcfFiles:
    i = i.strip()
    participants,locs = get_participantsAndLocs(i)
participants = list(set(participants))
locs = list(set(locs))

patientLocDict = {}
for i in participants:
    patientLocDict[i] = {}
    for j in locs:
        patientLocDict[i][j] =  
'''

def vcf2Dataframe(fileAbsPath,patientLocDict):
    fileHandle = open(fileAbsPath)
    for i in fileHandle:
        if i[0:2] != '##':
            i = i.strip().split('\t')
            if i[0][0] == '#':
                head = i
            else:
                chrom = i[0]
                pos = i[1]
                loc = chrom + '_' + pos
                ref = i[3]
                alt = i[4]
                #annoteInfo = annoteInfo2Dict(i[7])
                patientInfo = i[9:]
                if '1|1' in patientInfo or '1|0' in patientInfo or '0|1' in patientInfo:
                    cursor = 0
                    for j in patientInfo:
                        if j == '0|0':
                            base = ref
                        elif j == '1|0' or j == '0|1':
                            base = ref + '/' + alt
                        elif j == '1|1':
                            base = alt
                        patientID = head[cursor+9]
                        if patientID not in patientLocDict:
                            patientLocDict[patientID] = {}
                            patientLocDict[patientID][loc] = base
                        else:
                            patientLocDict[patientID][loc] = base
                        cursor += 1
    fileHandle.close()
    return patientLocDict

patientLocDict = {}
vcfFiles = os.popen('ls /home/lideyang/public_data/1000Genomes_vcf/HanWomenData_vcf/HanWomenData_annotated/425CancerPanelAnnotated/*.vcf')

for i in vcfFiles:
    i = i.strip()
    patientLocDict = vcf2Dataframe(i,patientLocDict)

vcfFiles = os.popen('ls /home/lideyang/public_data/1000Genomes_vcf/HanWomenData_vcf/HanWomenData_annotated/*invalid_input')
for i in vcfFiles:
    i = i.strip()
    patientLocDict = vcf2Dataframe(i,patientLocDict)

print 'Participant\t',
for k,v in patientLocDict.items():
    Loc = v.keys()
    print '\t'.join(Loc)
    break

for k,v in patientLocDict.items():
    print k,'\t',
    for i in Loc:
        print v[i],'\t',
    print

'''
def exchangeVcfRanks(fileAbsPath) :
    temp = []
    fileHandle = open(fileAbsPath)
    for i in fileHandle:
        if i[0:2] != '##':
            i = i.strip().split('\t')
            temp.append(i)
    fileHandle.close()
    return map(list,zip(*temp))

a = exchangeVcfRanks('/home/lideyang/public_data/1000Genomes_vcf/HanWomenData_vcf/HanWomenData_annotated/425CancerPanelAnnotated/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.HanWomenhg19_multianno.425Genes.vcf')
for i in a:
    for j in i:
        print j,'\t',
    print
'''

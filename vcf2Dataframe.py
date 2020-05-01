#!/usr/bin/python

import os

variants = []
f = open('/home/lideyang/public_data/1000Genomes_vcf/Dataframe_from_1000GenomeVCF/important_var.txt')
for i in f:
    var = i.strip().split(' ')[0]
    if var[1] != '_':
        var = var[1:]
        variants.append(var)
    else:
        variants.append(var)
f.close()

RefDict = {}
f = open('/home/lideyang/Reference/processedRefGenome_GRCh37_p13.txt')
for i in f:
    if len(i) > 1:
        i = i.strip()
        if i[0] == '>':
            i = i.split(',')[0]
            chrom = i.split(' ')[-1]
        else:
            RefDict[chrom] = i
f.close()

KY073PatientLocDict = {}
vcfs = os.popen('ls /home/lideyang/projects/Liaoning_cancer_hospital/breast_cancer/KY073/*vcf')
for i in vcfs:
    i = i.strip()
    sampleID = i.split('/')[7].split('-')[0]
    KY073PatientLocDict[sampleID] = {}
    f = open(i)
    for j in f:
        j = j.strip()
        if j[0] != '#':
            j = j.split('\t')
            chrom = j[0][3:]
            pos = j[1]
            loc = chrom + '_' + pos
            base = j[4]
            if loc in variants:
                KY073PatientLocDict[sampleID][loc] = base

for k,v in KY073PatientLocDict.items():
    for i in variants:
        if i not in v:
            chrom,pos = i.split('_')
            v[i] = RefDict[chrom][int(pos)-1].upper()

print 'Participant\t','category\t',
for k,v in KY073PatientLocDict.items():
    Loc = v.keys()
    print '\t'.join(Loc)
    break

for k,v in KY073PatientLocDict.items():
    print k,'\t','BreastCancerPatient\t',
    for i in Loc:
        print v[i],'\t',
    print

#!/usr/bin/python

import os
import gzip

sample_info = open('/home/lideyang/public_data/1000Genomes_vcf/igsr_samples.tsv')
expected_samples = []
for i in sample_info:
    i = i.strip().split('\t')
    if i[1] == 'female' and (i[3] == 'CHS' or i[3] == 'CHB'):
        expected_samples.append(i[0])
sample_info.close()

def extract_HanWomenData_from_1000GenomesVcf(i):
    f = gzip.GzipFile(i,'r')
    #outputfile = open('.'.join(i.strip().split('.')[0:-2]) + '.HanWomen' + '.vcf','w')
    #outputlist = []
    for j in f:
        j = j.strip()
        if j[0:2] == '##':
            print j
        else:
            j = j.split('\t')
            if j[0] == '#CHROM':
                header = j[0:9]
                expected_index = []
                for k in expected_samples:
                    if k in j:
                        expected_index.append(j.index(k))
                        header.append(k)
                #outputlist.append('\t'.join(header))
                print '\t'.join(header)
            else:
                line = j[0:9]
                for a in expected_index:
                    line.append(j[a])
                #outputlist.append('\t'.join(line))
                print '\t'.join(line)
    #outputfile.write('\n'.join(outputlist))
    f.close()

'''
chrs = os.popen('ls /home/lideyang/public_data/1000Genomes_vcf/*.vcf.gz')
for i in chrs:
    i = i.strip()
    extract_HanWomenData_from_1000GenomesVcf(i)
'''


#extract_HanWomenData_from_1000GenomesVcf('/home/lideyang/public_data/1000Genomes_vcf/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')

#extract_HanWomenData_from_1000GenomesVcf('/home/lideyang/public_data/1000Genomes_vcf/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz')

extract_HanWomenData_from_1000GenomesVcf('/home/lideyang/public_data/1000Genomes_vcf/ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz')



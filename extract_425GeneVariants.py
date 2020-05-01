#!/usr/bin/python

import os

f = open('/home/lideyang/projects/Liaoning_cancer_hospital/breast_cancer/KY073/ShiHe_GenePanel.txt')
ShiHe_GenePanel = []
for i in f:
    i = i.strip()
    if '(' in i:
        i = i.split('(')[0]
        ShiHe_GenePanel.append(i)
    elif '*' in i:
        i = i.split('*')[0]
        ShiHe_GenePanel.append(i)
    else:
        ShiHe_GenePanel.append(i)
f.close()

vcfFiles = os.popen('ls /home/lideyang/public_data/1000Genomes_vcf/HanWomenData_vcf/HanWomenData_annotated/*hg19_multianno.vcf')
vcfFiles_list = []
for i in vcfFiles:
    vcfFiles_list.append(i.strip())

def extract_425GenesVariants(vcfFile):
    f = open(vcfFile)
    outputFile = open('.'.join(vcfFile.split('.')[0:6]) + vcfFile.split('.')[8] + '.425Genes.vcf','w')
    for i in f:
        if i[0] == '#':
            outputFile.write(i)
        else:
            i_list = i.strip().split('\t')
            annotated_info = i_list[7].split(';')
            annotated_info_dict = {}
            for j in annotated_info:
                if '=' in j:
                    j = j.split('=')
                    annotated_info_dict[j[0]] = j[1]
            gene_name = annotated_info_dict['Gene.refGene'].split('\x3b')
            for k in gene_name:
                if k in ShiHe_GenePanel:
                    outputFile.write(i)
    f.close()
    outputFile.close()


#extract_425GenesVariants('/home/lideyang/public_data/1000Genomes_vcf/HanWomenData_vcf/HanWomenData_annotated/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.HanWomen.annoted.vcf.hg19_multianno.vcf')

for i in vcfFiles_list:
    extract_425GenesVariants(i)

#!/usr/bin/python

import random

f=open('/home/lideyang/cancer_somatic_mutation/analysis_needed_info/all_possible_base_substitution_annote.txt')
all_mutations = []
for i in f:
    i = i.strip().split('\t')
    info = (i[2],i[3],i[4].upper(),i[5].upper(),i[8],i[9])
    all_mutations.append(info)
f.close()

def selected_mutations(a,b,c):
    selected_base_sub = []
    while(1):
        index = random.randint(0,34184)
        selected = all_mutations[index]
        if selected[0] == a and selected[1] == b:
            selected_base_sub.append(all_mutations[index])
        if len(selected_base_sub) == c:
            break
    return selected_base_sub
mut_aas = {}
for i in range(1,10000):
    selected_AC = selected_mutations('A','C',57)
    selected_AG = selected_mutations('A','G',269)
    selected_AT = selected_mutations('A','T',34)
    selected_CA = selected_mutations('C','A',83)
    selected_CG = selected_mutations('C','G',11)
    selected_CT = selected_mutations('C','T',313)
    selected_GA = selected_mutations('G','A',1866)
    selected_GC = selected_mutations('G','C',65)
    selected_GT = selected_mutations('G','T',28)
    selected_TA = selected_mutations('T','A',20)
    selected_TC = selected_mutations('T','C',1210)
    selected_TG = selected_mutations('T','G',18)

    selected = [selected_AC,selected_AG,selected_AT,selected_CA,selected_CG,selected_CT,selected_GA,selected_GC,selected_GT,selected_TA,selected_TC,selected_TG]

    for i in selected:
        for j in i:
            mut_aa = j[4]
            if mut_aa not in mut_aas:
                mut_aas[mut_aa] = 1
            else:
                mut_aas[mut_aa] += 1

for k,v in mut_aas.items():
    print k,'\t',v

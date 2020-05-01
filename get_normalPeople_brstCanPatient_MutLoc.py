#!/usr/bin/python

def get_partLocDict(absPath):
    part_loc = {}
    f = open(absPath)
    for i in f:
        i = i.strip().split('\t')
        if i[0] == 'participant' or i[0] == 'Participant':
            part = i[1:]
        else:
            part_loc[i[0]] = {}
            cursor = 0
            for j in i[1:]:
                part_loc[i[0]][part[cursor]] = j
                cursor += 1
    f.close()
    return part_loc

part_loc_canPatient = get_partLocDict('/home/lideyang/projects/Liaoning_cancer_hospital/breast_cancer/KY073/annotatedMutation/LNCanHsp_mutData_patientLoc.txt')

part_loc_normalPeople = get_partLocDict('/home/lideyang/public_data/1000Genomes_vcf/HanWomenData_vcf/HanWomenData_annotated/425CancerPanelAnnotated/participantID_VarLoc.txt')

locs = open('/home/lideyang/projects/Liaoning_cancer_hospital/breast_cancer/KY073/annotatedMutation/LNCanHsp_mutData_patientLoc.txt').readline().strip().split('\t')[1:]

f = open('/home/lideyang/projects/Liaoning_cancer_hospital/breast_cancer/KY073/annotatedMutation/processedRefGenome.txt')
WGS_Dict = {}
for i in f:
    i = i.strip()
    if len(i) > 1:
        if i[0] == '>':
            chrom = i.split(',')[0].split(' ')[-1]
            WGS_Dict[chrom] =  ''
        else:
            WGS_Dict[chrom] += i.upper()
f.close()

normalPeople_Dict2 = {}
for k,v in part_loc_normalPeople.items():
    normalPeople_Dict2[k] = {}
    for i in locs:
        if i in v:
            normalPeople_Dict2[k][i] = v[i]
        else:
            chrom = i.split('_')[0]
            pos = int(i.split('_')[1])
            base = WGS_Dict[chrom][pos-1]
            normalPeople_Dict2[k][i] = base

normal_patient_Dict = part_loc_canPatient.copy()
normal_patient_Dict.update(normalPeople_Dict2)

print 'PartID','\t',
print '\t'.join(locs)
for k,v in normal_patient_Dict.items():
    print k,'\t',
    for i in locs:
        print v[i],'\t',
    print


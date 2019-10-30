import random
import numpy
import difflib
import fileinput
import argparse
import pandas
from pandas import DataFrame
import os
import time

os.chdir("/Users/maxime/Desktop/MutPred_Analysis_Astyanax/Modele_Neutre_Pygocentrus/")


n_mutation_CF = 58
TsTvratio = 2.28

#parser = argparse.ArgumentParser(description="Neutral sequence evolution for mutpred")
#parser.add_argument('-f', '--file', type=str, required=True, help="File name with sequences")        
        
#args = parser.parse_args()

#args.file

madf = pandas.read_csv("Pygocentrus_vision_CDS.tsv", sep="\t",header=None) #fichier avec les sequences
madf.columns = ["Gene", "cds_seq"] 

madf['cds_length']  = madf["cds_seq"].str.len() #ajouter la longeur de la CDS
madf['prot_length']  = madf["cds_length"]//3   #ajouter la longeur de la proteine

madf = madf.append(madf.agg(['sum']))


def translate(sequence):
    codon2aa = {"AAA":"K", "AAC":"N", "AAG":"K", "AAT":"N", 
                "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T", 
                "AGA":"R", "AGC":"S", "AGG":"R", "AGT":"S", 
                "ATA":"I", "ATC":"I", "ATG":"M", "ATT":"I", 

                "CAA":"Q", "CAC":"H", "CAG":"Q", "CAT":"H", 
                "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P", 
                "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", 
                "CTA":"L", "CTC":"L", "CTG":"L", "CTT":"L", 

                "GAA":"E", "GAC":"D", "GAG":"E", "GAT":"D", 
                "GCA":"A", "GCC":"A", "GCG":"A", "GCT":"A", 
                "GGA":"G", "GGC":"G", "GGG":"G", "GGT":"G", 
                "GTA":"V", "GTC":"V", "GTG":"V", "GTT":"V", 

                "TAA":"*", "TAC":"Y", "TAG":"*", "TAT":"T", 
                "TCA":"S", "TCC":"S", "TCG":"S", "TCT":"S", 
                "TGA":"*", "TGC":"C", "TGG":"W", "TGT":"C", 
                "TTA":"L", "TTC":"F", "TTG":"L", "TTT":"F"}

    protein_seq = ''
    for n in range(0, len(sequence), 3):
        if sequence[n:n+3] in codon2aa:
            protein_seq += codon2aa[sequence[n:n+3]]
        else:
            protein_seq += "X"
    return protein_seq


dna_bases = ["A", "T", "G", "C"]


p_transition = (1/(1+TsTvratio)) * TsTvratio
p_transversion = (1/(1+TsTvratio))/2


def mutate(ancestral_base):
    if ancestral_base == "A":
        new_base = "".join(numpy.random.choice(dna_bases, 1, p=[0,p_transversion,p_transition,p_transversion]))
    elif ancestral_base == "T":
      new_base = "".join(numpy.random.choice(dna_bases, 1, p=[p_transversion,0,p_transversion,p_transition]))
    elif ancestral_base == "G":
     new_base = "".join(numpy.random.choice(dna_bases, 1, p=[p_transition,p_transversion,0,p_transversion]))
    elif ancestral_base == "C":
      new_base = "".join(numpy.random.choice(dna_bases, 1, p=[p_transversion,p_transition,p_transversion,0]))
    return new_base

#Pygocentrus nattereri sequences


prot_seq=[]
for row in madf.itertuples(index=True, name='Pandas'):
    prot_seq.append(translate(row[2]))


madf['prot_seq'] = prot_seq


concatenated_CDS = madf.loc["sum","cds_seq"]  #sequence concatenee de tous les genes
concatenated_evoluated_CDS = concatenated_CDS #On prepare notre sequence evoluee, pr linstant aucun changement

original_prot_sequence = madf.loc["sum","prot_seq"]


list_concatenated_evoluated_CDS=list(concatenated_evoluated_CDS)

count_aa = 0

while count_aa != n_mutation_CF :
    random_position = random.randint(0, len(concatenated_CDS)-1) #Prendre un chiffre aleatoire 
    my_base = concatenated_CDS[random_position] #On prend la base du chiffre aleatoire
    new_base = mutate(my_base) #On mute la base en fonction du ratio Ts/Tv

    list_concatenated_evoluated_CDS[random_position] = new_base #On remplace la base d'origine par la nouvelle
    concatenated_evoluated_CDS = "".join(list_concatenated_evoluated_CDS) 

    evoluated_prot_sequence = translate(concatenated_evoluated_CDS) #Traduire la CDS evoluee
    liste_to_check_evolution = list(evoluated_prot_sequence) #transformation en liste

#si on a mit un stop, alors remettre le nucleotide de depart

    if any("*" in s for s in liste_to_check_evolution):
        list_concatenated_evoluated_CDS[random_position] = my_base
        
    concatenated_evoluated_CDS = "".join(list_concatenated_evoluated_CDS)
    evoluated_prot_sequence = translate(concatenated_evoluated_CDS)

    count_nt = sum(1 for a, b in zip(concatenated_CDS, concatenated_evoluated_CDS) if a != b) #Compte du nombre de mutations sur la CDS
    count_aa = sum(1 for a, b in zip(original_prot_sequence, evoluated_prot_sequence) if a != b) #Compte du nombre de mutations sur la prot



#print(original_prot_sequence)


#evoluated_prot_sequence[0:355] #(prot_length - 1)
#voluated_prot_sequence[355:370]


#evoluated_prot_sequence[:355]
#evoluated_prot_sequence[355:]


#Pour chaque ligne de madf, on va recuperer le gene depuis la sequence concatenee grace a la colonne prot_lengthUnicodeWarning



evoluated_genes=[]
for row in madf.itertuples(index=True, name='Pandas'):
    malen = row[4]
    print(malen)
    evoluated_genes.append(evoluated_prot_sequence[:malen])
    evoluated_prot_sequence = evoluated_prot_sequence[malen:]


madf['new_prot_seq'] = evoluated_genes
madf.drop(madf.tail(1).index,inplace=True)



timestr = time.strftime("%Y%m%d-%H%M%S")
myfile = open("neutral_pygocentrus"+timestr, 'w')

for row in madf.itertuples(index=True, name='Pandas'):
    gene_name=row[1]
    gene_seq=row[5]
    
    new_gene_seq=row[6]
    

    mutations_record = []

    for i in range(0, len(gene_seq), 1):
        aa_old = gene_seq[i]
        aa_new = new_gene_seq[i]
    
        if aa_old != aa_new:
            mutations_record.append(aa_old+str(i+1)+aa_new)


#    print(">"+gene_name+" "+" ".join(mutations_record))
#    print("".join(gene_seq))



    if len(mutations_record) != 0:
        myfile.write(">"+gene_name+" "+" ".join(mutations_record))
        myfile.write("\n")
        myfile.write("".join(gene_seq))
        myfile.write("\n")


myfile.close()




#grep "Pygocentrus" Supplemental_Table1\ -\ Zebrafish+Astyanax+Sinocyclocheilus\ eye\ genes.tsv | cut -f2,10 | sed 's/...$//g'
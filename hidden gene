rna_seq = open('hidden_gene.fasta','r')
gencode = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

def translate_rna(sequence):
    protein_seq = ''
    for n in range(0, len(sequence), 3):
        if sequence[n:n+3] in gencode:
            protein_seq += gencode[sequence[n:n+3]]
    return protein_seq


ids, sequences = [], []
n = -1
with open('hidden_gene.fasta') as fh:
    for line in fh:
        line = line.strip()
        if line[0] == '>':
            id = line.split()[0][1:]
            ids.append(id)
            sequences.append('')
            n += 1
        else:
            sequences[n] += line


for id, seq in zip(ids, sequences):
    print('>'+id)
    protein = translate_rna(seq)
    print(protein)

proteinseq_dict = {'>sequence_1':'ALAHKTLVLLLGVDPSRQLDHPLPTVHPQVTYAYMKNMWKSARKVCIPMHPEQTTFSRMLA*SASLCLFKARKETVFKPPTHNVRESLPVTKLFGLIS*SES*LRRWCVYH*SNQ*M*SQLYVFVSRKCWDIERAHFLPLREDGLFGTKKSRFVNTQTRT','>sequence_4':'QWDSMEEYTCMIPRDTHDGAFYRAVLALHQDLFSLAQQ','>sequence_9':'SLIGVEGGNATRIGRFANYLRNLLPSNDPVVMEMASKAIGRLAMAGDTFTAEYVEFEVKRALEWLGADRNEGRRHAA','>sequence_10':'MSQEESTRFYDQLNHHIFELVSSSDANERKGGILAI','>sequence_11':'GPLPLRDDNGIVLLGERAAKCRAYAKALHYKELEFQKGPTPAILESLI','>sequence_12':'SCMD*LLDQLFVMSFWLWVNPVSSCQLIPYIEKHLRLVLAPKLGSFIFFKAIITNMFFVCLFVFLRQ'}
import json
result = json.dumps(proteinseq_dict)
      

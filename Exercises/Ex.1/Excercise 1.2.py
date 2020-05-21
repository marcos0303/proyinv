#esta es de Javi, el código me lo explicó poco a poco. Todo quedó entendido.
s = 'GATGTAGGTCTGAAAGCAAAAGAAATTTAGGTAGCGGGCAAGTTGGTGACTGATGAGTTAGTTATCGCATTACTCAAAGAACGTATCACACAGGAAGATTGCCGCGTAGGTTTTCTGTTAGACGGGTTCCCGCGTACCATTCCTCAGGCAGATTGCCAGTGAAAGAAGCCGGTATCAAAGTTGATTATTGTGCTGGAGTTTGACTGTTCAACGAGCTGATTGTTGAGCGCATTGTCGGCCGTCGGGTACATCGCTGCTTCAGGCCGTTTATCACGTTAGAATTCAACCCACCAGTTGAAGATAAAGATCGAGTGTTACCGGTGAAGAGCTGACCTATTGTGAATAGAAACTGATCAGGAAGCGACTGTCCGTAAGCGTCTTATCCGAATATCATCAACAAACTGCACATTGGTTTCTTACTATCATAGAAGAGCGGATTGCAGGTAATCCAATATTTAAACTGGACGGAACCCGTAATTGTAGCAAAGTCAGTGCTGAACGGCGACTATTCTCGGTTAATTCTGGATCGGCCTTATAGCTAAGGCGGTTTAAGGCCGCCTTAGCTATTCAAGTAAGAAGGGCGTAGTACCTACAAAAGGAGATTTGGCATAGATAGCAAAGCAAACCCGGCGATTAATTGGTTAATTTGGGGACACCAGTAGCTC'
#s = s.replace('T','U')



# DEfine your trinucleotide to amino acid translation dictionary:
gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
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

gencode_U = dict()
for i in gencode.keys():
    gencode_U[i.replace('T','U')] = gencode[i]



# Finding indexes where ATG anr TAA start in the sequence:
TAAlocs = [x for x in range(0,len(s)-2) if s[x:x+3] == 'TAA']
ATGlocs = [x for x in range(0,len(s)-2) if s[x:x+3] == 'ATG']



# Iterate over all ATG and TAA indexes combinations and test if the sequence fraction they form:
#1) if multiple of 3, to create Amino_acids
#2) ATG index must be always lower than the TAA index
# If so, store each aminoacid into a protein with a uniq ID


proteins_dict = dict()

for ATG_index in ATGlocs:
    count = 0
    for TAA_index in TAAlocs:
        if (ATG_index < TAA_index) and ((TAA_index - ATG_index) % 3 == 0):
            protein_ID = '_'.join(['protein', str(ATG_index), str(TAA_index), str(count)])
            proteins_dict[protein_ID] = ''
            for trint_index in range(ATG_index, TAA_index,3):
                nt = s[trint_index:trint_index+3]
                proteins_dict[protein_ID] = proteins_dict[protein_ID] + gencode[nt]
            count = count + 1

# Find out the longest protein
length = 0

for i in proteins_dict.keys():
    test_length = len(proteins_dict[i])
    if test_length > length:
        ID = i
        length = test_length

print('The longest protein is ' + ID + ', ' + str(length) + ' amino acids long')

'''

This script assesses the impact of amino acid substitutions based upon three different matrixes:
	BLOSUM100, PAM30, GONNET
		These matrices can be adjusted to be more or less strigent depending on your preference

Adds additional weight to changes if there are charge changes

If you provide a file which predicts surface residues, you can also indicate which changes are surface exposed
	if you do not include a surface file, remove all the variables associated with surface file in order to execute script w/o errors


### In order to run this script you will need clustal w installed on your machine ###

'''

import os
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio import Phylo
from Bio.SubsMat.MatrixInfo import blosum90 as blosum
from Bio.SubsMat.MatrixInfo import pam30 as pam
from Bio.SubsMat.MatrixInfo import gonnet
from Bio.SeqUtils.ProtParam import ProteinAnalysis

groupings = {
	'G': 'Aliphatic',
	'A': 'Aliphatic',
	'V': 'Aliphatic',
	'C': 'Sulfur-containing/Aliphatic',
	'P': 'Aliphatic',
	'L': 'Aliphatic',
	'I': 'Aliphatic',
	'M': 'Sulfur-containing/Aliphatic',
	'T': 'Hydroxylic/polar',
	'F': 'Aromatic',
	'W': 'Aromatic',
	'Y': 'Aromatic',
	'D': 'Acidic',
	'E': 'Acidic',
	'R': 'Basic',
	'H': 'Basic',
	'K': 'Basic',
	'S': 'Hydroxylic/polar',
	'N': 'Amidic; polar',
	'Q': 'Amidic; polar',
	'-': 'Gap'
}

acidic = ['D', 'E']
basic = ['K', 'R', 'H']

# setting path and aligning sequences
clustalw_exe = r"path_to_clustal"
clustalw_cline = ClustalwCommandline(clustalw_exe, infile="sequences_in_FASTA_format.fasta")
assert os.path.isfile(clustalw_exe), "Clustal W exectuable missing"
stdout, stderr = clustalw_cline()

# reading alingment file and converting to FASTA format (from clustal)
align_fasta = AlignIO.convert("sequences_in_FASTA_format.aln", "clustal", "converting_alignment_to_FASTA.fasta", "fasta")

# opening FASTA file and creating dictionary of sequences (key=name: value=sequence)
my_file = 'converting_alignment_to_FASTA.fasta'

# opening surface resiude file (created with Surface Racer - http://biskit.pasteur.fr/install/applications/surfrace) 
# and creating dictionary 
# (key=site #: value=surface value)
# surface value arbitrary, only interested in 0 v not 0 (0 means not exposed)
my_surface = 'surface_racer_file.txt'

class ProteinManager():

	def __init__(self, my_file, my_surface):
	# def __init__(self, myfile):
		self.my_file = my_file
		self.my_surface = my_surface
		self.names = dict()
		self.names_complete = dict()
		self.headers = list()
		self.diffs = list()
		self.aminos = list()
		self.indi = list()
		self.same = list()
		self.consensus = ''
		self.numbers = list()
		self.surface = dict()
		self.exposed = ''
		with open(self.my_file, 'r') as file:
			lines = [line.strip() for line in file.readlines()]
			for line in lines:
				if '>' in line:
					header = line
					self.headers.append(header)
					self.names[header] = ''
					self.names_complete[header] = ''
					seq = ''
				elif '>' not in line:
					seq = seq + line
					self.names_complete[header] = seq
					self.names[header] = list(seq) # turing each residue into a seperate string character. 
									               # now can call on them by position (subracting 1 since starts at 0 and not 1)

		with open(self.my_surface, 'r') as file2:
			lines = [line.strip() for line in file2.readlines()]
			for line in lines:
				site = line[23:26].lstrip()
				if site in self.numbers:
					if self.surface[site] != 0:
						pass
					elif self.surface[site] == 0 and float(line[63:69].rstrip()) == 0:
						pass
					elif self.surface[site] == 0 and float(line[63:69].rstrip()) != 0:
						self.surface[site] = float(line[63:69].rstrip())
				elif site not in self.numbers:
					self.numbers.append(site)
					self.surface[site] =  float(line[63:69].rstrip())


# identifying which residues are different.
# searches by residue for adjacent sequences within dictionary
	def differences(self):
		for i in range(0, (len(self.headers)-1)):
			for n in range(len(self.names[self.headers[0]])):
				if self.names[self.headers[i]][n] != self.names[self.headers[i+1]][n]:
					if n+1 not in self.diffs:
						self.diffs.append(n+1)
		self.diffs = sorted(self.diffs)
		print(f'\nList of residue differences: {self.diffs}')
		print(f'\nTotal: {len(self.diffs)}')
		return self.diffs

	# generating the output
	def output(self):
		# calculating percent similarity (overall)
		for n in range(1, len(self.names[self.headers[0]])):
			if n not in self.diffs:
				self.same.append(n)
		percentage = (len(self.same) / len(self.names[self.headers[0]])) * 100

		# printing identity
		print(f'\n\nPercent identity')
		print('---------------') 
		print(f'{round(percentage, 1)}%')
		print(f'\n\nList of shared residues: {self.same}')
		print(f'\nTotal: {len(self.same)}')

		# uncomment to print out surface residues
		# surface_total = list()
		# for s in self.surface:
		# 	if self.surface[s] == 0:
		# 		surface_total.append(s)
		# print(f'\n\nList of shared surface exposed residues: {surface_total}')
		# print(f'\nTotal: {len(surface_total)}')
		

		# finding the pI for each protein
		print('\n\nIsoelectric Point')
		print('----------------------')
		print('Gene\t\tpI')
		for value in self.names_complete:
			analysis = ProteinAnalysis(self.names_complete[value])
			pI = analysis.isoelectric_point()
			print(f'{value[1:]}:\t{round(pI, 2)}')



		# identifying the mutations
		print('\n\n\nMutations')
		print('---------------') 
		print('---------------') 

		# nesting lists within each other
		# each list (self.indi) represents a different residue which is altered in one or more sequences
		for val in self.diffs:
			self.indi = []
			self.aminos = []
			self.aminos.append(self.indi)
			print()

			# creating table showing differences for each residue
			print(f'''\t\t{val}\nGene''')
			for i in range(0, (len(self.headers))):
				print(f'''{self.headers[i][1:]}\t{self.names[self.headers[i]][val-1]}''')
				self.indi.extend(self.names[self.headers[i]][val-1])
			for x in self.aminos:
				amino = x[i]
				# identifying consensus residue (needs fixing still!)
				self.consensus = self.counter(x)
				percent_consensus = (self.consensus[1]/len(self.headers))*100
				print(f'\nConsensus residue: {self.consensus[0]}')
				print(f'Consensus Usage: {percent_consensus}% of sequences')
			self.impact()

	# creating individual assessments (for each residue that's mutated)
	def impact(self):
		for i in range(0, (len(self.headers))):
			for x in self.aminos:
				amino = x[i]
				if amino == self.consensus[0]:
					pass
				elif amino != self.consensus[0]:
					indi_change = 0
					print(f'\n\n~{self.headers[i][1:]}~')
					print('Type of change')
					print('---------------')
					print('Consensus: {}, {} amino acid ---> Mutation: {}, {} amino acid'.format(self.consensus[0], groupings[self.consensus[0]], amino, groupings[amino]))
					print('\n\nMatrix scores')
					print('--------------')
					for key in blosum:
						if (self.consensus[0], amino) == key or (amino, self.consensus[0]) == key:
							print(f'BLOSUM90: {blosum[key]}')
							indi_change = float(indi_change) + float(blosum[key])
					for key in pam:
						if (self.consensus[0], amino) == key or (amino, self.consensus[0]) == key:
							print(f'PAM30: {pam[key]}')
							indi_change = float(indi_change) + float(pam[key])
					for key in gonnet:
						if (self.consensus[0], amino) == key or (amino, self.consensus[0]) == key:
							print(f'Gonnet: {gonnet[key]}\n')
							indi_change = float(indi_change) + float(gonnet[key])

					print()

					print('Charge Change')
					print('-------------')
					for x in self.aminos:
						if amino in acidic and self.consensus[0] in basic:
							print('Yes, acidic to basic*')
							charge_score = 3
						elif amino in basic and self.consensus[0] in acidic:
							print('Yes, basic to acidic*')
							charge_score = 3
						elif amino in acidic and self.consensus[0] not in acidic:
							print('Yes, neutral to acidic')
							charge_score = 2
						elif amino in basic and self.consensus[0] not in basic:
							print('Yes, neutral to basic')
							charge_score = 2
						elif self.consensus[0] in acidic and amino not in acidic:
							print('No, but neutralizing (acidic to neutral)*')
							charge_score = 1
						elif self.consensus[0] in basic and amino not in basic:
							print('No, but neutralizing (basic to neutral)*')
							charge_score = 1
						else:
							print('No')
							charge_score = 0
						indi_change_average = (float(indi_change) - float(charge_score)) / 4

						# print('\n\nSurface exposure')
						# print('----------')
						# print(f'Exposed: {self.exposed}')

						print('\n\nImpact')
						print('----------')
						if indi_change <= 0:
							print(f'Average score: {round(indi_change_average, 1)}; Could be impactful (score <= 0)\n\n')
						elif indi_change > 0:
							print(f'Average score: {round(indi_change_average, 1)}; May not be impactful (score > 0)\n\n')
						print('----------------------------------------------------------------------')
						print('----------------------------------------------------------------------')

	# identifying the consensus amino acid
	# right now only reports if 1. need to figure out how to evaluate for multiple....
	def counter(self, x):
		counter_dict = dict()
		counter_list = list()
		for a in x:
			if a not in counter_list:
				counter_dict[a] = 1
				counter_list.append(a)
			elif a in counter_list:
				counter_dict[a] = counter_dict[a] + 1
		counter_dict = sorted(counter_dict.items(), key=lambda x: x[1], reverse=True)
		dict_list = list(counter_dict)
		return dict_list[0]


test = ProteinManager(my_file, my_surface)
test.differences()
test.output()

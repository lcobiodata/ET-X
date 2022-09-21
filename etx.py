###====================================================================================================
### License
###====================================================================================================
# X-ET is a program that implements alphabet expansion to increase Evolutionary Trace sensitivity to marginally conserved protein sites.

# Copyright (C) 2018,   Lucas Carrijo de Oliveira (lucas.carrijodeoliveira@gmail.com)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
###====================================================================================================
### Libraries
###====================================================================================================
try:
	import argparse
	from Bio import SeqIO
	from Bio import Phylo
	from Bio.SubsMat import MatrixInfo
	import numpy as np
	import networkx as nx
	import matplotlib.pyplot as plt
	import pandas as pd
	from Bio.SeqUtils import seq3
except:
	raise ImportError('Need the following packages installed: argparse, numpy, scipy, matplotlib, networkx, Bio, pyclustering, sklearn, pandas and seaborn.')
###====================================================================================================
### Parameters
###====================================================================================================
parser = argparse.ArgumentParser(description = "Density Based Residue Clustering by Dissimilarity Between Sequence SubSets (DBRC/DBS3).")
parser.add_argument("msa_file", help = "Path to multiple sequence alignment file in FASTA format.", type = str)
parser.add_argument("tree_file", help = "Path to tree file in NEWICK format.", type = str)
parser.add_argument("-o", "--out", help = "Label to be used for naming output files (default: stdout)", type = str, default = None, required = False)
parser.add_argument("-x", "--plus_aa", help = "True for expanding alphabet (default: False).", type = bool, default = False, required = False)
args = parser.parse_args()	# returns data from the options specified (echo)
if args.out == None:
	args.out = args.msa_file.split('.')[0]
###====================================================================================================
### Classes
###====================================================================================================
class MSA(object):
	"""A class representing a Multiple Sequence Alignment"""
	def __init__(self, MSA_file):
		super(MSA, self).__init__()
		self.alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
		self.sthereochemistry = {
			'Aliphatic':['G','A','V','L','I'],
			'Amide':['N','Q'],
			'Aromatic':['F','Y','W'],
			'Basic':['H','K','R'],
			'Big':['M','I','L','K','R'],
			'Hydrophilic':['R','K','N','Q','P','D'],
			'Median':['E','V','Q','H'],
			'Negatively charged':['D','E'],
			'Non-polar':['F','G','V','L','A','I','P','M','W'],
			'Polar':['Y','S','N','T','Q','C'],
			'Positively charged':['K','R'],
			'Similar (Asn or Asp)':['N','D'],
			'Similar (Gln or Glu)':['Q','E'],
			'Small':['C','D','P','N','T'],
			'Tiny':['G','A','S'],
			'Very hydrophobic':['L','I','F','W','V','M'],
			'With hydroxyl':['S','T','Y'],
			'With sulfur':['C','M']
		}
		self.data = self.parse(MSA_file)
		self.headers, self.sequences = self.read()
		self.size, self.length = np.array(self.sequences).shape
		self.weights = self.Henikoff()
		self.sequence_indices = {x:i for i,x in enumerate(self.headers)}
		self.collection = self.collect()
		self.gap_content = self.count_gaps()
	def parse(self, MSA_file):
		try:
			return list(SeqIO.parse(MSA_file, "fasta"))
		except FileNotFoundError:
			raise SystemExit("No such file or directory: '%s'" % args.msa)
	def read(self):
		headers, sequences = [], []
		for seq in self.data:
			headers.append(seq.id)
			sequences.append(list(seq.seq))
		return (headers, np.array(sequences))
	def Henikoff(self):
		weights = []
		for i in range(self.size):
			row = []
			for j in range(self.length):
				x = self.sequences[i][j]
				k = float(len(set(self.sequences[:, j])))
				n = 0.
				for y in self.sequences[:, j]:
					if y == x:
						n += 1.
				row.append(1. / (k * n)) 
			weights.append(sum(row) / float(self.length))
		return weights
	def collect(self):
		collection = {}
		for l in range(self.length):
			collection[l]=[]
			for a in self.alphabet:
				if a in self.sequences[:, l]:
					sequence_indices = []
					for n in range(self.size):
						if self.sequences[n][l]==a:
							sequence_indices.append(n)
					collection[l].append(Residue(self, [a], l, sequence_indices))
		if args.plus_aa:
			for l in range(self.length): # For each column of the alignment, it looks for all possible subsets of similar amino acids.
				temp = {}
				for k, v in list(self.sthereochemistry.items()):
					if len(set(self.sequences[:, l]) & set(v)) > 0:
						temp[k] = ([],[])
				for n in range(self.size):
					for k in list(temp.keys()):
						if self.sequences[n][l] in self.sthereochemistry[k]:
							if self.sequences[n][l] not in temp[k][0]:
								temp[k][0].append(self.sequences[n][l])
							temp[k][1].append(n)
				for k, (x, y) in list(temp.items()):
					temp[k] = (tuple(x), tuple(y))
				aux={x: [] for x in set(temp.values())}
				for k, v in list(temp.items()):
					aux[v].append(k)
				for (aa, idx), ftr in list(aux.items()):
					if len(aa) > 1:
						for i in range(len(ftr)):
							if 'Similar' in ftr[i]:
								ftr[i] = 'Similar'
						label = ', '.join(ftr)
						collection[l].append(Residue(self, list(aa), l, list(idx), ftr))
		return collection
	def count_gaps(self):
		gap_content = []
		for l in range(self.length):
			temp = []
			for n in range(self.size):
				if self.sequences[n][l] == '-':
					temp.append(n)
			gap_content.append(Residue(self, ['-'], l, temp, None, '-%d'%(l+1)))
			# gap_content.append(sum([self.weights[n] for n in range(self.size) if self.sequences[n][l] == '-']))
		return gap_content
###====================================================================================================
class Subset(object):
	"""A class for subsets of sequences from the MSA"""
	def __init__(self, MSA_object, sequence_indices, label = None):
		super(Subset, self).__init__()
		self.msa = MSA_object
		self.sequence_indices = set(sequence_indices)
		self.label = label
		self.p = self.probability(self)
	def __repr__(self):
		return self.label
	class probability:
		def __init__(self, subset):
			self.subset = subset
			self.result = float(sum(map(lambda x: self.subset.msa.weights[x], self.subset.sequence_indices)))
		def __call__(self):
			return self.result
		def given(self, other_subset):
			return float(sum(map(lambda x: self.subset.msa.weights[x], self.subset.sequence_indices & other_subset.sequence_indices))) / other_subset.p()
	def get_consensus(self):
		consensus = ''
		for l in range(self.msa.length):
			K = [x for x in self.msa.collection[l] if len(x.amino_acids) == 1 and not x.sequence_indices.isdisjoint(self.sequence_indices)]
			if len(K)>0:
				consensus += max(K, key=lambda x: x.p.given(self)).amino_acids[0]
			else:
				consensus += '-'
		return consensus
	def Shannon(self, position):
		K = [x for x in self.msa.collection[position] if not x.sequence_indices.isdisjoint(self.sequence_indices)]
		if len(K) == 0:
			return np.log(20)
		if args.plus_aa:
			content = set([self.msa.sequences[i][position] for i in self.sequence_indices])
			aux, K = K, []
			aux = sorted(aux, key=lambda x: x.p.given(self), reverse=True)
			done, i = False, 0
			while not done and i < len(aux):
				x = aux[i]
				if i == 0:
					K.append(x)
					i += 1
				else:
					temp = set([self.msa.sequences[j][position] for j in set.union(*map(lambda k: k.sequence_indices, K))])
					if temp == content:
						done = True
					else:
						if set(x.amino_acids).isdisjoint(temp):
							K.append(x)
						i += 1
		g = self.msa.gap_content[position].p.given(self)
		S = ((1-g)*np.log((1-g)/20))+np.log(20)
		for x in K:
			p = x.p.given(self)
			S -= p*np.log(p)
		# print(position, self, g, K, S)
		return S
###====================================================================================================
class Residue(Subset): # A special kind of Subset, inherited from it, for sets of sequences defined by having a specific residue (or sthereochemistry) at a given position.
	def __init__(self, MSA_object, amino_acid, position, sequence_indices, sthereochemistry=None, label = None):
		Subset.__init__(self, MSA_object, sequence_indices, label)
		self.sthereochemistry = sthereochemistry
		self.amino_acids = amino_acid
		self.position = position
		self.label = self.labeling(label)
	def labeling(self, label):
		if label == None:
			if len(self.amino_acids) == 1:
				return seq3(self.amino_acids[0])+str(self.position+1)
			else:
				self.amino_acids = sorted(self.amino_acids)
				sets = []
				for x in self.sthereochemistry:
					if 'Similar' in x:
						sets.append(x.split('(')[0])
					else:
						sets.append(x)
				sets = sorted(list(set(sets)))
				if len(sets) > 1:
					feature = ', '.join(sets)
				else:
					feature = sets[0]
				return '%s (%s or %s) %d' % (feature, ', '.join(map(seq3, self.amino_acids[:-1])), seq3(self.amino_acids[-1]), self.position+1)
		return label
###====================================================================================================
class Clade(Subset):
	def __init__(self, MSA_object, sequence_indices, branch_length=None, label=None, clades=None):
		Subset.__init__(self, MSA_object, sequence_indices, label)
		self.branch_length = branch_length
		self.clades = clades
###====================================================================================================
### FUNCTIONS
###====================================================================================================
def similarity(aa1, aa2):
	alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	pair = (aa1, aa2)
	if aa1 == aa2:
		return 1
	if aa1 == '-' or aa2 == '-':
		return 0
	if pair not in MatrixInfo.blosum62.keys():
		pair = (aa2, aa1)
	if MatrixInfo.blosum62[pair]>=np.mean(list(map(lambda x: MatrixInfo.blosum62[(x,x)], alphabet))):
		return 1
	return 0
###====================================================================================================
def ET(MSA_object, DAG, groups, root, leaf):
		path = nx.shortest_path(DAG, root, leaf)
		ranking = {}
		for l in range(MSA_object.length):
			score = []
			for n in path:
				if len(n.clades)>0:
					host, i = None, 0
					while host == None and i < len(groups[n]):
						if MSA_object.sequence_indices[leaf.label] in groups[n][i].sequence_indices:
							host = groups[n][i]
						else:
							i += 1
					d0, entropy = .05,  []
					for g in groups[n]:
						seq1, seq2 = host.get_consensus(), g.get_consensus()
						temp = []
						for a,b in list(zip(seq1,seq2)):
							temp.append(similarity(a,b))
						d = 1.-np.mean(temp)
						Wg = np.exp(-((d**2)/(d0**2)))
						s = g.Shannon(l)
						if s != None:
							entropy.append(Wg*s)
					score.append(sum(entropy))
			ranking[l]=1.+float(sum(score))
		return ranking
###====================================================================================================
### MAIN
###====================================================================================================
if __name__=="__main__":
	msa = MSA(args.msa_file)
	tree = Phylo.read(args.tree_file, 'newick')
	tree.ladderize()   # Flip branches so deeper clades are displayed at top
	clades = list(tree.find_clades(order='level'))
	subfamily = {}
	leaves = []
	for i, clade in enumerate(clades):
		leaf = False
		if clade.is_terminal():
			leaf = True
		if not leaf:
			clade.name = 'N%d'%i
		subfamily[clade] = Clade(msa, [msa.sequence_indices[x.name] for x in list(clade.get_terminals())], clade.branch_length, clade.name)
		if leaf:
			leaves.append(subfamily[clade])
	net = nx.DiGraph()
	for clade in clades:
		for child in clade.clades:
			net.add_edge(subfamily[clade], subfamily[child])
	nodes = list(net.nodes())
	for node in nodes:
		node.clades = list(net.successors(node))
	root, i = None, 0
	while root == None and i<len(nodes):
		if net.in_degree(nodes[i]) == 0:
			root = nodes[i]
		else:
			i+=1
	D = {}
	for node in nodes:
		d = 0.
		path = list(nx.shortest_path(net, root, node))
		for step in path[1:]:
			d += step.branch_length
		D[node]=d
	aux = nx.DiGraph()
	for u,v in list(net.edges()):
		aux.add_edge(u,v)
	groups = {}
	for n, d in sorted(D.items(), key=lambda x: x[1]):
		for x in list(aux.nodes()):
			if D[x]<d:
				aux.remove_node(x)
		if len(n.clades)>0:
			groups[n] = [x for x in list(aux.nodes()) if aux.in_degree(x)==0 and len(x.clades)>0]

	

	# for leaf in leaves:
	# 	print(leaf)
	# 	ranking = ET(msa, net, groups, root, leaf)
	# 	count = 1
	# 	for k,v in sorted(ranking.items(), key=lambda x: x[1]):
	# 		aa = msa.sequences[msa.sequence_indices[leaf.label]][k]
	# 		if aa != '-':
	# 			print(count, seq3(msa.sequences[msa.sequence_indices[leaf.label]][k])+str(k+1), v)
	# 			count+=1
	# 	print('\n')


	leaf = leaves[0]
	print(leaf)
	ranking = ET(msa, net, groups, root, leaf)
	count = 1
	for k,v in sorted(ranking.items(), key=lambda x: x[1]):
		aa = msa.sequences[msa.sequence_indices[leaf.label]][k]
		if aa != '-':
			print('%d\t%s\t%d\t%f'%(count, seq3(msa.sequences[msa.sequence_indices[leaf.label]][k]), k+1, v))
			count+=1
	print('\n')



				







	
"""rosalind.py by https://github.com/KEYS248

Solves bioinformatic problems from the website https://rosalind.info

Takes the filename for the rosalind problem as input as such:
	$ python3 rosalind.py rosalind_orf.txt

Some problems require a codon table, which should be included in the directory with this file
"""

from Bio import SeqIO
import sys


def main():
	if len(sys.argv) == 1:
		print("Please run again with a text file specified as an argument")
	else:
		function = sys.argv[1].split(".")[0]
		eval(function + '(sys.argv[1])')


def rosalind_cons(filename):
	records = list(SeqIO.parse(filename, "fasta"))
	master = [ [0] * len(records[0].seq), [0] * len(records[0].seq), [0] * len(records[0].seq), [0] * len(records[0].seq) ]
	for record in records:
		for i in range(len(record.seq)):
			if record.seq[i] == "A":
				master[0][i] += 1
			elif record.seq[i] == "C":
				master[1][i] += 1
			elif record.seq[i] == "G":
				master[2][i] += 1
			else:
				master[3][i] += 1
	consensus = ""
	for i in range(len(master[0])):
		total = master[0][i] + master[1][i] + master[2][i] + master[3][i]
		if master[0][i]*4 > total:
			consensus += "A"
		elif master[1][i]*4 > total:
			consensus += "C"
		elif master[2][i]*4 > total:
			consensus += "G"
		else:
			consensus += "T"
	print(len(consensus) == len(records[0].seq))
	print(len(consensus) == len(records[-1].seq))
	f = open("results.txt", "w")
	f.write(consensus)
	f.write("\nA: {}".format(" ".join(str(x) for x in master[0])))
	f.write("\nC: {}".format(" ".join(str(x) for x in master[1])))
	f.write("\nG: {}".format(" ".join(str(x) for x in master[2])))
	f.write("\nT: {}".format(" ".join(str(x) for x in master[3])))
	f.close()


def rosalind_fib(filename):
	file = open(filename, "r")
	data = file.read().split(" ")
	file.close()
	repro = 1
	growing = 0
	for _ in range(int(data[0])-2):
		temp = int(repro)
		repro += growing
		growing = temp*int(data[1])
	print(repro+growing)


def rosalind_iprb(filename):
	file = open(filename, "r")
	data = file.read().split(" ")
	file.close()
	k = int(data[0])
	m = int(data[1])
	n = int(data[2])
	t = k + m + n
	num = (t-1)*k+m*(k+3/4*(m-1)+1/2*n)+n*(k+1/2*m)
	print(num/(t*(t-1)))


def rosalind_lcsm(filename):
	records = list(SeqIO.parse(filename, "fasta"))
	base_seq = records[0].seq
	for record in records:
		if len(record.seq) < len(base_seq):
			base_seq = record.seq
	subseq = "", 0
	for length in range(2, len(base_seq)):
		for start in range(0, len(base_seq)):
			if start+length >= len(base_seq) or len(subseq) == length:
				break
			target = base_seq[start:start+length]
			shared = [False] * len(records)
			for i in range(0, len(records)):
				compare = records[i].seq
				index = compare.find(target)
				if index >= 0:
					shared[i] = True
			if False not in shared:
				subseq = target
				print(target)
				continue
	print(subseq)


def rosalind_sseq(filename):
	records = list(SeqIO.parse(filename, "fasta"))
	seq = records[0].seq
	match, curr_base = records[1].seq, 0
	indexes = []
	for start in range(len(seq)):
		if seq[start:start+1] == match[curr_base]:
			indexes.append(str(start+1))
			curr_base += 1
		if curr_base >= len(match):
			break
	print(" ".join(indexes))


def rosalind_orf(filename):
	record = list(SeqIO.parse(filename, "fasta"))
	fwd = str(record[0].seq).replace("T", "U")
	rev = str(fwd[::-1])

	import csv
	f = open("rna_codon_table.csv", 'r')
	codons = list(csv.reader(f))
	f.close()

	prots = []
	for frame in range(0, 3):
		for start in range(frame, len(fwd), 3):
			if fwd[start:start+3] == "AUG":
				prot = translate(fwd[start:], codons)
				if prot not in prots:
					prots.append(prot)
		for start in range(frame, len(rev), 3):
			if rev[start:start + 3] == "AUG":
				prot = translate(rev[start:], codons)
				if prot not in prots:
					prots.append(prot)
	for prot in prots:
		print(prot)


def rosalind_grph(filename):
	records = list(SeqIO.parse(filename, "fasta"))
	node_list = []
	for record in records:
		temp = Node(record.id, record.seq[-3:])
		for rec in records:
			if rec.seq[:3] == temp.suffix and rec.id != temp.id:
				temp.add_child(rec.id)
		node_list.append(temp)
	for node in node_list:
		node.print_children()


def translate(rna, codons):
	output = ""
	for i in range(0, len(rna), 3):
		for j in range(len(codons)):
			if rna[i:i+3] == codons[j][0]:
				if codons[j][1] == "Stop":
					return output
				output += codons[j][1]
	return output


class Node:
	def __init__(self, id, suffix):
		self.id = id
		self.suffix = suffix
		self.children = []

	def add_child(self, id):
		self.children.append(id)

	def print_children(self):
		for child in self.children:
			print(self.id, child)


def rosalind_perm(data):
	import random
	import math
	num = int(data)
	perms = []
	output = str(math.factorial(num)) + "\n"
	for _ in range(math.factorial(num)):
		add = [str(random.randint(1, num)) for i in range(num)]
		while add in perms:
			add = [str(random.randint(1, num)) for i in range(num)]
		perms.append(add)
		output += " ".join(add) + "\n"
	return output

def rosalind_gc(data):
	data = data.replace("\n", '')
	lines = data.split(">")
	highest_gc, highest_index = 0.0, 0

	for i in range(len(lines)):
		if len(lines[i]) == 0:
			continue
		seq = lines[i][13:]
		gc_count = 0
		for letter in seq:
			if letter == "G" or letter == "C":
				gc_count += 1
		gc_percent = gc_count / len(seq)
		if gc_percent > highest_gc:
			highest_gc = gc_percent
			highest_index = i
	return lines[highest_index][:13] + "\n" + str(highest_gc)

def rosalind_dna(data):
	output = [0, 0, 0, 0]

	for i in data:
		if i == 'A':
			output[0] += 1
		elif i == 'C':
			output[1] += 1
		elif i == 'G':
			output[2] += 1
		elif i == 'T':
			output[3] +=1
	return output

def rosalind_rna(data):
	output = data.replace('T', 'U')
	return output

def rosalind_revc(data):
	rev = "".join(reversed(data))
	output = ""
	for i in rev:
		if i == 'A':
			output += 'T'
		elif i == 'C':
			output += 'G'
		elif i == 'G':
			output += 'C'
		elif i == 'T':
			output += 'A'
	return output

def rosalind_hamm(data):
	output = 0
	seqs = data.split("\n")
	for i in range(len(seqs[0])):
		if seqs[0][i] != seqs[1][i]:
			output += 1
	return output

def rosalind_subs(data):
	locs = []
	seq = data.split("\n")[0]
	tar = data.split("\n")[1]
	for i in range(len(seq)):
		if i + len(tar) >= len(seq):
			continue
		if seq[i:i+len(tar)] == tar:
			locs.append(str(i+1))
	return " ".join(locs)


main()

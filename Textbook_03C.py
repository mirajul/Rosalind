#!/usr/bin/env python
'''
A solution to a code challenges that accompanies Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.
The textbook is hosted on Stepic and the problem is listed on ROSALIND under the Textbook track.

Problem Title: Profile-most Probable k-mer Problem
Chapter #: 03
Problem ID: C
URL: https://beta.stepic.org/Bioinformatics-Algorithms-2/Greedy-Motif-Search-159/#step-3
'''

with open('data/textbook/rosalind_3c.txt') as input_data:
	dna = input_data.readline().strip()
	k = int(input_data.readline())
	profile = [map(float,line.strip().split()) if i!=0 else line.strip().split() for i,line in enumerate(input_data.readlines())]

# A dictionary relating nucleotides to their position within the profile.
nuc_loc = {nucleotide:index for index,nucleotide in enumerate(profile[0])}

# Initialize the maximum probabily.
max_prob = [-1, None]

# Compute the probability of the each k-mer, store it if it's currently a maximum.
for i in xrange(len(dna)-k+1):
	current_prob = 1
	for j, nucleotide in enumerate(dna[i:i+k]):
		current_prob *= profile[j+1][nuc_loc[nucleotide]]
	if current_prob > max_prob[0]:
		max_prob = [current_prob, dna[i:i+k]]

# Print and save the answer.
print max_prob[1]
with open('output/textbook/Textbook_03C.txt', 'w') as output_data:
	output_data.write(max_prob[1])

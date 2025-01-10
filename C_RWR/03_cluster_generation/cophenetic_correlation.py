import numpy as np
import nimfa
import sys
import os

#loading the signature matrix generated with the virus_signature.py script
f1=open("/results/signature_matrix.txt")
seq=f1.readline()
mat=[]
while(seq!=""):
	seq=seq.strip().split("\t")
	mat.append(np.array(seq[1:],float))
	seq=f1.readline()

mat=np.array(mat)

#performing the NMF with random initialization nd writing the cophenetic correlation
rank_cands = range(2,11, 1)
nmf = nimfa.Nmf(mat, seed='random', max_iter=100)
summary = nmf.estimate_rank(rank_range=rank_cands, n_run=100, what=["cophenetic"])
coph = [summary[rank]['cophenetic'] for rank in rank_cands]

output_dir = "/results/cophenetic_scores_100/"
os.makedirs(output_dir, exist_ok=True)

f1=open("/results/cophenetic_scores_100/"+sys.argv[1]+".txt","w")
f1.write("coph\t"+"\t".join([str(i) for i in coph])+"\n")
f1.close()

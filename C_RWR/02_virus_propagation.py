'''
This script is a modified version of the propagation_virus.py file from 
the paper 'Large-scale phage-based screening reveals extensive pan-viral 
mimicry of host short linear motifs' by Mihalič et al. Nature 2023, in 
order to run the RWR algorithm for each virus family instead of virus taxa.

Modified by: Minna Sayehban, student at Uppsala University
E-mail: minna.sayehban.1224@student.uu.se or minna.sayehban@gmail.com
Date: 26-11-2024
Last modified: 19-12-2024
'''

############### IMPORT & LOAD NECESSARY MODULES ################
from igraph import *
import numpy as np
from scipy.stats import *
import os,sys,scipy,math
from scipy.spatial.distance import jensenshannon
import fisher_test_virus
from win32 import win32file
win32file._setmaxstdio(2048)


############# GET STARTING DATA POINTS (SEED NODES) ##############
def load_seeds():

	f1=open(in_data_path+virus_family)
	seq=f1.readline()
	seeds={}
	while(seq!=""):
		seq= seq.strip().split("\t")
		seq[0]=float(seq[0])
		seeds[seq[1]]=seq[0]
		seq=f1.readline()

	uniprot_to_gene={}
	f1=open(uniprot,"r")
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split("\t")
		uniprot_to_gene[seq[0]]=seq[1].split(";")[0].strip()
		seq=f1.readline()
	return seeds,uniprot_to_gene


########### PERFORM RWR ON THE NETWORKS ############
def perform_rwr(network, results_path):
	number_of_nodes=network.vcount()

	graph_nodes=network.vs["name"]
	number_of_nodes=network.vcount()

	# CREATE VECTORS FOR EMPIRICAL RWR-SCORES AND P-VALUES
	empirical_values={}
	pvalues={}
	for i in network.vs:
		pvalues[i["name"]]=0.0
		empirical_values[i["name"]]=0.0

	# ASSIGN VALUES FOR SEED NODES
	reset_vertex=np.zeros(number_of_nodes)
	for j in network.vs.select(name_in=seeds.keys()):
		reset_vertex[j.index]=seeds[j["name"]]
	
	# RUN EMPIRICAL RWR
	pagerank=np.array(network.personalized_pagerank(reset=reset_vertex,directed=False, damping=damping, weights='weight'))

	# STORE VALUES FOR EMPIRICAL RUN
	for i in enumerate(graph_nodes):
		empirical_values[graph_nodes[i[0]]]=pagerank[i[0]]

	# RUN RANDOMISATION WITH 1000 NETWORK
	for ii in range(1000):
		file_path = os.path.join(f"{network_path}_random", f"{ii}.txt")
		print(f"Trying to read file: {file_path}")
		try:
			print(f"Processing file: {file_path}")
			network_random = Graph.Read_Ncol(file_path, weights=True, directed=False)
			print(f"Graph processed: {network_random.summary()}")
			del network_random
		except Exception as e:
			print(f"Error processing file: {file_path}. Error: {e}")
			continue		
		
		network_random = Graph.Read_Ncol(os.path.join(f"{network_path}_random", f"{ii}.txt"), weights=True, directed=False)
		random_nodes=network_random.vs["name"]
		reset_vertex=[0.0]*number_of_nodes
		for j in network_random.vs.select(name_in=seeds.keys()):
			reset_vertex[j.index]=seeds[j["name"]]

		# RUN RWR ON RANDOM NETWORK
		prandom=np.array(network_random.personalized_pagerank(reset=reset_vertex,directed=False, damping=damping, weights='weight'))

		# COUNT WHEN EMPIRICAL RWR-SCORE > RANDOM RWR-SCORE
		for i in enumerate(random_nodes):
			if empirical_values[i[1]]>prandom[i[0]]:
				pvalues[i[1]]+=1

	# SAVE RESULTS
	result_folder = os.path.join(results_path, os.path.splitext(virus_family)[0])
	os.makedirs(result_folder, exist_ok=True)

	f1=open(result_folder+"/rwr_scores.txt","w")
	for i in empirical_values:
		f1.write(i+"\t"+str(empirical_values[i])+"\n")
	f1.close()


	# ANNOTATE SIGNIFICANCE HITS 
	f1=open(result_folder+"/significance_hits.txt","w")
	enriched_proteins=[]
	for i in pvalues:
		if pvalues[i]>990:
			enriched_proteins.append(i)
		f1.write(i+"\t"+str(pvalues[i])+"\n")
	f1.close()

	# PERFORM ENRICHMENT ANALYSIS
	enriched_proteins=enriched_proteins+list(seeds.keys())
	fisher_folder= result_folder+"/fisher/"
	if not os.path.exists(fisher_folder):
		os.makedirs(fisher_folder)
	fisher_test_virus.load(list(set(enriched_proteins)),0.05, ["RT"],fisher_folder,list(seeds.keys()),uniprot_to_gene,"proteome")

	f1=open(result_folder+"/start_seeds.txt","w")
	for i in seeds:
		f1.write(i+"\t"+str(seeds[i])+"\n")
	f1.close()

if __name__ == '__main__':
	# SET PARAMETERS
	damping=0.7
	
	# SET VARIABLES
	virus_family=sys.argv[1]
	in_data_path=sys.argv[2]
	uniprot=sys.argv[3]
	network_path=sys.argv[4]
	results_path=sys.argv[5]

	# LOAD EMPIRICAL NETWORK
	network = Graph.Read_Ncol(network_path+".txt", weights=True, directed=False)
	graph_nodes=network.vs["name"]
	number_of_nodes=network.vcount()

	seeds,uniprot_to_gene=load_seeds()
	perform_rwr(network, results_path)
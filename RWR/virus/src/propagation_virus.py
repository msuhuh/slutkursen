'''
This script is a modified version of the propagation_virus.py file from 
the paper 'Large-scale phage-based screening reveals extensive pan-viral 
mimicry of host short linear motifs' by Mihalič et al. Nature 2023, in 
order to run the algorithm for each virus family instead.

Modified by: Minna Sayehban, student at Uppsala University
E-mail: minna.sayehban.1224@student.uu.se or minna.sayehban@gmail.com
Date: 26-11-2024
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

	f1=open(data_folder+test)
	seq=f1.readline()
	seeds={}
	while(seq!=""):
		seq= seq.strip().split("\t")
		seq[0]=float(seq[0])
		seeds[seq[1]]=seq[0]
		seq=f1.readline()

	uniprot_to_gene={}
	f1=open("../data/uniprot_to_gene.tab","r")
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split("\t")
		uniprot_to_gene[seq[0]]=seq[1].split(";")[0].strip()
		seq=f1.readline()
	return seeds,uniprot_to_gene


########### PERFORM RANDOM NETWORK ############
def perform_rwr(network):
	number_of_nodes=network.vcount()

	graph_nodes=network.vs["name"]
	number_of_nodes=network.vcount()

	# CREATE VECTORS FOR EMPIRICAL- AND P-VALUES
	empirical_values={}
	pvalues={}
	for i in network.vs:
		pvalues[i["name"]]=0.0
		empirical_values[i["name"]]=0.0

	# ASSIGN VALUES FOR SEED NODES
	reset_vertex=np.zeros(number_of_nodes)
	for j in network.vs.select(name_in=seeds.keys()):
		reset_vertex[j.index]=seeds[j["name"]]
	
	# FUNCTION FOR RWR
	pagerank=np.array(network.personalized_pagerank(reset=reset_vertex,directed=False, damping=damping, weights='weight'))

	# STORE VALUES FOR EMPIRICAL RUN
	for i in enumerate(graph_nodes):
		empirical_values[graph_nodes[i[0]]]=pagerank[i[0]]

	# RUN RANDOMISATION WITH 1000 NETWORK
	for ii in range(1000):
		file_path = os.path.join("..", "networks", f"{sim_type}_random", f"{ii}.txt")
		print(f"Trying to read file: {file_path}")
		try:
			print(f"Processing file: {file_path}")
			network_random = Graph.Read_Ncol(file_path, weights=True, directed=False)
			print(f"Graph processed: {network_random.summary()}")
			del network_random
		except Exception as e:
			print(f"Error processing file: {file_path}. Error: {e}")
			continue

		network_random = Graph.Read_Ncol(os.path.join("..", "networks", f"{sim_type}_random", f"{ii}.txt"), weights=True, directed=False)

		random_nodes=network_random.vs["name"]
		reset_vertex=[0.0]*number_of_nodes
		for j in network_random.vs.select(name_in=seeds.keys()):
			reset_vertex[j.index]=seeds[j["name"]]

		# RUN RWR ON RANDOM NETWORK
		prandom=np.array(network_random.personalized_pagerank(reset=reset_vertex,directed=False, damping=damping, weights='weight'))

		# COUNT WHEN EMPIRICAL RWR DISTR. > RANDOM RWR DISTR.
		for i in enumerate(random_nodes):
			if empirical_values[i[1]]>prandom[i[0]]:
				pvalues[i[1]]+=1

	os.makedirs(os.path.join(res_folder, test),exist_ok=True)

	# SAVE RESULTS
	f1=open(res_folder+test+"/rwr.txt","w")
	for i in empirical_values:
		f1.write(i+"\t"+str(empirical_values[i])+"\n")
	f1.close()


	# 990/1000 = 0.01
	f1=open(res_folder+test+"/pvalues.txt","w")
	enriched_proteins=[]
	for i in pvalues:
		if pvalues[i]>990:
			enriched_proteins.append(i)
		f1.write(i+"\t"+str(pvalues[i])+"\n")
	f1.close()

	# PERFORM ENRICHMENT ANALYSIS
	enriched_proteins=enriched_proteins+list(seeds.keys())
	folder= res_folder+test+"/fisher/"
	if not os.path.exists(folder):
		os.makedirs(folder)
	fisher_test_virus.load(list(set(enriched_proteins)),0.05, ["RT"],folder,list(seeds.keys()),uniprot_to_gene,"proteome")

	f1=open(res_folder+"/"+test+"/start_seeds.txt","w")
	for i in seeds:
		f1.write(i+"\t"+str(seeds[i])+"\n")
	f1.close()

if __name__ == '__main__':
	# SET PARAMETERS
	sim_type="TCSS"
	damping=0.7
	
	# SET VARIABLES
	data_folder=sys.argv[1]
	test=sys.argv[2]
	res_folder=sys.argv[3]

	# LOAD EMPIRICAL NETWORK
	network = Graph.Read_Ncol("../networks/"+sim_type+".txt", weights=True, directed=False)
	graph_nodes=network.vs["name"]
	number_of_nodes=network.vcount()

	seeds,uniprot_to_gene=load_seeds()
	perform_rwr(network)
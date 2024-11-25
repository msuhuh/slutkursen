'''
Thic code uses the .... 
Originally developed by: .... 
Modified by: Minna Sayehban 25-11-2024
'''


############ IMPORT & LOAD NECESSARY MODULES ################
from igraph import *
import numpy as np
from scipy.stats import *
import os,sys,scipy,math
from scipy.spatial.distance import jensenshannon
import fisher_test_virus
from win32 import win32file
win32file._setmaxstdio(2048)


############# GET THE VIRUS FAMILIES ############
def load_family_to_taxa_mapping(family_file):
    family_to_taxa = {}
    with open(family_file, 'r') as f:
        for line in f:
            family, taxa = line.strip().split("\t")
            if family not in family_to_taxa:
                family_to_taxa[family] = []
            family_to_taxa[family].append(taxa)
    return family_to_taxa


#load data
def load_seeds():
	data_path="../data/datasets/"+data_folder+"/"
	#load the seeds
	f1=open(data_path+test)
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

def perform_rwr(network):
	number_of_nodes=network.vcount()

	graph_nodes=network.vs["name"]
	number_of_nodes=network.vcount()

	#here I initialize the vectors where the empirical_values and pvalues are stores
	empirical_values={}
	pvalues={}
	for i in network.vs:
		pvalues[i["name"]]=0.0
		empirical_values[i["name"]]=0.0

	#here the values are assigned to the seed nodes
	reset_vertex=np.zeros(number_of_nodes)
	for j in network.vs.select(name_in=seeds.keys()):
		reset_vertex[j.index]=seeds[j["name"]]
	#this is the function for the RWR
	pagerank=np.array(network.personalized_pagerank(reset=reset_vertex,directed=False, damping=damping, weights='weight'))

	#storing the values of the empirical
	for i in enumerate(graph_nodes):
		empirical_values[graph_nodes[i[0]]]=pagerank[i[0]]

	#here we start the randomization (1000 network with the same seed nodes (is the same function above but looping around 1000 random network)
	for ii in range(1000):
		file_path = os.path.join("..", "networks", f"{sim_type}_random", f"{ii}.txt")
		print(f"Trying to read file: {file_path}")
		try:
			print(f"Processing file: {file_path}")
			network_random = Graph.Read_Ncol(file_path, weights=True, directed=False)
			print(f"Graph processed: {network_random.summary()}")
			del network_random  # Free memory immediately
		except Exception as e:
			print(f"Error processing file: {file_path}. Error: {e}")
			continue

		network_random = Graph.Read_Ncol(os.path.join("..", "networks", f"{sim_type}_random", f"{ii}.txt"), weights=True, directed=False)

		random_nodes=network_random.vs["name"]
		reset_vertex=[0.0]*number_of_nodes
		for j in network_random.vs.select(name_in=seeds.keys()):
			reset_vertex[j.index]=seeds[j["name"]]
		#RWR on random network
		prandom=np.array(network_random.personalized_pagerank(reset=reset_vertex,directed=False, damping=damping, weights='weight'))

		#here I count how many times the empirical rwr distribution is greater than the random
		for i in enumerate(random_nodes):
			if empirical_values[i[1]]>prandom[i[0]]:
				pvalues[i[1]]+=1

	os.makedirs(os.path.join(res_folder, test),exist_ok=True)

	#saving the results
	f1=open(res_folder+"/"+test+"/rwr.txt","w")
	for i in empirical_values:
		f1.write(i+"\t"+str(empirical_values[i])+"\n")
	f1.close()


	f1=open(res_folder+"/"+test+"/pvalues.txt","w")
	enriched_proteins=[]
	for i in pvalues:
		#990 corresponds to a pvalue<0.01
		if pvalues[i]>990:
			enriched_proteins.append(i)
		f1.write(i+"\t"+str(pvalues[i])+"\n")
	f1.close()

	#performing enrichement analysis
	enriched_proteins=enriched_proteins+list(seeds.keys())
	folder= res_folder+"/"+test+"/fisher/"
	if not os.path.exists(folder):
		os.makedirs(folder)
	fisher_test_virus.load(list(set(enriched_proteins)),0.05, ["RT"],folder,list(seeds.keys()),uniprot_to_gene,"proteome")

	f1=open(res_folder+"/"+test+"/start_seeds.txt","w")
	for i in seeds:
		f1.write(i+"\t"+str(seeds[i])+"\n")
	f1.close()

if __name__ == '__main__':
	#parameters
	sim_type="TCSS"
	damping=0.7
	#this is the file containing the uniprot id involved in a side effect
	test=sys.argv[1]
	res_folder=sys.argv[2]
	data_folder=sys.argv[3]
	family_file = sys.argv[4]

	family_to_taxa = load_family_to_taxa_mapping(family_file)

	# loading empirical network
	network = Graph.Read_Ncol("../networks/"+sim_type+".txt", weights=True, directed=False)
	graph_nodes=network.vs["name"]
	number_of_nodes=network.vcount()

	seeds,uniprot_to_gene=load_seeds()
	perform_rwr(network)

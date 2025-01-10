from scipy.stats import fisher_exact
import sys
def load(protein_list,threshold,component,path_def,starting_proteins,uniprot_to_gene,background):
	protein=protein_list
	path="../data/"+background+"/"
	temp={}
	for ii in component:
		descr={}
		temp[ii]=[]
		f1=open(path+ii+"_descr.txt","r")
		seq=f1.readline()
		while (seq!=""):
			seq=seq.strip().split("\t")
			descr[seq[0]]=seq[1]
			seq=f1.readline()
		f1.close()
		f1=open(path+ii+".txt","r")
		seq=f1.readline()
		fisher={}
		fisher_count={}
		while(seq!=""):
			seq=seq.strip().split("\t")
			fisher[seq[0]]=seq[1:]
			seq=f1.readline()
		f1.close()
		f1=open(path+ii+"_count.txt","r")
		seq=f1.readline()
		while(seq!=""):
			seq=seq.strip().split("\t")
			fisher_count[seq[0]]=int(seq[1])
			seq=f1.readline()
		f1.close()
		fisherset={}
		fisherset_count=[]
		protein_annotation={}
		for i in protein:
			if i in fisher:
				for j in fisher[i]:
					if j in protein_annotation:
						protein_annotation[j].append(i)
					else:
						protein_annotation[j]=[]
						protein_annotation[j].append(i)

					if j in fisherset:
						fisherset[j]=fisherset[j]+1
					else:
						fisherset[j]=1

			else:
				fisherset_count.append(i)


		totalfisher=len(fisher)
		numberofproteins=len(protein)-len(list(set(fisherset_count)))
		fisher={}
		fisher_value_ic={}
		fisher_value={}
		fisher_value_no={}
		lenfisherset=len(fisherset)
		for i in fisherset:

			a=fisherset[i]
			b=numberofproteins-a
			c=fisher_count[i]-a
			d=totalfisher-a-b-c
			table=[[a,b],[c,d]]
			fisher[i]=fisher_exact(table,alternative ="greater")[1]
			if fisher[i]<(threshold/lenfisherset):
				if fisher[i] in fisher_value:
					fisher_value[fisher[i]].append(i)
				else:
					fisher_value[fisher[i]]=[]
					fisher_value[fisher[i]].append(i)

		f2=open(path_def+ii+"fisher.txt","w")
		for i in sorted(fisher_value):
			for j in fisher_value[i]:
				#temp_gene=[]
				temp=list(set(protein_list).intersection(set(protein_annotation[j])))
				#	temp_gene.append(uniprot_to_gene.get(k))
				f2.write(j+"\t"+str(i)+"\t"+descr[j]+"\t"+"\t".join(temp)+"\n")

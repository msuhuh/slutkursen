import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from os.path import exists

#code used to generate supp Figure 3
values={}
val=[]
for i in range(0,100,1):

	if exists("/results/cophenetic_scores_100/"+str(i)+".txt"):
		f1=open("/results/cophenetic_scores_100/"+str(i)+".txt")
		seq=f1.readline()

		while(seq!=""):
			seq=seq.strip().split("\t")
			if seq[0]=="coph":
				temp=np.array(seq[1:],float)
				for j in enumerate(temp):
					val.append([j[0]+2,j[1]])
			seq=f1.readline()


df = pd.DataFrame(val, columns = ['number of clusters', 'cophenetic correlation coefficient'])
plt.rcParams['font.family'] = 'Arial'
sns.set_theme(style="whitegrid")
sns.set_context("paper", font_scale=1.5, rc={ 'lines.markersize': 4.0,'grid.linewidth':0.5,"lines.linewidth": 0.5, 'lines.linewidth': 0.25,'ytick.major.size': 1.5})
sns.violinplot(x="number of clusters", y="cophenetic correlation coefficient", data=df)
plt.savefig("/results/cophenetic_100.pdf")
plt.show()

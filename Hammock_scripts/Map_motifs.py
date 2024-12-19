import pandas as pd 

# Script to map expected human motifs to viral hits for each human bait id 



def map_help(hit, motif):
    if len(hit) == 0 or len(motif) == 0:
        return ""  
    elif motif[0] == "[":
        i = 1
        options = ""
        while motif[i] != "]":
            options += motif[i]
            i += 1  
        if hit[0] in options:
            return hit[0] + map_help(hit[1:], motif[(i+1):]) 
        else:
            return ""
    elif motif[0] == ".":
        return hit[0] + map_help(hit[1:], motif[1:])
    elif hit[0] == motif[0]:
        return hit[0] + map_help(hit[1:], motif[1:])
    else: 
        return ""
        
def score_map(hit, motif, match, score):
    if len(motif) == 0:
        return([match, score]) 
    elif len(hit) == 0: 
        score += 1
        return([match, score]) 
    else: 
        match += hit[0]
    if motif[0] == "[":
        i = 1
        options = ""
        while motif[i] != "]":
            options += motif[i]
            i += 1  
        if hit[0] not in options:
            score += 0.5
        return score_map(hit[1:], motif[(i+1):], match, score) 
    elif motif[0] != "." and hit[0] != motif[0]:
         score += 1 
    return score_map(hit[1:], motif[1:], match, score) 
       

def len_mot(word):
    if word == '*':
        return 0
    elif len(word) == 0:
        return 0
    elif word[0] == "[":
        a = 1
        while word[a] != "]":
                a += 1 
        return 1 + len_mot(word[(a+1):])
    else:
        return 1 + len_mot(word[1:])


def n_set(motif):
    if motif == '*':
        return 0
    elif len(motif) == 0:
        return 0
    elif motif[0] == "[":
        a = 1
        while motif[a] != "]":
                a += 1 
        return 1 + n_set(motif[(a+1):])
    elif motif[0] == '.':
        return n_set(motif[1:]) 
    else:
        return 1 + n_set(motif[1:])


def hits(hit, mot):
    mot_l = len_mot(mot)
    hit_list = []
    for i in range (0, (len(hit) - (mot_l) + 1)):
        hit_list.append(hit[i:(mot_l+i)])
    return hit_list

def slims(hit_long, mot):
    lst = hits(hit_long, mot)
    match_dict = {}
    for h in lst: 
        match = score_map(h, mot, "", 0)
        hit = match[0]
        score = match[1]
        n = n_set(mot)
        if n == 0:
            rel_err = 1
        else: 
            rel_err = score/n
            if rel_err <= 0.3:
                match_dict[hit] = rel_err 
    if match_dict: 
        min_key = min(match_dict, key=match_dict.get)
        min_score = match_dict[min_key]
        res = [min_key, min_score]
        return res
    else:
        return ['*','*']
    

 
bait_file = '20241119_human_baits_metadata.xlsx' # Replace with your path to human baits metadata, file containing the expected human motifs
viral_file = '20241119_Compiled_viral_results_compiled.xlsx' # Relace with your path to viral results, file containing the viral - human bait bindings and viral peptide sequence hit

# Create dataframes 

df1 = pd.read_excel(bait_file, usecols=['bait_id', 'slimfinder_motif_regex']) # Create dataframe with the human bait metadata, only need collumns bait id that contains the bait if and slimfinder_motif_regex that contains the expected motif
df2 = pd.read_excel(viral_file, usecols=['human_bait_id', 'virus_collapsed_hit']) # Create dataframe with human bait id and cirus collapsed hit from virus data. 
df3 = pd.read_excel(viral_file)  # Create data with entire viral file

df1 = df1[df1['slimfinder_motif_regex'].notna() & (df1['slimfinder_motif_regex'] != '')] # Removes the columns that lack motifs from the metadata
bait_to_motif = df1.set_index('bait_id')['slimfinder_motif_regex'].to_dict() # Create a dictionary with the human bait ids as keys and the motifs to be mapped as values


motifs = [] # Create empty motif list
scores = [] # Create empty score list
slimfinder_motif_regexs =[] # Create empty motif list

for _, row in df2.iterrows():
    bait_id = row['human_bait_id']
    virus_hit = row['virus_collapsed_hit']
    
    # Find the corresponding motif regex for the same bait_id
    slimfinder_motif = bait_to_motif.get(bait_id, '*') 
    res = slims(virus_hit, slimfinder_motif)# Default to None if bait_id not found
    motif = res[0]
    score = res[1]
    # Append the result for further processing
    motifs.append({'new_motif': motif})
    scores.append({'new_scores': score})
    slimfinder_motif_regexs.append({'slimfinder_motif_regex': slimfinder_motif})


slimfinder_motif_regexs_df = pd.DataFrame(slimfinder_motif_regexs)
motifs_df = pd.DataFrame(motifs)
scores_df = pd.DataFrame(scores)

# Add the new_motif column to df2
df3['new_motif'] = motifs_df['new_motif']
df3['new_scores'] = scores_df['new_scores']
df3['slimfinder_motif_regex'] = slimfinder_motif_regexs_df['slimfinder_motif_regex']

# Save the updated DataFrame to a new Excel file
df3.to_excel('20241122_motif_mapped_data.xlsx', index=False)

# Print the result to check
print(df3.head())

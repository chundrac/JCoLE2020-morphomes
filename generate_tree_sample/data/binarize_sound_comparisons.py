from collections import defaultdict
import numpy as np

lang_replace = [l.strip('\n').split('\t') for l in open('ASJP_merged.csv','r')]
lang_replace = {l[1]:l[3] for l in lang_replace[1:]}

text = [l.strip('\n').split('\t') for l in open('cognates_wide.tsv','r')]

clts = [l.strip('\n').split('\t') for l in open('clts_graphemes.csv','r')]
clts_values = {}
for l in clts[1:]:
    clts_values[l[0]]=l[1].split()
    if 'sibilant' in clts_values[l[0]]:
        clts_values[l[0]].pop(clts_values[l[0]].index('sibilant'))
    if 'lateral' in clts_values[l[0]]:
        clts_values[l[0]].pop(clts_values[l[0]].index('lateral'))

shortconcepts = ['I', 'blood', 'bone', 'breast', 'come', 'die', 'dog', 'drink', 'ear', 'eye', 'fire', 'fish', 'full', 'hand', 'hear', 'horn', 'knee', 'leaf', 'liver', 
            'louse', 'mountain', 'name', 'new', 'night', 'nose', 'one', 'path', 'person', 'see', 'skin', 'star', 'stone', 'sun', 'tongue', 'tooth', 'tree', 'two', 'water', 'we', 'you']

conceptinds = [text[0].index(w) for w in shortconcepts if w in text[0]]

text = [[l[0]]+[l[i] for i in conceptinds] for l in text]

ngrams = defaultdict(list)

for l in text[1:]:
    for i,w in enumerate(l):
        if i > 0:
            concept = text[0][i]
            #w_ = ['#']+list(w)+['$']
            w_ = list(w)
            #for j in range(len(w_)-2):
            #    ngrams[l[0]].append(concept+':'+''.join(w_[j:j+2])
            #for j in range(len(w_)-1):
            #    ngrams[l[0]].append(concept+':'+''.join(w_[j:j+1]))
            for s in w_:
                if s in clts_values.keys():
                    ngrams[l[0]].append(concept+':'+'|'.join(clts_values[s][1:]))
                else:
                    ngrams[l[0]].append(concept+':'+s)

characters = sorted(set([v for k in ngrams.keys() for v in ngrams[k]]))

languages = list(ngrams.keys())
if 'LG_FILEPATHPART' in languages:
    languages.pop(languages.index('LG_FILEPATHPART'))

bin_data = np.zeros([len(languages),len(characters)])

for i,l in enumerate(languages):
    for j,w in enumerate(characters):
        if w in ngrams[l]:
            bin_data[i,j] = 1

f = open('soundcomp_romance_binarized.csv','w')
print('\t'.join(['language']+characters),file=f)
for k in lang_replace.keys():
    l = lang_replace[k]
    if l in languages:
        print('\t'.join([k]+[str(s) for s in list(bin_data[languages.index(l),])]),file=f)
    else:
        print(k,l)

f.close()
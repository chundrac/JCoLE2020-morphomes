import pandas as pd
import numpy as np
from ete3 import Tree
from collections import defaultdict

def nexCharOutput(chMtx, outfile, datatype='STANDARD'):
    names = chMtx.index
    with open(outfile, 'w') as f:
        f.write('#NEXUS\n\n')
        f.write('BEGIN DATA;\n')
        f.write('DIMENSIONS ntax=' + str(len(chMtx)) +
                ' NCHAR='+str(len(chMtx.T))+';\n')
        f.write('FORMAT DATATYPE=' + datatype +
                ' GAP=? MISSING=- interleave=yes;\n')
        f.write('MATRIX\n\n')
        txLgth = max(map(len, names))
        for i in range(len(chMtx)):
            f.write(names[i].ljust(txLgth+2))
            for ch in chMtx.values[i]:
                if ch == -1:
                    ch = '-'
                else:
                    ch = str(ch)
                f.write(ch)
            f.write('\n')
        f.write('\n;\n\nEND;\n')
        f.close()


nexus_data = pd.read_csv('../data/all_data_binarized.csv',
                     na_filter=False,index_col=0)
                     
part_ind = [i for i,s in enumerate(list(nexus_data.columns)) if '.' in s][0]

ccMtx = nexus_data.iloc[:,:part_ind]
#ccMtx.loc[ccMtx.index != 'Italic_Latino-Faliscan_Latin',]
ccMtx.drop('Italic_Latino-Faliscan_Latin',inplace=True)
scMtx = nexus_data.iloc[:,part_ind:]
#scMtx.loc[scMtx.index != 'Italic_Latino-Faliscan_Latin',]
scMtx.drop('Italic_Latino-Faliscan_Latin',inplace=True)

nexCharOutput(ccMtx, '../data/cc.nex', datatype="Standard")
nexCharOutput(scMtx, '../data/sc.nex', datatype="Standard")

def nname(x):
    if x.is_leaf():
        return '"'+x.name+'"'
    else:
        return x.name


#ctree = Tree('../data/atlantic-constraint-tree-relaxed.newick',8)
#Italic_Latino-Faliscan_Latin

taxa = list(scMtx.index)
clades = defaultdict(list)
for t in taxa:
    clades[t.split('_')[0]].append(t)

skeleton_tree = '(\n'+'\n,\n'.join(['('+','.join([t for t in clades[k]])+')' for k in clades.keys()])+'\n);'

f = open('../data/skeleton_tree.txt','w')
print(skeleton_tree,file=f)
f.close()

ctree = Tree('../data/constraint_tree.newick',8)
ctree.ladderize()

#ctree.prune(scMtx.index)

assert(sorted(ctree.get_leaf_names())==sorted(scMtx.index))

nodes = np.array([nd for nd in ctree.get_descendants()
                  if not nd.is_root() and not nd.is_leaf()])

#for i, nd in enumerate(nodes):
#    nd.name = ''

for i, nd in enumerate(nodes):
    if nd.name == '':
        nd.name = 'clade'+str(i+1).rjust(2, '0')


#print([nn.name for nn in nodes])

rev = ""
for nd in reversed(nodes):
    rev += nd.name + " = clade("
    rev += ', '.join([nname(x) for x in nd.get_children()])
    rev += ")\n"

rev += 'constraints = ['+', '.join([nname(nd) for nd in reversed(nodes)])
rev += ']\n'

#print(rev)

with open('constraints.Rev', 'w') as f:
    f.write(rev)

dating = [l.strip('\n').split('\t') for l in open('../data/romance_node_date_constraints.tsv','r')]
dates = {}
for l in dating[1:]:
    dates[l[0]] = (l[1],l[2])

rev = ""
for nd in nodes:
    if nd.name in dates.keys():
        rev += "{}_age := tmrca(tree, {})\n".format(nd.name,nd.name)
        rev += "obs_{}_age ~ dnUnif({}, {})\n".format(nd.name,dates[nd.name][1],dates[nd.name][0])
        rev += "obs_{}_age.clamp({}_age)\n".format(nd.name,nd.name)
        rev += '\n'

for nd in ctree.get_leaf_names():
    rev += '#{}_age := tmrca(tree, clade("{}"))\n'.format(nd.replace('-',''),nd)
    rev += '#obs_{}_age ~ dnUnif({}, {})\n'.format(nd.replace('-',''),300,800)
    rev += '#obs_{}_age.clamp({}_age)\n'.format(nd.replace('-',''),nd.replace('-',''))
    rev += '\n'

with open('calibrations.Rev', 'w') as f:
    f.write(rev)
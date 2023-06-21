
import pandas as pd
import collections
import re
import sys
import networkx as nx
import mygene
from networkx.algorithms.dag import transitive_reduction    
from collections import Iterable
from six import string_types
from functools import reduce
sys.path.append("/booleanfs2/sahoo/Hegemon/")
sys.path = ["/booleanfs2/sahoo/BoNE/"] + sys.path
import HegemonUtil as hu
import bone

try:
    reload  # Python 2.7
except NameError:
    try:
        from importlib import reload  # Python 3.4+
    except ImportError:
        from imp import reload  # Python 3.0 - 3.3



reload(hu)

reload(bone)

db = hu.Database("/booleanfs2/sahoo/Hegemon/explore.conf")
h = hu.Hegemon(db.getDataset("GL1"))
h.init()
h.initPlatform()


def init():
    global NAG #
    NAG = []
    global saved_shared

    global saved_dys

    global edge_write

    global nodes_in_path

    


    global g_right

    
    



# find the best probeID for each gene
def BESTID(LIST): 

    IDS = {}
    global NAG #not available genes, try to see if any of its aliases is available
    
    for i in LIST:
        
        if  i=='SMO':
            IDS[i+'#29'] = '218629_at'
        elif i=='BMPR1A':
            IDS[i + '#32'] = '204832_s_at'
        elif i=='ACVRL1':
            IDS[i + '#950'] = '226950_at'
                       
        elif i== 'FN1':
            IDS[i + '#199'] = '1558199_at'             

        elif i=='IHH':
            IDS[i + '#358'] = '229358_at'

        elif i=='FBXO7':
            IDS[i + '#178'] = '201178_at'

        elif i=='ANGPT1':
            IDS[i + '#939'] = '1552939_at'
        elif i=='ANXA5':
            IDS[i + '#782'] = '200782_at'
  
        elif i=='RUNX1':
            IDS[i + '#129'] = '208129_x_at'

        else:
            j=h.getBestID(h.getIDs(i).keys())
            if j is not None:

                IDS[i+'#'+j.split('_')[0][-3:]]=j
            else:
                #print("This gene name is not available:",i)
                NAG.append(i)

    return IDS


# find Boolean relationships between genes and store the matrix as a file
def relmatrix(list1,list2,net,n):

    tbM = collections.defaultdict(dict)
    I=0
    for k_i, i in list1.items():
        ru = bone.BIGraph.getNCodes(net, [i])

        for k_j, j in list2.items():
            rel=[k for k in set(ru.keys()) if j in ru[k].keys()]
            if len(rel)==0:
                tbM[k_i][k_j]=0

            elif len(rel)==1:
                tbM[k_i][k_j]=rel[0]

                
        tbM[k_i][k_i]=5

    rel_matrix = pd.DataFrame(tbM)

    rel_matrix.to_csv('./'+n, header=True, sep='\t')
    return  rel_matrix



# sort the genes within a cluster alphabetically.
def sortGroup(str): 
    if '_' in str:
        arr = str.split('_')
        arr = sorted(arr)
    else:
        return str
    return '_'.join(arr)




# Search for the first duplicate value
def find_first_duplicate(nums):
    num_set = []
    no_duplicate = -1

    for i in range(len(nums)):

        if nums[i].name.split('#')[0] in num_set:
            return i
        else:
            num_set.append(nums[i].name.split('#')[0])

    return no_duplicate




#take a list of gene symbols as input and retrieve their aliases using the MyGene
def add_aliases_set(L):

    mg = mygene.MyGeneInfo()

    geneSymbol = mg.querymany(L, scopes='symbol', fields='alias', species='human')
    geneSymbol

    result = {}

    for d in geneSymbol:
        result[d["query"]]=[d["query"]]

        if 'alias' in d:
            i=d['alias'] 
            if isinstance(i, string_types):
                result[d["query"]].append(i)
            else:
                
                result[d["query"]]+=d['alias'] 
     
    return result




class Gene:
    def __init__(self, name, cell_type,color):
        self.name = name #gene symbol+#+ the last three digits of probeID
        self.ct = cell_type # functional, dysfunctional or shared
        self.edges = collections.defaultdict(list) #the relationships between the current gene and other genes
        self.color=color # Represents the color associated with the gen

    def relate_to(self, gene, rel): #establishes a relationship between the two genes 
        if gene not in self.edges[(rel, gene.ct)]:
            self.edges[(rel, gene.ct)].append(gene)

    def path(self, rel, type, i, pt, prev):

        for gn in self.edges[(rel, type)]:
            if gn not in prev and gn not in pt:
                pt.append(gn)
                gn.path(rel, type, i + 1, pt, prev)
        return pt

    def dfs(self, rel, type, i, p, prev, f,f_n='nlolo'):

        for gn in self.edges[(rel, type)]:
            if gn.name not in [k.name for k in p] and gn.name not in [k.name for k in prev]:
                if f_n =='lolo':
                    f(p, gn)
                    gn.dfs(rel, type, i + 1, p + [gn], prev, f, f_n)
                else:
                    if gn.name.split('#')[0]  not in [k.name.split('#')[0] for k in p]:
                        f(p, gn)
                        gn.dfs(rel, type, i + 1, p + [gn], prev, f, f_n)


    def dfs_eq(self, rel, type, i, p, prev, l_eqv):

        for gn in self.edges[(rel, type)]:
            if gn not in p and gn.name not in [k.name for k in prev]:
                info = (p[-1].name, gn.name)
                if info not in l_eqv:
                    l_eqv.append(info)
                gn.dfs_eq(rel, type, i + 1, p + [gn], prev + [gn], l_eqv)

    def all_eq_dif(self, rel, type_1, type_2, i, pt, prev):
        for gn in self.edges[(rel, type_2)]:
            if gn not in prev and gn not in pt:
                pt.append(gn)
        return pt

    def pretty(self):
        return '(' + self.name + ', ' + self.ct + ')'


def create(table, gr, per): #gr:functional, per:dysfunctional
    genes = {}
    for f, t, r in table: #f:from, t:to, r:Boolean relationship
        if f in genes:
            g1 = genes[f]
        else:

            if f in gr: #gr or functional
                type = 'func'
                color='lightgreen'
            elif f in per: # per or dysfunctional
                type = 'dysfunc'
                color='skyblue' 
            else:
                type = 'shared'
                color='darksalmon'   
        
            g1 = Gene(f, type,color)
            genes[f] = g1

        if t in genes:
            g2 = genes[t]
        else:
            type = 'func' if t in gr else 'dysfunc' if t in per else 'shared'
            color = 'lightgreen' if t in gr else 'skyblue' if t in per else 'darksalmon'

            g2 = Gene(t, type, color)
            genes[t] = g2

        g1.relate_to(g2, r)
        if r == '2':
            g2.relate_to(g1, '3')
        elif r == '3':
            g2.relate_to(g1, '2')
        else:
            g2.relate_to(g1, r)
    return genes
    
    
    
    
    
    
    


def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)
        


def find_final_clusters(genes):  
    #gr:functional, per:dysfunctional
    eqv_edges_shared = []
    eqv_edges_gr = []
    eqv_edges_per = []
    eqv_sh_gr = []
    eqv_gr_per = []
    eqv_sh_per = []

    for gg in genes.keys():
        g1 = genes[gg]
        print(g1.pretty())

        if g1.ct == 'func':
            g1.dfs_eq('5', g1.ct, 1, [g1], [], eqv_edges_gr)
            gr_sh = g1.all_eq_dif('5', g1.ct, 'shared', 0, [], [])
            if len(gr_sh) != 0:
                for k in gr_sh:
                    eqv_sh_gr.append((g1.name, k.name))
            gr_p = g1.all_eq_dif('5', g1.ct, 'dysfunc', 0, [], [])
            if len(gr_p) != 0:
                for k in gr_p:
                    eqv_gr_per.append((g1.name, k.name))
        if g1.ct == 'shared':
            g1.dfs_eq('5', g1.ct, 1, [g1], [], eqv_edges_shared)
            sh_gr = g1.all_eq_dif('5', g1.ct, 'func', 0, [], [])
            if len(sh_gr) != 0:
                for k in sh_gr:
                    eqv_sh_gr.append((g1.name, k.name))

            sh_pr = g1.all_eq_dif('5', g1.ct, 'dysfunc', 0, [], [])
            if len(sh_pr) != 0:
                for k in sh_pr:
                    eqv_sh_per.append((g1.name, k.name))
        if g1.ct == 'dysfunc':

            g1.dfs_eq('5', g1.ct, 1, [g1], [], eqv_edges_per)
            sh_pr = g1.all_eq_dif('5', g1.ct, 'shared', 0, [], [])
            if len(sh_pr) != 0:
                for k in sh_pr:
                    eqv_sh_per.append((g1.name, k.name))

            gr_p = g1.all_eq_dif('5', g1.ct, 'func', 0, [], [])
            if len(gr_p) != 0:
                for k in gr_p:
                    eqv_gr_per.append((g1.name, k.name))

    g = nx.Graph(eqv_gr_per)
    rm_eqv_gr_per = []
    for subgraph in connected_component_subgraphs(g):
        rm_eqv_gr_per.append(list(subgraph.nodes()))

    g = nx.Graph(eqv_sh_per)
    rm_eqv_sh_per = []
    for subgraph in connected_component_subgraphs(g):
        rm_eqv_sh_per.append(list(subgraph.nodes()))

    g = nx.Graph(eqv_sh_gr)
    rm_eqv_sh_gr = []
    for subgraph in connected_component_subgraphs(g):
        rm_eqv_sh_gr.append(list(subgraph.nodes()))

    g = nx.Graph(eqv_sh_gr + eqv_gr_per + eqv_sh_per)
    rm_eqv_sh_gr_per = []
    for subgraph in connected_component_subgraphs(g):
        rm_eqv_sh_gr_per.append(list(subgraph.nodes()))
    # ----------------------
    g = nx.Graph(eqv_edges_per)
    mg_eqv_per = []
    for subgraph in connected_component_subgraphs(g):
        mg_eqv_per.append(list(subgraph.nodes()))

    g = nx.Graph(eqv_edges_shared)
    mg_eqv_sh = []
    for subgraph in connected_component_subgraphs(g):
        mg_eqv_sh.append(list(subgraph.nodes()))

    g = nx.Graph(eqv_edges_gr)
    mg_eqv_gr = []
    for subgraph in connected_component_subgraphs(g):
        mg_eqv_gr.append(list(subgraph.nodes()))

    print('done mg')
    print('mg_eqv_gr=', mg_eqv_gr)
    print('mg_eqv_sh=', mg_eqv_sh)
    print('mg_eqv_per=', mg_eqv_per)

    res_sh_gr = []
    res_sh_per = []
    for cluster in rm_eqv_sh_gr:
        gr_part = []
        sh_part = []
        for c in cluster:
            info1 = [mg_eqv_gr[i] for i, gene_n in enumerate(mg_eqv_gr) if c in gene_n]
            if len(info1) != 0:
                gr_part.append(info1[-1])
            info2 = [mg_eqv_sh[i] for i, gene_n in enumerate(mg_eqv_sh) if c in gene_n]
            if len(info2) != 0:
                sh_part.append(info2[-1])
        l = gr_part + sh_part + [cluster]
        res_sh_gr.append(list(set(reduce(lambda x, y: x + y, l))))

    for cluster in rm_eqv_sh_per:
        per_part = []
        sh_part2 = []
        for c in cluster:
            info1 = [mg_eqv_per[i] for i, gene_n in enumerate(mg_eqv_per) if c in gene_n]
            if len(info1) != 0:
                per_part.append(info1[-1])
            info2 = [mg_eqv_sh[i] for i, gene_n in enumerate(mg_eqv_sh) if c in gene_n]
            if len(info2) != 0:
                sh_part2.append(info2[-1])
        l = per_part + sh_part2 + [cluster]
        res_sh_per.append(list(set(reduce(lambda x, y: x + y, l))))

    res_sh_per_gr = []
    res_per_gr = []
    for cluster in rm_eqv_sh_gr_per:
        per_part2 = []
        gr_part2 = []
        sh_part3 = []
        for c in cluster:
            info1 = [mg_eqv_per[i] for i, gene_n in enumerate(mg_eqv_per) if c in gene_n]
            if len(info1) != 0:
                per_part2.append(info1[-1])
            info2 = [mg_eqv_gr[i] for i, gene_n in enumerate(mg_eqv_gr) if c in gene_n]
            if len(info2) != 0:
                gr_part2.append(info2[-1])
            info3 = [mg_eqv_sh[i] for i, gene_n in enumerate(mg_eqv_sh) if c in gene_n]
            if len(info3) != 0:
                sh_part3.append(info3[-1])
        l = per_part2 + gr_part2 + [cluster]
        l2 = per_part2 + gr_part2 + sh_part3 + [cluster]
        res_sh_per_gr.append(list(set(reduce(lambda x, y: x + y, l2))))

        res_per_gr.append(list(set(reduce(lambda x, y: x + y, l))))



    return res_sh_per, res_sh_gr, res_sh_per_gr, res_per_gr
    
#-----------------
# update the edges based on the res_sh_dys_fn (genes clustered together) and ensure that the clustered nodes inherit all relations that their members have to other genes.
#-----------------    
    
def Update_edges_after_clustering(moved_and_shared_genes,res_sh_dys_fn, edges):
    
    
    edges_new = []
    for f,t,r in edges:
        if r!='5':
            if f not in moved_and_shared_genes and t not in moved_and_shared_genes:
                edges_new.append((f,t,r))


    #-----------------------
    genes_cluster = collections.defaultdict()   #it creates a new name for a clustered node by joining the individual genes inside the cluster with underscores and sorting the genes' symbols within a cluster alphabetically.

    for i in res_sh_dys_fn:
        new_name = "_".join(i)
        for k in i:
            genes_cluster[k] = sortGroup(new_name)


    for i in res_sh_dys_fn:
        new_name = "_".join(i)
        new_name =sortGroup(new_name)
        for k in i:
            for f,t,r in edges:
                if k == f:
                    if t in genes_cluster.keys():
                        if genes_cluster[t] != new_name:
                            edges_new.append((new_name, genes_cluster[t], r))
                    else:
                        edges_new.append((new_name, t, r))
                elif k == t:
                    if f in genes_cluster.keys():
                        if genes_cluster[f] != new_name:
                            edges_new.append((genes_cluster[f], new_name, r))
                    else:
                        edges_new.append((f, new_name, r))
                        
                        
    return edges_new,genes_cluster







def lolo_path(p_hh_fn, fn_hl_shrd_p_ll_dysfn, fn_ll_dysfn, t_m):
    global edge_write
    global nodes_in_path

    global saved_dys#=saved_per
    global saved_shared

    if t_m == 'hilo_shared':
        last_shrd = fn_hl_shrd_p_ll_dysfn[-1]

        last_shrd.dfs('2', 'dysfunc', 1, fn_hl_shrd_p_ll_dysfn, [], lambda p, g: lolo_path([], p, g, t_m='hilo'), 'lolo')

        for fn_ll_shrd in last_shrd.edges[('2', 'shared')]:
            if fn_ll_shrd not in fn_hl_shrd_p_ll_dysfn:
                lolo_path(p_hh_fn, fn_hl_shrd_p_ll_dysfn + [fn_ll_shrd], [], t_m)
    if t_m == 'hilo':
        if len(fn_ll_dysfn.edges[('2', 'dysfunc')]) == 0 :

            jj = [gl for gl in fn_hl_shrd_p_ll_dysfn + [fn_ll_dysfn]]  
            ind=find_first_duplicate(jj)
            if ind != -1:
                jj=jj[0:ind]
            if len(jj)<2:# 
                return False


            list_path_edges = p_hh_fn + fn_hl_shrd_p_ll_dysfn + [fn_ll_dysfn]
            list_path_edges = p_hh_fn + jj
            if fn_hl_shrd_p_ll_dysfn[0].ct=='dysfunc':
                
                saved_dys.append(fn_hl_shrd_p_ll_dysfn[0].name)
            else:
                saved_shared.append(fn_hl_shrd_p_ll_dysfn[0].name)
            

            # print("------------EDGES-------------------")
            for i in range(len(list_path_edges) - 1):
                r = ('3' if i < len(p_hh_fn) - 1 else '4' if i == len(p_hh_fn) - 1 else '2')
                info = ','.join([list_path_edges[i].name, list_path_edges[i + 1].name, r])
                if info not in edge_write:
                    edge_write.append(info)
                if list_path_edges[i].name not in nodes_in_path:
                    nodes_in_path.append(list_path_edges[i].name)
                if list_path_edges[i + 1].name not in nodes_in_path:
                    nodes_in_path.append(list_path_edges[i + 1].name)

    return True

    
########################APPLY transitive reduction on lolo pathways #########################    
def transitive_reduction_part(res_sh_dys_fn,genes_cluster):
    global edge_write
    global nodes_in_path
    eqv_edges_cluster = []
    node_ev = []

    
    ###############START
    for i in res_sh_dys_fn:
        # if (any(i in x for x in edges_new)):
        for j in i:
            if (any(genes_cluster[j] in x for x in edge_write)):
                gg = genes_cluster[j]
                eqv_edges_cluster.append(','.join([gg,j,'5']))
                node_ev.append(gg)
                
                
    e_t = list(edge_write) + eqv_edges_cluster
    edge_write = set(e_t)
    nodes_in_path = set(list(nodes_in_path) + node_ev)
    edge_write_tuple = [(x.split(",")[0], x.split(",")[1], x.split(",")[2]) for x in edge_write]
    edge_wo_tr = list(edge_write)                
                



    
    G = nx.DiGraph()
    edgesR23_3 = [(f, t, r) for f, t, r in edge_write_tuple if r == '2' or r == '3']
    edgesR23_2 = [(f + r, t + r) for f, t, r in edge_write_tuple if r == '2' or r == '3']
    G.add_edges_from(edgesR23_2)


#    G2 = nx.DiGraph(edgesR23_2)
#    try:
#        nx.find_cycle(G2, orientation='original')
#    except:
#        pass
#    dag = list(nx.find_cycle(G2, orientation='ignore'))

    
#    with open('cyclefile.txt', 'w') as filehandle:
#        for listitem in dag:
#            d = str(listitem)
#            filehandle.write('%s\n' % d)


# G.add_nodes_from(nodes)
    TR_edges = transitive_reduction(G).edges()
    TR_edges_3 = [(f[:-1], t[:-1], t[-1]) for f, t in TR_edges]
    removed_edg2 = [sub for sub in edgesR23_3 if sub not in TR_edges_3]
    #
    print('step3 done')
    remaining_edges = []
    for item in edge_write_tuple:
        if item not in removed_edg2:
            remaining_edges.append(item)
    print('step4 done')
    return remaining_edges
    
                
def build_net(Graph,R_edges,genes,genes_old):

    for f, t, r in R_edges:  # f: from, t:to,r: relation or Boolean relationship

        if r=='2':  # low-low
            Graph.add_edge(f, t, rel='2', color='navy',arrows={'to':{'enabled':True, 'type':'arrow'}})  
            Graph.add_node(f, label=genes[f].pretty(), size=5*(f.count("#")), color=genes[f].color)
            Graph.add_node(t, label=genes[t].pretty(), size=5*(t.count("#")),color=genes[t].color)

        elif r=='5': # equivalent
            Graph.add_edge(f, t, rel='5', color='cyan',arrows={'to':{'enabled':False, 'type':'-'}})  
            Graph.add_node(f, label=genes[f].pretty(), size=5*(f.count("#")),title='eqv', color="darksalmon")
            Graph.add_node(t, label=t, size=5*(t.count("#")), color=genes_old[t].color)

        elif r=='3': #high-high
            Graph.add_edge(f, t, rel='3', color='green',arrows={'to':{'enabled':True, 'type':'arrow'}})  
            Graph.add_node(f, label=genes[f].pretty(), size=5*(f.count("#")), color=genes[f].color)
            Graph.add_node(t, label=genes[t].pretty(), size=5*(t.count("#")), color=genes[t].color)      

        elif r=='4': # high-low
            Graph.add_edge(f, t, rel='4', color='red')  
            Graph.add_node(f, label=genes[f].pretty(), size=5*(f.count("#")),color=genes[f].color)
            Graph.add_node(t, label=genes[t].pretty(), size=5*(t.count("#")),color=genes[t].color)    
        else: # opposite
            Graph.add_edge(f, t, rel='6', color='black')  
            Graph.add_node(f, label=genes[f].pretty(), size=5*(f.count("#")),color=genes[f].color)
            Graph.add_node(t, label=genes[t].pretty(), size=5*(t.count("#")),color=genes[t].color) 
    return Graph
                








def add_single_hilo_op(add_hilo_op, genes):

    global g_right
    

    
    ending_points=set([])
    edge_write_p=set([])
    
    
    


    g_right_plus_hilo= g_right.copy()
    for g in add_hilo_op:

        if len(genes[g].edges[('4', 'shared')]) != 0: # add hilo
            for hl in genes[g].edges[('4', 'shared')]:


                info = ','.join([hl.name, g, '4'])

                edge_write_p.add(info)

                ending_points.add( hl.name)
                f=hl.name
                t=g
                g_right_plus_hilo.add_edge(f, t, rel='4', color='red')  
                g_right_plus_hilo.add_node(f, label=genes[f].pretty(), size=5*(f.count("#")),color=genes[f].color)
                g_right_plus_hilo.add_node(t, label=genes[t].pretty(), size=5*(t.count("#")),color=genes[t].color)                       
                    
                    

        if len(genes[g].edges[('6', 'shared')]) != 0: # add opo
            for hl in genes[g].edges[('6', 'shared')]:


                info = ','.join([hl.name, g, '6'])

                edge_write_p.add(info)

                ending_points.add( hl.name)
                f=hl.name
                t=g
                g_right_plus_hilo.add_edge(f, t, rel='6', color='black')  
                g_right_plus_hilo.add_node(f, label=genes[f].pretty(), size=5*(f.count("#")),color=genes[f].color)
                g_right_plus_hilo.add_node(t, label=genes[t].pretty(), size=5*(t.count("#")),color=genes[t].color)                       
                                
        
    return ending_points,edge_write_p,g_right_plus_hilo
            
            
            

def hihi_path(p_hh_func, g_hl_shrd_p_ll_per, g_ll_per, t_m):

    global edge_write_F
    global nodes_in_path_F
    global success
    global saved_func
    global saved_shared
    global edge_write_p
    
    saved_shared=set(saved_shared)


    if t_m == 'hilo_shared':
        
        last_shrd = g_hl_shrd_p_ll_per[-1]

        

        last_shrd.dfs('2', 'func', 1, g_hl_shrd_p_ll_per, [], lambda p, g: hihi_path([], p, g, t_m='hilo'), 'lolo')

        for func_ll_shrd in last_shrd.edges[('2', 'shared')]:
            if func_ll_shrd not in g_hl_shrd_p_ll_per:
                hihi_path(p_hh_func, g_hl_shrd_p_ll_per + [func_ll_shrd], [], t_m)
    if t_m == 'hilo':
        if len(g_ll_per.edges[('2', 'func')]) == 0 :


            jj = [gl for gl in g_hl_shrd_p_ll_per + [g_ll_per]]  
            ind=find_first_duplicate(jj)
            if ind != -1:
                jj=jj[0:ind]
            if len(jj)<2:
                return False

            
         
            

            list_path_edges = p_hh_func + g_hl_shrd_p_ll_per + [g_ll_per]
            list_path_edges = p_hh_func + jj
            if g_hl_shrd_p_ll_per[0].ct=='func':
                
                saved_func.append(g_hl_shrd_p_ll_per[0].name)
            else:
                saved_shared.add(g_hl_shrd_p_ll_per[0].name)
        
            for s in edge_write_p:
                if g_hl_shrd_p_ll_per[-1].name in s:
                    edge_write_F.add(s)
                    success.add(g_hl_shrd_p_ll_per[-1].name)
                    
                    
                    
            
            for i in range(len(list_path_edges) - 1):

                r = "3" 
                info = ','.join([ list_path_edges[i + 1].name,list_path_edges[i].name, r])

                if info not in edge_write_F:
                    edge_write_F.add(info)
                if list_path_edges[i].name not in nodes_in_path_F:
                    nodes_in_path_F.append(list_path_edges[i].name)
                if list_path_edges[i + 1].name not in nodes_in_path_F:
                    nodes_in_path_F.append(list_path_edges[i + 1].name)
                

    return True

            
            
                 
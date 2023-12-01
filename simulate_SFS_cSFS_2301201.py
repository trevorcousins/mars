import numpy as np
import msprime
import tskit
import pickle
import matplotlib.pyplot as plt
import pdb
import argparse
from datetime import datetime
import os
import time
from configure_model import *


# example usage 
# python /home/trc468/mars/simulate_SFS_cSFS_2301201.py -model C -num_MH 12 -T_super_archaic 1e+06 -T_modern_archaic 8e+05 -T_den_nea 4e+05 -T_pulse_ghost_to_MH 1e+05 -T_pulse_ghost_to_NEA 3e+05 -T_pulse_NEA_to_DEN 8e+04 -p_pulse_ghost_to_MH 0.3 -p_pulse_ghost_to_NEA 0.01 -p_pulse_NEA_to_DEN 0.03 -mu 1.25e-08 -N_super_ancestral 10000 -N_ancestral 10000 -N_ghost 10000 -N_modern 20000 -N_archaic 10000 -N_Neanderthal 1000 -N_Denisovan 1000 -out_SFS /tmp/deleteme -replicates 1e+04 
# python /home/trc468/mars/simulate_SFS_cSFS_2301201.py -model A -T_supersuper_archaic 1.2e+06 -T_super_archaic 1e+06 -T_modern_archaic 8e+05 -T_den_nea 4e+05 -T_pulse_ghost_to_MH 1e+05 -T_pulse_MH_to_NEA 3e+05 -T_pulse_superghost_to_DEN 1e+05 -T_pulse_NEA_to_DEN 8e+04 -p_pulse_ghost_to_MH 0.1 -p_pulse_ghost_to_NEA 0.4 -p_pulse_superghost_to_DEN 0.05 -p_pulse_MH_to_NEA 0.06 -p_pulse_NEA_to_DEN 0.03 -mu 1.25e-08 -N_supersuper_ancestral 10000 -N_superghost 10000 -N_super_ancestral 10000 -N_ancestral 10000 -N_ghost 10000 -N_modern 20000 -N_archaic 10000 -N_Neanderthal 1000 -N_Denisovan 1000 -out_SFS /tmp/deletemetest -num_MH 50 -mu 1.25e-08 -gen 29 -replicates 1e+04

def get_all_modern_coaltimes(tree,MH_indices):
    coaltimes = []
    for i in MH_indices:
        for j in MH_indices:
            if i!=j:
                coaltimes.append(tree.time(tree.mrca(i,j)))
    return coaltimes

def get_all_modern_archaic_coaltimes(tree,MH_indices,archaic):
    coaltimes = []
    for i in MH_indices:
        coaltimes.append(tree.time(tree.mrca(i,archaic)))
    return coaltimes

def index_columns_val_to_key(index_columns,value):
    return [key for key in index_columns.keys() if index_columns[key]==value ][0]

def clean_genotype_matrix(gen_mat):
    newgen_mat = np.copy(gen_mat)
    newgen_mat[newgen_mat>1]=0
    return newgen_mat

def new_clean_genotype_matrix(msim):
    gtmat=np.zeros(num_derived_alleles+2,dtype=int)
    for var in msim.variants():
        zedge = var.site.mutations[0].edge # only consider first mutation
        znodes = [node for node in tree.nodes()]
        zmutnode  = [i for i in znodes if tree.edge(i)==zedge][0]
        descs = [i for i in tree.get_leaves(zmutnode)]
        # znumMHderived = sum([tree.is_descendant(i,zmutnode) for i in MH_indices])
        gtmat[descs] = 1
    # return znumMHderived
    return gtmat

parser = argparse.ArgumentParser(description="Set parameters for simulation")
parser.add_argument('-num_MH','--num_MH',help='Number of modern human chromosomes to simulate',required=True,type=int,default=12)
parser.add_argument('-gen','--generation_time',help='Years per generation',required=False,type=int,default=29)
parser.add_argument('-NEA_age','--NEA_age',help='age of sampling for Neanderthal',required=False,type=int,default=3000)
parser.add_argument('-DEN_age','--DEN_age',help='age of sampling for Neanderthal',required=False,type=int,default=2500)
parser.add_argument('-model','--model',help='Simulate from model A or model C',required=True,type=str)
parser.add_argument('-T_super_archaic','--T_super_archaic',help='split time in years of root of (ancestral_humans,ghost)',required=False,type=float)
parser.add_argument('-T_modern_archaic','--T_modern_archaic',help='split time in years of main human lineage and lineage leading to (DEN,NEA)',required=False,type=float)
parser.add_argument('-T_den_nea','--T_den_nea',help='split time in years between NEA and DEN',required=False,type=float)
parser.add_argument('-T_pulse_ghost_to_MH','--T_pulse_ghost_to_MH',help='time in years at which super archaic ghost introgresses into modern humans',required=False,type=float)
parser.add_argument('-T_pulse_ghost_to_NEA','--T_pulse_ghost_to_NEA',help='time in years at which super archaic ghost introgresses into Neanderthals',required=False,type=float)
parser.add_argument('-T_pulse_NEA_to_DEN','--T_pulse_NEA_to_DEN',help='time in years at which which Neanderthals introgress into Denisovans',required=False,type=float)
parser.add_argument('-p_pulse_ghost_to_MH','--p_pulse_ghost_to_MH',help='pulse fraction for archaic ghost introgression into modern humans',required=False,type=float)
parser.add_argument('-p_pulse_ghost_to_NEA','--p_pulse_ghost_to_NEA',help='pulse fraction for archaic ghost introgression into Neanderthals',required=False,type=float)
parser.add_argument('-p_pulse_NEA_to_DEN','--p_pulse_NEA_to_DEN',help='pulse fraction for introgression into Denisovans',required=False,type=float)
parser.add_argument('-N_super_ancestral','--N_super_ancestral',help='Population size of super ancestral branch',required=False,type=int)
parser.add_argument('-N_ancestral','--N_ancestral',help='Population size of ancestral branch',required=False,type=int)
parser.add_argument('-N_ghost','--N_ghost',help='Population size of ghost lineage',required=False,type=int)
parser.add_argument('-N_modern','--N_modern',help='Population size of modern human lineage',required=False,type=int)
parser.add_argument('-N_archaic','--N_archaic',help='Population size of archaic branch',required=False,type=int)
parser.add_argument('-N_Neanderthal','--N_Neanderthal',help='Population size of Neanderthal lineage',required=False,type=int)
parser.add_argument('-N_Denisovan','--N_Denisovan',help='Population size of Denisovan lineage',required=False,type=int)
parser.add_argument('-T_supersuper_archaic','--T_supersuper_archaic',help='(model A) split time in years of root of (T_super_archaic,super_ghost)',required=False,type=float)
parser.add_argument('-T_pulse_superghost_to_DEN','--T_pulse_superghost_to_DEN',help='(model A) time in years at which superghost lineage introgresses into Denisovans',required=False,type=float)
parser.add_argument('-T_pulse_MH_to_NEA','--T_pulse_MH_to_NEA',help='(model A) time in years at which super archaic ghost introgresses into Neanderthals',required=False,type=float)
parser.add_argument('-p_pulse_MH_to_NEA','--p_pulse_MH_to_NEA',help='(model A) pulse fraction for MH introgression into Neanderthals',required=False,type=float)
parser.add_argument('-p_pulse_superghost_to_DEN','--p_pulse_superghost_to_DEN',help='(model A) pulse fraction for introgression from superghost into Denisovans',required=False,type=float)
parser.add_argument('-N_supersuper_ancestral','--N_supersuper_ancestral',help='(model A) Population size of super super ancestral branch',required=False,type=int)
parser.add_argument('-N_superghost','--N_superghost',help='(model A) Population size of superghost lineage (the one that introgresses into Neanderthal)',required=False,type=int)
parser.add_argument('-mu','--mu',help='mutation rate per generation per base pair',required=False,type=float)
# parser.add_argument('-r','--r',help='recombination rate per generation per base pair',required=True,type=float)
parser.add_argument('-out_SFS','--out_SFS',help='output path for SFS',required=True,type=str)
parser.add_argument('-replicates','--replicates',help='Number of replicates for drawing trees',required=True,type=float)
parser.add_argument('-mutationmodel','--mutationmodel',help='Model to use for mutation',required=False,type=str,default="jc69")
parser.add_argument('-record_muts','--record_muts',help='Binary for recording mutations or not',action='store_true')




args = parser.parse_args()
zargs = dir(args)
zargs = [zarg for zarg in zargs if zarg[0]!='_']
for zarg in zargs:
    print(f'{zarg} is ',end='',flush=True)
    exec(f'{zarg}=args.{zarg}')
    exec(f'print(args.{zarg},flush=True)')

if num_MH<5:
    print(f'ERROR; num_MH must be greater than 4. Aborting')
    sys.exit()

if os.path.isdir(os.path.dirname(out_SFS)) is False:
    os.makedirs(os.path.dirname(out_SFS))

zdemography = configure_demography(
    model,
    generation_time,
    T_super_archaic,
    T_modern_archaic, 
    T_den_nea,
    T_pulse_ghost_to_MH,
    T_pulse_ghost_to_NEA, 
    T_pulse_NEA_to_DEN,
    N_super_ancestral,
    N_ancestral,
    N_ghost,
    N_modern,
    N_archaic,
    N_Neanderthal,
    N_Denisovan,
    T_supersuper_archaic,
    p_pulse_ghost_to_MH,
    p_pulse_ghost_to_NEA,
    p_pulse_NEA_to_DEN,
    N_supersuper_ancestral,
    N_superghost,
    T_pulse_superghost_to_DEN,
    T_pulse_MH_to_NEA,
    p_pulse_superghost_to_DEN,
    p_pulse_MH_to_NEA)

MH_indices = [j for j in range(0,num_MH)]
NEA_index = num_MH # Neanderthal chromosome is the (num_MH+1)th chromosome, which in python indexing is num_MH
DEN_index = num_MH+1 # Neanderthal chromosome is the (num_MH+2)th chromosome, which in python indexing is num_MH+1
indices_without_DEN = MH_indices + [NEA_index]
indices_without_NEA = MH_indices + [DEN_index]

num_derived_alleles = num_MH+1 # include monomorphic and fixed derived
SFS_YRI_ndxx_mut = np.zeros(num_derived_alleles) # SFS for YRI with no conditions on site of Neanderthal or Denisovan. Include monomorphic and fixed derived
SFS_YRI_ndxx_mut_test = np.zeros(num_derived_alleles) # SFS for YRI with no conditions on site of Neanderthal or Denisovan. Include monomorphic and fixed derived
SFS_YRI_nd1x_mut = np.zeros(num_derived_alleles) # SFS for YRI with condition that Nea is derived and no conditions on Den. Include monomorphic and fixed derived
SFS_YRI_ndx1_mut = np.zeros(num_derived_alleles) # SFS for YRI with condition that Den is derived and no conditions on Nea. Include monomorphic and fixed derived
SFS_YRI_nd10_mut = np.zeros(num_derived_alleles) # SFS for YRI with condition that Nea is derived and Den is ancestral. Include monomorphic and fixed derived
SFS_YRI_nd01_mut = np.zeros(num_derived_alleles) # SFS for YRI with condition that Den is derived and Nea is ancestral. Include monomorphic and fixed derived
SFS_YRI_nd11_mut = np.zeros(num_derived_alleles) # SFS for YRI with condition that Den is derived and Nea is derived. Include monomorphic and fixed derived
SFS_YRI_nd00_mut = np.zeros(num_derived_alleles) # SFS for YRI with condition that Den is ancestral and Nea is ancestral. Include monomorphic and fixed derived

SFS_YRI_ndxx_branch = np.zeros(num_derived_alleles) # branch SFS for YRI with no conditions on site of Neanderthal or Denisovan. Include monomorphic and fixed derived
SFS_YRI_nd1x_branch = np.zeros(num_derived_alleles) # branch SFS for YRI with condition that Nea is derived and no conditions on Den. Include monomorphic and fixed derived
SFS_YRI_ndx1_branch = np.zeros(num_derived_alleles) # branch SFS for YRI with condition that Den is derived and no conditions on Nea. Include monomorphic and fixed derived
SFS_YRI_nd10_branch = np.zeros(num_derived_alleles) # branch SFS for YRI with condition that Nea is derived and Den is ancestral. Include monomorphic and fixed derived
SFS_YRI_nd01_branch = np.zeros(num_derived_alleles) # branch SFS for YRI with condition that Den is derived and Nea is ancestral. Include monomorphic and fixed derived
SFS_YRI_nd11_branch = np.zeros(num_derived_alleles) # branch SFS for YRI with condition that Nea is derived and Den is derived. Include monomorphic and fixed derived
SFS_YRI_nd00_branch = np.zeros(num_derived_alleles) # branch SFS for YRI with condition that Den is ancestral and Nea is ancestral. Include monomorphic and fixed derived


start = time.time()
for i in range(0,int(replicates)):
    sim = msprime.sim_ancestry( # simulate ancestry
        samples=[msprime.SampleSet(num_samples=num_MH, ploidy=1,population='modern',time=0), \
                msprime.SampleSet(num_samples=1, ploidy=1,population='neanderthal',time=NEA_age), \
                msprime.SampleSet(num_samples=1, ploidy=1,population='denisovan',time=DEN_age)],\
                demography=zdemography,sequence_length=1)
    if record_muts:
        msim = msprime.sim_mutations(sim, rate=mu,model=mutationmodel) # add mutations  
    tree = sim.at(0)
    # record branchs of YRI - don't want monomorphic
    # SFS_YRI_ndxx_branch[1:]+=msim.allele_frequency_spectrum(sample_sets=[MH_indices],mode="branch",polarised=True)[1:]
    SFS_YRI_ndxx_branch[1:]+=sim.allele_frequency_spectrum(sample_sets=[MH_indices],mode="branch",polarised=True)[1:]



    NEA_DEN_mrca = tree.get_mrca(NEA_index,DEN_index) # ancestral node to NEA and DEN

    # record branchs of YRI conditional on Nea being derived (nd1x)
    parents_NEA = [] # get parents on NEA lineage
    j = NEA_index
    while j!=-1:    
        parents_NEA.append(tree.parent(j))
        j = tree.parent(j)
    parents_NEA = np.array(parents_NEA[0:-1])
    parents_NEA_times = np.array([tree.time(j) for j in parents_NEA]) # get coal time of NEA parents
    zz = np.array([parents_NEA_times[0]] + [parents_NEA_times[i]- parents_NEA_times[i-1] for i in range(1,len(parents_NEA_times))]) # difference in tree heights for parents of NEA including singleton NEA branch
    zz[0] = zz[0] - NEA_age # adjust for sampling age
    z = [sum([tree.is_descendant(i,j) for i in MH_indices]) for j in parents_NEA] # number of YRI descendants, for the nodes that are the parents of NEA
    zprob_k_YRI_der_given_der_NEA = np.zeros(num_derived_alleles) # array for "given that NEA is derived, what is the probability that there are k derived alleles in YRI?", include 0 der in YRI
    zprob_k_YRI_der_given_der_NEA[0] = zz[0]
    for k,j in enumerate(z[0:-1]):
        zprob_k_YRI_der_given_der_NEA[j] += zz[k+1]
    SFS_YRI_nd1x_branch+=zprob_k_YRI_der_given_der_NEA

    # record branchs of YRI conditional on Nea being derived AND Den ancestral (nd10)
    if NEA_DEN_mrca==parents_NEA[-1]: # if DEN is outgroup wrt to the rest of the tree (i.e. NEA is derived and DEN is ancestral)
        SFS_YRI_nd10_branch+=zprob_k_YRI_der_given_der_NEA
    if NEA_DEN_mrca!=parents_NEA[-1] and NEA_DEN_mrca!=parents_NEA[0]: # if DEN join the (NEA,YRI) at some non trivial location (i.e. not the first or last coalescence event wrt Nea)
        # zq=zprob_k_YRI_der_given_der_NEA[0:z[np.where(parents_NEA==NEA_DEN_mrca)[0][0]]]
        # SFS_YRI_nd10_branch[0:zq.shape[0]]+=zq
        zq = z[0:np.where(parents_NEA==NEA_DEN_mrca)[0][0]]
        zprob_k_YRI_der_given_der_NEA_anc_DEN = np.zeros(num_derived_alleles) # array for "given that NEA is derived and DEN is ancestral, what is the probability that there are k derived alleles in YRI?"
        zprob_k_YRI_der_given_der_NEA_anc_DEN[0] = zz[0]
        for k,j in enumerate(zq):
            zprob_k_YRI_der_given_der_NEA_anc_DEN[j] += zz[k+1]
        SFS_YRI_nd10_branch+=zprob_k_YRI_der_given_der_NEA_anc_DEN

    # record banchs of YRI conditional NEA and DEN being ancestral (nd00), i.e. mutation arises more recently than mi(Nea_joins_tree,DEN_joins_tree)
    coal_nodes = np.array([i for i in tree.nodes() if i>num_derived_alleles+1]) # all the coalescent nodes in the tree
    coal_nodes_times = np.array([tree.time(i) for i in coal_nodes]) # times for all the coalescent nodes in the tree
    coal_nodes_woac = [i for i in coal_nodes if not tree.is_descendant(NEA_index,i) and not tree.is_descendant(DEN_index,i)] # coalescent nodes without archaic descendants
    coal_nodes_woac_nhd = [sum([tree.is_descendant(i,j) for i in MH_indices]) for j in coal_nodes_woac] # number of human descendants for each coalescent node without archaic children
    singletons = 0
    for qq in MH_indices:
        singletons += tree.time(tree.parent(qq)) # count each singleton 
    post_NEA_DEN = np.zeros(num_derived_alleles)
    post_NEA_DEN[1]+=singletons
    for qq,j in enumerate(coal_nodes_woac): #  for every coalescent node without archaic children
        num_leaves = coal_nodes_woac_nhd[qq] # how many MH children does it have
        length_branch = tree.time(tree.parent(j)) - tree.time(j) # how long is this subtending branch 
        post_NEA_DEN[num_leaves]+=length_branch 
    SFS_YRI_nd00_branch+=post_NEA_DEN

    # record banchs of YRI conditional NEA and DEN being dervied (nd11), i.e. mutation arises more anciently than coal_time(Nea,Den)
    zprob_k_YRI_der_given_der_NEADEN = np.zeros(num_derived_alleles) # array for "given that NEA and DEN is derived, what is the probability that there are k derived alleles in YRI?", include 0 der in YRI
    coal_nodes_wac = np.sort(np.array([i for i in coal_nodes if tree.is_descendant(NEA_index,i) and tree.is_descendant(DEN_index,i)])) # coalescent nodes with both archaic as descendants
    coal_nodes_wac_times = np.array([tree.time(i) for i in coal_nodes_wac]) # times of these nodes
    zaq = np.array([coal_nodes_wac_times[i+1]-coal_nodes_wac_times[i] for i in range(0,len(coal_nodes_wac_times)-1)]) # differences in times
    coal_nodes_wac_nhd = [sum([tree.is_descendant(i,j) for i in MH_indices]) for j in coal_nodes_wac] # number of human descendants for each coalescent node with both archaic children
    for k,j in enumerate(coal_nodes_wac_nhd[0:-1]):
        zprob_k_YRI_der_given_der_NEADEN[j] += zaq[k]
    SFS_YRI_nd11_branch+=zprob_k_YRI_der_given_der_NEADEN    

    # record branches of YRI conditional on Denisova (ndx1)
    parents_DEN = [] # get parents on DEN lineage
    j = DEN_index
    while j!=-1:    
        parents_DEN.append(tree.parent(j))
        j = tree.parent(j)
    parents_DEN = np.array(parents_DEN[0:-1])
    parents_DEN_times = np.array([tree.time(j) for j in parents_DEN]) # get coal time of DEN parents
    zz = np.array([parents_DEN_times[0]] + [parents_DEN_times[i]- parents_DEN_times[i-1] for i in range(1,len(parents_DEN_times))]) # difference in tree heights for parents of DEN including singleton DEN branch
    zz[0] = zz[0] - DEN_age # adjust for sampling age
    z = [sum([tree.is_descendant(i,j) for i in MH_indices]) for j in parents_DEN] # number of YRI descendants, for the nodes that are the parents of DEN
    zprob_k_YRI_der_given_der_DEN = np.zeros(num_derived_alleles) # array for "given that DEN is derived, what is the probability that there are k derived alleles in YRI?", include 0 der in YRI
    zprob_k_YRI_der_given_der_DEN[0] = zz[0]
    for k,j in enumerate(z[0:-1]):
        zprob_k_YRI_der_given_der_DEN[j] += zz[k+1]
    SFS_YRI_ndx1_branch+=zprob_k_YRI_der_given_der_DEN

    # record branchs of YRI conditional on Den being derived AND Nea ancestral (nd01)
    if NEA_DEN_mrca==parents_DEN[-1]: # if NEA is outgroup wrt to the rest of the tree (i.e. DEN is derived and NEA is ancestral)
        SFS_YRI_nd01_branch+=zprob_k_YRI_der_given_der_DEN
    if NEA_DEN_mrca!=parents_DEN[-1] and NEA_DEN_mrca!=parents_DEN[0]: # if NEA join the (DEN,YRI) at some non trivial location (i.e. not the first or last coalescence event wrt Den)
        # zq=zprob_k_YRI_der_given_der_DEN[0:z[np.where(parents_DEN==NEA_DEN_mrca)[0][0]]]
        # SFS_YRI_nd01_branch[0:zq.shape[0]]+=zq
        zq = z[0:np.where(parents_DEN==NEA_DEN_mrca)[0][0]]
        zprob_k_YRI_der_given_der_DEN_anc_NEA = np.zeros(num_derived_alleles) # array for "given that DEN is derived and NEA is ancestral, what is the probability that there are k derived alleles in YRI?"
        zprob_k_YRI_der_given_der_DEN_anc_NEA[0] = zz[0]
        for k,j in enumerate(zq):
            zprob_k_YRI_der_given_der_DEN_anc_NEA[j] += zz[k+1]
        SFS_YRI_nd01_branch+=zprob_k_YRI_der_given_der_DEN_anc_NEA
        


    # record muts
    if record_muts:
        if len(msim.genotype_matrix())>0:
            if max(msim.genotype_matrix()[0])>1: # if there is more than one type of mutation, skip
                continue
            # record YRI muts
            SFS_YRI_ndxx_mut[1:]+=msim.allele_frequency_spectrum(sample_sets=[MH_indices],mode="site",polarised=True)[1:]
            numMHderived = np.sum(msim.genotype_matrix()[0][0:num_MH])  
            SFS_YRI_ndxx_mut_test[numMHderived]+=1
            
            # record YRI conditioned Neanderthal muts
            if msim.genotype_matrix()[0,NEA_index]==1: # if Nea derived
                SFS_YRI_nd1x_mut[numMHderived]+=1
                if msim.genotype_matrix()[0,DEN_index]==0: # if Nea derived and Den ancestral
                    SFS_YRI_nd10_mut[numMHderived]+=1
                elif msim.genotype_matrix()[0,DEN_index]==1: # if Nea derived and Den derived
                    SFS_YRI_nd11_mut[numMHderived]+=1

            # record YRI conditioned Denisovan muts
            if msim.genotype_matrix()[0,DEN_index]==1: # if Den derived
                SFS_YRI_ndx1_mut[numMHderived]+=1
                if msim.genotype_matrix()[0,NEA_index]==0: # if Den derived and Nea ancestral
                    SFS_YRI_nd01_mut[numMHderived]+=1
            
            if msim.genotype_matrix()[0,DEN_index]==0 and msim.genotype_matrix()[0,NEA_index]==0: # if Den ancestral and Nea ancestral
                SFS_YRI_nd00_mut[numMHderived]+=1

 
    if i%100==0:
        print(f'on replicate={i} out of {replicates}',flush=True)
        if record_muts:    
            np.savetxt(out_SFS+'YRIndxx_muts.txt',SFS_YRI_ndxx_mut)
            np.savetxt(out_SFS+'YRIndxx_muts_test.txt',SFS_YRI_ndxx_mut_test)
            np.savetxt(out_SFS+'YRInd1x_muts.txt',SFS_YRI_nd1x_mut)
            np.savetxt(out_SFS+'YRIndx1_muts.txt',SFS_YRI_ndx1_mut)
            np.savetxt(out_SFS+'YRInd10_muts.txt',SFS_YRI_nd10_mut)
            np.savetxt(out_SFS+'YRInd01_muts.txt',SFS_YRI_nd01_mut)
            np.savetxt(out_SFS+'YRInd11_muts.txt',SFS_YRI_nd11_mut)
            np.savetxt(out_SFS+'YRInd00_muts.txt',SFS_YRI_nd00_mut)
        np.savetxt(out_SFS+'YRIndxx_branch.txt',SFS_YRI_ndxx_branch)
        np.savetxt(out_SFS+'YRInd1x_branch.txt',SFS_YRI_nd1x_branch)
        np.savetxt(out_SFS+'YRIndx1_branch.txt',SFS_YRI_ndx1_branch)
        np.savetxt(out_SFS+'YRInd10_branch.txt',SFS_YRI_nd10_branch)
        np.savetxt(out_SFS+'YRInd01_branch.txt',SFS_YRI_nd01_branch)
        np.savetxt(out_SFS+'YRInd11_branch.txt',SFS_YRI_nd11_branch)
        np.savetxt(out_SFS+'YRInd00_branch.txt',SFS_YRI_nd00_branch)

end = time.time()
print(f'time taken for {replicates} replicates = {end-start} seconds')

if record_muts:    
    np.savetxt(out_SFS+'YRIndxx_muts.txt',SFS_YRI_ndxx_mut)
    np.savetxt(out_SFS+'YRIndxx_muts_test.txt',SFS_YRI_ndxx_mut_test)
    np.savetxt(out_SFS+'YRInd1x_muts.txt',SFS_YRI_nd1x_mut)
    np.savetxt(out_SFS+'YRIndx1_muts.txt',SFS_YRI_ndx1_mut)
    np.savetxt(out_SFS+'YRInd10_muts.txt',SFS_YRI_nd10_mut)
    np.savetxt(out_SFS+'YRInd01_muts.txt',SFS_YRI_nd01_mut)
    np.savetxt(out_SFS+'YRInd11_muts.txt',SFS_YRI_nd11_mut)
    np.savetxt(out_SFS+'YRInd00_muts.txt',SFS_YRI_nd00_mut)
np.savetxt(out_SFS+'YRIndxx_branch.txt',SFS_YRI_ndxx_branch)
np.savetxt(out_SFS+'YRInd1x_branch.txt',SFS_YRI_nd1x_branch)
np.savetxt(out_SFS+'YRIndx1_branch.txt',SFS_YRI_ndx1_branch)
np.savetxt(out_SFS+'YRInd10_branch.txt',SFS_YRI_nd10_branch)
np.savetxt(out_SFS+'YRInd01_branch.txt',SFS_YRI_nd01_branch)
np.savetxt(out_SFS+'YRInd11_branch.txt',SFS_YRI_nd11_branch)
np.savetxt(out_SFS+'YRInd00_branch.txt',SFS_YRI_nd00_branch)
qwes = ['xx','1x','x1','10','01','11','00']
if record_muts:
    qwers = ['muts','branch']
else:
    qwers = ['branch']
for qwer in qwers:
    for qwe in qwes:
        print(f'saved to {out_SFS}YRInd{qwe}_{qwer}.txt',flush=True)
        if qwe=='xx' and qwer=='muts':
            print(f'saved to {out_SFS}YRInd{qwe}_{qwer}_test.txt',flush=True)

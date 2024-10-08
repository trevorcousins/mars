"""
Copied from /home/trc468/mars/simulate_archaic_model_231116.py

Simulate from model A or C, record sfs, genotypes, coal data
Usage 
    python /home/trc468/mars/simulate_archaic_model_231201.py -model C -T_super_archaic 1e+06 -T_modern_archaic 8e+05 -T_den_nea 4e+05 -T_pulse_ghost_to_MH 1e+05 -T_pulse_ghost_to_NEA 3e+05 -T_pulse_NEA_to_DEN 8e+04 -p_pulse_ghost_to_MH 0.1 -p_pulse_ghost_to_NEA 0.4 -p_pulse_NEA_to_DEN 0.03 -mu 1.25e-08 -N_super_ancestral 10000 -N_ancestral 10000 -N_ghost 10000 -N_modern 20000 -N_archaic 10000 -N_Neanderthal 1000 -N_Denisovan 1000 -out_prefix /tmp/deletemetest -num_MH 10 -mu 1.25e-08 -r 1e-08 -L 5e+06 -gen 29
    python /home/trc468/mars/simulate_archaic_model_231201.py -model A -T_supersuper_archaic 1.2e+06 -T_super_archaic 1e+06 -T_modern_archaic 8e+05 -T_den_nea 4e+05 -T_pulse_ghost_to_MH 1e+05 -T_pulse_MH_to_NEA 3e+05 -T_pulse_superghost_to_DEN 1e+05 -T_pulse_NEA_to_DEN 8e+04 -p_pulse_ghost_to_MH 0.1 -p_pulse_ghost_to_NEA 0.4 -p_pulse_superghost_to_DEN 0.05 -p_pulse_MH_to_NEA 0.06 -p_pulse_NEA_to_DEN 0.03 -mu 1.25e-08 -N_supersuper_ancestral 10000 -N_superghost 10000 -N_super_ancestral 10000 -N_ancestral 10000 -N_ghost 10000 -N_modern 20000 -N_archaic 10000 -N_Neanderthal 1000 -N_Denisovan 1000 -out_prefix /tmp/deletemetest -num_MH 10 -mu 1.25e-08 -r 1e-08 -L 1e+05 -gen 29

    python /home/trc468/mars/simulate_archaic_model_231201.py -model C -T_super_archaic 1e+06 -T_modern_archaic 8e+05 -T_den_nea 4e+05 -T_pulse_ghost_to_MH 1e+05 -T_pulse_ghost_to_NEA 3e+05 -T_pulse_NEA_to_DEN 8e+04 -p_pulse_ghost_to_MH 0.1 -p_pulse_ghost_to_NEA 0.4 -p_pulse_NEA_to_DEN 0.03 -mu 1.25e-08 -N_super_ancestral 10000 -N_ancestral 10000 -N_ghost 10000 -N_modern 20000 -N_archaic 10000 -N_Neanderthal 1000 -N_Denisovan 1000 -out_prefix /home/trc468/notebooks/archaics/simulations_testing_231118/231118_1132 -num_MH 50 -mu 1.25e-08 -r 1e-08 -L 1e+07 -gen 29
    python /home/trc468/mars/simulate_archaic_model_231201.py -model A -T_supersuper_archaic 1.2e+06 -T_super_archaic 1e+06 -T_modern_archaic 8e+05 -T_den_nea 4e+05 -T_pulse_ghost_to_MH 1e+05 -T_pulse_MH_to_NEA 3e+05 -T_pulse_superghost_to_DEN 1e+05 -T_pulse_NEA_to_DEN 8e+04 -p_pulse_ghost_to_MH 0.1 -p_pulse_ghost_to_NEA 0.4 -p_pulse_superghost_to_DEN 0.05 -p_pulse_MH_to_NEA 0.06 -p_pulse_NEA_to_DEN 0.03 -mu 1.25e-08 -N_supersuper_ancestral 10000 -N_superghost 10000 -N_super_ancestral 10000 -N_ancestral 10000 -N_ghost 10000 -N_modern 20000 -N_archaic 10000 -N_Neanderthal 1000 -N_Denisovan 1000 -out_prefix /tmp/deletemetest -num_MH 50 -mu 1.25e-08 -r 1e-08 -L 1e+07 -gen 29

    python /home/trc468/mars/simulate_archaic_model_231201.py -model C -T_super_archaic 1e+06 -T_modern_archaic 8e+05 -T_den_nea 4e+05 -T_pulse_ghost_to_MH 1e+05 -T_pulse_ghost_to_NEA 3e+05 -T_pulse_NEA_to_DEN 8e+04 -p_pulse_ghost_to_MH 0.1 -p_pulse_ghost_to_NEA 0.4 -p_pulse_NEA_to_DEN 0.03 -mu 1.25e-08 -N_super_ancestral 10000 -N_ancestral 10000 -N_ghost 10000 -N_modern 20000 -N_archaic 10000 -N_Neanderthal 1000 -N_Denisovan 1000 -out_prefix /home/trc468/notebooks/archaics/simulations_testing_231118/231118_1132 -num_MH 50 -mu 1.25e-08 -r 1e-08 -L 1e+06 -gen 29 -T_AMH_expand 20000 -N_MH_expand_rate 0.00001
    python /home/trc468/mars/simulate_archaic_model_231201.py -model C -T_super_archaic 7.5e+05 -T_modern_archaic 6e+05 -T_den_nea 3.5e+05 -T_pulse_ghost_to_MH 1e+05 -T_pulse_ghost_to_NEA 2.5e+05 -T_pulse_NEA_to_DEN 8e+04 -p_pulse_ghost_to_MH 0.9 -p_pulse_ghost_to_NEA 0.1 -p_pulse_NEA_to_DEN 0.05 -N_super_ancestral 20000 -N_ancestral 20000 -N_ghost 20000 -N_modern 20000 -N_AMH 30000 -N_archaic 15000 -N_Neanderthal 1000 -N_Denisovan 1000  -num_MH 42 -gen 29 -DEN_age 80000 -NEA_age 100000 -T_AMH_expand 1e+04 -N_AMH_present 1e+06 -mu 1.25e-08 -r 1e-08 -L 1e+07 -out_prefix /n/scratch3/users/t/trc468/mars_simulations_231204/231211_2200_modelC_numMH42_L1e07_7.5e+05_6e+05_3.5e+05_1e+05_2.5e+05_8e+04_0.9_0.1_0.05_20000_20000_20000_20000_15000_1000_1000_1e+04_1e+06

"""
from configure_model_231227 import *
import numpy as np
import msprime
import tskit
import pickle
import matplotlib.pyplot as plt
import pdb
import argparse
from datetime import datetime
import os
import sys



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



def get_sfs(gtmat,num_MH):
    ndxx = np.zeros(num_MH+1)
    nd00 = np.zeros(num_MH+1)
    nd01 = np.zeros(num_MH+1)
    nd10 = np.zeros(num_MH+1)
    nd11 = np.zeros(num_MH+1)

    # remove double mutations
    nodoublemutations = np.where(gtmat[:,1:].max(axis=1)==1)[0] 
    gtmat = gtmat[nodoublemutations,:]

    num_monomorphic_ancestral = L -len(nodoublemutations) - gtmat.shape[0]

    ndxx_counts = gtmat[:,1:-4].sum(axis=1)
    ndxx = np.histogram(ndxx_counts,bins=np.arange(0,num_MH+1))[0]
    ndxx[0]+=num_monomorphic_ancestral

    nd00_counts = gtmat[(gtmat[:,num_MH+1]==0) & (gtmat[:,num_MH+3]==0),1:-4].sum(axis=1) # condition on NEA and DEN hap being ancestral
    nd00 = np.histogram(nd00_counts,bins=np.arange(0,num_MH+1))[0]
    nd00[0]+=num_monomorphic_ancestral


    nd10_counts = gtmat[(gtmat[:,num_MH+1]==1) & (gtmat[:,num_MH+3]==0),1:-4].sum(axis=1) # condition on NEA being derived and DEN ancestral
    nd10 = np.histogram(nd10_counts,bins=np.arange(0,num_MH+1))[0]

    nd01_counts = gtmat[(gtmat[:,num_MH+1]==0) & (gtmat[:,num_MH+3]==1),1:-4].sum(axis=1) # condition on DEN being derived and NEA ancestral
    nd01 = np.histogram(nd01_counts,bins=np.arange(0,num_MH+1))[0]

    nd11_counts = gtmat[(gtmat[:,num_MH+1]==1) & (gtmat[:,num_MH+3]==1),1:-4].sum(axis=1) # condition on NEA and DEN being derived
    nd11 = np.histogram(nd11_counts,bins=np.arange(0,num_MH+1))[0]
    return ndxx,nd00,nd10,nd01,nd11

def combine_gtmat_coaldata(gtmat,coaldata):
    width_gtmat = gtmat.shape[1]
    width_coaldata = coaldata.shape[1]

    combined_gt_coal = np.zeros(shape=(gtmat.shape[0],4+width_coaldata))
    zstart = 0 
    zcount = 0
    for i in gtmat:
        zz = coal_data[(coal_data[zstart:,0]<=i[0]) & (coal_data[zstart:,1]>i[0]),:]
        combined_gt_coal[zcount,0] = i[0]
        combined_gt_coal[zcount,1] = i[MH_indices+1].sum()
        combined_gt_coal[zcount,2] = i[NEA_indices+1].sum()
        combined_gt_coal[zcount,3] = i[DEN_indices+1].sum()
        combined_gt_coal[zcount,4:] = zz
        zcount+=1
    return combined_gt_coal

parser = argparse.ArgumentParser(description="Set parameters for simulation")
parser.add_argument('-num_MH','--num_MH',help='Number of MH haplotypes to simulate',required=True,type=int)

parser.add_argument('-gen','--generation_time',help='Years per generation',required=False,type=int,default=29)
parser.add_argument('-NEA_age','--NEA_age',help='age of sampling for Neanderthal in years',required=False,type=int,default=100000)
parser.add_argument('-DEN_age','--DEN_age',help='age of sampling for Denisovans in years',required=False,type=int,default=80000)
parser.add_argument('-model','--model',help='Simulate from model A or model C',required=True,type=str)
parser.add_argument('-T_super_archaic','--T_super_archaic',help='split time in years of root of (ancestral_humans,ghost)',required=False,type=float)
parser.add_argument('-T_modern_archaic','--T_modern_archaic',help='split time in years of main human lineage and lineage leading to (DEN,NEA)',required=False,type=float)
parser.add_argument('-T_den_nea','--T_den_nea',help='split time in years between NEA and DEN',required=False,type=float)
parser.add_argument('-T_pulse_ghost_to_MH','--T_pulse_ghost_to_MH',help='time in years at which super archaic ghost introgresses into modern humans',required=False,type=float)
parser.add_argument('-T_pulse_ghost_to_NEA','--T_pulse_ghost_to_NEA',help='(model C) time in years at which super archaic ghost introgresses into Neanderthals',required=False,type=float)
parser.add_argument('-T_pulse_NEA_to_DEN','--T_pulse_NEA_to_DEN',help='time in years at which which Neanderthals introgress into Denisovans',required=False,type=float)
parser.add_argument('-p_pulse_ghost_to_MH','--p_pulse_ghost_to_MH',help='pulse fraction for archaic ghost introgression into modern humans',required=False,type=float)
parser.add_argument('-p_pulse_ghost_to_NEA','--p_pulse_ghost_to_NEA',help='(model C) pulse fraction for archaic ghost introgression into Neanderthals',required=False,type=float)
parser.add_argument('-p_pulse_NEA_to_DEN','--p_pulse_NEA_to_DEN',help='pulse fraction for introgression into Denisovans',required=False,type=float)
parser.add_argument('-N_super_ancestral','--N_super_ancestral',help='Population size of super ancestral branch',required=False,type=float)
parser.add_argument('-N_ancestral','--N_ancestral',help='Population size of ancestral branch',required=False,type=float)
parser.add_argument('-N_ghost','--N_ghost',help='Population size of ghost lineage',required=False,type=float)
parser.add_argument('-N_modern','--N_modern',help='Population size of modern human lineage',required=False,type=float)
parser.add_argument('-N_AMH','--N_AMH',help='Population size of modern human lineage after (forward in time) the ghost contributes to modern humans',required=False,type=float)
parser.add_argument('-N_archaic','--N_archaic',help='Population size of archaic branch',required=False,type=float)
parser.add_argument('-N_Neanderthal','--N_Neanderthal',help='Population size of Neanderthal lineage',required=False,type=float)
parser.add_argument('-N_Denisovan','--N_Denisovan',help='Population size of Denisovan lineage',required=False,type=float)

parser.add_argument('-T_supersuper_archaic','--T_supersuper_archaic',help='(model A) split time in years of root of (T_super_archaic,super_ghost)',required=False,type=float)
parser.add_argument('-T_pulse_superghost_to_DEN','--T_pulse_superghost_to_DEN',help='(model A) time in years at which superghost lineage introgresses into Denisovans',required=False,type=float)
parser.add_argument('-T_pulse_MH_to_NEA','--T_pulse_MH_to_NEA',help='(model A) time in years at which super archaic ghost introgresses into Neanderthals',required=False,type=float)
parser.add_argument('-p_pulse_MH_to_NEA','--p_pulse_MH_to_NEA',help='(model A) pulse fraction for MH introgression into Neanderthals',required=False,type=float)
parser.add_argument('-p_pulse_superghost_to_DEN','--p_pulse_superghost_to_DEN',help='(model A) pulse fraction for introgression from superghost into Denisovans',required=False,type=float)
parser.add_argument('-N_supersuper_ancestral','--N_supersuper_ancestral',help='(model A) Population size of super super ancestral branch',required=False,type=int)
parser.add_argument('-N_superghost','--N_superghost',help='(model A) Population size of superghost lineage (the one that introgresses into Neanderthal)',required=False,type=float)

parser.add_argument('-T_AMH_expand','--T_AMH_expand',help='Time at which MH begin to exponentially expand',required=False,type=float)
parser.add_argument('-N_AMH_present','--N_AMH_present',help='Population size of modern human lineage at present; at time T_AMH_expand start growing exponentially (forwards in time) to reach size N_AMH_present',required=False,type=float)
parser.add_argument('-N_ghost_recent','--N_ghost_recent',help='Population size of ghost lineage after (forwards in time) it introgresses into Neanderthal but before it introgresses into MH',required=False,type=float)
parser.add_argument('-T_ghost_BN_start','--T_ghost_BN_start',help='Start time in years (going forwards in time) of bottleneck in ghost lineage',required=False,type=float)
parser.add_argument('-T_ghost_BN_end','--T_ghost_BN_end',help='End time in years (going forwards in time) of bottleneck in ghost lineage',required=False,type=float)
parser.add_argument('-N_ghost_BN_intensity','--N_ghost_BN_intensity',help='Intensity of bottleneck in ghost lineage; population size contracts to N_ghost_BN_intensity*N_ghost at time T_ghost_BN_start and remains until T_ghost_BN_end; after which it returns to N_ghost',required=False,type=float)
parser.add_argument('-Cprime_introgression','--Cprime_introgression',help='If model is Cprime, which population does the superghost introgress in to? If "a", then into Denisovans; if "b" then into the ancestors of Neanderthals and Denisovans ("archaic"); if "c" then the ancestors of modern humans, Neanderthals and Denisovans ("ancestral"). T_pulse_superghost_to_DEN must reflect this.',required=False,type=str,default="a")

parser.add_argument('-N_Neanderthal_arch','--N_Neanderthal_arch',help='(model Cprime) The size of the Neanderthal lineage more anciently than the admixture event "T_Nea_Admix" . If not given defaults to "N_Neanderthal" ',required=False,type=float)
parser.add_argument('-T_NEA_admix','--T_NEA_admix',help='The time at which the ghost (ghost_nea) introgresses into the Neanderthal lineage',required=False,type=float)
parser.add_argument('-N_ghost_nea','--N_ghost_nea',help='(model Cprime) The size of the ghost lineage that contributes to Neanderthals',required=False,type=float)
parser.add_argument('-N_ghost_human','--N_ghost_human',help='(model Cprime) The size of the ghost lineage that contributes to huamns',required=False,type=float)

parser.add_argument('-L','--L',help='Length of haplotypes',required=True,type=float)
parser.add_argument('-mu','--mu',help='mutation rate per generation per base pair',required=True,type=float)
parser.add_argument('-r','--r',help='recombination rate per generation per base pair',required=True,type=float)
parser.add_argument('-out_prefix','--out_prefix',help='output path (prefix) for data ',required=True,type=str)

arguments = sys.argv[1:]
command_line = 'python ' + ' '.join(['"{}"'.format(arg) if ' ' in arg else arg for arg in [sys.argv[0]] + arguments])
print(f'Command line: {command_line}')
print()

args = parser.parse_args()
zargs = dir(args)
zargs = [zarg for zarg in zargs if zarg[0]!='_']
for zarg in zargs:
    print(f'{zarg} is ',end='')
    exec(f'{zarg}=args.{zarg}')
    exec(f'print(args.{zarg})')

if num_MH<4:
    print(f'ERROR: num_MH must be bigger than or equal to 4; currently it is {num_MH}. Aborting')
    sys.exit()
if model!='A' and model!='C' and model!='Cprime':
    print(f'ERROR: model must be A or C; currently it is {model}. Aborint')
    sys.exit()

if out_prefix==None:
    out_prefix = currentdir

if os.path.isdir(os.path.dirname(out_prefix)) is False:
    os.makedirs(os.path.dirname(out_prefix))

if T_AMH_expand=='None':
    T_AMH_expand = T_pulse_ghost_to_MH-10
    N_AMH_present=N_AMH
else:
    T_AMH_expand = int(T_AMH_expand)
    if N_AMH_present==None:
        print(f'T_AMH_expand is given but N_AMH_present is not given; this is contraditcory')
        N_AMH_present=N_AMH
N_MH_expand_rate = -np.log(N_AMH/N_AMH_present)/(T_AMH_expand/generation_time)
print(f'N_MH_expand_rate={N_MH_expand_rate}')


if T_ghost_BN_start==None:
    T_ghost_BN_start = T_super_archaic-2
    T_ghost_BN_end = T_super_archaic-3
    N_ghost_BN_intensity = N_ghost

if (T_super_archaic >= T_supersuper_archaic):
    print(f'ERROR: T_super_archaic={T_super_archaic} must be bigger than or equal to T_supersuper_archaic={T_supersuper_archaic}. Aborting.')
    sys.exit()


if T_ghost_BN_start!=None:
    if T_ghost_BN_start<=T_ghost_BN_end:
        print(f'T_ghost_BN_start must be greater than T_ghost_BN_end. Aborting')
        sys.exit()

if N_Neanderthal_arch==None:
    N_Neanderthal_arch = N_Neanderthal
if T_NEA_admix>T_pulse_ghost_to_NEA:
    print(f'ERROR: T_NEA_admix = {T_NEA_admix} must be smaller than T_pulse_ghost_to_NEA={T_pulse_ghost_to_NEA}')
if N_ghost_nea==None:
    N_ghost_nea = N_ghost
if N_ghost_human==None:
    N_ghost_human = N_ghost

if N_ghost_recent==None:
    N_ghost_recent = N_ghost

MH_indices = np.array([i for i in range(0,num_MH)]) # indices of MH haplotypes
NEA_indices = np.array([num_MH,num_MH+1]) # indices of NEA haplotypes
DEN_indices = np.array([num_MH+2,num_MH+3]) # indices of DEN haplotypes

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
    p_pulse_MH_to_NEA,
    T_AMH_expand,
    N_MH_expand_rate,
    N_AMH_present,
    N_ghost_recent,
    T_ghost_BN_start,
    T_ghost_BN_end,
    N_ghost_BN_intensity,
    N_AMH,
    Cprime_introgression,
    N_Neanderthal_arch,
    T_NEA_admix,
    N_ghost_nea,
    N_ghost_human
    )


NEA_age = NEA_age/generation_time
DEN_age = DEN_age/generation_time

sim = msprime.sim_ancestry( # simulate ancestry
    samples=[msprime.SampleSet(num_samples=num_MH, ploidy=1,population='modern',time=0), 
            msprime.SampleSet(num_samples=1, ploidy=2,population='neanderthal',time=NEA_age), \
            msprime.SampleSet(num_samples=1, ploidy=2,population='denisovan',time=DEN_age)],\
            demography=zdemography,sequence_length=L,recombination_rate=r)
            
msim = msprime.sim_mutations(sim, rate=mu) # add mutations
pos = np.array([int(i.position) for i in msim.variants()])
zgtmat = msim.genotype_matrix()
gtmat = np.zeros(shape=(zgtmat.shape[0],zgtmat.shape[1]+1),dtype=int)
gtmat[:,0] = pos
gtmat[:,1:] = zgtmat
zb = gtmat.shape[0]
ze = len(np.where(gtmat[:,1:].sum(axis=1)==0)[0])
zd = len(np.where(gtmat[:,1:].max(axis=1)==2)[0])
gtmat = gtmat[np.where(gtmat[:,1:].max(axis=1)==1)[0],:] # remove double mutations
za = gtmat.shape[0]
print(f'\tOut of {zb-ze} sites, removed {zd} sites where there are double mutations')


zstring = [f'MHhap{i}' for i in MH_indices]
gtmat_columns = ['Position'] + zstring + ['NEAhap0','NEAhap1','DENhap0','DENhap1']
gtmat_columns_str = " ".join(gtmat_columns)
ndxx, nd00,nd10,nd01,nd11 = get_sfs(gtmat,num_MH)
# want left right TMRCAmax TMRCAmin NEA1_NEA2 DEN1_DEN2 NEA1_DEN1 NEA1_DEN2 NEA2_DEN1 NEA2_DEN2 MH0_MH1 MH2_MH3 MH0_NEA0 MH0_NEA1 MH0_DEN0 MH0_DEN1 MH1_NEA0 MH1_NEA1 MH1_DEN0 MH1_DEN1
coaldata_columns = ['left','right','TMRCAmax','TMRCAmin','NEA1_NEA2','DEN1_DEN2','NEA1_DEN1','NEA1_DEN2','NEA2_DEN1','NEA2_DEN2','MH0_MH1','MH2_MH3','MH0_NEA0','MH0_NEA1','MH0_DEN0','MH0_DEN1','MH1_NEA0','MH1_NEA1','MH1_DEN0','MH1_DEN1']
coaldata_columns_str = " ".join(coaldata_columns)
gtmat_coal_data_columns = ['Position','allele_count_HM','allele_count_Neanderthal','allele_count_Denisovan'] + coaldata_columns
gtmat_coal_data_columns_str = " ".join(gtmat_coal_data_columns)
coal_data = np.zeros(shape=(sim.get_num_trees(),20))
for tree in sim.trees():
    coal_data[tree.index,0] = tree.interval[0] # left interval
    coal_data[tree.index,1] = tree.interval[1] # right interval
    allMH_TMRCA = np.array([tree.time(tree.mrca(i,j)) for i in range(0,len(MH_indices)) for j in range(i+1,len(MH_indices)) ])
    coal_data[tree.index,2] = np.max(allMH_TMRCA) # TMRCAmax
    coal_data[tree.index,3] = np.min(allMH_TMRCA) # TMRCAmin

    coal_data[tree.index,4] = tree.time(tree.mrca(NEA_indices[0],NEA_indices[1])) # within NEA coal time
    coal_data[tree.index,5] = tree.time(tree.mrca(DEN_indices[0],DEN_indices[1])) # within DEN coal time
    zz = 6
    # write cross NEA_DEN coal times
    for jj in [0,1]: # NEA hap
        for kk in [0,1]: # DEN hap
            coal_data[tree.index,zz] = tree.time(tree.mrca(NEA_indices[jj],DEN_indices[kk])) 
            zz+=1
    coal_data[tree.index,zz] = tree.time(tree.mrca(MH_indices[0],MH_indices[1])) # within MH0 MH1 coal time
    coal_data[tree.index,zz+1] = tree.time(tree.mrca(MH_indices[2],MH_indices[3])) # within MH2 MH3 coal time
    zz+=2
    for jj in [0,3]: # MH haps
        for arch in ["NEA","DEN"]: # arch type
            for kk in [0,1]: # arch haps
                exec(f'arch_indices={arch}_indices')    
                coal_data[tree.index,zz] = tree.time(tree.mrca(MH_indices[jj],arch_indices[kk]))
                zz+=1


for sfstype in ['ndxx','nd00','nd10','nd01','nd11']:
    exec(f'zarray={sfstype}') 
    fname = f'{out_prefix}_{sfstype}.txt'
    np.savetxt(fname,zarray)
    print(f'Saved {sfstype} to {fname}')
print()

fname = f'{out_prefix}_gtmat.txt.gz'
np.savetxt(fname,gtmat,header=gtmat_columns_str, fmt='%i')
print(f'Saved haplotype data to {fname}')

fname = f'{out_prefix}_coaldata.txt.gz'
np.savetxt(fname,coal_data,header=coaldata_columns_str)
print(f'Saved coalescent data to {fname}')

fname = f'{out_prefix}_gtmat_coaldata.txt.gz'
combined_gtmat_coaldata = combine_gtmat_coaldata(gtmat,coal_data)
np.savetxt(fname,combined_gtmat_coaldata,header=gtmat_coal_data_columns_str)
print(f'Saved combined allele count and coalescent data to {fname}')


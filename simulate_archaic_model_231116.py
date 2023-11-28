"""
Simulate from model A or C, record sfs, genotypes, coal data
Usage 
    python /home/trc468/mars/simulate_archaic_model_231116.py -model C -T_super_archaic 1e+06 -T_modern_archaic 8e+05 -T_den_nea 4e+05 -T_pulse_ghost_to_MH 1e+05 -T_pulse_ghost_to_NEA 3e+05 -T_pulse_NEA_to_DEN 8e+04 -p_pulse_ghost_to_MH 0.1 -p_pulse_ghost_to_NEA 0.4 -p_pulse_NEA_to_DEN 0.03 -mu 1.25e-08 -N_super_ancestral 10000 -N_ancestral 10000 -N_ghost 10000 -N_modern 20000 -N_archaic 10000 -N_Neanderthal 1000 -N_Denisovan 1000 -out_prefix /tmp/deletemetest -num_MH 10 -mu 1.25e-08 -r 1e-08 -L 5e+06 -gen 29
    python /home/trc468/mars/simulate_archaic_model_231116.py -model A -T_supersuper_archaic 1.2e+06 -T_super_archaic 1e+06 -T_modern_archaic 8e+05 -T_den_nea 4e+05 -T_pulse_ghost_to_MH 1e+05 -T_pulse_MH_to_NEA 3e+05 -T_pulse_superghost_to_DEN 1e+05 -T_pulse_NEA_to_DEN 8e+04 -p_pulse_ghost_to_MH 0.1 -p_pulse_ghost_to_NEA 0.4 -p_pulse_superghost_to_DEN 0.05 -p_pulse_MH_to_NEA 0.06 -p_pulse_NEA_to_DEN 0.03 -mu 1.25e-08 -N_supersuper_ancestral 10000 -N_superghost 10000 -N_super_ancestral 10000 -N_ancestral 10000 -N_ghost 10000 -N_modern 20000 -N_archaic 10000 -N_Neanderthal 1000 -N_Denisovan 1000 -out_prefix /tmp/deletemetest -num_MH 10 -mu 1.25e-08 -r 1e-08 -L 1e+05 -gen 29

    python /home/trc468/mars/simulate_archaic_model_231116.py -model C -T_super_archaic 1e+06 -T_modern_archaic 8e+05 -T_den_nea 4e+05 -T_pulse_ghost_to_MH 1e+05 -T_pulse_ghost_to_NEA 3e+05 -T_pulse_NEA_to_DEN 8e+04 -p_pulse_ghost_to_MH 0.1 -p_pulse_ghost_to_NEA 0.4 -p_pulse_NEA_to_DEN 0.03 -mu 1.25e-08 -N_super_ancestral 10000 -N_ancestral 10000 -N_ghost 10000 -N_modern 20000 -N_archaic 10000 -N_Neanderthal 1000 -N_Denisovan 1000 -out_prefix /home/trc468/notebooks/archaics/simulations_testing_231118/231118_1132 -num_MH 50 -mu 1.25e-08 -r 1e-08 -L 1e+07 -gen 29
    python /home/trc468/mars/simulate_archaic_model_231116.py -model A -T_supersuper_archaic 1.2e+06 -T_super_archaic 1e+06 -T_modern_archaic 8e+05 -T_den_nea 4e+05 -T_pulse_ghost_to_MH 1e+05 -T_pulse_MH_to_NEA 3e+05 -T_pulse_superghost_to_DEN 1e+05 -T_pulse_NEA_to_DEN 8e+04 -p_pulse_ghost_to_MH 0.1 -p_pulse_ghost_to_NEA 0.4 -p_pulse_superghost_to_DEN 0.05 -p_pulse_MH_to_NEA 0.06 -p_pulse_NEA_to_DEN 0.03 -mu 1.25e-08 -N_supersuper_ancestral 10000 -N_superghost 10000 -N_super_ancestral 10000 -N_ancestral 10000 -N_ghost 10000 -N_modern 20000 -N_archaic 10000 -N_Neanderthal 1000 -N_Denisovan 1000 -out_prefix /tmp/deletemetest -num_MH 50 -mu 1.25e-08 -r 1e-08 -L 1e+07 -gen 29


"""

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
parser.add_argument('-gen','--generation_time',help='Years per generation',required=True,type=int,default=29)
parser.add_argument('-NEA_age','--NEA_age',help='age of sampling for Neanderthal',required=False,type=int,default=3000)
parser.add_argument('-DEN_age','--DEN_age',help='age of sampling for Neanderthal',required=False,type=int,default=2500)
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

parser.add_argument('-num_MH','--num_MH',help='Number of MH haplotypes to simulate',required=True,type=int)
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
if model!='A' and model!='C':
    print(f'ERROR: model must be A or C; currently it is {model}. Aborint')
    sys.exit()

if out_prefix==None:
    out_prefix = currentdir

if os.path.isdir(os.path.dirname(out_prefix)) is False:
    os.makedirs(os.path.dirname(out_prefix))


MH_indices = np.array([i for i in range(0,num_MH)]) # indices of MH haplotypes
NEA_indices = np.array([num_MH,num_MH+1]) # indices of NEA haplotypes
DEN_indices = np.array([num_MH+2,num_MH+3]) # indices of DEN haplotypes


if model=='C':
    variable_names = [
        "N_super_ancestral",
        "N_ancestral",
        "N_ghost",
        "N_modern",
        "N_archaic",
        "N_Neanderthal",
        "N_Denisovan",
        "T_super_archaic",
        "T_modern_archaic",
        "T_den_nea",
        "T_pulse_ghost_to_MH",
        "T_pulse_ghost_to_NEA",
        "T_pulse_NEA_to_DEN"
        "p_pulse_ghost_to_MH",
        "p_pulse_ghost_to_NEA",
        "p_pulse_NEA_to_DEN"
        ]
    try:
        T_super_archaic = T_super_archaic/generation_time # split time in generations of root of (ancestral_humans,ghost)
        T_modern_archaic = T_modern_archaic/generation_time # split time in generations of main human lineage and lineage leading to (DEN,NEA)
        T_den_nea = T_den_nea/generation_time # split time in generations between NEA and DEN
        T_pulse_ghost_to_MH = T_pulse_ghost_to_MH/generation_time  # time in generations at which super archaic ghost introgresses into modern humans
        T_pulse_ghost_to_NEA = T_pulse_ghost_to_NEA/generation_time # time in generations at which super archaic ghost introgresses into Neanderthals
        T_pulse_NEA_to_DEN = T_pulse_NEA_to_DEN/generation_time # time in generations at which Neanderthals introgress into Denisovans   

        # configure demography
        demography = msprime.Demography()
        demography.add_population(name="super_ancestral", description="super Ancestral population, including ghost", initial_size=N_super_ancestral)
        demography.add_population(name="ancestral", description="Ancestral population", initial_size=N_ancestral)
        demography.add_population(name="ghost", description="ghost population", initial_size=N_ghost)
        demography.add_population(name="modern", description="modern human lineage", initial_size=N_modern)
        demography.add_population(name="archaic", description="lineage ancestral to Neanderthals and Denisovans, but derived wrt 'ancestral' ", initial_size=N_archaic)
        demography.add_population(name="neanderthal", description="neanderthals' ", initial_size=N_Neanderthal)
        demography.add_population(name="denisovan", description="denisovans", initial_size=N_Denisovan)

        # add split events
        demography.add_population_split(time=T_super_archaic, ancestral="super_ancestral", derived=["ancestral","ghost"])
        demography.add_population_split(time=T_modern_archaic, ancestral="ancestral", derived=["modern","archaic"])
        demography.add_population_split(time=T_den_nea, ancestral="archaic", derived=["neanderthal","denisovan"])

            
        # add pulse introgression events
        demography.add_mass_migration(time=T_pulse_ghost_to_MH, source='modern', dest='ghost', proportion=p_pulse_ghost_to_MH)
        demography.add_mass_migration(time=T_pulse_ghost_to_NEA, source='neanderthal', dest='ghost', proportion=p_pulse_ghost_to_NEA)    
        demography.add_mass_migration(time=T_pulse_NEA_to_DEN, source='denisovan', dest='neanderthal', proportion=p_pulse_NEA_to_DEN)    
        demography.sort_events()
    except:
        print(f'ERROR: missing input parameters for model C. The required parameters are')
        for qq in variable_names:
            print(f'\t{qq}')
        print(f'Aborting.')
        sys.exit()

elif model=='A':
    # configure demography
    variable_names = [
    "N_supersuper_ancestral",
    "N_superghost",
    "N_super_ancestral",
    "N_ancestral",
    "N_ghost",
    "N_modern",
    "N_archaic",
    "N_Neanderthal",
    "N_Denisovan",
    "T_supersuper_archaic",
    "T_super_archaic",
    "T_modern_archaic",
    "T_den_nea",
    "T_pulse_superghost_to_DEN",
    "T_pulse_ghost_to_MH",
    "T_pulse_MH_to_NEA",
    "T_pulse_NEA_to_DEN"
    "p_pulse_superghost_to_DEN",
    "p_pulse_ghost_to_MH",
    "p_pulse_MH_to_NEA",
    "p_pulse_NEA_to_DEN"    
    ]
    try:
        demography = msprime.Demography()
        demography.add_population(name="supersuper_ancestral", description="supersuper ancestral population, including 'super ancestral' and 'super ghost", initial_size=N_supersuper_ancestral)
        demography.add_population(name="super_ghost", description="super ghost which introgresses into Denisovan", initial_size=N_superghost)
        demography.add_population(name="super_ancestral", description="super Ancestral population, including ghost", initial_size=N_super_ancestral)
        demography.add_population(name="ancestral", description="Ancestral population", initial_size=N_ancestral)
        demography.add_population(name="ghost", description="ghost population", initial_size=N_ghost)
        demography.add_population(name="modern", description="modern human lineage", initial_size=N_modern)
        demography.add_population(name="archaic", description="lineage ancestral to Neanderthals and Denisovans, but derived wrt 'ancestral' ", initial_size=N_archaic)
        demography.add_population(name="neanderthal", description="neanderthals' ", initial_size=N_Neanderthal)
        demography.add_population(name="denisovan", description="denisovans", initial_size=N_Denisovan)

        # add split events
        demography.add_population_split(time=T_supersuper_archaic, ancestral="supersuper_ancestral", derived=["super_ancestral","super_ghost"])
        demography.add_population_split(time=T_super_archaic, ancestral="super_ancestral", derived=["ancestral","ghost"])
        demography.add_population_split(time=T_modern_archaic, ancestral="ancestral", derived=["modern","archaic"])
        demography.add_population_split(time=T_den_nea, ancestral="archaic", derived=["neanderthal","denisovan"])

        # add pulse introgression events
        demography.add_mass_migration(time=T_pulse_superghost_to_DEN, source='denisovan', dest='super_ghost', proportion=p_pulse_superghost_to_DEN)
        demography.add_mass_migration(time=T_pulse_ghost_to_MH, source='modern', dest='ghost', proportion=p_pulse_ghost_to_MH)
        demography.add_mass_migration(time=T_pulse_MH_to_NEA, source='neanderthal', dest='modern', proportion=p_pulse_MH_to_NEA)    
        demography.add_mass_migration(time=T_pulse_NEA_to_DEN, source='denisovan', dest='neanderthal', proportion=p_pulse_NEA_to_DEN)    
        demography.sort_events()
    except:
        print(f'ERROR: missing input parameters for model A. The required parameters are')
        for qq in variable_names:
            print(f'\t{qq}')
        print(f'Aborting.')
        sys.exit()



sim = msprime.sim_ancestry( # simulate ancestry
    samples=[msprime.SampleSet(num_samples=num_MH, ploidy=1,population='modern',time=0), 
            msprime.SampleSet(num_samples=1, ploidy=2,population='neanderthal',time=NEA_age), \
            msprime.SampleSet(num_samples=1, ploidy=2,population='denisovan',time=DEN_age)],\
            demography=demography,sequence_length=L,recombination_rate=r)
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


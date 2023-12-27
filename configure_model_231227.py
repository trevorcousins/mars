import msprime
import sys
import pdb 

def configure_demography(
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
    N_modern_present,
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
    ):
    if model=='C':
        variable_names = [ # required
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
            "T_pulse_NEA_to_DEN",
            "p_pulse_ghost_to_MH",
            "p_pulse_ghost_to_NEA",
            "p_pulse_NEA_to_DEN",
            "N_AMH"
            ]
        for i in variable_names:
            if eval(i)==None:
                print(f"Parameter {i} is required for model C, but it is not given. Aborting.")
                sys.exit()

        try:
            T_super_archaic = T_super_archaic/generation_time # split time in generations of root of (ancestral_humans,ghost)
            T_modern_archaic = T_modern_archaic/generation_time # split time in generations of main human lineage and lineage leading to (DEN,NEA)
            T_den_nea = T_den_nea/generation_time # split time in generations between NEA and DEN
            T_pulse_ghost_to_MH = T_pulse_ghost_to_MH/generation_time  # time in generations at which super archaic ghost introgresses into modern humans
            T_pulse_ghost_to_NEA = T_pulse_ghost_to_NEA/generation_time # time in generations at which super archaic ghost introgresses into Neanderthals
            T_pulse_NEA_to_DEN = T_pulse_NEA_to_DEN/generation_time # time in generations at which Neanderthals introgress into Denisovans   
            T_AMH_expand = T_AMH_expand/generation_time
            T_ghost_BN_start = T_ghost_BN_start/generation_time
            T_ghost_BN_end = T_ghost_BN_end/generation_time


            # configure demography
            demography = msprime.Demography()
            demography.add_population(name="super_ancestral", description="super Ancestral population, including ghost", initial_size=N_super_ancestral)
            demography.add_population(name="ancestral", description="Ancestral population", initial_size=N_ancestral)
            demography.add_population(name="ghost", description="ghost population", initial_size=N_ghost_recent)
            demography.add_population(name="modern", description="modern human lineage", initial_size=N_modern)
            demography.add_population(name="archaic", description="lineage ancestral to Neanderthals and Denisovans, but derived wrt 'ancestral' ", initial_size=N_archaic)
            demography.add_population(name="neanderthal", description="neanderthals' ", initial_size=N_Neanderthal)
            demography.add_population(name="denisovan", description="denisovans", initial_size=N_Denisovan)
            demography.add_population_parameters_change(time=0, population="modern", initial_size=N_modern_present,growth_rate=N_MH_expand_rate)
            demography.add_population_parameters_change(time=T_AMH_expand, population="modern", initial_size=N_AMH,growth_rate=0)
            demography.add_population_parameters_change(time=T_pulse_ghost_to_MH, population="modern", initial_size=N_modern,growth_rate=0)
            
            demography.add_population_parameters_change(time=T_ghost_BN_start, population="ghost", initial_size=N_ghost,growth_rate=0)
            demography.add_population_parameters_change(time=T_ghost_BN_end, population="ghost", initial_size=N_ghost*N_ghost_BN_intensity,growth_rate=0)
            demography.add_population_parameters_change(time=T_pulse_ghost_to_NEA, population="ghost", initial_size=N_ghost,growth_rate=0)

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
    elif model=='Cprime':
        variable_names = [ # required
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
            "T_pulse_NEA_to_DEN",
            "p_pulse_ghost_to_MH",
            "p_pulse_ghost_to_NEA",
            "p_pulse_NEA_to_DEN",
            "N_AMH",
            "T_supersuper_archaic",
            "T_pulse_superghost_to_DEN",
            "N_supersuper_ancestral",
            "N_superghost",
            "p_pulse_superghost_to_DEN",
            "N_Neanderthal_arch",
            "T_NEA_admix",
            "N_ghost_nea",
            "N_ghost_human",
            ]
        for i in variable_names:
            if eval(i)==None:
                print(f"Parameter {i} is required for model C, but it is not given. Aborting.")
                sys.exit()

        try:
            T_supersuper_archaic = T_supersuper_archaic/generation_time # split time in generations of root of (superghost,(ancestral_humans,ghost))
            T_super_archaic = T_super_archaic/generation_time # split time in generations of root of (ancestral_humans,ghost)
            T_modern_archaic = T_modern_archaic/generation_time # split time in generations of main human lineage and lineage leading to (DEN,NEA)
            T_den_nea = T_den_nea/generation_time # split time in generations between NEA and DEN
            T_pulse_ghost_to_MH = T_pulse_ghost_to_MH/generation_time  # time in generations at which super archaic ghost introgresses into modern humans
            T_pulse_ghost_to_NEA = T_pulse_ghost_to_NEA/generation_time # time in generations at which the ghost population branches into two other ghosts, one that contributes to humans and the other to Neanderthals
            T_pulse_NEA_to_DEN = T_pulse_NEA_to_DEN/generation_time # time in generations at which Neanderthals introgress into Denisovans  
            T_pulse_superghost_to_DEN = T_pulse_superghost_to_DEN/generation_time # time in generations at which superghost introgress into Denisovans   
            T_NEA_admix = T_NEA_admix/generation_time # time in generations at which the ghost population (ghost_nea) contributes to Neanderthals
            
            T_AMH_expand = T_AMH_expand/generation_time
            T_ghost_BN_start = T_ghost_BN_start/generation_time
            T_ghost_BN_end = T_ghost_BN_end/generation_time


            # configure demography
            demography = msprime.Demography()
            demography.add_population(name="supersuper_ancestral", description="supersuper ancestral population, ancestor of (superghost,(ancestral_humans,ghost))", initial_size=N_supersuper_ancestral)
            demography.add_population(name="super_ancestral", description="super Ancestral population, including ghost", initial_size=N_super_ancestral)
            demography.add_population(name="ancestral", description="Ancestral population", initial_size=N_ancestral)
            demography.add_population(name="ghost", description="ghost population that is the ancestor of the ghost populations that introgress into neanderthals and humans", initial_size=N_ghost)
            demography.add_population(name="ghost_nea", description="ghost population that contributes to neanderthals", initial_size=N_ghost_nea)
            demography.add_population(name="ghost_human", description="ghost population that contributes to modern humans", initial_size=N_ghost_human)
            demography.add_population(name="superghost", description="super ghost population (that introgresses into Neanderthal)", initial_size=N_superghost)
            demography.add_population(name="modern", description="modern human lineage", initial_size=N_modern_present)
            demography.add_population(name="archaic", description="lineage ancestral to Neanderthals and Denisovans, but derived wrt 'ancestral' ", initial_size=N_archaic)
            demography.add_population(name="neanderthal", description="neanderthals' ", initial_size=N_Neanderthal)
            demography.add_population(name="denisovan", description="denisovans", initial_size=N_Denisovan)
            demography.add_population_parameters_change(time=0, population="modern", initial_size=N_modern_present,growth_rate=N_MH_expand_rate)
            demography.add_population_parameters_change(time=T_AMH_expand, population="modern", initial_size=N_AMH,growth_rate=0)
            demography.add_population_parameters_change(time=T_pulse_ghost_to_MH, population="modern", initial_size=N_modern,growth_rate=0)
            demography.add_population_parameters_change(time=T_ghost_BN_start, population="ghost", initial_size=N_ghost,growth_rate=0)
            demography.add_population_parameters_change(time=T_ghost_BN_end, population="ghost", initial_size=N_ghost*N_ghost_BN_intensity,growth_rate=0)
            demography.add_population_parameters_change(time=T_pulse_ghost_to_NEA, population="ghost", initial_size=N_ghost,growth_rate=0)
            demography.add_population_parameters_change(time=T_NEA_admix, population="neanderthal", initial_size=N_Neanderthal_arch,growth_rate=0)


            # add split events
            demography.add_population_split(time=T_supersuper_archaic, ancestral="supersuper_ancestral", derived=["super_ancestral","superghost"])
            demography.add_population_split(time=T_super_archaic, ancestral="super_ancestral", derived=["ancestral","ghost"])
            demography.add_population_split(time=T_modern_archaic, ancestral="ancestral", derived=["modern","archaic"])
            demography.add_population_split(time=T_den_nea, ancestral="archaic", derived=["neanderthal","denisovan"])
            demography.add_population_split(time=T_pulse_ghost_to_NEA, ancestral="ghost", derived=["ghost_nea","ghost_human"]) # ghost splits to two derived ghosts

                
            # pdb.set_trace()
            # add pulse introgression events
            demography.add_mass_migration(time=T_pulse_ghost_to_MH, source='modern', dest='ghost_human', proportion=p_pulse_ghost_to_MH)
            demography.add_mass_migration(time=T_NEA_admix, source='neanderthal', dest='ghost_nea', proportion=p_pulse_ghost_to_NEA)    
            demography.add_mass_migration(time=T_pulse_NEA_to_DEN, source='denisovan', dest='neanderthal', proportion=p_pulse_NEA_to_DEN)
            if Cprime_introgression=='a': # superghost introgresses into Denisovans
                if not T_pulse_superghost_to_DEN<T_den_nea:
                    print(f'ERROR: Introgression from superghost into Denisovans must be less than the time at which Denisovans and Neanderthals split, but T_pulse_superghost_to_DEN={T_pulse_superghost_to_DEN} and  T_den_nea={T_den_nea}. Aborting.')
                    sys.exit()
                demography.add_mass_migration(time=T_pulse_superghost_to_DEN, source='denisovan', dest='superghost', proportion=p_pulse_superghost_to_DEN)
            elif Cprime_introgression=='b': # superghost introgresses into the ancestors of (Denisovans,Neanderthals)
                if not (T_pulse_superghost_to_DEN>T_den_nea and T_pulse_superghost_to_DEN<T_modern_archaic) :
                    print(f'ERROR: Introgression from superghost into archaic (ancestors of Neanderthals and Denisovans) must be more than the time at which Denisovans and Neanderthals split, and less than time at which (Neanderthals,Denisovans) and humans split, but T_pulse_superghost_to_DEN={T_pulse_superghost_to_DEN} and T_den_nea={T_den_nea} and T_modern_archaic={T_modern_archaic}.Aborting.')
                    sys.exit()    
                demography.add_mass_migration(time=T_pulse_superghost_to_DEN, source='archaic', dest='superghost', proportion=p_pulse_superghost_to_DEN)
            elif Cprime_introgression=='c': # superghost introgresses into the ancestors of (modern humans, ((Denisovans,Neanderthals)))
                if not (T_pulse_superghost_to_DEN<T_super_archaic and T_pulse_superghost_to_DEN>T_modern_archaic) :
                    print(f'ERROR: Introgression from superghost into ancestral (ancestors of (Neanderthals,Denisovans) and modern humans) must be more than the time at which (Neanderthals,Denisovans) and modern humans split, and less than time at which ((Neanderthals,Denisovans),modern_humans) split from the ghost that introgresses into humans and Neanderthals, but T_pulse_superghost_to_DEN={T_pulse_superghost_to_DEN} and T_modern_archaic={T_modern_archaic} and T_super_archaic={T_super_archaic}.Aborting.')
                    sys.exit()
                demography.add_mass_migration(time=T_pulse_superghost_to_DEN, source='ancestral', dest='superghost', proportion=p_pulse_superghost_to_DEN)
            else:
                print(f'ERROR: Cprime_introgression={Cprime_introgression} is not a or b or c. Aborting')
                sys.exit()
            demography.sort_events()
        except:
            print(f'ERROR: missing input parameters for model C. The required parameters are')
            for qq in variable_names:
                print(f'\t{qq}')
            print(f'Aborting.')
            sys.exit()
    elif model=='A':
        # configure demography
        variable_names = [ # required
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
        "p_pulse_NEA_to_DEN",    
        "N_AMH"
        ]
        for i in variable_names:
            if eval(i)==None:
                print(f"Parameter {i} is required for model A, but it is not given. Aborting.")
                sys.exit()

        try:
            T_supersuper_archaic = T_supersuper_archaic/generation_time # split time in generations of root of (ancestral_humans,ghost)
            T_super_archaic = T_super_archaic/generation_time # split time in generations of root of (ancestral_humans,ghost)
            T_modern_archaic = T_modern_archaic/generation_time # split time in generations of main human lineage and lineage leading to (DEN,NEA)
            T_den_nea = T_den_nea/generation_time # split time in generations between NEA and DEN
            T_pulse_superghost_to_DEN = T_pulse_superghost_to_DEN/generation_time # time in generations at which superghost lineage introgresses into Denisovans
            T_pulse_ghost_to_MH = T_pulse_ghost_to_MH/generation_time  # time in generations at which super archaic ghost introgresses into modern humans
            T_pulse_MH_to_NEA = T_pulse_MH_to_NEA/generation_time # time in generations at which modern humans introgresses into Neanderthals
            T_pulse_NEA_to_DEN = T_pulse_NEA_to_DEN/generation_time # time in generations at which Neanderthals introgress into Denisovans   
            T_AMH_expand = T_AMH_expand/generation_time
            
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
            demography.add_population_parameters_change(time=0, population="modern", initial_size=N_modern_present,growth_rate=N_MH_expand_rate)
            demography.add_population_parameters_change(time=T_AMH_expand, population="modern", initial_size=N_AMH,growth_rate=0)
            demography.add_population_parameters_change(time=T_pulse_ghost_to_MH, population="modern", initial_size=N_modern,growth_rate=0)
            
            demography.add_population_parameters_change(time=T_ghost_BN_start, population="ghost", initial_size=N_ghost,growth_rate=0)
            demography.add_population_parameters_change(time=T_ghost_BN_end, population="ghost", initial_size=N_ghost*N_ghost_BN_intensity,growth_rate=0)


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
    print(demography)
    print(demography.debug())
    return demography
    

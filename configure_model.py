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
    

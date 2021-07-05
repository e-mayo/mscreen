from pathlib import Path
from shutil import which

# import os
# for r,d,f in os.walk("c:\\"):
#     for files in f:
#          if files == "something.exe":
#               print os.path.join(r,files)
class DockInFile:
    def __init__(self, DOCK_HOME=None, **kwargs):
        #Rigid Body and Flexible Ligand Docking Input Parameters
        default_parameter={"conformer_search_type":"flex", # Choose the type of docking simulation: Rigid body docking (rigid):Flexible ligand docking (flex) 	 #  de novo ligand design (denovo) 
                            "user_specified_anchor":"no",   # Will the user specify an anchor file?
                            "atom_in_anchor":'"C1",1',   # If the user specifies an anchor which atom label in the anchor? ["C1",1]
                            "limit_max_anchors":"no",   # Will the user limit the maximum number of anchors docked?
                            "max_anchor_num":1, 	 # If the user limits the number - maximum number of anchors allowed
                            "min_anchor_size":5, 	 # Minimum number of atoms in the anchor
                            "pruning_use_clustering":"yes", 	 # Will pruning the conformers use a clustering algorithm?
                            "pruning_max_orients":1000, 	 # How many orients will be generated prior to pruning?
                            "pruning_clustering_cutoff":100, 	 # Maximum number of clusterheads retained from pruning
                            "pruning_conformer_score_cutoff":100, 	 # Maximum score allowed for conformers (kcal/mol)
                            "pruning_conformer_score_scaling_factor":1, 	 # Score cutoff scaling factor to increase of reduce the score cutoff as molecules rebuild
                            "use_clash_overlap":"no", 	 # Flag to check for overlapping atomic volumes during anchor and grow
                            "clash_overlap":0.5, 	 # A clash exists id the distance between a pair of atoms is less than the clash overlap times the sum of their atom type radii
                            "write_growth_trees":"no",  # Generate large growth tree files (increases memory usage - recommended to concatenate and compress growth tree branches),
                            #Ligand RMSD Parameters
                            "calculate_rmsd":"no",  # Calculate root mean square deviation?	
                            "use_rmsd_reference_mol":"no", #Does the user want to use a reference molecule to calculate rmsd?
                            "rmsd_reference_filename":None,  #The path to the rmsd reference molecule	
                            # Internal Energy Parameters
                            "use_interal_energy":"yes",    # Does the user want to use internal energy for growth and or minimization (only repulsive VDW)
                            "internal_energy_rep_exp":12,  # The VDW exponent only when use internal energy is turned on(DOCK is optimized for default value)
                            "internal_energy_cutoff":100.0,  # All conformers with an internal energy value above this cutoff are pruned(only turned on use internal energy is used)	
                            # Bump Filter Parameters
                            "bump_filter":"no", # Does the user want to perform bump filter?	no
                            "bump_grid_prefix":"grid", # The prefix to the grid file containing the desired bump grid (only turned on when bump filter is used)	grid
                            "max_bumps_anchor":12, #	The maximum allowed number of bumps for an anchor to pass the filter	12
                            "max_bumps_growth":12,  #	The maximum allowed number of bumps for a molecule to pass the filter	12
                            # Contact Score Parameters
                            "contact_score_primary":"no",#	Does the user want to perform contact scoring as primary scoring function	no
                            "contact_score_secondary":"no",	# Does the user want to perform contact scoring as a secondary scoring function?	no
                            "contact_score_cutoff_distance":4.5,	# The distance threshold defining a contact when contact scoring is turned on	4.5
                            "contact_score_clash_overlap":0.75, # Contact definition for use with intramolecular scoring when contact scoring is turned on	0.75
                            "contact_score_clash_penalty":50, # The penalty for each contact overlap made when contact score is turned on	50
                            "contact_score_grid_prefix":"grid",	# The prefix to the grid files containing the desired contact when contact score is turned on	grid
                            # Grid Score Parameters
                            "grid_score_primary":"yes", 	      # Does the user want to perform grid-based energy scoring as the primary scoring function?	yes
                            "grid_score_secondary":"yes",      # Does the user want to perform grid-based energy scoring as the secondary scoring function?	yes
                            "grid_score_rep_rad_scale":1.0,    # Scalar multiplier of the radii for the repulsive portion of the VDW energy component only when grid score is turned on	1.0
                            "grid_score_vdw_scale":1,          # Scalar multiplier of the VDW energy component	1
                            "grid_score_turn_off_vdw":"yes",   # A flag to turn off vdw portion of scoring function when grid score vdw scale:0	yes
                            "grid_score_es_scale":1,           # Flag to scale up or down the es portion of the scoring function when es scale is turned on	1
                            "grid_score_turn_off_es":"yes",   # A flag to turn off es portion of scoring function when grid score es scale:0	yes
                            "grid_score_grid_prefix":"grid", # The prefix to the grid files containing the desired nrg/bmp grid	grid
                            # Structure Input
                            "ligand_atom_file": "lig.mol2", # 1HW9.lig.min_scored.mol2
                            "ligand_outfile_prefix":"lig",  #,
                            # Others
                            "use_internal_energy":"yes",
                            "limit_max_ligands":"no",
                            "skip_molecule":"no",
                            "read_mol_solvation":"no",
                            "use_database_filter":"no",
                            "orient_ligand":"yes",
                            "score_molecules":"yes",
                            "multigrid_score_secondary":"no",
                            "dock35_score_secondary":"no",
                            "continuous_score_secondary":"no",
                            "footprint_similarity_score_secondary":"no",
                            "pharmacophore_score_secondary":"no",
                            "descriptor_score_secondary":"no",
                            "gbsa_zou_score_secondary":"no",
                            "gbsa_hawkins_score_secondary":"no",
                            "SASA_score_secondary":"no",
                            "amber_score_secondary":"no",
                            "minimize_ligand":"yes",
                            "simplex_max_iterations":1000,
                            "simplex_tors_premin_iterations":0,
                            "simplex_max_cycles":1,
                            "simplex_score_converge":0.1,
                            "simplex_cycle_converge":1,
                            "simplex_trans_step":1,
                            "simplex_rot_step":0.1,
                            "simplex_tors_step":10,
                            "simplex_random_seed":0,
                            "simplex_restraint_min":"no",
                            "atom_model":"all",
                            "vdw_defn_file":None,
                            "flex_defn_file":None,
                            "flex_drive_file":None,
                            "write_orientations":"yes",
                            "write_conformations":"yes",
                            "num_scored_conformers":100,
                            "rank_ligands":"no",
                            "minimize_anchor":"yes"}
        
        self.default_parameter = {**default_parameter,**kwargs}
        if not DOCK_HOME:
            DOCK_HOME = which("dock6")
        if not DOCK_HOME:
            print('no DOCK exectable found')
        else:
            self.DOCK_HOME = Path(DOCK_HOME)
            self.PARAMETERS_HOME = self.DOCK_HOME / 'parameters'
    
    def write_energy_mininization_input(self,file_name,**kwargs):
        minim_parameters = {"conformer_search_type":"rigid",
                            "calculate_rmsd":"yes",
                            "use_rmsd_reference_mol":"yes",
                            "grid_score_secondary":"no",
                            "grid_score_grid_prefix":2,# grid_score_grid_prefix
                            "ligand_atom_file":"lig",# ligand_atom_file
                            "ligand_outfile_prefix":"lig.min",
                            "orient_ligand":"no",
                            "simplex_restraint_min":"yes",
                            "simplex_coefficient_min":10.0}
        minim_parameters = {**minim_parameters, **kwargs}   
        self.__write__dockin_file__(file_name,minim_parameters)

        
            
    def write_flexible_docking(self,file_name,**kwargs):
        flex_parameters = {"conformer_search_type":"flex",
                           "atom_in_anchor":"no",
                           "grid_score_secondary":"no",
                           "grid_score_grid_prefix":"grid_score_grid_prefix",
                           "ligand_atom_file":"lig.mol2",
                           "ligand_outfile_prefix":"lig.flex",
                           "num_scored_conformers":20}
        self.__write__dockin_file__(file_name,flex_parameters)

            
    def write_fixed_anchor_docking(self,file_name,**kwargs):
        fixed_parameters = {"conformer_search_type":"flex",
                           "atom_in_anchor":"no",
                           "grid_score_secondary":"no",
                           "grid_score_grid_prefix":"grid_score_grid_prefix",
                           "ligand_atom_file":"lig.mol2",
                           "ligand_outfile_prefix":"lig.fixed",
                           "num_scored_conformers":20,
                           "orient_ligand":"no"}
        self.__write__dockin_file__(file_name,fixed_parameters)

    def write_virtual_screening(self,file_name,**kwargs):
        vs_parameters = {"conformer_search_type":"flex",
                        "atom_in_anchor":"no",
                        "grid_score_secondary":"no",
                        "grid_score_grid_prefix":"grid_score_grid_prefix",
                        "ligand_atom_file":"lig.mol2",
                        "ligand_outfile_prefix":"lig.fixed",
                        "write_orientations":"no",
                        "num_scored_conformers":10,
                        "rank_ligands":"no"}
        self.__write__dockin_file__(file_name,vs_parameters)
                         

    def write_cartesian_minimization(self,file_name,**kwargs):
        cmin_parameters = {"conformer_search_type":"rigid",
                            "calculate_rmsd":"yes",
                            "use_rmsd_reference_mol":"no",
                            "grid_score_secondary":"no",
                            "grid_score_grid_prefix":2,# grid_score_grid_prefix
                            "continuous_score_primary":"yes",
                            "continuous_score_secondary":"no",
                            "cont_score_rec_filename":"rec_file",
                            "minimize_ligand":"yes",
                            "ligand_atom_file":"lig",# ligand_atom_file
                            "ligand_outfile_prefix":"lig.min",
                            "orient_ligand":"no",
                            "simplex_restraint_min":"no"}    
        self.__write__dockin_file__(file_name,cmin_parameters)
        
        
    def __write__dockin_file__(self,file_name,parameters):
        with open(file_name, 'w') as f:
            for key,val in parameters.items():
                f.write(f'{key:<20}{val}\n')
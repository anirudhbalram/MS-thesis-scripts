import numpy as np
import pandas as pd
import caesar
import yt
from tqdm import tqdm
import random
random.seed(10)

#this function makes a file "Clouds_per_Core_m{config["boxsize"]}_z_2_MS.txt" and writes into it the gas particle indicies of all clouds that we need to run. Each cloud is separated by a SPACE and "n_clouds_per_line" is the number of clouds in each line of the txt file aka number of clouds that will run in each core

def create_cloudspercore_list(config):
    
    obj = caesar.load(config["caesarfilename"])

    if config["mode"] == 'randomize':
        galaxies_list = [(gal.GroupID,gal.glist) for gal in obj.galaxies if (gal.masses['total']>config["min_mass"] and gal.masses['total']<=config["max_mass"])]
        galaxies_list = random.sample(galaxies_list, config["n_galaxies_sample"])
    else:    #although we set mode=total, this will only take into account the selected galaxies as saved in snap_z_2_MS.csv
        selected_ids_df = pd.read_csv("snap_z_2_MS.csv")
        selected_ids_array = np.array(selected_ids_df['GroupID'])

        galaxies_list = [(gal.GroupID,gal.glist) for gal in [obj.galaxies[int(GroupID)] for GroupID in selected_ids_array]]    #making a list of galaxy ids and their gas particle indices

    n_of_gals = len(galaxies_list)
    n_of_clouds = sum([len(galaxies_list[i][1]) for i in np.arange(n_of_gals)])
    n_clouds_per_line = 900    #setting the number of clouds each core will run. This can be increased or decreased accordingly to use more/less cores. Accordingly however, you must change the number of cores used in the slick_run_jobscript.sh script for running Despotic!!!

    param_filename = f'{config["output_dir"]}/Clouds_per_Core_m{config["boxsize"]}_z_2_MS.txt'
    with open(param_filename, 'w') as f:
        num_lines = 1
        gal_clouds = np.concatenate(np.array([gal[1] for gal in galaxies_list],dtype=object))
        for (i,), cloud in np.ndenumerate(gal_clouds.flatten()):
            if i % n_clouds_per_line == 0:
                if i != 0:
                    f.write("\n")
                    num_lines += 1
            f.write(f"{cloud} ")
    return param_filename, num_lines

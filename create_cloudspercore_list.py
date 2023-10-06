import numpy as np
import pandas as pd
import caesar
import yt
from tqdm import tqdm
import random
random.seed(10)

def create_cloudspercore_list(config):
    
    #yt_snap = yt.load(config["ytfilename"])
    #yt_snap_data = yt_snap305.all_data()
    obj = caesar.load(config["caesarfilename"])

    if config["mode"] == 'randomize':
        galaxies_list = [(gal.GroupID,gal.glist) for gal in obj.galaxies if (gal.masses['total']>config["min_mass"] and gal.masses['total']<=config["max_mass"])]
        galaxies_list = random.sample(galaxies_list, config["n_galaxies_sample"])
    else:
        #gal = obj.galaxies[11]
        selected_ids_df = pd.read_csv("snap_z_0.8_MS.csv")
        selected_ids_array = np.array(selected_ids_df['GroupID'])

        galaxies_list = [(gal.GroupID,gal.glist) for gal in [obj.galaxies[int(GroupID)] for GroupID in selected_ids_array]]
        #galaxies_list = [(gal.GroupID,gal.glist) for gal in obj.galaxies[1:12]]
        #galaxies_list = [(gal.GroupID,gal.glist) for gal in obj.galaxies]

    n_of_gals = len(galaxies_list)
    n_of_clouds = sum([len(galaxies_list[i][1]) for i in np.arange(n_of_gals)])
    n_clouds_per_line = 900
    #n_clouds_per_line = n_of_clouds/900

    param_filename = f'{config["output_dir"]}/Clouds_per_Core_m{config["boxsize"]}_z_0.8_MS.txt'
    with open(param_filename, 'w') as f:
        num_lines = 1
        gal_clouds = np.concatenate(np.array([gal[1] for gal in galaxies_list],dtype=object))
            #sfr_gal = np.sum(yt_snap_data["PartType0", "StarFormationRate"][gal_clouds].value)
            #if sfr_gal <= 0:
            #    continue
        for (i,), cloud in np.ndenumerate(gal_clouds.flatten()):
            if i % n_clouds_per_line == 0:
                if i != 0:
                    f.write("\n")
                    num_lines += 1
                    #f.write(f"{gal[0]} ")
            f.write(f"{cloud} ")
        #f.write("\n")
        #num_lines += 1
    return param_filename, num_lines
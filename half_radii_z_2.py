'''
This script is to compute the half light and half mass radii of galaxies based on the output stored in lim_df. CO sleds are also computed and stored. The half radii for each (selected) galaxy at a given snapshot is computed in individual cores.

!!!Run the following few lines of code in this comment individually in a login node!!! This ensures that an empty half-radii file and a sled file exist. 

radii_df = pd.DataFrame({"Galaxy_ID":[],"g_Metallicity":[],"g_SFR":[],"g_stellar_mass":[],"g_gas_mass":[],"g_dust_mass":[],"g_Mol_gas_mass":[],"g_Metal_mass":[],"alpha_CO":[],"clean_pct":[],"half_CO10_radius":[],"half_CO21_radius":[],"half_CO32_radius":[],"half_CO43_radius":[],"half_CO54_radius":[],"half_CO65_radius":[],"half_CO76_radius":[],"half_CO87_radius":[],"half_CO98_radius":[],"half_CI10_radius":[],"half_CI21_radius":[],"half_CII_radius":[],"half_Mol_gas_mass_radius":[],"half_total_gas_mass_radius":[],"half_SFR_radius":[],"half_Metal_mass_radius":[],"r25_CO10_radius":[],"r25_CO21_radius":[],"r25_CO32_radius":[],"r25_CO43_radius":[],"r25_CO54_radius":[],"r25_CO65_radius":[],"r25_CO76_radius":[],"r25_CO87_radius":[],"r25_CO98_radius":[],"r25_CI10_radius":[],"r25_CI21_radius":[],"r25_CII_radius":[],"r25_Mol_gas_mass_radius":[],"r25_total_gas_mass_radius":[],"r25_SFR_radius":[],"r25_Metal_mass_radius":[],"r75_CO10_radius":[],"r75_CO21_radius":[],"r75_CO32_radius":[],"r75_CO43_radius":[],"r75_CO54_radius":[],"r75_CO65_radius":[],"r75_CO76_radius":[],"r75_CO87_radius":[],"r75_CO98_radius":[],"r75_CI10_radius":[],"r75_CI21_radius":[],"r75_CII_radius":[],"r75_Mol_gas_mass_radius":[],"r75_total_gas_mass_radius":[],"r75_SFR_radius":[],"r75_Metal_mass_radius":[]})
radii_df.to_csv("results/snapshot_2_radii_projected.csv",index=False)

sled_df = pd.DataFrame({"Galaxy_ID":[],"sled_CO10":[],"sled_CO21":[],"sled_CO32":[],"sled_CO43":[],"sled_CO54":[],"sled_CO65":[],"sled_CO76":[],"sled_CO87":[],"sled_CO98":[],})
sled_df.to_csv("results/sled_snapshot_2.csv",index=False)

'''

import numpy as np
import pandas as pd
import caesar
import yt
import argparse

df = pd.read_csv("/blue/narayanan/a.ravishankar/slick/repo/results/lim_df_z_2_MS.csv")    #reading the saved luminosity datafile given as output by slick for a particular redshift
obj = caesar.load("/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0160_z2.000.hdf5")
selected_ids_df = pd.read_csv("/blue/narayanan/a.ravishankar/slick/repo/snap_z_2_MS.csv")    #loads a pandas dataframe showing the list of galaxies that have been selected based on the SFR-M* Main Sequence fit (see code in separate folder)
selected_ids_array = np.array(selected_ids_df['GroupID'])    #converts this dataframe to an array for ease of use

clouds_in_each_galaxy = [(gal.GroupID,gal.glist,gal.masses['gas'].in_units('Msun').value,gal.masses['H2'].in_units('Msun').value,gal.radii['baryon_half_mass'].in_units('pc').value,gal.metallicities['mass_weighted'].value,gal.masses['stellar'].in_units('Msun').value,gal.sfr.in_units("Msun/yr").value,gal.masses["dust"].in_units("Msun").value) for gal in [obj.galaxies[int(GroupID)] for GroupID in selected_ids_array]]    #creates a list containing properties of selected galaxies and also their gas particle list indices

yt_snap = yt.load("/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_160.hdf5")
yt_data = yt_snap.all_data()

#simple function to compute euclidean distance aka the radius in 2-D (i.e., projected radius)

def radii(coords_array,pos,proj):
    
    if proj==0:
        a = 1
        b = 2
    elif proj==1:
        a = 0
        b = 2
    else:
        a = 0
        b = 1
    
    radius = np.sqrt((coords_array[a]-pos[a])**2+(coords_array[b]-pos[b])**2)
    return float(radius)

#initialising a number of lists corresponding to sleds, half radii, etc.

sled_CO10_array = []
sled_CO21_array = []
sled_CO32_array = []
sled_CO43_array = []
sled_CO54_array = []
sled_CO65_array = []
sled_CO76_array = []
sled_CO87_array = []
sled_CO98_array = []

half_CO10_radius_array = []
half_CO21_radius_array = []
half_CO32_radius_array = []
half_CO43_radius_array = []
half_CO54_radius_array = []
half_CO65_radius_array = []
half_CO76_radius_array = []
half_CO87_radius_array = []
half_CO98_radius_array = []
half_CI10_radius_array = []
half_CI21_radius_array = []
half_CII_radius_array = []
half_Mol_gas_mass_radius_array = []
half_total_gas_mass_radius_array = []
half_SFR_radius_array = []
half_Metal_mass_radius_array = []

r25_CO10_radius_array = []
r25_CO21_radius_array = []
r25_CO32_radius_array = []
r25_CO43_radius_array = []
r25_CO54_radius_array = []
r25_CO65_radius_array = []
r25_CO76_radius_array = []
r25_CO87_radius_array = []
r25_CO98_radius_array = []
r25_CI10_radius_array = []
r25_CI21_radius_array = []
r25_CII_radius_array = []
r25_Mol_gas_mass_radius_array = []
r25_total_gas_mass_radius_array = []
r25_SFR_radius_array = []
r25_Metal_mass_radius_array = []

r75_CO10_radius_array = []
r75_CO21_radius_array = []
r75_CO32_radius_array = []
r75_CO43_radius_array = []
r75_CO54_radius_array = []
r75_CO65_radius_array = []
r75_CO76_radius_array = []
r75_CO87_radius_array = []
r75_CO98_radius_array = []
r75_CI10_radius_array = []
r75_CI21_radius_array = []
r75_CII_radius_array = []
r75_Mol_gas_mass_radius_array = []
r75_total_gas_mass_radius_array = []
r75_SFR_radius_array = []
r75_Metal_mass_radius_array = []

clean_pct_array = []
total_Mol_array = []
total_Metal_array = []
alpha_CO_array = []

parser = argparse.ArgumentParser(prog='halfradii')
parser.add_argument("--galinfoline", type=int)    #galinfoline stores the SLURM_ARRAY_TASK_ID which is used to pick a particular galaxy for computing the half radii

args = parser.parse_args()

for cloud in clouds_in_each_galaxy[args.galinfoline-1:args.galinfoline]:
        
    clouds_in_this_galaxy = cloud[1]
    gal_id = cloud[0]
    boxsize = obj.simulation.boxsize.in_units("kpc").value
    coords_gas = yt_data["PartType0","Coordinates"][clouds_in_this_galaxy].in_units("kpc").value    #gas particle coordinates (in kpc) for a given galaxy
    vel_gas = yt_data["PartType0","velocity"][clouds_in_this_galaxy].in_units("km/s").value    #gas particle velocity (in km/s) for a given galaxy
    mass_gas = yt_data["PartType0","Masses"][clouds_in_this_galaxy].in_units("Msun").value    #gas particle masses (in Msun) for a given galaxy
    gal_gas_mass = np.sum(mass_gas)    #galaxy gas mass
    
    #the following lines of code are used to ensure that box periodicity is taken into account
    
    if max(coords_gas[:,0])-min(coords_gas[:,0])>boxsize/2:
        coords_gas[:,0][np.where(coords_gas[:,0]<boxsize/2)[0]]+=boxsize
    if max(coords_gas[:,1])-min(coords_gas[:,1])>boxsize/2:
        coords_gas[:,1][np.where(coords_gas[:,1]<boxsize/2)[0]]+=boxsize
    if max(coords_gas[:,2])-min(coords_gas[:,2])>boxsize/2:
        coords_gas[:,2][np.where(coords_gas[:,2]<boxsize/2)[0]]+=boxsize
    
    #This uses the caesar COM position (which is the COM position of all baryons in the galaxy) for radii calculations
    
    pos=obj.galaxies[gal_id].pos.in_units("kpc").value
    gal_list = np.where(df["Galaxy_ID"]==gal_id)[0]
    x_df = np.array(df["x"][gal_list])
    y_df = np.array(df["y"][gal_list])
    z_df = np.array(df["z"][gal_list])
    
    #the following lines of code are used to ensure that box periodicity is taken into account
    
    if max(x_df)-min(x_df)>boxsize/2:
        x_df[np.where(x_df<boxsize/2)[0]]+=boxsize
        if pos[0]<boxsize/2:
            pos[0]+=boxsize
    if max(y_df)-min(y_df)>boxsize/2:
        y_df[np.where(y_df<boxsize/2)[0]]+=boxsize
        if pos[1]<boxsize/2:
            pos[1]+=boxsize
    if max(z_df)-min(z_df)>boxsize/2:
        z_df[np.where(z_df<boxsize/2)[0]]+=boxsize
        if pos[2]<boxsize/2:
            pos[2]+=boxsize
    
    #computing the angular momenta of each gas particle using L vector = r vector cross p vector (p vector = m * v vector)
    
    angmom_gas = np.zeros((len(clouds_in_this_galaxy),3))
    for i in range(len(clouds_in_this_galaxy)):
        angmom_gas[i][0]+=(coords_gas[i][1]-pos[1])*mass_gas[i]*(vel_gas[i][2])-(coords_gas[i][2]-pos[2])*mass_gas[i]*(vel_gas[i][1])
        angmom_gas[i][1]+=(coords_gas[i][2]-pos[2])*mass_gas[i]*(vel_gas[i][0])-(coords_gas[i][0]-pos[0])*mass_gas[i]*(vel_gas[i][2])
        angmom_gas[i][2]+=(coords_gas[i][0]-pos[0])*mass_gas[i]*(vel_gas[i][1])-(coords_gas[i][1]-pos[1])*mass_gas[i]*(vel_gas[i][0])
    
    #computing the resultant angular momentum vector based on the gas bound to the galaxy. This is used to project the galaxy accordingly
    
    L_x = np.sum(angmom_gas[:,0])
    L_y = np.sum(angmom_gas[:,1])
    L_z = np.sum(angmom_gas[:,2])
    vector_L = np.array([L_x,L_y,L_z])
    norm_L = np.linalg.norm(vector_L)
    L_yz_norm = np.sqrt(L_y**2+L_z**2)
    rot_matrix = np.array([[L_yz_norm/norm_L,-1*L_x*L_y/(L_yz_norm*norm_L),-1*L_x*L_z/(L_yz_norm*norm_L)],[L_x/norm_L,L_y/norm_L,L_z/norm_L],[0,-1*L_z/L_yz_norm,L_y/L_yz_norm]])    #this rotation matrix is computed by hand based on the matrix required for rotating a galaxy in alpha, beta directions so as to align it with the resultant angular momentum of the galaxy. Here, I fix the rotation matrix such that the new y-axis aligns with the resultant angular momenta, and therefore, the galaxy is projected along the new y axis.

    clean_list = np.where((df["conv_flag"]==16) & (df["Galaxy_ID"]==gal_id))[0]
    clean_pct = float(len(clean_list)/len(obj.galaxies[gal_id].glist)*100)    #percentage of fully converged clouds in the galaxy
    
    #now considering only the converged clouds for half-radii calculations
    
    x_df = np.array(df["x"][clean_list])
    y_df = np.array(df["y"][clean_list])
    z_df = np.array(df["z"][clean_list])
    
    #the following lines of code are used to ensure that box periodicity is taken into account
    
    if max(x_df)-min(x_df)>boxsize/2:
        x_df[np.where(x_df<boxsize/2)[0]]+=boxsize
    if max(y_df)-min(y_df)>boxsize/2:
        y_df[np.where(y_df<boxsize/2)[0]]+=boxsize
    if max(z_df)-min(z_df)>boxsize/2:
        z_df[np.where(z_df<boxsize/2)[0]]+=boxsize
    
    #new x,y,z coordinates found using the rotation matrix
    
    new_x,new_y,new_z = np.array([]),np.array([]),np.array([])

    for ind in range(0,len(x_df)):
    
        gas_vector = np.array([x_df[ind],y_df[ind],z_df[ind]])
        temp_x,temp_y,temp_z = np.dot(rot_matrix,gas_vector)
        
        new_x = np.append(new_x,temp_x)
        new_y = np.append(new_y,temp_y)
        new_z = np.append(new_z,temp_z)
        
    new_pos = np.dot(rot_matrix,pos)    #also rotating the COM position of the galaxy using the rotation matrix
    
    coords = np.array([new_x,new_y,new_z]).transpose()
        
    radius_array = []
    
    for i in range(len(clean_list)):
    
        radius_array.append(radii(coords[i],new_pos,1))    #finding the radii to each gas particle in the y axis-projected image
        
    df_aux = pd.DataFrame({"Radius":radius_array,"CO10":df["CO10"][clean_list],"CO21":df["CO21"][clean_list],"CO32":df["CO32"][clean_list],"CO43":df["CO43"][clean_list],"CO54":df["CO54"][clean_list],"CO65":df["CO65"][clean_list],"CO76":df["CO76"][clean_list],"CO87":df["CO87"][clean_list],"CO98":df["CO98"][clean_list],"CI10":df["CI10"][clean_list],"CI21":df["CI21"][clean_list],"H2_lcii":df["H2_lcii"][clean_list],"Mol_gas":df["Mol_gas"][clean_list],"Mcloud":df["Mcloud"][clean_list],"SFR":df["c_SFR"][clean_list],"Metallicity":df["Metallicity"][clean_list],"CO10_areal_TB":df["CO10_areal_TB"][clean_list]})
    df_aux = df_aux.sort_values(by='Radius')    #making a dataframe of the cloud properties index by their radial distance from the center
    
    #computing total luminosities by using a cumulative summation    
    
    total_CO10 = np.cumsum(np.array(df_aux["CO10"]))
    total_CO21 = np.cumsum(np.array(df_aux["CO21"]))
    total_CO32 = np.cumsum(np.array(df_aux["CO32"]))
    total_CO43 = np.cumsum(np.array(df_aux["CO43"]))
    total_CO54 = np.cumsum(np.array(df_aux["CO54"]))
    total_CO65 = np.cumsum(np.array(df_aux["CO65"]))
    total_CO76 = np.cumsum(np.array(df_aux["CO76"]))
    total_CO87 = np.cumsum(np.array(df_aux["CO87"]))
    total_CO98 = np.cumsum(np.array(df_aux["CO98"]))
    total_CI10 = np.cumsum(np.array(df_aux["CI10"]))
    total_CI21 = np.cumsum(np.array(df_aux["CI21"]))
    total_CII = np.cumsum(np.array(df_aux["H2_lcii"]))
    total_Mol = np.cumsum(np.array(df_aux["Mol_gas"]))
    total_mass = np.cumsum(np.array(df_aux["Mcloud"]))
    total_SFR = np.cumsum(np.array(df_aux["SFR"]))
    total_Metal = np.cumsum(np.array(df_aux["Mcloud"])*np.array(df_aux["Metallicity"]))
    total_CO10_areal_TB = np.cumsum(np.array(df_aux["CO10_areal_TB"]))
    
    #computing sleds
    
    sled_CO10_array.append(1)
    sled_CO21_array.append(total_CO21[-1]/total_CO10[-1]*115.27/230.54)
    sled_CO32_array.append(total_CO32[-1]/total_CO10[-1]*115.27/345.80)
    sled_CO43_array.append(total_CO43[-1]/total_CO10[-1]*115.27/461.04)
    sled_CO54_array.append(total_CO54[-1]/total_CO10[-1]*115.27/576.27)
    sled_CO65_array.append(total_CO65[-1]/total_CO10[-1]*115.27/691.47)
    sled_CO76_array.append(total_CO76[-1]/total_CO10[-1]*115.27/806.65)
    sled_CO87_array.append(total_CO87[-1]/total_CO10[-1]*115.27/921.80)
    sled_CO98_array.append(total_CO98[-1]/total_CO10[-1]*115.27/1036.9)
    
    #computing half radii (in principle, the variable "xx" can be modified to get r0.9, r0.99, and similar quantities)
    
    xx = 0.5
    
    half_CO10_radius = np.array(df_aux["Radius"])[np.where(total_CO10>=float(total_CO10[-1]*xx))[0][0]]
    half_CO21_radius = np.array(df_aux["Radius"])[np.where(total_CO21>=float(total_CO21[-1]*xx))[0][0]]
    half_CO32_radius = np.array(df_aux["Radius"])[np.where(total_CO32>=float(total_CO32[-1]*xx))[0][0]]
    half_CO43_radius = np.array(df_aux["Radius"])[np.where(total_CO43>=float(total_CO43[-1]*xx))[0][0]]
    half_CO54_radius = np.array(df_aux["Radius"])[np.where(total_CO54>=float(total_CO54[-1]*xx))[0][0]]
    half_CO65_radius = np.array(df_aux["Radius"])[np.where(total_CO65>=float(total_CO65[-1]*xx))[0][0]]
    half_CO76_radius = np.array(df_aux["Radius"])[np.where(total_CO76>=float(total_CO76[-1]*xx))[0][0]]
    half_CO87_radius = np.array(df_aux["Radius"])[np.where(total_CO87>=float(total_CO87[-1]*xx))[0][0]]
    half_CO98_radius = np.array(df_aux["Radius"])[np.where(total_CO98>=float(total_CO98[-1]*xx))[0][0]]
    half_CI10_radius = np.array(df_aux["Radius"])[np.where(total_CI10>=float(total_CI10[-1]*xx))[0][0]]
    half_CI21_radius = np.array(df_aux["Radius"])[np.where(total_CI21>=float(total_CI21[-1]*xx))[0][0]]
    half_CII_radius = np.array(df_aux["Radius"])[np.where(total_CII>=float(total_CII[-1]*xx))[0][0]]
    half_Mol_radius = np.array(df_aux["Radius"])[np.where(total_Mol>=float(total_Mol[-1]*xx))[0][0]]
    half_mass_radius = np.array(df_aux["Radius"])[np.where(total_mass>=float(total_mass[-1]*xx))[0][0]]
    half_SFR_radius = np.array(df_aux["Radius"])[np.where(total_SFR>=float(total_SFR[-1]*xx))[0][0]]
    half_Metal_mass_radius = np.array(df_aux["Radius"])[np.where(total_Metal>=float(total_Metal[-1]*xx))[0][0]]
    
    half_CO10_radius_array.append(half_CO10_radius)
    half_CO21_radius_array.append(half_CO21_radius)
    half_CO32_radius_array.append(half_CO32_radius)
    half_CO43_radius_array.append(half_CO43_radius)
    half_CO54_radius_array.append(half_CO54_radius)
    half_CO65_radius_array.append(half_CO65_radius)
    half_CO76_radius_array.append(half_CO76_radius)
    half_CO87_radius_array.append(half_CO87_radius)
    half_CO98_radius_array.append(half_CO98_radius)
    half_CI10_radius_array.append(half_CI10_radius)
    half_CI21_radius_array.append(half_CI21_radius)
    half_CII_radius_array.append(half_CII_radius)
    half_Mol_gas_mass_radius_array.append(half_Mol_radius)
    half_total_gas_mass_radius_array.append(half_mass_radius) #obj.galaxies[gal_id].radii["gas_half_mass"].in_units("kpc").value
    half_SFR_radius_array.append(half_SFR_radius)
    half_Metal_mass_radius_array.append(half_Metal_mass_radius)
    
    #copmuting quarter radii
    
    xx = 0.25
    
    r25_CO10_radius = np.array(df_aux["Radius"])[np.where(total_CO10>=float(total_CO10[-1]*xx))[0][0]]
    r25_CO21_radius = np.array(df_aux["Radius"])[np.where(total_CO21>=float(total_CO21[-1]*xx))[0][0]]
    r25_CO32_radius = np.array(df_aux["Radius"])[np.where(total_CO32>=float(total_CO32[-1]*xx))[0][0]]
    r25_CO43_radius = np.array(df_aux["Radius"])[np.where(total_CO43>=float(total_CO43[-1]*xx))[0][0]]
    r25_CO54_radius = np.array(df_aux["Radius"])[np.where(total_CO54>=float(total_CO54[-1]*xx))[0][0]]
    r25_CO65_radius = np.array(df_aux["Radius"])[np.where(total_CO65>=float(total_CO65[-1]*xx))[0][0]]
    r25_CO76_radius = np.array(df_aux["Radius"])[np.where(total_CO76>=float(total_CO76[-1]*xx))[0][0]]
    r25_CO87_radius = np.array(df_aux["Radius"])[np.where(total_CO87>=float(total_CO87[-1]*xx))[0][0]]
    r25_CO98_radius = np.array(df_aux["Radius"])[np.where(total_CO98>=float(total_CO98[-1]*xx))[0][0]]
    r25_CI10_radius = np.array(df_aux["Radius"])[np.where(total_CI10>=float(total_CI10[-1]*xx))[0][0]]
    r25_CI21_radius = np.array(df_aux["Radius"])[np.where(total_CI21>=float(total_CI21[-1]*xx))[0][0]]
    r25_CII_radius = np.array(df_aux["Radius"])[np.where(total_CII>=float(total_CII[-1]*xx))[0][0]]
    r25_Mol_radius = np.array(df_aux["Radius"])[np.where(total_Mol>=float(total_Mol[-1]*xx))[0][0]]
    r25_mass_radius = np.array(df_aux["Radius"])[np.where(total_mass>=float(total_mass[-1]*xx))[0][0]]
    r25_SFR_radius = np.array(df_aux["Radius"])[np.where(total_SFR>=float(total_SFR[-1]*xx))[0][0]]
    r25_Metal_mass_radius = np.array(df_aux["Radius"])[np.where(total_Metal>=float(total_Metal[-1]*xx))[0][0]]
    
    r25_CO10_radius_array.append(r25_CO10_radius)
    r25_CO21_radius_array.append(r25_CO21_radius)
    r25_CO32_radius_array.append(r25_CO32_radius)
    r25_CO43_radius_array.append(r25_CO43_radius)
    r25_CO54_radius_array.append(r25_CO54_radius)
    r25_CO65_radius_array.append(r25_CO65_radius)
    r25_CO76_radius_array.append(r25_CO76_radius)
    r25_CO87_radius_array.append(r25_CO87_radius)
    r25_CO98_radius_array.append(r25_CO98_radius)
    r25_CI10_radius_array.append(r25_CI10_radius)
    r25_CI21_radius_array.append(r25_CI21_radius)
    r25_CII_radius_array.append(r25_CII_radius)
    r25_Mol_gas_mass_radius_array.append(r25_Mol_radius)
    r25_total_gas_mass_radius_array.append(r25_mass_radius) #obj.galaxies[gal_id].radii["gas_r25_mass"].in_units("kpc").value
    r25_SFR_radius_array.append(r25_SFR_radius)
    r25_Metal_mass_radius_array.append(r25_Metal_mass_radius)
    
    #computing three-fourths radius
    
    xx = 0.75
    
    r75_CO10_radius = np.array(df_aux["Radius"])[np.where(total_CO10>=float(total_CO10[-1]*xx))[0][0]]
    r75_CO21_radius = np.array(df_aux["Radius"])[np.where(total_CO21>=float(total_CO21[-1]*xx))[0][0]]
    r75_CO32_radius = np.array(df_aux["Radius"])[np.where(total_CO32>=float(total_CO32[-1]*xx))[0][0]]
    r75_CO43_radius = np.array(df_aux["Radius"])[np.where(total_CO43>=float(total_CO43[-1]*xx))[0][0]]
    r75_CO54_radius = np.array(df_aux["Radius"])[np.where(total_CO54>=float(total_CO54[-1]*xx))[0][0]]
    r75_CO65_radius = np.array(df_aux["Radius"])[np.where(total_CO65>=float(total_CO65[-1]*xx))[0][0]]
    r75_CO76_radius = np.array(df_aux["Radius"])[np.where(total_CO76>=float(total_CO76[-1]*xx))[0][0]]
    r75_CO87_radius = np.array(df_aux["Radius"])[np.where(total_CO87>=float(total_CO87[-1]*xx))[0][0]]
    r75_CO98_radius = np.array(df_aux["Radius"])[np.where(total_CO98>=float(total_CO98[-1]*xx))[0][0]]
    r75_CI10_radius = np.array(df_aux["Radius"])[np.where(total_CI10>=float(total_CI10[-1]*xx))[0][0]]
    r75_CI21_radius = np.array(df_aux["Radius"])[np.where(total_CI21>=float(total_CI21[-1]*xx))[0][0]]
    r75_CII_radius = np.array(df_aux["Radius"])[np.where(total_CII>=float(total_CII[-1]*xx))[0][0]]
    r75_Mol_radius = np.array(df_aux["Radius"])[np.where(total_Mol>=float(total_Mol[-1]*xx))[0][0]]
    r75_mass_radius = np.array(df_aux["Radius"])[np.where(total_mass>=float(total_mass[-1]*xx))[0][0]]
    r75_SFR_radius = np.array(df_aux["Radius"])[np.where(total_SFR>=float(total_SFR[-1]*xx))[0][0]]
    r75_Metal_mass_radius = np.array(df_aux["Radius"])[np.where(total_Metal>=float(total_Metal[-1]*xx))[0][0]]
    
    r75_CO10_radius_array.append(r75_CO10_radius)
    r75_CO21_radius_array.append(r75_CO21_radius)
    r75_CO32_radius_array.append(r75_CO32_radius)
    r75_CO43_radius_array.append(r75_CO43_radius)
    r75_CO54_radius_array.append(r75_CO54_radius)
    r75_CO65_radius_array.append(r75_CO65_radius)
    r75_CO76_radius_array.append(r75_CO76_radius)
    r75_CO87_radius_array.append(r75_CO87_radius)
    r75_CO98_radius_array.append(r75_CO98_radius)
    r75_CI10_radius_array.append(r75_CI10_radius)
    r75_CI21_radius_array.append(r75_CI21_radius)
    r75_CII_radius_array.append(r75_CII_radius)
    r75_Mol_gas_mass_radius_array.append(r75_Mol_radius)
    r75_total_gas_mass_radius_array.append(r75_mass_radius) #obj.galaxies[gal_id].radii["gas_r75_mass"].in_units("kpc").value
    r75_SFR_radius_array.append(r75_SFR_radius)
    r75_Metal_mass_radius_array.append(r75_Metal_mass_radius)
    
    clean_pct_array.append(clean_pct)
    total_Mol_array.append(total_Mol[-1])
    total_Metal_array.append(total_Metal[-1])
    alpha_CO_array.append(total_Mol[-1]/total_CO10_areal_TB[-1])
    
#appending a specific galaxy's half-radii and sleds to the datafiles initialised earlier

radii_df = pd.DataFrame({"Galaxy_ID":cloud[0],"g_Metallicity":cloud[5],"g_SFR":cloud[7],"g_stellar_mass":cloud[6],"g_gas_mass":cloud[2],"g_dust_mass":cloud[8],"g_Mol_gas_mass":total_Mol_array,"g_Metal_mass":total_Metal_array,"alpha_CO":alpha_CO_array,"clean_pct":clean_pct_array,"half_CO10_radius":half_CO10_radius_array,"half_CO21_radius":half_CO21_radius_array,"half_CO32_radius":half_CO32_radius_array,"half_CO43_radius":half_CO43_radius_array,"half_CO54_radius":half_CO54_radius_array,"half_CO65_radius":half_CO65_radius_array,"half_CO76_radius":half_CO76_radius_array,"half_CO87_radius":half_CO87_radius_array,"half_CO98_radius":half_CO98_radius_array,"half_CI10_radius":half_CI10_radius_array,"half_CI21_radius":half_CI21_radius_array,"half_CII_radius":half_CII_radius_array,"half_Mol_gas_mass_radius":half_Mol_gas_mass_radius_array,"half_total_gas_mass_radius":half_total_gas_mass_radius_array,"half_SFR_radius":half_SFR_radius_array,"half_Metal_mass_radius":half_Metal_mass_radius_array,"r25_CO10_radius":r25_CO10_radius_array,"r25_CO21_radius":r25_CO21_radius_array,"r25_CO32_radius":r25_CO32_radius_array,"r25_CO43_radius":r25_CO43_radius_array,"r25_CO54_radius":r25_CO54_radius_array,"r25_CO65_radius":r25_CO65_radius_array,"r25_CO76_radius":r25_CO76_radius_array,"r25_CO87_radius":r25_CO87_radius_array,"r25_CO98_radius":r25_CO98_radius_array,"r25_CI10_radius":r25_CI10_radius_array,"r25_CI21_radius":r25_CI21_radius_array,"r25_CII_radius":r25_CII_radius_array,"r25_Mol_gas_mass_radius":r25_Mol_gas_mass_radius_array,"r25_total_gas_mass_radius":r25_total_gas_mass_radius_array,"r25_SFR_radius":r25_SFR_radius_array,"r25_Metal_mass_radius":r25_Metal_mass_radius_array,"r75_CO10_radius":r75_CO10_radius_array,"r75_CO21_radius":r75_CO21_radius_array,"r75_CO32_radius":r75_CO32_radius_array,"r75_CO43_radius":r75_CO43_radius_array,"r75_CO54_radius":r75_CO54_radius_array,"r75_CO65_radius":r75_CO65_radius_array,"r75_CO76_radius":r75_CO76_radius_array,"r75_CO87_radius":r75_CO87_radius_array,"r75_CO98_radius":r75_CO98_radius_array,"r75_CI10_radius":r75_CI10_radius_array,"r75_CI21_radius":r75_CI21_radius_array,"r75_CII_radius":r75_CII_radius_array,"r75_Mol_gas_mass_radius":r75_Mol_gas_mass_radius_array,"r75_total_gas_mass_radius":r75_total_gas_mass_radius_array,"r75_SFR_radius":r75_SFR_radius_array,"r75_Metal_mass_radius":r75_Metal_mass_radius_array})
radii_df.to_csv("results/snapshot_2_radii_projected.csv",index=False,mode='a',header=False)

sled_df = pd.DataFrame({"Galaxy_ID":cloud[0],"sled_CO10":sled_CO10_array,"sled_CO21":sled_CO21_array,"sled_CO32":sled_CO32_array,"sled_CO43":sled_CO43_array,"sled_CO54":sled_CO54_array,"sled_CO65":sled_CO65_array,"sled_CO76":sled_CO76_array,"sled_CO87":sled_CO87_array,"sled_CO98":sled_CO98_array,})
sled_df.to_csv("results/sled_snapshot_2.csv",index=False,mode='a',header=False)

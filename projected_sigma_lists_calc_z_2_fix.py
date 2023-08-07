import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import caesar
import yt
from scipy.spatial import KDTree
import argparse
import scipy.constants as physcons

kB = physcons.k*1e7
mH = physcons.m_p*1e3
G = physcons.G*1e3

#RUN THESE TWO COMMENTED LINES OF CODE IN A PYTHON TERMINAL BEFORE EXECUTING THE ENTIRE THING!! This ensures that the output csv file is created first. All the parallelized codes then append the output to it.
#sigma_df = pd.DataFrame({"Galaxy_ID":[],"CO10_sigma":[],"CO21_sigma":[],"CO32_sigma":[],"CO43_sigma":[],"CO54_sigma":[],"CI10_sigma":[],"CI21_sigma":[],"CII_sigma":[],"Mol_gas_sigma":[],"HI_gas_sigma":[],"Metal_sigma":[],"gas_sigma":[],"sfr_sigma":[],"gas_temp":[],"dust_temp":[],"n_dens":[],"col_dens":[],"RadField":[],"gal_metal":[],"gal_stellar_mass":[],"dist":[],"half_radius":[],"x":[],"z":[],"x_pos":[],"z_pos":[]})
#sigma_df.to_csv("projected_sigma_lists_z_2_parr_rot.csv",index=False)

df = pd.read_csv("/blue/narayanan/a.ravishankar/slick/repo/results/lim_df_z_2_MS_coords.csv")
selected_ids_df = pd.read_csv("/blue/narayanan/a.ravishankar/slick/repo/snap_z_2_MS.csv")
selected_ids_array = np.array(selected_ids_df["GroupID"])

obj = caesar.load("/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0160_z2.000.hdf5")
yt_snap = yt.load("/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_160.hdf5")
yt_data = yt_snap.all_data()

clouds_in_each_galaxy = [(gal.GroupID,gal.glist,gal.masses['gas'].in_units('Msun').value,gal.masses['H2'].in_units('Msun').value,gal.radii['baryon_half_mass'].in_units('pc').value,gal.metallicities['mass_weighted'].value,gal.masses['stellar'].in_units('Msun').value) for gal in [obj.galaxies[int(GroupID)] for GroupID in selected_ids_array]]

dist_list = np.array([])
half_radius_list = np.array([])
x_coord_list = np.array([])
z_coord_list = np.array([])
x_pos_list = np.array([])
z_pos_list = np.array([])
CO10_sigma_list = np.array([])
CO21_sigma_list = np.array([])
CO32_sigma_list = np.array([])
CO43_sigma_list = np.array([])
CO54_sigma_list = np.array([])
CI10_sigma_list = np.array([])
CI21_sigma_list = np.array([])
CII_sigma_list = np.array([])
Mol_gas_sigma_list = np.array([])
HI_gas_sigma_list = np.array([])
Metal_sigma_list = np.array([])
sfr_sigma_list = np.array([])
gas_sigma_list = np.array([])
gas_temp_list = np.array([])
dust_temp_list = np.array([])
n_dens_list = np.array([])
col_dens_list = np.array([])
chi_list = np.array([])
sigmaNT_list = np.array([])
gal_metal_list = np.array([])
stellar_mass_list = np.array([])

parser = argparse.ArgumentParser(prog='surfdens')
parser.add_argument("--galinfoline", type=int)

args = parser.parse_args()

for cloud in clouds_in_each_galaxy[args.galinfoline-1:args.galinfoline]:
    
    clouds_in_this_galaxy = cloud[1]
    gal_id = cloud[0]
    boxsize = obj.simulation.boxsize.in_units("kpc").value
    coords_gas = yt_data["PartType0","Coordinates"][clouds_in_this_galaxy].in_units("kpc").value
    vel_gas = yt_data["PartType0","velocity"][clouds_in_this_galaxy].in_units("km/s").value
    mass_gas = yt_data["PartType0","Masses"][clouds_in_this_galaxy].in_units("Msun").value
    gal_gas_mass = np.sum(mass_gas)
    #angmom_gas = yt_data["PartType0","angular_momentum"][clouds_in_this_galaxy]
    
    if max(coords_gas[:,0])-min(coords_gas[:,0])>boxsize/2:
        coords_gas[:,0][np.where(coords_gas[:,0]<boxsize/2)[0]]+=boxsize
    if max(coords_gas[:,1])-min(coords_gas[:,1])>boxsize/2:
        coords_gas[:,1][np.where(coords_gas[:,1]<boxsize/2)[0]]+=boxsize
    if max(coords_gas[:,2])-min(coords_gas[:,2])>boxsize/2:
        coords_gas[:,2][np.where(coords_gas[:,2]<boxsize/2)[0]]+=boxsize
    
    pos = np.array([np.sum(coords_gas[:,0]*mass_gas)/gal_gas_mass,np.sum(coords_gas[:,1]*mass_gas)/gal_gas_mass,np.sum(coords_gas[:,2]*mass_gas)/gal_gas_mass])
    #vel = np.array([np.sum(vel_gas[:,0]*mass_gas)/gal_gas_mass,np.sum(vel_gas[:,1]*mass_gas)/gal_gas_mass,np.sum(vel_gas[:,2]*mass_gas)/gal_gas_mass])
    
    angmom_gas = np.zeros((len(clouds_in_this_galaxy),3))
    for i in range(len(clouds_in_this_galaxy)):
        angmom_gas[i][0]+=(coords_gas[i][1]-pos[1])*mass_gas[i]*(vel_gas[i][2])-(coords_gas[i][2]-pos[2])*mass_gas[i]*(vel_gas[i][1])
        angmom_gas[i][1]+=(coords_gas[i][2]-pos[2])*mass_gas[i]*(vel_gas[i][0])-(coords_gas[i][0]-pos[0])*mass_gas[i]*(vel_gas[i][2])
        angmom_gas[i][2]+=(coords_gas[i][0]-pos[0])*mass_gas[i]*(vel_gas[i][1])-(coords_gas[i][1]-pos[1])*mass_gas[i]*(vel_gas[i][0])
    
    L_x = np.sum(angmom_gas[:,0])
    L_y = np.sum(angmom_gas[:,1])
    L_z = np.sum(angmom_gas[:,2])
    vector_L = np.array([L_x,L_y,L_z])
    norm_L = np.linalg.norm(vector_L)
    L_yz_norm = np.sqrt(L_y**2+L_z**2)
    rot_matrix = np.array([[L_yz_norm/norm_L,-1*L_x*L_y/(L_yz_norm*norm_L),-1*L_x*L_z/(L_yz_norm*norm_L)],[L_x/norm_L,L_y/norm_L,L_z/norm_L],[0,-1*L_z/L_yz_norm,L_y/L_yz_norm]])
    
    gal_list = np.where(df["Galaxy_ID"]==gal_id)[0]
    x_df = np.array(df["x"][gal_list])
    y_df = np.array(df["y"][gal_list])
    z_df = np.array(df["z"][gal_list])
    
    """
    This uses the caesar COM position (which is the COM position of all baryons in the galaxy) for radii calculations
    
    pos=obj.galaxies[gal_id].pos.in_units("kpc").value
    
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
    """
    
    #The following uses the COM position of JUST the gas particles in the galaxy
    
    if max(x_df)-min(x_df)>boxsize/2:
        x_df[np.where(x_df<boxsize/2)[0]]+=boxsize
    if max(y_df)-min(y_df)>boxsize/2:
        y_df[np.where(y_df<boxsize/2)[0]]+=boxsize
    if max(z_df)-min(z_df)>boxsize/2:
        z_df[np.where(z_df<boxsize/2)[0]]+=boxsize
    
    #pos = np.array([np.sum(x_df*df["Mcloud"][gal_list])/np.sum(df["Mcloud"][gal_list]),np.sum(y_df*df["Mcloud"][gal_list])/np.sum(df["Mcloud"][gal_list]),np.sum(z_df*df["Mcloud"][gal_list])/np.sum(df["Mcloud"][gal_list])])
    
    new_x,new_y,new_z = np.array([]),np.array([]),np.array([])

    for ind in range(0,len(x_df)):
    
        gas_vector = np.array([x_df[ind],y_df[ind],z_df[ind]])
        temp_x,temp_y,temp_z = np.dot(rot_matrix,gas_vector)
        
        new_x = np.append(new_x,temp_x)
        new_y = np.append(new_y,temp_y)
        new_z = np.append(new_z,temp_z)
        
    new_pos = np.dot(rot_matrix,pos)
        
    x_axis = np.arange(min(new_x),max(new_x)+1,0.5)
    z_axis = np.arange(min(new_z),max(new_z)+1,0.5)

    for i in range(0,len(x_axis)-1):
    
        for j in range(0,len(z_axis)-1):
        
            unclean_list = np.where((new_x>=x_axis[i]) & (new_x<x_axis[i+1]) & (new_z>=z_axis[j]) & (new_z<z_axis[j+1]))[0]
            
            CO10_sigma,CO21_sigma,CO32_sigma,CO43_sigma,CO54_sigma,CI10_sigma,CI21_sigma,CII_sigma,Mol_gas_sigma,HI_gas_sigma,gas_temp,dust_temp = np.zeros(12)
            
            clean_list = np.where((np.array(df["conv_flag"][gal_list])==16) & (new_x>=x_axis[i]) & (new_x<x_axis[i+1]) & (new_z>=z_axis[j]) & (new_z<z_axis[j+1]))[0]
            
            if len(clean_list)>0:
            
                CO10_sigma = np.sum(np.array(df["CO10_areal_TB"][gal_list])[clean_list])/0.25
                CO21_sigma = np.sum(np.array(df["CO21"][gal_list])[clean_list])*1e11/(3.826e33*0.25*3*(230.538**3))
                CO32_sigma = np.sum(np.array(df["CO32"][gal_list])[clean_list])*1e11/(3.826e33*0.25*3*(345.7959899**3))
                CO43_sigma = np.sum(np.array(df["CO43"][gal_list])[clean_list])*1e11/(3.826e33*0.25*3*(461.0407682**3))
                CO54_sigma = np.sum(np.array(df["CO54"][gal_list])[clean_list])*1e11/(3.826e33*0.25*3*(576.2679305**3))
                CI10_sigma = np.sum(np.array(df["CI10"][gal_list])[clean_list])*1e11/(3.826e33*0.25*3*(492.160651**3))
                CI21_sigma = np.sum(np.array(df["CI21"][gal_list])[clean_list])*1e11/(3.826e33*0.25*3*(809.34197**3))
                CII_sigma = np.sum(np.array(df["H2_lcii"][gal_list])[clean_list])*1e11/(3.826e33*0.25*3*(1900.5369**3))
                Mol_gas_sigma = np.sum(np.array(df["Mol_gas"][gal_list])[clean_list])/0.25
                HI_gas_sigma = np.sum(np.array(df["fH"][gal_list])[clean_list]*np.array(df["Mcloud"][gal_list])[clean_list])/0.25
                gas_temp = np.sum(np.array(df["gas_temp"][gal_list])[clean_list]*np.array(df["Mcloud"][gal_list])[clean_list])/np.sum(np.array(df["Mcloud"][gal_list])[clean_list])
                dust_temp = np.sum(np.array(df["dust_temp"][gal_list])[clean_list]*np.array(df["Mcloud"][gal_list])[clean_list])/np.sum(np.array(df["Mcloud"][gal_list])[clean_list])
                
            if len(unclean_list)>0:
                
                dist = np.sqrt(((x_axis[i]+x_axis[i+1])/2-new_pos[0])**2+((z_axis[j]+z_axis[j+1])/2-new_pos[2])**2)
                dist_list = np.append(dist_list,dist)
                x_coord_list = np.append(x_coord_list,(x_axis[i]+x_axis[i+1])/2)
                z_coord_list = np.append(z_coord_list,(z_axis[j]+z_axis[j+1])/2)
                x_pos_list = np.append(x_pos_list,new_pos[0])
                z_pos_list = np.append(z_pos_list,new_pos[2])
                half_radius_list = np.append(half_radius_list,obj.galaxies[gal_id].radii["gas_half_mass"].in_units("kpc"))
                CO10_sigma_list = np.append(CO10_sigma_list,CO10_sigma)
                CO21_sigma_list = np.append(CO21_sigma_list,CO21_sigma)
                CO32_sigma_list = np.append(CO32_sigma_list,CO32_sigma)
                CO43_sigma_list = np.append(CO43_sigma_list,CO43_sigma)
                CO54_sigma_list = np.append(CO54_sigma_list,CO54_sigma)
                CI10_sigma_list = np.append(CI10_sigma_list,CI10_sigma)
                CI21_sigma_list = np.append(CI21_sigma_list,CI21_sigma)
                CII_sigma_list = np.append(CII_sigma_list,CII_sigma)
                Mol_gas_sigma_list = np.append(Mol_gas_sigma_list,Mol_gas_sigma)
                HI_gas_sigma_list = np.append(HI_gas_sigma_list,HI_gas_sigma)
                Metal_sigma_list = np.append(Metal_sigma_list,np.sum(np.array(df["Metallicity"][gal_list])[unclean_list]*np.array(df["Mcloud"][gal_list])[unclean_list])*0.0196/0.25)
                sfr_sigma_list = np.append(sfr_sigma_list,np.sum(np.array(df["c_SFR"][gal_list])[unclean_list])/0.25)
                gas_sigma_list = np.append(gas_sigma_list,np.sum(np.array(df["Mcloud"][gal_list])[unclean_list])/0.25)
                gas_temp_list = np.append(gas_temp_list,gas_temp)
                dust_temp_list = np.append(dust_temp_list,dust_temp)
                n_dens = np.sum(np.array(df["n_dens"][gal_list])[unclean_list]*np.array(df["Mcloud"][gal_list])[unclean_list])/np.sum(np.array(df["Mcloud"][gal_list])[unclean_list])
                n_dens_list = np.append(n_dens_list,n_dens)
                col_dens = np.sum(np.array(df["col_dens"][gal_list])[unclean_list]*np.array(df["Mcloud"][gal_list])[unclean_list])/np.sum(np.array(df["Mcloud"][gal_list])[unclean_list])
                col_dens_list = np.append(col_dens_list,col_dens)
                chi_list = np.append(chi_list,np.mean(np.array(df["RadField"][gal_list])[unclean_list]))
                #sigmaTotSqr = (3.0*np.pi*G/20.0)*col_dens**2/n_dens*(1.4)*mH
                #sigmaThSqr = kB*gas_temp/(self.comp.mu*mH)
                #sigma_NT = np.sqrt(sigmaTotSqr - sigmaThSqr)
                gal_metal_list = np.append(gal_metal_list,cloud[5])
                stellar_mass_list = np.append(stellar_mass_list,cloud[6])
        
sigma_df = pd.DataFrame({"Galaxy_ID":[gal_id]*len(gas_sigma_list),"CO10_sigma":CO10_sigma_list,"CO21_sigma":CO21_sigma_list,"CO32_sigma":CO32_sigma_list,"CO43_sigma":CO43_sigma_list,"CO54_sigma":CO54_sigma_list,"CI10_sigma":CI10_sigma_list,"CI21_sigma":CI21_sigma_list,"CII_sigma":CII_sigma_list,"Mol_gas_sigma":Mol_gas_sigma_list,"HI_gas_sigma":HI_gas_sigma_list,"Metal_sigma":Metal_sigma_list,"gas_sigma":gas_sigma_list,"sfr_sigma":sfr_sigma_list,"gas_temp":gas_temp_list,"dust_temp":dust_temp_list,"n_dens":n_dens_list,"col_dens":col_dens_list,"RadField":chi_list,"gal_metal":gal_metal_list,"gal_stellar_mass":stellar_mass_list,"dist":dist_list,"half_radius":half_radius_list,"x":x_coord_list,"z":z_coord_list,"x_pos":x_pos_list,"z_pos":z_pos_list},index=[0]*len(gas_sigma_list))
sigma_df.to_csv("projected_sigma_lists_z_2_parr_rot.csv",index=False,mode='a',header=False)

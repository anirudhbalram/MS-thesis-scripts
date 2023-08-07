import numpy as np
from hyperion.model import ModelOutput
from astropy.cosmology import Planck13
from astropy import units as u
from astropy import constants
import yt
import caesar as cs
from scipy import interpolate
import pickle
import scipy.stats
import tqdm
import pandas as pd


# Filepaths and parameters to adjust
snapid = '078'
csfile = f'/orange/narayanan/desika.narayanan/gizmo_runs/simba/m100n1024/Groups/m100n1024_{snapid}.hdf5'
outfile = '/blue/narayanan/hollis.akins/data/UVJ_m100n1024_ap3_z2.pickle'
filters = ['U','V','J_wfcam','u_SDSS','r_SDSS','g_SDSS','F300W_WFPC2','galex2500','galex1500','H1','H2','K_uv','irac_ch4','mips24', 'W3']

pdpath1 = f'/blue/narayanan/hollis.akins/m100n1024/pd_runs_ap/snap{snapid}/'
pdpath2 = f'/blue/narayanan/hollis.akins/m100n1024/pd_runs_ap_nodust/snap{snapid}/'

timescale = 200
MS_type = 'bending'

# Set generate_id_file to True if you need to re-generate the .npz file of relevant IDs that you want to get properties for
# (it became necessary to have this stored separately so that we can ensure that the dusty and dust-free pd runs are matched)
generate_id_file = True


def get_ids(pdpath,snapid):
    import glob
    filenames = glob.glob(pdpath+'*.rtout.sed')
    prefix, suffix = f'{pdpath}snap{snapid}.galaxy', '.rtout.sed'    
    ids = np.sort(np.array([int(filename[len(prefix):len(filename)-len(suffix)]) for filename in filenames]))
    indices_to_remove = []
    for i in tqdm.tqdm(ids):
        run = pdpath + f'snap{snapid}.galaxy{i}.rtout.sed'
        m = ModelOutput(run)
        try:
            wav,flux = m.get_sed()
        except:
            #print(f'Ignoring galaxy{i}')
            indices_to_remove.append(i)
            continue
    ids = np.array([i for i in ids if i not in indices_to_remove])
    return ids

def get_pd_ids(pdpath, snapid):
    import glob
    filenames = glob.glob(pdpath+'*.rtout.sed')
    prefix, suffix = f'{pdpath}snap{snapid}.galaxy', '.rtout.sed'    
    ids = np.sort(np.array([int(filename[len(prefix):len(filename)-len(suffix)]) for filename in filenames]))
    return ids

def get_filtered_mag(flux, wav, filterpath):
    transmission_curve = np.loadtxt(filterpath)
    filt_wav = transmission_curve[:,0]*u.micron
    filt_pass = transmission_curve[:,1]
    
    
    s = interpolate.interp1d(filt_wav,filt_pass, bounds_error=False, fill_value=0.)
    pass_int = s(wav)
    
    dlamb = np.abs(wav[1:]-wav[:-1]).to(u.cm)
    filtered_flux = np.sum(flux[:-1] * pass_int[:-1] * dlamb) / np.sum(pass_int[:-1]*dlamb)
    filtered_mag = -2.5 * np.log10(filtered_flux * u.cm * u.cm / u.erg) - 48.6 # from a ppt https://www.uio.no/studier/emner/matnat/astro/AST4320/h14/timeplan/lecture1.pdf
    return filtered_mag

def get_mags_from_SEDs(filters,pdpath, snapid,ids):
    print(f'Getting UVJ magnitudes from SEDs for {len(ids)} galaxies...')
        
    distance = 10*u.pc
    distance = distance.to(u.cm)
    
    magnitudes = dict()
    for key in filters:
        magnitudes[key] = np.array([]) 
    
    fluxes, V_band_fluxes, indices_to_remove = [], [], []
    
    for i in tqdm.tqdm(ids):
        run = pdpath + f'snap{snapid}.galaxy{i}.rtout.sed'
        m = ModelOutput(run)
        try:
            wav,flux = m.get_sed(inclination='all',aperture=-1, distance=distance.value, units='ergs/cm^2/s/Hz') # what do these mean? look up in Hyperion documentaiton
        except ValueError:
            #print(f'ValueError for galaxy{i}')
            indices_to_remove.append(i)
            continue
        except KeyError: 
            #print(f'KeyError for galaxy{i}')
            indices_to_remove.append(i)
            continue

        wav  = np.asarray(wav)*u.micron #wav is in micron
        flux = np.asarray(flux[0])*u.erg/(u.cm*u.cm*u.s*u.Hz)
        
        V_band_flux = flux[np.argmin(np.abs(wav/u.micron - 0.5729700000000001))] 
        fluxes.append(flux)
        V_band_fluxes.append(V_band_flux)
        
        for filt in filters:
            mags = magnitudes[filt]
            filtpath = f'/home/hollis.akins/Code/powderday/filters/{filt}.filter'
            filtered_mag = get_filtered_mag(flux, wav, filtpath)
            magnitudes[filt] = np.append(mags, float(filtered_mag))
                
    return magnitudes, wav, fluxes, V_band_fluxes
            
def Speagle_MS_SFR(mstar,t):
    t = t.to(u.Gyr)
    slope = 0.84 - 0.026 * t / u.Gyr
    yint = -(6.51 - 0.11 * t / u.Gyr)
    logSFR = slope*np.log10(mstar/u.Msun) + yint
    return 10**logSFR * u.Msun / u.yr
    
    
id_file_path = outfile.replace('.pickle','_ids.npz')
if generate_id_file:
    print('Generating ID file (may take a while)...')
    ids1 = get_ids(pdpath1, snapid)
    ids2 = get_ids(pdpath2, snapid)
    ids = np.intersect1d(ids1,ids2)
    np.savez(id_file_path, ids)
else:
    print(f'Grabbing ID file from {id_file_path}')
    ids = np.load(id_file_path)['arr_0']
    
print(f'Found {len(ids)} galaxy ids')


print('Loading caesar file...')
csobj = cs.load(csfile)
t = float(csobj.simulation.time)*u.s
t = t.to(u.Gyr)
z = float(csobj.simulation.redshift)


print(f'Getting magnitudes from {pdpath1}')
magnitudes, wav, fluxes, V_band_fluxes = get_mags_from_SEDs(filters,pdpath1,snapid,ids)


print(f'Getting magnitudes from {pdpath2}')
magnitudes_nodust, wav_nodust, fluxes_nodust, V_band_fluxes_nodust = get_mags_from_SEDs(filters,pdpath2,snapid,ids)


assert len(fluxes) == len(fluxes_nodust)

print('Getting MS from full caesar file...')
### compute the main sequence for the z=2 snapshot
all_caesar_ids = [gal.GroupID for gal in csobj.galaxies]
mstar = np.array([gal.masses['stellar'] for gal in csobj.galaxies])
t_obs = float(t.value)

with open('/orange/narayanan/hollis.akins/data/sfh_m100n1024_z2.pickle', 'rb') as f:
    data = pickle.load(f)
    data = pd.DataFrame(data)

sfr_200 = np.array([])
for i in tqdm.tqdm(all_caesar_ids):
    massform = data['massform'][data.id==i].tolist()[0]
    tform = data['tform'][data.id==i].tolist()[0]
    
    age = t_obs - tform
    age = age*u.Gyr
    age = age.to(u.Myr)
    age = age.value
        
    w = age < timescale
    if len(w[w])==0:
        sfr_200 = np.append(sfr_200, 1e-2)
    else:
        sfr_200 = np.append(sfr_200, np.sum(massform[w])/(timescale*1e6))

sfr = sfr_200
sSFR = sfr/mstar
sfr[sfr == 1e-2] = 1e-13*mstar[sfr==1e-2] 

bins = np.arange(np.min(np.log10(mstar)), np.max(np.log10(mstar))+0.5, 0.2)
bincenters = 0.5*(bins[1:]+bins[:-1])

t_obs = (t_obs*u.Gyr).to(u.yr).value

# iteratively fit polynomial to the MS
star_forming_galaxies = np.ones(shape=mstar.shape, dtype=bool)
for i in range(5):
    median, sigma = np.array([]), np.array([])
    for left_edge, right_edge in zip(bins[:-1], bins[1:]):
        left_edge, right_edge = np.round(left_edge, 1), np.round(right_edge,1)
        s = sfr[(np.log10(mstar) > left_edge)&(np.log10(mstar)<=right_edge)&(star_forming_galaxies)]
        s = np.log10(s)
        m = np.median(s)
        sig = np.sqrt(np.sum((s-m)**2)/len(s)) # basically the std deviation but with median instead of mean
        median = np.append(median,m)
        sigma = np.append(sigma, sig)
    
    if MS_type == 'linear': # if you want a linear MS for SIMBA
        slope, yint, r, p, sig = scipy.stats.linregress(bincenters[(bincenters < 11)], median[(bincenters < 11)]) # is 11.2 a good number? check later
        ms_sfr = 10**(slope*np.log10(mstar)+yint)
    elif MS_type == 'bending': # if you want a bending MS for SIMBA
        p = np.polyfit(bincenters[~np.isnan(median)], median[~np.isnan(median)], 2)
        ms_sfr = 10**(np.log10(mstar)**2 * p[0] + np.log10(mstar)*p[1] + p[2])
    else: raise Exception('Unrecognized MS_type')
    star_forming_galaxies = (np.log10(sfr/ms_sfr) > -0.5) 

# now that we have the main sequence SFR, we restrict our result to just the POWDERDAY galaxies and continue
ms_sfr = ms_sfr[np.isin(all_caesar_ids, ids)]
sfr = sfr[np.isin(all_caesar_ids, ids)]
mstar = mstar[np.isin(all_caesar_ids, ids)]
sSFR = sfr/mstar
ms_sSFR = ms_sfr / mstar
delta_SFR = np.log10(sfr/ms_sfr)


print('Getting remaining caesar properties...')

R_half = np.array([gal.radii['stellar_half_mass'].to('kpc') for gal in csobj.galaxies if gal.GroupID in ids])
mH2 = np.array([gal.masses['H2'] for gal in csobj.galaxies if gal.GroupID in ids])
mHI = np.array([gal.masses['HI'] for gal in csobj.galaxies if gal.GroupID in ids])
mBH = np.array([gal.masses['bh'] for gal in csobj.galaxies if gal.GroupID in ids])
mdust = np.array([gal.masses['dust'] for gal in csobj.galaxies if gal.GroupID in ids])
mgas = np.array([gal.masses['gas'] for gal in csobj.galaxies if gal.GroupID in ids])
sat = ~np.array([gal.central for gal in csobj.galaxies if gal.GroupID in ids],dtype=bool)
local_mass_density = np.array([gal.local_mass_density['1000'] for gal in csobj.galaxies if gal.GroupID in ids])
Zstar = np.array([gal.metallicities['stellar'] for gal in csobj.galaxies if gal.GroupID in ids])
caesar_U = np.array([gal.absmag['sdss_u'] for gal in csobj.galaxies if gal.GroupID in ids])                      # should change these to match the filters used for pd
caesar_V = np.array([gal.absmag['v'] for gal in csobj.galaxies if gal.GroupID in ids])
caesar_J = np.array([gal.absmag['wfcam_j'] for gal in csobj.galaxies if gal.GroupID in ids])
caesar_U_nodust = np.array([gal.absmag_nodust['sdss_u'] for gal in csobj.galaxies if gal.GroupID in ids])
caesar_V_nodust = np.array([gal.absmag_nodust['v'] for gal in csobj.galaxies if gal.GroupID in ids])
caesar_J_nodust = np.array([gal.absmag_nodust['wfcam_j'] for gal in csobj.galaxies if gal.GroupID in ids])


assert len(ids) == len(mstar) == len(fluxes)


print(f'Dumping data to {outfile}')

out_dict = {
    'pd_id':ids,
    'caesar_U':caesar_U,
    'caesar_V':caesar_V,
    'caesar_J':caesar_J,
    'caesar_U_nodust':caesar_U_nodust,
    'caesar_V_nodust':caesar_V_nodust,
    'caesar_J_nodust':caesar_J_nodust,
    'sSFR':sSFR,
    'sfr':sfr,
    'mstar':mstar,
    'ms_sfr':ms_sfr,
    'ms_sSFR':ms_sSFR,
    'delta_SFR':delta_SFR,
    'wav':[wav]*len(fluxes),
    'flux':fluxes,
    'flux_nodust':fluxes_nodust, 
    'V_band_flux':V_band_fluxes,
    'V_band_flux_nodust':V_band_fluxes_nodust,
    'R_half':R_half, 
    'mH2':mH2,
    'mHI':mHI,
    'mBH':mBH,
    'mdust':mdust,
    'mgas':mgas,
    'sat':sat,
    'local_mass_density':local_mass_density,
    'Zstar':Zstar
}

for key in magnitudes.keys():
    out_dict[key] = magnitudes[key]

for key in magnitudes_nodust.keys():
    out_dict[key+'_nodust'] = magnitudes_nodust[key]

with open(outfile, 'wb') as f:
    pickle.dump(out_dict, f)

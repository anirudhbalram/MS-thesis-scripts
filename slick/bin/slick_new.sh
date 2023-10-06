cat > parameters.ini << ENDOFFILE
[sample]
mode=total

[snap]
boxsize=400
ytfilename=/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_160.hdf5
caesarfilename=/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/Groups/caesar_0160_z2.000.hdf5

[run]
output_dir=/blue/narayanan/a.ravishankar/slick/results    #insert alternate output directory here
skip_run=True
ENDOFFILE

cat > run.sh << ENDOFFILE
#!/bin/bash
#SBATCH --qos=narayanan-b
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=1:00:00

slick_frontend.sh parameters.ini
ENDOFFILE

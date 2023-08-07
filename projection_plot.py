import yt

snap = '/orange/narayanan/desika.narayanan/gizmo_runs/simba/m25n512/output/snapshot_150.hdf5'
ds  = yt.load(snap)

p = yt.ProjectionPlot(ds, "z", ("gas", "temperature"))
p.save("plot_snapshot_150_gas_temperature.jpg")

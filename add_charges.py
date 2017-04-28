#!/usr/bin/env python
from read_write import topos,plumed
from modify_structure import man,plane,com
from read_write import topos 
import numpy as np

# Converting reference file to pdb file
fn_in ='data.water'
npoints = 10
r = 5.25

data,index, molecule_tag, atom_type, charge, coordinates,velocities = topos.read_lammps(fn_in)
symbols = topos.ids_symbols(atom_type)
xhi = data['xhi']
yhi = data['yhi']
zhi = max(coordinates[:,2])
zlo = min(coordinates[:,2])
xcenter = xhi/2.0
ycenter = yhi/2.0

inds_Os = [ind+1 for ind,val in enumerate(symbols) if val == 'O']
inds_Hs = [ind+1 for ind,val in enumerate(symbols) if val == 'H']
zs = np.linspace(zlo,zhi,npoints)
points_coords = np.zeros((npoints,3))
for n in range(npoints):
    x = (r - -r)*np.random.random() - r
    y = (1 - -1)*np.random.random() - 1
    if y < 0:
        y = -1
    else:
        y = 1
    points_coords[n,0] += x + xcenter
    points_coords[n,1] += y*np.sqrt(r**2-x**2) + ycenter
    points_coords[n,2] += zs[n]

new_coords = np.vstack((coordinates,points_coords))
new_atom_type = np.hstack((atom_type,np.ones(npoints)*3))
new_molecule_tag = np.hstack((molecule_tag,np.ones(npoints)*2.0))
new_charge = np.hstack((charge,np.ones(npoints)*-1.0))
new_index = np.hstack((index,np.arange(1.0,npoints+1.0)+len(coordinates)))
natoms = len(coordinates) + npoints

data['symbols'] = ['H','O','H']
data['masses'] = [1.007940,15.999400,1.007940]
data['ids'] = [1,2,3]
data['atom_types'] = 3

topos.write_lammps(fn_out = "new_data.water",
             index = new_index,
             atom_type = new_atom_type, 
             charge = new_charge,
             coordinates = new_coords,
             molecule_tag = new_molecule_tag,
             in_opts = data,
            )


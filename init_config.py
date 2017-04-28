#!/usr/bin/env python -u 
from read_write import topos,plumed, plumed2
from modify_structure import man,plane,com
import math
import numpy as np
import mdtraj as mdt

# Constructs initial geometry by cutting from bulk system

def build_system(
                in_opts = {},
                ):

    # => Parameters for Constructing Initial Geometry from bulk trajectory <= #
    options = {    

    'cut_sphere' : True,
    'center_sphere' : True,
    'nwaters' : 300,
    
    'cut_plane' : False,
    'plane_gap' : 0.2,
    
    'expand_system' : False,
    'perp_expand' : False,
    'axis' : [1,2],
    'times' : [2,2],
    'exp_gap' : 0.1,
    
    'scale_system' : False,
    'axis_scale' : 0.8,
    'scale_axis' : [1,2],
    
    'ostar_id' : 25,
    
    'periodic' : False,
    
    # plumed parameters
    'sigma' : 4.0,
    'epsilon' : 0.500,
    'alpha' : 8.0,
    'kappa' : 250.0,
    'exp' : 2.0,
    'edge_gap' : 0.1,
    
    # Filenames 
    'fn_in' : "data.water-out",
    'fn_lammps' : "data.water",
    'fn_pdb' : "water.pdb",
    'fn_int' : "ref.pdb",

    }

    # => Override Default Options <= #

    for key,val in in_opts.items():
       if key not in options.keys():
          raise ValueError('%s option unavailable, RTFM' % (key))
       if options[key] != val:
          print "Default Option Overridden: %s = %s" % (key,str(val))
       options[key] = val
    print '\n'
    
    
    # => Declaring Variables from Options <= #
    print "Build: Declaring variables from options"

    cut_sphere     = options['cut_sphere']
    center_sphere  = options['center_sphere']
    nwaters        = options['nwaters']
    cut_plane      = options['cut_plane'] 
    plane_gap      = options['plane_gap'] 
    expand_system  = options['expand_system']
    perp_expand  = options['perp_expand']
    axis           = options['axis']
    times          = options['times']
    exp_gap        = options['exp_gap']
    scale_system   = options['scale_system']
    axis_scale     = options['axis_scale']
    scale_axis     = options['scale_axis']
    ostar_id       = options['ostar_id']
    periodic       = options['periodic']

    kappa          = options['kappa']
    exp            = options['exp']
    edge_gap       = options['edge_gap']
    epsilon        = options['epsilon']
    sigma          = options['sigma']
    alpha          = options['alpha']
    
    fn_in          = options['fn_in']
    fn_lammps      = options['fn_lammps']
    fn_pdb         = options['fn_pdb']
    fn_int         = options['fn_int']
    
    # => Loading Data <= #
    
    # Converting reference file to pdb file
    print "Build: Reading data from lammps file"
    data,index, molecule_tag, atom_type, charge, coordinates,velocities = topos.read_lammps(fn_in)
    symbols = topos.ids_symbols(atom_type)
    
    print "Build: Writing new reference pdb file"
    topos.write_pdb(
        fn_out=fn_int,
        index=index,
        symbols=symbols,
        coordinates=coordinates,
        )

    print "Build: Eliminating incomplete water molecules due to pbc conditions"
    # cutting incomplete molecules from system
    new_atoms, indices = man.cull_incomplete(fn_int)
    
    # assigning new indices and updating reference file
    new_index = np.arange(1,len(new_atoms)+1)
    new_symbols = symbols[indices]
    new_charge = charge[indices]
    new_atoms *= 10.0
    print "Build: Updating reference pdb file"
    topos.write_pdb(
        fn_out=fn_int,
        index=new_index,
        symbols=new_symbols,
        coordinates=new_atoms,
        )
    
    if cut_sphere == True and cut_plane == True:
        print "Warning: Plane cut from sphere!"
    
    # => Cutting Sphere out of Bulk <= #
    if cut_sphere == True:
        if 2.0*nwaters > len(coordinates)/3:
            cube_water = math.ceil((2.0*(3.0*nwaters / (4.0*np.pi))**(1.0/3.0))**3)
            factor = int(math.ceil((cube_water/(len(coordinates)/3))**(1.0/3.0)))
           
            print "Build: Expanding base system to accomodate new sphere size" 
            print "factor: %d" % (factor)
            
            old_natoms = len(new_atoms)
            new_atoms = man.expand(new_atoms,
                                   axis=[0,1,2],
                                   times=[factor,factor,factor],
                                   gap=exp_gap,
                                  )

            # assigning indices, charges and symbols to new atoms
            new_index = np.arange(1,len(new_atoms)+1)        
            new_charge = charge[indices]
            new_symbols = symbols[indices]
            
            add_charge = new_charge
            add_symbols = new_symbols
            natoms = len(new_atoms)
            
            for n in range(natoms/old_natoms-1):
                new_charge = np.hstack((new_charge,add_charge))
                new_symbols = np.hstack((new_symbols,add_symbols))
    
            new_index = np.arange(1,len(new_atoms)+1)
            print "Build: Updating reference pdb file"
            topos.write_pdb(fn_out=fn_int,
                               index=new_index,
                               symbols=new_symbols,
                               coordinates=new_atoms,
                               )
            print "Build: Calculating new center of mass"
            center = com.com(new_atoms,new_symbols)
            center = [val/10.0 for val in center]
        else:
            center = com.com(new_atoms,new_symbols)
            center = [val/10.0 for val in center]
        print "Build: Cutting sphere from bulk"
        new_atoms,indices = man.n_sphere(
            fn_int,
            nwaters=nwaters,
            centerX=center[0],
            centerY=center[1],
            centerZ=center[2],
            )

        new_charge = new_charge[indices]
        new_symbols = new_symbols[indices]
    
    # => Cutting system into 2D Plane <= #
    
    if cut_plane == True:
        if perp_expand == True:
            old_natoms = len(new_atoms)
            new_atoms = man.expand(new_atoms,
                                   axis=[0],
                                   times=[3],
                                   gap=exp_gap,
                                  )
            # assigning indices, charges and symbols to new atoms
            new_index = np.arange(1,len(new_atoms)+1)        
            new_charge = charge[indices]
            new_symbols = symbols[indices]
            
            add_charge = new_charge
            add_symbols = new_symbols
            natoms = len(new_atoms)
            
            for n in range(natoms/old_natoms-1):
                new_charge = np.hstack((new_charge,add_charge))
                new_symbols = np.hstack((new_symbols,add_symbols))
    
            new_index = np.arange(1,len(new_atoms)+1)
            topos.write_pdb(fn_out=fn_int,
                               index=new_index,
                               symbols=new_symbols,
                               coordinates=new_atoms,
                               )
            center = com.com(new_atoms,new_symbols)
            center = [val/10.0 for val in center]
        new_atoms,indices = plane.plane(fn_int,
                                        gap=plane_gap,
                                        ax_perp=0,
                                        )
        new_charge = new_charge[indices]
        new_symbols = new_symbols[indices]

    # => Enlarging System <= #
    
    if expand_system == True:
        old_natoms = len(new_atoms)
        new_atoms = man.expand(new_atoms,
                               axis=axis,
                               times=times,
                               gap=exp_gap,
                               )
    
        # assigning indices, charges and symbols to new atoms
        new_index = np.arange(1,len(new_atoms)+1)        
        
        add_charge = new_charge
        add_symbols = new_symbols
        natoms = len(new_atoms)
        
        for n in range(natoms/old_natoms-1):
            new_charge = np.hstack((new_charge,add_charge))
            new_symbols = np.hstack((new_symbols,add_symbols))
    
    
    # writing new reference pdb file
    new_atoms *= 10.0
    topos.write_pdb(fn_out=fn_int,
                       index=new_index,
                       symbols=new_symbols,
                       coordinates=new_atoms,
                       )
    
    # => Scaling density of system <= #
    
    if scale_system == True:
        new_atoms =  man.scale_density(
            fn_int,
            axis=scale_axis,
            axis_scale=axis_scale,
            )
        new_atoms *= 10.0
        topos.write_pdb(fn_out=fn_int,
                           index=new_index,
                           symbols=new_symbols,
                           coordinates=new_atoms,
                           )
    
    # => Adding Proton Defect <= #
    print "Build: Adding proton defect"
    final_atoms = man.proton(fn_int,
                            ostar_id,
                           )
    
    # updating indices,symbols and charges
    new_index = np.arange(1,len(final_atoms)+1)        
    new_symbols = np.hstack((new_symbols,['H']))
    new_charge = np.hstack((new_charge,[0.3]))
   
    #EDIT HERE: TODO 
    # writing final pdb file
    final_atoms *= 10.0 
    if cut_sphere == True:
        actual_center = com.com(final_atoms,new_symbols)
        xhi = max(final_atoms[:,0])+10.0
        yhi = max(final_atoms[:,1])+10.0
        zhi = max(final_atoms[:,2])+10.0
        x_disp = xhi/2.0-actual_center[0]
        y_disp = yhi/2.0-actual_center[1]
        z_disp = zhi/2.0-actual_center[2]
        
        print "Build: Centering system in unit cell"
        final_atoms = man.displace(final_atoms,disp=[x_disp,y_disp,z_disp])
        final_opts = {
                     'xhi' : xhi,
                     'yhi' : yhi,
                     'zhi' : zhi,
                     }
    else:
        final_opts = {}
   
    print "Build: Writing final pdb file" 
    topos.write_pdb(fn_out=fn_pdb,
                       index=new_index,
                       symbols=new_symbols,
                       coordinates=final_atoms,     
                       in_opts = final_opts
                       )
    
    # => Writing Lammps data file <= #
    if center_sphere == True:
        actual_center = com.com(final_atoms,new_symbols)
        x_disp = -actual_center[0]
        y_disp = -actual_center[1]
        z_disp = -actual_center[2]
        
        print "Build: Centering system in unit cell"
        final_atoms = man.displace(final_atoms,disp=[x_disp,y_disp,z_disp])
        xhi = max(final_atoms[:,0])+10.0
        yhi = max(final_atoms[:,1])+10.0
        zhi = max(final_atoms[:,2])+10.0
        xlo = min(final_atoms[:,0])-10.0
        ylo = min(final_atoms[:,1])-10.0
        zlo = min(final_atoms[:,2])-10.0
        final_opts = {
                     'xhi' : xhi,
                     'yhi' : yhi,
                     'zhi' : zhi,
                     'xlo' : xlo,
                     'ylo' : ylo,
                     'zlo' : zlo,
                     }
    
    # converting symbols to lammps id numbers
    atom_type = topos.symbols_ids(new_symbols)
    if cut_plane == True:
        final_atoms = man.displace(final_atoms,disp=[10.0,10.0,10.0]) 
    print "Build: Writing final lammps file" 
    
    topos.write_lammps(fn_out=fn_lammps,
                       index=new_index,
                       atom_type=atom_type,
                       charge=new_charge,
                       coordinates=final_atoms,
                       in_opts = final_opts
                       ) 
    
    # => Writing Plumed File <= #
   
    print "Build: Writing plumed file" 

    natoms = len(final_atoms)
    print "Number of Water Molecules: %d" % ((natoms-1)/3)
    
    max_X = max(final_atoms[:,0])
    max_Y = max(final_atoms[:,1])
    max_Z = max(final_atoms[:,2])
    
    min_X = min(final_atoms[:,0])
    min_Y = min(final_atoms[:,1])
    min_Z = min(final_atoms[:,2])
    
   # centerX = (max_X+min_X)/2.0
   # centerY = max_Y/2.0
   # centerZ = max_Z/2.0
    centerXYZ = com.com(final_atoms,new_symbols)
    centerX = centerXYZ[0]
    centerY = centerXYZ[1]
    centerZ = centerXYZ[2]
    
    xuat = centerX + plane_gap/2
    yuat = max_Y + edge_gap
    zuat = max_Z + edge_gap
    
    xlat = centerX - plane_gap/2
    ylat = min_Y - edge_gap
    zlat = min_Z - edge_gap
    
    inds_atoms = [ind+1 for ind,val in enumerate(new_symbols) if val == 'O']

    if cut_plane == True:    
        if periodic == False:
            plumed.plane(inds_atoms,
                     xuat = xuat/10.0,
                     xlat = xlat/10.0,
                     yuat = yuat/10.0,
                     zuat = zuat/10.0,
                     ylat = ylat/10.0,
                     zlat = zlat/10.0,
                     kappa=kappa,
                     exp=exp)
        
        else:
            plumed.plane(inds_atoms,
                     xuat = xuat,
                     xlat = xlat,
                     kappa=kappa,
                     exp=exp)

    #if cut_sphere == True:

    #    at = (max_X-min_X)/2.0
    #    plumed2.micelle(
    #        atomIDs=inds_atoms,
    #        centerX=0.0,
    #        centerY=0.0,
    #        centerZ=0.0,
    #        Rmax=Rmax,
    #        sigma=sigma,
    #        epsilon=epsilon,
    #        alpha=alpha,
    #        arg='rad',
    #        name='bias',
    #        )

    return

    
def test():

    in_opts = {    

    'cut_sphere' : True,
    'nwaters' : 500,
    
    'cut_plane' : False,
    'plane_gap' : 0.2,
    
    'expand_system' : False,
    'axis' : [1,2],
    'times' : [2,2],
    'exp_gap' : 0.1,
    
    'scale_system' : False,
    'axis_scale' : 0.8,
    'scale_axis' : [1,2],
    
    'ostar_id' : 25,
    
    'periodic' : False,
    
    # plumed parameters
    'kappa' : 250.0,
    'exp' : 2.0,
    'edge_gap' : 0.1,
    
    # Filenames 
    'fn_in' : "data.water-out",
    'fn_lammps' : "data.water",
    'fn_pdb' : "water.pdb",
    'fn_int' : "ref.pdb",

    }

    build_system(in_opts)    
    
    return

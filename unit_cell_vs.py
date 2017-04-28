#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as mdt

def vectors(
        fn_top='water.pdb',
	fn_traj = 'positions.xtc',
	fn_out = "unit_cell.npy",
	len_chunk = 1,
        cut=1E23,
        ):
	
	# => Load Info <= #
	
	top = mdt.load(fn_top).topology
	traj = mdt.iterload(fn_traj, top=top, chunk=len_chunk)

        step = 0
	for chunk in traj:
            if step == 0:
                vectors = chunk.unitcell_lengths
            else:
                vectors = np.vstack((vectors,chunk.unitcell_lengths))
            step += len_chunk
            if step == cut:
                np.array(vectors).reshape(-1,3).dump(fn_out)
                break
            print 'step: %d\r' % (step),

        np.array(vectors).reshape(-1,3).dump(fn_out)
        
       


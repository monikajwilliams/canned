#!/usr/bin/env python
import os, sys, re, math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def unwrap_pos(
                fn_in,
                fn_out='unwrapped.npy',
                fn_vectors='unit_cell.npy',
                ):

    Ostar = np.load(fn_in)
    vectors = np.load(fn_vectors)
    new_Ostar = np.zeros_like(Ostar)


    for coord in range(3):    
        
        # original coordinates
        x_wrapped = Ostar[:,coord]
        x_cell = vectors[:,coord]

        # array for new coordinates
        x_unwrapped = np.copy(x_wrapped)
        # cell lengths
        x_cell = vectors[:,coord]

        dx = np.array([(x_wrapped[ind]-x_wrapped[ind-1]) for ind in range(len(x_wrapped)-1)])
        ind_jumps = np.array([ind for ind,val in enumerate(dx) if abs(dx[ind]) > (x_cell[ind]/2.0)])
        val_jumps = np.array([-x_cell[ind] if val > 0.0 else x_cell[ind] for ind,val in enumerate(dx)])

        for ind,val in enumerate(ind_jumps):
            x_unwrapped[val:] += val_jumps[val]
        
        new_Ostar[:,coord] = x_unwrapped
    np.array(new_Ostar).dump(fn_out)


    return

# Takes unwrapped coordinates and wraps back into periodic box
def wrap_pos(
    fn_in,
    fn_out='wrapped.npy',
    fn_vectors='unit_cell.npy',
    ):

    Ostar = np.load(fn_in)
    vectors = np.load(fn_vectors)

    new_Ostar = np.mod(Ostar / vectors,1) * vectors

    np.array(new_Ostar).dump(fn_out)


    return





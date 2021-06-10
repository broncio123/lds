#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 12:36:43 2020

@author: vk17590
"""
import os
import pathlib
from os.path import abspath
import numpy as np
import matplotlib.pyplot as plt
import h5py
from ldds.base import fit_pes

def discretise_potential(coords, potential):
    """
    Returns a 1- or 2-dimensional array of function (potential energy) values on a grid of points.

    Parameters
    ----------
    coords : list of ndarrays,
        [x] or [x,y] contain coordinates.

    potential : function,
        Function/potential energy to be discretised.

    Returns
    -------
    pes_data : ndarray, shape(len(x)) or shape(len(y),len(x)),
        Array of function/potential values.
    """

    if len(coords) == 1:
        pes_data = potential(coords)
        return pes_data

    elif len(coords) == 2:
        x, y = coords
        X,Y = np.meshgrid(x,y,indexing='xy')
        points = np.column_stack((X.flatten(),Y.flatten()))
        pes_data = potential(points)
        return pes_data.reshape(len(x),len(y))

    else:
        print('splines in +3D are not implemented yet')

def generate_pes_data(coords, potential, filename):
    """
    Saves 1- or 2-dimensional array of function (potential energy) values on a grid of points to 
    ldds/pes_files/filename.hdf5. File format fixed to HDF5 by default.

    Parameters
    ----------
    coords : list (or ndarray) of ndarrays
        [x] or [x,y] contain coordinates.

    potential : function
        Potential to be discretised.

    filename : string
    Name of hdf5 file to be saved.
    """

    pes_data = discretise_potential(coords, potential)
    
    dirname = "pes_files"
    dirpath = os.path.join(pathlib.Path(__file__).parent.absolute(), dirname)
    filepath = os.path.join(dirpath, filename+'.hdf5')
    
    hf = h5py.File(filepath,'w')
    hf.create_dataset('coords', data=np.array(coords).astype('float64'))
    hf.create_dataset('pes_data', data=pes_data.astype('float64'))
    hf.close()

def discretise_vector_field(sample_time_points, sample_coords, vector_field):
    """
    Returns a 1-dimensional array of 2-D vector field values evaluated on a 2D grid of points 
    for a 1-dimensional array of time-points. 

    Parameters
    ----------
    sample_time_points: 1d array,
        time-points in a sample time-interval.
    
    sample_coords : list of ndarrays,
        [x,y] contain coordinates.

    vector_field : function,
        Function/vector field to be discretised.

    Returns
    -------
    vector_field_data : ndarray, len(t),
        Array of function/vector field values.
    """
    x, y = sample_coords
    X,Y = np.meshgrid(x,y,indexing='xy')
    sample_xy_points = np.column_stack([X.ravel(),Y.ravel()])
    
    vector_field_data = [vector_field(t, sample_xy_points) for t in sample_time_points]

    return np.array(vector_field_data)

def generate_vector_field_data(sample_time_points, sample_coords, vector_field, filename):
    """
    Saves 2-dimensional array of function (vector field) values on a grid of points to 
    ldds/vector_field_files/filename.hdf5. File format fixed to HDF5 by default.

    Note that storage fails if x_axis and y_axis in sample_coords are not the same length.

    Parameters
    ----------
    sample_time_points: 1d array,
        time-points in a sample time-interval.
    
    sample_coords : list of ndarrays,
        [x,y] contain coordinates.

    vector_field : function,
        Function/vector field to be discretised.

    filename : string
        Name of hdf5 file to be saved.
    """
    
    vector_field_data = discretise_vector_field(sample_time_points, sample_coords, vector_field)
    
    dirname = "vector_field_files"
    dirpath = os.path.join(pathlib.Path(__file__).parent.absolute(), dirname)
    filepath = os.path.join(dirpath, filename+'.hdf5')
    
    hf = h5py.File(filepath,'w')
    
    hf.create_dataset('sample_time', data=np.array(sample_time_points).astype('float64'))
    hf.create_dataset('sample_coords', data=np.array(sample_coords).astype('float64'))
    hf.create_dataset('vector_field_data', data=vector_field_data.astype('float64'))
    
    hf.close()

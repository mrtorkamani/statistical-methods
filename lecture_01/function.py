#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 19:56:30 2023

@author: mr
"""
import numpy as np

def sample_sphere(n):
    # Generate random points in spherical coordinates
    phi = 2 * np.pi * np.random.rand(n)  # Azimuthal angle (longitude)
    costheta = 2 * np.random.rand(n) - 1  # Cosine of polar angle (colatitude)
    theta = np.arccos(costheta)  # Polar angle (latitude)

    # Return the spherical coordinates
    return theta, phi

def get_cartesian_coords(theta,phi,r=1):
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    return x, y, z
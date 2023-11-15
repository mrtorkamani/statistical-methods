#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 14:29:23 2023

@author: mr
"""
import numpy as np
from scipy.optimize import minimize
from astropy.io import fits

# Load data
path = "data/k_nmf_derived.newdefault.fits"
temps = fits.open(path)
tspec = temps[1].data
lam = temps[11].data

filters = np.loadtxt('data/bandfilters.txt')
mufilters = (filters[1:] + filters[:-1]) * 0.5

# Load catalog and true redshifts
catalog = np.loadtxt('data/catalog1.txt')
true_redshifts = np.loadtxt('data/redshifts1.txt')

# Define the likelihood function
def log_likelihood(amplitude, redshift, template=tspec[0], wavelengths=lam, observed_fluxes=None, errors=None):
    if observed_fluxes is None:
        raise ValueError("Observed fluxes must be provided.")
    
    if errors is None:
        errors = np.ones_like(observed_fluxes)  # Assuming equal errors for simplicity
    
    # Apply redshift
    shifted_wavelengths = wavelengths * (1 + redshift)
    
    # Interpolate template onto observed wavelengths
    template_flux = np.interp(shifted_wavelengths, wavelengths, template)
    
    # Calculate residuals for each band
    residuals = (template_flux * amplitude - observed_fluxes) / errors
    
    # Calculate log-likelihood
    log_likelihood_value = -0.5 * np.sum(residuals**2 + np.log(2 * np.pi * errors**2))
    
    # Print diagnostic information
    print(f"Amplitude = {amplitude}, Redshift = {redshift}, Likelihood = {log_likelihood_value}")
    
    return -log_likelihood_value  # Minimize the negative log-likelihood

# Iterate over galaxies in the catalog
estimated_redshifts = []

# Redshift range
z_range = np.linspace(0, 1.5, 100)

for i, galaxy_data in enumerate(catalog):
    # Define the likelihood function for a given redshift
    likelihood_for_redshift = lambda z: log_likelihood(np.max(galaxy_data), z, template=tspec[0], wavelengths=lam, observed_fluxes=galaxy_data)
    
    # Find the redshift that maximizes the likelihood
    best_redshift = minimize(likelihood_for_redshift, x0=0.0, bounds=[(0, 1.5)])
    
    # Print diagnostic information
    print(f"Galaxy {i + 1}: Initial Redshift = 0.0, Best Redshift = {best_redshift.x[0]}, Likelihood = {-best_redshift.fun}")
    
    # Store the estimated redshift
    estimated_redshifts.append(best_redshift.x[0])

# Convert the result to a numpy array
estimated_redshifts = np.array(estimated_redshifts)

# Print or analyze the estimated redshifts
print("Estimated Redshifts:", estimated_redshifts)
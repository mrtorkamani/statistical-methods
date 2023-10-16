#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 20:12:53 2023

@author: mr
"""


import matplotlib.pyplot as plt

def plot_scatter(x,y,z):
    fig, ax = plt.subplots(3, 1,figsize=(5,15))
    ax[0].scatter(x,y,s=1)
    ax[0].set_xlabel('x')
    ax[0].set_ylabel('y')
    ax[1].scatter(x,z,s=1)
    ax[1].set_xlabel('x')
    ax[1].set_ylabel('z')
    ax[2].scatter(y,z,s=1)
    ax[2].set_xlabel('y')
    ax[2].set_ylabel('z')
    
def plot_histogram(x,y,z):
    fig, axs = plt.subplots(3, 1,figsize=(5,15))
    axs[0].hist(x,bins=30)   
    axs[1].hist(y,bins=30)   
    axs[2].hist(z,bins=30)
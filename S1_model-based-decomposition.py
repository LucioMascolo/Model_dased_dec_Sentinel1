#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pathlib

import numpy as np
from spectral import envi

#import psp_tools as psp 

import matplotlib.pyplot as plt

from PIL import Image

###################################################################################
##                                  DESCRIPTION 
## Decompose the Sentinel-1 Stokes vector with the method published in the article:
## L. Mascolo, S.R. Cloude, J.M. Lopez-Sanchez, "Model-based decomposition of 
## dual-pol SAR data: application to Sentinel-1", Trans. and Geosci. Remote Sens.
##
## Inputs:
##
## - Elements of the wave coherency matrix C2 (either from HH-HV or VH-VV data)
## Outputs:
## - mv, mp, alpha and delta
## - RGB, alpha and HSV images
## 
## Authors:
## - Lucio Mascolo, formerly at the Institute for Computer Research (IUII), University of Alicante
## - Shane R. Cloude, formerly at Applied Electromagnetics Consultants (AELC)
## - Juan Manuel Lopez-Sanchez, Institute for Computer Research (IUII), University of Alicante
## Date: 20/01/2025
################################################################################


## Path to the main C2 folder (time series of C2 matrices)
c2_folder  = '/Volumes/..../C2/'

## Path to ouptup folder 
folder_out = '/Volumes/..../C2/'

## Sentinel-1 acqusition dates: these are the final folders containing the C2 at those dates
dates     = ['20211226', '20220107'] 


## Loop for all the images in the time series
ndates    = len(dates)
for date in dates:
       
    print('date = ' + str(date))

    # Folder       
    folder_date_in = c2_folder + '/' + date + '/'
    
    folder_date_out = folder_out + '/' + date 
    pathlib.Path(folder_date_out).mkdir(parents=True, exist_ok=True)

    ## Load the C2 elements
    
    data_hdr = envi.open(folder_date_in + '/C11.hdr', folder_date_in + '/C11.bin')
    C11 = data_hdr.read_band(0)
    
    data_hdr = envi.open(folder_date_in + '/C12_real.hdr', folder_date_in + '/C12_real.bin')
    C12_real = data_hdr.read_band(0)
    
    data_hdr = envi.open(folder_date_in + '/C12_imag.hdr', folder_date_in + '/C12_imag.bin')
    C12_imag = data_hdr.read_band(0)
    
    data_hdr = envi.open(folder_date_in + '/C22.hdr', folder_date_in + '/C22.bin')
    C22 = data_hdr.read_band(0)
    
    ## Apply model-based decomposition
    
    ## - Build Stokes vector from C2 
    ## ** Updated Stokes convention: C11 = VV backscattering ; C22 = VH backscattering
    
    s1               = C11 + C22
    s2               = C22 - C11
    s3               = 2*C12_real
    s4               = -2*C12_imag
    

    ## - Solve the quadratic (solution for the volume power)
    a                = 0.75
    b                = -2*( s1 - (0.5*s2) )
    c                = s1*s1 - s2*s2 - s3*s3 - s4*s4
    delta1           = b**2
    delta2           = 4*a*c
    delta            = delta1 - delta2
    
    mv1              = -b + np.sqrt(delta)
    mv1              = mv1 / (2*a)
    mv2              = (-b - np.sqrt(delta))
    mv2              = mv2 / (2*a)
    
    ## - check which solution satisfes s1>mv
    
    ind1   =  (s1[:] > 0) # exclude pixels outside the imaged scene
    cond_1 = s1[ind1] > mv1[ind1] 
    
    ## - volume power solution
    
    if cond_1.size == 0:
        mv = mv1
    else:
        mv = mv2
        
    ## - Obtain mp (polarized power)
    mp = s1 - mv
    
    ## - obtain alpha and delta
    alpha = 0.5 * np.arccos( (s2-0.5*mv) / mp )
    
    delta = np.angle( s3 + 1j*s4)

    ## Convert powers to dB
    
    mpdb  = 10*np.log10(mp + 1e-12)
    mvdb  = 10*np.log10(mv + 1e-12)
    ratdb = mpdb-mvdb    
    
    ## - Save outputs
    
    envi.save_image(folder_date_out + '/' + 'alpha' + '.hdr', alpha*180/np.pi, \
                    force = True, interleave = 'BSQ', ext = '.img', metadata=data_hdr.metadata) 

    envi.save_image(folder_date_out + '/' + 'mp_dB' + '.hdr', mpdb, \
                   force = True, interleave = 'BSQ', ext = '.img', metadata=data_hdr.metadata) 

    envi.save_image(folder_date_out + '/' + 'mv_dB' + '.hdr', mvdb, \
                   force = True, interleave = 'BSQ', ext = '.img', metadata=data_hdr.metadata) 

    ##- Save RGB image
    
    min_m = -20
    max_m = -2
    min_r = -10
    max_r = 10

    mpdb [ mpdb < min_m ] = min_m
    mpdb [ mpdb > max_m ] = max_m

    mvdb [ mvdb < min_m ] = min_m
    mvdb [ mvdb > max_m ] = max_m

    ratdb [ ratdb < min_r ] = min_r
    ratdb [ ratdb > max_r ] = max_r
    
    Nrow, Ncol = C11.shape
    RGB_image = np.zeros([Nrow, Ncol, 3], dtype=np.uint8)

    RGB_image[:,:,0] = np.array(np.ceil(255 * (mpdb - min_m) / (max_m - min_m)), dtype=np.uint8)
    RGB_image[:,:,1] = np.array(np.ceil(255 * (mvdb - min_m) / (max_m - min_m)), dtype=np.uint8)
    RGB_image[:,:,2] = np.array(np.ceil(255 * (ratdb - min_r) / (max_r - min_r)), dtype=np.uint8)
   
    pIm = Image.fromarray(RGB_image).save(folder_date_out + '/RGB_stokes_decomposition.png') 

   
    ## Show and save alpha image
    
    fig, ax = plt.subplots()

    im = ax.imshow(alpha*180/np.pi, vmin = 0, vmax = 90, cmap='jet', interpolation='none')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    cbar = fig.colorbar(im, ax = ax, shrink = 0.7)
    cbar.set_label(''r'$\alpha$ (degrees)')

    fig.savefig(folder_date_out + '/' + 'alpha.pdf', dpi=150)
    
 
    ## Show and save HSV image 
        
    HSV_image = np.zeros([Nrow, Ncol, 3], dtype=np.uint8)

    # Hue (delta)
    pscale        = (delta+np.pi)/(2*np.pi) # scale the phase in the range 0 to 1
    HSV_image[:,:,0] = np.array(np.ceil(255 * pscale), dtype=np.uint8)
    
    # Saturation (cross-polarized coherence)
    ro            = np.sqrt( C12_real*C12_real + C12_imag*C12_imag )
    ro            = ro / np.sqrt( C11*C22 )
    HSV_image[:,:,1] = np.array(np.ceil(255 * ro), dtype=np.uint8)
    
    # Value (span)
    spandB        = 10*np.log10(C11+C22 + 1e-12)
    
    spandB [ spandB < -25 ] = -25
    spandB [ spandB > 5 ] = 5

    HSV_image[:,:,2] = np.array(np.ceil(255 * (spandB + 25) / 30), dtype=np.uint8)

    pIm = Image.fromarray(HSV_image, 'HSV').convert("RGB").save(folder_date_out + '/HSV_stokes_decomposition.png') 
    


    
    
    
    
    
    
    
    
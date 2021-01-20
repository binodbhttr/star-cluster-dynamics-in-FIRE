import numpy as np
import pdb
from importlib import reload

#pdb.set_trace()  #<--in case need to troubleshoot

#center of mass of an array (e.g. x, y, z) 
def cm(x, mass):   
    total_mass=sum(mass)  
    cm = sum(mass[i]*x[i] for i in range(len(x)))/ total_mass
    return cm

#array of distances between 1d positions and center of mass on that axis
def diff_1d(x, xcm):
    diff_arr = []
    x = np.array(x)
    xcm = np.array(xcm) 
    
    for xpos in x:
        if xpos > 0 and xcm > 0:
            if xcm > xpos:
                xdiff = xcm - xpos
            else:
                xdiff = xpos - xcm
        if xpos < 0 and xcm > 0:
            xdiff = xcm + np.abs(xpos)
        if xpos > 0 and xcm < 0:
            xdiff = xpos + np.abs(xcm)
        if xpos < 0 and xcm < 0:
            if xcm > xpos: #xcm=-.1 xpos=-.2 xdiff=.1
                xdiff = np.abs(xpos) - np.abs(xcm)
            else: #xpos=-.1 xcm=-.3 xdiff = .2
                xdiff = np.abs(xcm) - np.abs(xpos)
        diff_arr.append(xdiff)
    diff_arr = np.array(diff_arr)
        
    return diff_arr

#spherical radial distance from center of mass to each object
def dr(x,y,z,mass):  
    xcm = cm(x, mass)
    ycm = cm(y, mass)
    zcm = cm(z, mass)

    xdiff = diff_1d(x, xcm) #generates a np.array of differences
    ydiff = diff_1d(y, ycm)
    zdiff = diff_1d(z, zcm)
            
    dr = np.sqrt(xdiff**2. + ydiff**2. + zdiff**2.)
    dr = dr.flatten() #flattened (not sure need this now)
    
    return dr

#max distance from center of mass
def drmax(x,y,z,mass):
    deltar = dr(x,y,z,mass)
    drmax = max(deltar)
    
    return drmax

#to sort dr in order of increasing size
def dr_sort(dr):
    print('watch out order matters here: returns sorted array and then indices of input array')
    dr = np.array(dr) #just in case make sure it is a numpy array
    sortind = np.argsort(dr) 
    dr_sorted = dr[sortind]
    
    return(dr_sorted, sortind) 

#sort any array (x, y, z, what have you) according to index order you give it
def arr_sort(arr, sortind): #sort according to index order passed
    arr_sort = arr[sortind]
    arr_sort = np.array(arr_sort)
    return arr_sort

#pass indices sorted by a particular order (say increasing dr) and then figure out what fractional mass is enclosed (mass sorted for you)
def frac_enc_mass(mass, sortind):
    mass = np.array(mass) #just in case make sure it is a numpy array
    total_mass = np.sum(mass)
    mass_sort = arr_sort(mass, sortind) #mass in order of increasing dr
    mass_sum = 0
    mass_arr = []
    for m in mass_sort: #does this in order
        mass_sum =  mass_sum + m
        mass_arr.append(mass_sum)
        
    mass_arr = np.array(mass_arr)
    frac_enc_mass_arr = mass_arr / (1.*total_mass)
    return frac_enc_mass_arr

#take fraction enclosed mass array and find the location relative to dr_sort that it is greater than or equal to fraction value (for example 90% of the mass)
def dr_mass_frac(frac_enc_mass, fraction, dr_sort):
    indfrac = np.where(frac_enc_mass >= fraction)
    if len(indfrac[0]) == 1:   #only 1 entry >= 90% --> 100%. apd nearest entry
        first = indfrac[0][0] - 1
        second = indfrac[0][0]          
        temp_dr_gt_frac = dr_sort[[first, second]]
    else:
        temp_dr_gt_frac = dr_sort[indfrac] #array of dr
        
    findexfrac = np.where(temp_dr_gt_frac == min(temp_dr_gt_frac))
    drfrac = temp_dr_gt_frac[findexfrac][0]
    
    return(drfrac)

#individual embedded check: is the small object embedded fullying inside the large object
#assumes pass it the relevant individual individual positions
def individual_embedded_check(xcm_small, ycm_small, zcm_small, rmax_small, xcm_big, ycm_big, zcm_big, rmax_big):  

    xmax_big = xcm_big + rmax_big
    xmin_big = xcm_big - rmax_big
    ymax_big = ycm_big + rmax_big
    ymin_big = ycm_big - rmax_big
    zmax_big = zcm_big + rmax_big
    zmin_big = zcm_big - rmax_big
    
    xmax_small = xcm_small + rmax_small
    xmin_small = xcm_small - rmax_small
    ymax_small = ycm_small + rmax_small
    ymin_small = ycm_small - rmax_small
    zmax_small = zcm_small + rmax_small
    zmin_small = zcm_small - rmax_small

    if (xmax_big > xmax_small) & (xmin_big < xmin_small):
        xcheck = True
    else:
        xcheck = False
        
    if (ymax_big > ymax_small) & (ymin_big < ymin_small):
        ycheck = True
    else:
        ycheck = False
        
    if (zmax_big > zmax_small) & (zmin_big < zmin_small):
        zcheck = True
    else:
        zcheck = False

    if (xcheck == True) & (ycheck == True) & (zcheck == True):
        boolean_result = True
    else:
        boolean_result = False
    
    return boolean_result

#find closest gmc to individual star cluster (xcm, ycm, zcm)
#this can match multiple clusters to the same gmc
def closest_cluster_index(xcm1, ycm1, zcm1, xcm_arr2, ycm_arr2, zcm_arr2):

    xdiff = diff_1d(xcm_arr2, xcm1) #generates a np.array of differences
    ydiff = diff_1d(ycm_arr2, ycm1)
    zdiff = diff_1d(zcm_arr2, zcm1)
            
    dr = np.sqrt(xdiff**2. + ydiff**2. + zdiff**2.)
    dr = dr.flatten() #flattened (not sure need this now)

    drmin = min(dr)
    index = np.where(dr == drmin)
    index = index[0] #select the first one if there are multiple?
    
    return (index, drmin)

#xcm_arr1 = star cluster arr
def array_embedded_check(xcm_arr1, ycm_arr1, zcm_arr1, rmax_arr1, xcm_arr2, ycm_arr2, zcm_arr2, rmax_arr2):

    bool_arr = np.zeros(len(xcm_arr1), dtype=bool)
    
    for i in range(len(xcm_arr1)):
        xcm1  = xcm_arr1[i]
        ycm1  = ycm_arr1[i]
        zcm1  = zcm_arr1[i]
        rmax1 = rmax_arr1[i]

        #get the index of the closest gmc from gmc arr
        ind2, drmin = closest_cluster_index(xcm1, ycm1, zcm1, xcm_arr2, ycm_arr2, zcm_arr2)

        #check if it is embedded
        boolean_check = individual_embedded_check(xcm1, ycm1, zcm1, rmax1, xcm_arr2[ind2], ycm_arr2[ind2], zcm_arr2[ind2], rmax_arr2[ind2])
        bool_arr[i] = boolean_check

    return(bool_arr)


def aperture_selection_ind(x, y, xcm, ycm, aperture_parsec):

    aperture_kpc        = aperture_parsec / 1000.
    aperture_radius_kpc = aperture_kpc / 2.0
    
    xdiff = diff_1d(x, xcm) #generates a np.array of differences
    ydiff = diff_1d(y, ycm)
            
    dr = np.sqrt(xdiff**2. + ydiff**2.)

    ind = np.where(dr <= aperture_radius_kpc)
    ind = ind[0]
    return(ind)

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


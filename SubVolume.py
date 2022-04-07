##################################################################################################################################################### 
# This script discards positions that are outside the investigation volume or that belong to short trajectories. These short trajectories typically
# result from impurities in the water, bright elements in the image background, or reflections on the surface of the aquarium. 
   
# The code flags blocks of coordinates and identifies the starting and ending indices of each trajectory within these blocks. This approach 
# is much faster than looping over trajectories.

# Input and output ----------------------------------------------------------------------------------------------------------------------------------
# File (txt) containing [cx, cy, cz, nb, fg, ts, rk] 

# where: cx, cy and cz are the coordinates (mm)
#        nb is the trajectory number
#        fg is the flag for added positions (1 if added and 0 otherwise)
#        ts is the time step, starting at one and without padding
#        rk is the row index of the particle in the corresponding ptv file (only for particles that are not added)
#####################################################################################################################################################

#####################################################################################################################################################
import numpy as np
# import pandas as pd

import os
# import re
#####################################################################################################################################################

#####################################################################################################################################################
inpt = ["D:/Input folder for Recording 01/",
        "D:/Input folder for Recording 02/",
        "D:/Input folder for Recording 03/"
        ]
    
oupt = ["D:/Output folder for Recording 01/",
        "D:/Output folder for Recording 02/",
        "D:/Output folder for Recording 03/"
        ] 

nmat = ["Recording_01_100001_101000_Buildingptv", # Input file name
        "Recording_02_100001_101000_Buildingptv", # Input file name
        "Recording_03_100001_101000_Buildingptv"  # Input file name
        ]  
#####################################################################################################################################################

#####################################################################################################################################################
leng = 5  # Minimal length (in mm) 
dura = 20 # Minimal duration (in frame(s)) 
#####################################################################################################################################################

#####################################################################################################################################################
# Limits of the investigation volume (in mm) --------------------------------------------------------------------------------------------------------
volu = {'X':[10, 90], 
        'Y':[10, 90], 
        'Z':[10, 90]} # OpenPTV reference frame
#####################################################################################################################################################

#####################################################################################################################################################
# File loop -----------------------------------------------------------------------------------------------------------------------------------------
for index, entry in enumerate(nmat):
      
    name = os.path.basename(entry)
    
    print(f"Processing file {name}")
  
    data = np.loadtxt(os.path.join(inpt[index], "".join((name, '.txt'))))
                        
    # Flag coordinates outside the investigation volume (in mm) -------------------------------------------------------------------------------------    
    disc = np.logical_or.reduce((data[:,0] < volu['X'][0], data[:,0] > volu['X'][1],    
                                 data[:,1] < volu['Y'][0], data[:,1] > volu['Y'][1], 
                                 data[:,2] < volu['Z'][0], data[:,2] > volu['Z'][1]))                   
         
    # Extract indices of coordinates inside the investigation volume --------------------------------------------------------------------------------   
    disc = np.invert(disc) 
    disc = disc.astype(int) 
        
    ivec = np.diff(np.concatenate(([0], disc, [0])))
            
    stti = np.where(ivec == 1)[0] 
    stpi = np.where(ivec == -1)[0]
    
    keep = []
    
    # Length and duration filtering and numbering ---------------------------------------------------------------------------------------------------
    for segi, vala in enumerate(stti):
   
        tvec = np.add(np.where(np.diff(data[stti[segi]:stpi[segi], 3]))[0],1) 
        
        bidx = np.hstack([stti[segi], stti[segi] + tvec])
        eidx = np.hstack([stti[segi] + tvec, stpi[segi]])
        
        for posi, valb in enumerate(bidx):
            
            if (np.size(np.arange(bidx[posi], eidx[posi])) >= dura and 
                np.any(np.subtract(np.amax(data[bidx[posi]:eidx[posi], 0:3], axis = 0), 
                                   np.amin(data[bidx[posi]:eidx[posi], 0:3], axis = 0)) >= leng) == True):           
            
                # Fragment is longer than threshold in space in at least one dimension.            
                              
                keep.append(np.arange(bidx[posi], eidx[posi]))
            
    # New trajectory numbering ----------------------------------------------------------------------------------------------------------------------        
    for coun, valc in enumerate(keep, start = 0): 
        
        data[keep[coun], 3] = np.tile(np.add(coun,1), np.size(keep[coun]))
        
    # Store in array with single precision ----------------------------------------------------------------------------------------------------------
    temp = data[np.hstack(keep), :].astype('float32')   

    # Save trajectories -----------------------------------------------------------------------------------------------------------------------------
    ntmp = "".join((name, "_SubVolume.txt"))
    
    np.savetxt(os.path.join(oupt[index], ntmp), temp, delimiter = '\t')
#####################################################################################################################################################

#####################################################################################################################################################
# This function reconstructs trajectories without linking broken segments, using xuap files (txt) from the postprocessing code associated 
# with OpenPTV. It gives the same results as the MATLAB function.   

# The full sequence of input files is divided into several subsequences that reflect changes in experimental conditions during the recording. 
# To process the entire sequence without creating subsequences, set seqNumb to one. 

# Output --------------------------------------------------------------------------------------------------------------------------------------------
# [nb, ts, fg, cx, cy, cz] in file coor (txt)
# [nb, ts, fg, ux, uy, uz] in file velo (txt)
# [nb, ts, fg, ax, ay, az] in file acce (txt)

# where: nb is the trajectory number
#        ts is the time step, starting at one and without padding 
#        fg is the flag for filtered positions
#        cx, cy and cz are the coordinates (m)
#        u is the velocity (m/s)
#        a is the acceleration (m/s^2)

# Time indexing -------------------------------------------------------------------------------------------------------------------------------------
# Indicate the first and last frames of the sequence used to build trajectories in the name of the output file. The code will retrieve 
# this information to determine which files to process. This also allows to keep the correspondence between the time step in column ts,
# which always starts at one, and the number contained in the name of the xuap files, in case the first frame of the sequence used to 
# build the trajectories does not correspond to the first frame of the sequence processed with OpenPTV.   

# The time step ts does not reset between subsequences, but the trajectory number nb does reset between subsequences.
#####################################################################################################################################################

#####################################################################################################################################################
import numpy as np
# import pandas as pd

import os
import re
# import glob

from pathlib import Path
#####################################################################################################################################################

#####################################################################################################################################################
seqNumb = 3        # Number of subsequences
colNumb = 15       # Number of columns in xuap files
maxPart = 6000     # Max number of particles per frame
minSave = 2        # Min trajectory duration for saving (in frames) 

# The minimal duration is two frames no matter the value of minSave, because the initialization of a new trajectory requires an existing link
# to the next time step.
#####################################################################################################################################################

#####################################################################################################################################################
inpt = ["D:/Path to res folder for Recording 01/",
        "D:/Path to res folder for Recording 02/",
        "D:/Path to res folder for Recording 03/",
        ]
    
oupt = ["D:/Output folder for Recording 01/",
        "D:/Output folder for Recording 02/",
        "D:/Output folder for Recording 03/",
        ]  

nmat = ["Recording_01_100001_101000_Buildingxuap", # Output file name
        "Recording_02_100001_101000_Buildingxuap", # Output file name
        "Recording_03_100001_101000_Buildingxuap", # Output file name
        ]
#####################################################################################################################################################

#####################################################################################################################################################
def buildTraj(traj):
    
    global partMat
    
    ifra = traj[-1,0] # ifra is already integer
    ipar = traj[-1,1] # ipar is already integer
       
    while (partMat[ifra, ipar, 1] > 0) and (ifra < np.subtract(np.shape(partMat)[0], 1)): # The ranking of the particles starts at one in xuap files      
        
            ipar = int(np.subtract(partMat[ifra, ipar, 1],1)) 

            ifra = int(np.add(ifra, 1))      
        
            traj = np.append(traj, np.reshape(np.array([ifra, ipar]), (1, 2)), axis = 0)
    
    return traj
#####################################################################################################################################################

#####################################################################################################################################################
def formatTraj(traj, fbeg, svec):
    
    global partMat
    global trajInd 
    
    # Coordinates -----------------------------------------------------------------------------------------------------------------------------------
    cx = partMat[traj[:,0], traj[:,1], np.tile(5, (np.shape(traj)[0]))]
    cy = partMat[traj[:,0], traj[:,1], np.tile(6, (np.shape(traj)[0]))]    
    cz = partMat[traj[:,0], traj[:,1], np.tile(7, (np.shape(traj)[0]))]  
    
    # Velocity --------------------------------------------------------------------------------------------------------------------------------------
    ux = partMat[traj[:,0], traj[:,1], np.tile(8, (np.shape(traj)[0]))]
    uy = partMat[traj[:,0], traj[:,1], np.tile(9, (np.shape(traj)[0]))]    
    uz = partMat[traj[:,0], traj[:,1], np.tile(10, (np.shape(traj)[0]))]  

    # Acceleration ----------------------------------------------------------------------------------------------------------------------------------
    ax = partMat[traj[:,0], traj[:,1], np.tile(11, (np.shape(traj)[0]))]
    ay = partMat[traj[:,0], traj[:,1], np.tile(12, (np.shape(traj)[0]))]    
    az = partMat[traj[:,0], traj[:,1], np.tile(13, (np.shape(traj)[0]))] 

    # Flag for filtered positions -------------------------------------------------------------------------------------------------------------------     
    fg = partMat[traj[:,0], traj[:,1], np.tile(14, (np.shape(traj)[0]))]
              
    # Export data -----------------------------------------------------------------------------------------------------------------------------------
    temp = np.hstack((np.tile(trajInd, (np.shape(traj)[0])).reshape(-1,1), traj[:,0, np.newaxis] + fbeg - svec[0] + 1))
    
    data = dict([('coor', np.hstack((temp, fg[:, np.newaxis], cx[:, np.newaxis], cy[:, np.newaxis], cz[:, np.newaxis]))), 
                 ('velo', np.hstack((temp, fg[:, np.newaxis], ux[:, np.newaxis], uy[:, np.newaxis], uz[:, np.newaxis]))), 
                 ('acce', np.hstack((temp, fg[:, np.newaxis], ax[:, np.newaxis], ay[:, np.newaxis], az[:, np.newaxis])))])
        
    trajInd = np.add(trajInd, 1)   
      
    return data
#####################################################################################################################################################

#####################################################################################################################################################
# Directory loop ------------------------------------------------------------------------------------------------------------------------------------         
for idir, duma in enumerate(nmat):
    name = os.path.basename(duma)
            
    strp = [match.start() for match in re.finditer('_', name)]
    
    nbeg = int(duma[np.add(strp[1],1):strp[2]]) # Adjust here depending on the name of your output file 
    nend = int(duma[np.add(strp[2],1):strp[3]]) # Adjust here depending on the name of your output file      
    
    # Subsequence loop -----------------------------------------------------------------------------------------------------------------------------
    svec = np.linspace(nbeg, nend, num = seqNumb + 1, endpoint = True, dtype = int)  
    
    for iseq, dumb in enumerate(svec[:-1]):
         
        fbeg = svec[iseq] # fbeg is the first frame of the subsequence (padded format)
        fend = svec[np.add(iseq,1)] # fend is the last frame of the subsequence (padded format)
           
        sdur = fend - fbeg + 1 # Number of time steps in the subsequence
        
        trajInd = 1 # Initialize new trajectory number for each subsequence
               
        # Build partMat matrix ----------------------------------------------------------------------------------------------------------------------        
        partMat = np.zeros((sdur, maxPart, np.add(colNumb,1)))

        partMat[:,:,0] = -1 # This works both for ptv and xuap files
        partMat[:,:,1] = -2 # This works both for ptv and xuap files

        posiVec = np.zeros((np.shape(partMat)[0]))
        
        for indx in np.linspace(fbeg, fend, num = sdur, endpoint = True, dtype = int):
 
            temp = np.loadtxt(os.path.join(Path(inpt[idir]), "".join(('xuap.', str(indx))))) 
            
            posiVec[indx - fbeg] = np.shape(temp)[0]   
            partMat[indx - fbeg, 0:np.shape(temp)[0], 0:colNumb] = temp           

            # We can initialize all the trajectories of the current time step if it is the first time step of the subsequence

            if np.equal(indx, fbeg):  
                
               partMat[indx - fbeg, 0:np.shape(temp)[0], 0] = -1               
               # partMat[indx - fbeg, 0:maxPart, 0] = -1
               # This works both for ptv and xuap files                                                                      
               
        # Open output files -------------------------------------------------------------------------------------------------------------------------
        ouptCoor = open(os.path.join(Path(oupt[idir]), '{n}_{s:02d}_coor.txt'.format(n = nmat[idir], s = np.add(iseq,1))), 'a')    
        ouptVelo = open(os.path.join(Path(oupt[idir]), '{n}_{s:02d}_velo.txt'.format(n = nmat[idir], s = np.add(iseq,1))), 'a')
        ouptAcce = open(os.path.join(Path(oupt[idir]), '{n}_{s:02d}_acce.txt'.format(n = nmat[idir], s = np.add(iseq,1))), 'a')

        # Time step loop ----------------------------------------------------------------------------------------------------------------------------                                     
        for ifra in np.linspace(0, np.shape(partMat)[0], num = np.shape(partMat)[0], endpoint = False, dtype = int):
            
            print('Processing frame number {na:0{width}} out of {nb:0{width}}'.format(na = np.add(ifra,1), nb = sdur, width = len(str(sdur)) ))
            
            # Particle loop -------------------------------------------------------------------------------------------------------------------------
            for ipar in np.linspace(0, int(posiVec[ifra]), num = int(posiVec[ifra]), endpoint = False, dtype = int):
                
                # New trajectory --------------------------------------------------------------------------------------------------------------------                            
                if (partMat[ifra, ipar, colNumb] == 0) and (partMat[ifra, ipar, 1] > 0):
                    
                   # For xuap files it is partMat[ifra, ipar, 1] > 0 because the ranking of the particles starts at one and the flag for missing 
                   # link is (-1). For ptv files it is partMat(ifra, ipar, 1] >= 0 because the ranking of the particles starts at zero and the flag  
                   # for missing link is (-2).                    
                       
                   # If no initialization of all the trajectories at the first frame of each subsequence, do not add partMat[ifra, ipar, 0] == (-1) 
                   # as a condition for the if statement since it prevents the processing of the first time step. This is valid both for ptv and  
                   # xuap files because we have populated the first and second pages of partMat with (-1) and (-2), respectively.  
                
                   # Start and build trajectory -----------------------------------------------------------------------------------------------------
                   traj = np.reshape(np.array([ifra, ipar]), (1, 2)).astype(int)
                   
                   traj = buildTraj(traj)
                                      
                   # Update partMat -----------------------------------------------------------------------------------------------------------------                                                        
                   partMat[traj[:,0], traj[:,1], np.tile(colNumb, (np.shape(traj)[0]))] = 1 # We flag particles already in a trajectory                   
                 
                   # Format and save trajectory -----------------------------------------------------------------------------------------------------
                   if np.shape(traj)[0] >= minSave:
                       
                      data = formatTraj(traj, fbeg, svec)
                      
                      # Save [nb, ts, fg, cx, cy, cz] -----------------------------------------------------------------------------------------------
                      np.savetxt(ouptCoor, data['coor'], fmt = '%i %i %i %f %f %f')        
            
                      # Save [nb, ts, fg, ux, uy, uz] -----------------------------------------------------------------------------------------------
                      np.savetxt(ouptVelo, data['velo'], fmt = '%i %i %i %f %f %f')    
            
                      # Save [nb, ts, fg, ax, ay, az] -----------------------------------------------------------------------------------------------
                      np.savetxt(ouptAcce, data['acce'], fmt = '%i %i %i %f %f %f')                       

        # Close output files ------------------------------------------------------------------------------------------------------------------------
        ouptCoor.close()                 
        ouptVelo.close() 
        ouptAcce.close()      
#####################################################################################################################################################        
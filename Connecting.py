##################################################################################################################################################### 
# This code links trajectory segments that have been previously reconstructed with Buildingptv. 

# Input ---------------------------------------------------------------------------------------------------------------------------------------------
# File (txt) containing [cx, cy, cz, nb, fg, ts, rk] 

# where: cx, cy and cz are the coordinates (mm)
#        nb is the trajectory number
#        fg is the flag for added positions (always set to zero)
#        ts is the time step, starting at one and without padding
#        rk is the row index of the particle in the corresponding ptv file

# Output --------------------------------------------------------------------------------------------------------------------------------------------
# File (txt) containing [cx, cy, cz, nb, fg, ts, rk] 

# where: cx, cy and cz are the coordinates (mm)
#        nb is the trajectory number
#        fg is the flag for added positions (1 if added and 0 otherwise)
#        ts is the time step, starting at one and without padding
#        rk is the row index of the particle in the corresponding ptv file (only for particles that are not added)
#####################################################################################################################################################

#####################################################################################################################################################
import numpy as np
import pandas as pd

import os
# import glob

from pathlib import Path
#####################################################################################################################################################

#####################################################################################################################################################
inpt = Path("D:/Path to input folder/")      
oupt = Path("D:/Path to output folder/")

nmat = ["Recording_01_100001_101000_Buildingptv", # Input file name
        "Recording_02_100001_101000_Buildingptv", # Input file name
        "Recording_03_100001_101000_Buildingptv", # Input file name
        ]   
#####################################################################################################################################################

#####################################################################################################################################################
maxiDisp = 10    # Max jump displacement (mm)               
maxiDura = 10    # Max jump duration (frames) 
polySpan = 10    # Max duration of the regression window (frames) 

sepaRadi = (1,5) # Radius of the search volume around predicted positions (mm)  

radi = np.linspace(sepaRadi[0], sepaRadi[1], num = maxiDura, endpoint = True, dtype = float)
#####################################################################################################################################################

#####################################################################################################################################################
def addposi(ta, tb, jp) :
                                     
    dx = (tb[0,0] - ta[-1,0]) / jp 
    dy = (tb[0,1] - ta[-1,1]) / jp             
    dz = (tb[0,2] - ta[-1,2]) / jp             
           
    rg = np.linspace(1, np.subtract(jp,1), num = np.subtract(jp,1), endpoint = True, dtype = int)
            
    # Coordinates -----------------------------------------------------------------------------------------------------------------------------------
    co = np.hstack((ta[-1,0] + dx*rg[:, np.newaxis], 
                    ta[-1,1] + dy*rg[:, np.newaxis], 
                    ta[-1,2] + dz*rg[:, np.newaxis])) 
    
    # Trajectory number and row indice of the particle in ptv file ----------------------------------------------------------------------------------
    nb = np.full(np.size(rg), np.nan).reshape(-1,1)
    rk = np.full(np.size(rg), np.nan).reshape(-1,1)   
        
    # Flag for added positions ----------------------------------------------------------------------------------------------------------------------
    fg = np.ones(np.size(rg), dtype = int).reshape(-1,1)
    
    # Time step -------------------------------------------------------------------------------------------------------------------------------------
    ts = np.add(ta[-1,5], rg).reshape(-1,1)
    
    # Format [cx, cy, cz, nb, fg, ts, rk] -----------------------------------------------------------------------------------------------------------
    se = np.hstack((co, nb, fg, ts, rk))
                
    return se
#####################################################################################################################################################

#####################################################################################################################################################
# File loop -----------------------------------------------------------------------------------------------------------------------------------------
for fidx in nmat:
    
    name = os.path.basename(fidx)
    
    print('Processing file', name)
  
    data = np.loadtxt(os.path.join(inpt, "".join((name, '.txt'))))
        
    #################################################################################################################################################
    # Remove short fragments, renumber fragments, and retrieve starting frame and coordinates of each fragment --------------------------------------         
    stti = np.where(np.diff(np.vstack((0, data[:,[3]])), axis = 0) == 1)[0]    
    stpi = np.where(np.vstack((np.diff(data[:,[3]], axis = 0), 1)) == 1)[0]
        
    dura = np.add(np.subtract(stpi, stti),1)
            
    data = np.split(data, np.cumsum(dura[:-1], dtype = int), axis = 0) # Convert from one single numpy array to a list of numpy arrays
    
    data = np.delete(data, np.where(dura < polySpan)[0]) # Remove entries corresponding to short trajectories in time
    
    cbeg = [] # Coordinates of the first position
    cend = [] # Coordinates of the last position
    time = [] # Starting frame (without padding, starting at one)
    
    for itra, traj in enumerate(data): 
        
        traj[:,3] = np.tile(np.add(itra,1), np.shape(traj)[0]) # Renumber fragment
        
        cbeg.append(traj[0, 0:3])
        
        cend.append(traj[-1, 0:3])
        
        time.append(traj[0, 5])                         
       
    data = np.concatenate(data, axis = 0) # Convert back to one single numpy array
    #################################################################################################################################################

    #################################################################################################################################################
    # Predict positions for each segment and store in dataframe -------------------------------------------------------------------------------------
    fragStor = [] # Initialize list
    framStor = [] # Initialize list
    jumpStor = [] # Initialize list
    distStor = [] # Initialize list 
    candStor = [] # Initialize list
    
    for itra in np.unique(data[:,3]): 
        
        indx = np.where(data[:,3] == itra)[0]
        indx = indx[-polySpan:] # Restrict to last positions for polynomial fitting
          
        cx = np.squeeze(data[indx, 0])
        cy = np.squeeze(data[indx, 1])
        cz = np.squeeze(data[indx, 2])
                        
        a = np.linspace(0, np.size(indx), num = np.size(indx), endpoint = False, dtype = int)
        
        b = np.linspace(np.size(indx), np.add(np.size(indx), maxiDura), num = maxiDura, endpoint = False, dtype = int)
                
        bx = np.polyfit(a, cx, 2) 
        by = np.polyfit(a, cy, 2) 
        bz = np.polyfit(a, cz, 2) 
        
        ex = np.polyval(bx, b)
        ey = np.polyval(by, b)        
        ez = np.polyval(bz, b)
        
        tvec = np.linspace(np.add(data[indx[-1], 5], 1), np.add(data[indx[-1], 5], maxiDura), num = maxiDura, endpoint = True, dtype = int)
                
        for tidx, tval in enumerate(tvec):
                        
            cand = np.where(time == tval)[0] # List of candidates (starting at zero)
            
            temp = list(cbeg[i] for i in cand) # List of arrays containing the starting coordinates of the candidates
            
            funbeg = lambda vari : np.sqrt(np.sum(np.square(vari - [ex[tidx], ey[tidx], ez[tidx]])))
            
            funend = lambda vari : np.sqrt(np.sum(np.square(vari - cend[np.subtract(itra,1).astype(int)])))
            
            sepbeg = np.array(list(map(funbeg, temp))) # Distance between predicted position and starting position of candidates
            
            sepend = np.array(list(map(funend, temp))) # Distance between last position of fragment and starting position of candidates
          
            keep = np.where((sepbeg < radi[tidx]) & (sepend < maxiDisp))[0] # Reject separation distances that are too large
                                              
            # Store information ---------------------------------------------------------------------------------------------------------------------
            for ican in keep:
            
                fragStor.append(itra) # The fragment index starts at one
            
                framStor.append(tval)   
                
                distStor.append(sepbeg[ican])                                                              
                      
                jumpStor.append(np.add(tidx, 1)) # The jump index starts at one
                                       
                candStor.append(np.add(cand[ican], 1)) # The candidate index starts at one          
            
    tota = pd.DataFrame(list(zip(fragStor, candStor, framStor, jumpStor, distStor)),
                                 columns = ["Fragment", "Candidate", "Frame", "Jump", "Distance"]) 
    #################################################################################################################################################
       
    #################################################################################################################################################
    # Retrieve links between segments ---------------------------------------------------------------------------------------------------------------
    fragStor = [] # Initialize list
    jumpStor = [] # Initialize list
    candStor = [] # Initialize list
    
    for ifra in np.unique(tota["Fragment"]):
        
        mata = tota.query("Fragment == @ifra")
        # mata = tota.loc[tota["Fragment"] == ifra]
        
        mata.index = pd.RangeIndex(len(mata.index)) # Reset index to avoid wrong indexing
                        
        for ijmp in np.unique(mata["Jump"]):
                        
            matb = mata.query("Jump == @ijmp")
            
            matb.index = pd.RangeIndex(len(matb.index)) # Reset index to avoid wrong indexing
            
            vecb = np.ones(np.shape(matb)[0], dtype = bool)
            
            for rowb in matb.itertuples(): 
                                
                matc = tota.query("Candidate == @rowb.Candidate & Fragment != @ifra") # Do not include trajectory that is being processed
                
                matc.index = pd.RangeIndex(len(matc.index)) # Reset index to avoid wrong indexing
                                                              
                if matc.empty == False: # The candidate is shared between the fragment being processed and other fragments
                                    
                   for rowc in matc.itertuples(): 
                       
                       if rowc.Jump < ijmp: # The other fragment has time priority for this candidate over the fragment being processed
                           
                          vecb[rowb.Index] = False # Remove the candidate from the list of candidates
       
                       if (rowc.Jump == ijmp) and (rowc.Distance < rowb.Distance):
                                                                                                              
                          vecb[rowb.Index] = False # Remove the candidate from the list of candidates
                   
            matb = matb.loc[vecb == True]
            
            matb.index = pd.RangeIndex(len(matb.index)) # Reset index to avoid wrong indexing

            if matb.empty == False: # This fragment has priority over all the other fragments
                
               link = matb["Distance"].idxmin() # Returns the dataframe index, which here also corresponds to the row index
                              
               tota.loc[tota["Candidate"] == rowb.Candidate] = np.nan
               
               jumpStor.append(matb.loc[link, "Jump"])
               
               fragStor.append(matb.loc[link, "Fragment"])
               
               candStor.append(matb.loc[link, "Candidate"])
             
               break # Exit jump loop because we have found a link
    #################################################################################################################################################
               
    #################################################################################################################################################              
    # Connect couples -------------------------------------------------------------------------------------------------------------------------------                   
    jumpStor = np.array(jumpStor, dtype = int)  
    fragStor = np.array(fragStor, dtype = int)              
    candStor = np.array(candStor, dtype = int)   
    
    flag = np.zeros(np.size(fragStor), dtype = bool)
  
    mtot = [] # Initialize list to store trajectories
         
    for itra, duma in enumerate(fragStor):
    
        if flag[itra] == False:
            
           traj = [] # Initialize trajectory
                            
           fragTemp = fragStor[itra] # Define a temporary fragment index
           candTemp = candStor[itra] # Define a temporary candidate index
           
           jp = jumpStor[itra]
           
           ta = data[np.where(data[:,3] == fragTemp)[0],:]
           tb = data[np.where(data[:,3] == candTemp)[0],:]
                      
           se = addposi(ta, tb, jp) 
           
           traj.append(np.vstack((ta, se, tb)))
                                
           nfra = np.where(fragStor == candTemp)[0] # Check for next fragment   
        
           while np.size(nfra) != 0:
                        
                 fragTemp = candTemp # The former candidate is now the new fragment               
                 candTemp = candStor[nfra][0] 
                 
                 jp = jumpStor[nfra][0]
                 
                 ta = data[np.where(data[:,3] == fragTemp)[0],:]
                 tb = data[np.where(data[:,3] == candTemp)[0],:]
                      
                 se = addposi(ta, tb, jp)                
                 
                 traj.append(np.vstack((se, tb))) # Add next fragment to growing list
              
                 flag[nfra] = True # We will not process this couple anymore

                 nfra = np.where(fragStor == candTemp)[0]
                 
           traj = np.vstack((traj)) # Concatenate all fragments and added positions into one single trajectory
                                                              
           mtot.append(traj)
    #################################################################################################################################################    
           
    #################################################################################################################################################                  
    # Add back single fragments ---------------------------------------------------------------------------------------------------------------------
    disc = np.unique(np.concatenate((fragStor, candStor)))
    
    keep = np.setdiff1d(np.linspace(1, np.max(data[:,3]), num = int(np.max(data[:,3])), endpoint = True, dtype = int), disc)
    
    for itra in keep:
        
        traj = data[np.where(data[:,3] == itra)[0],:]
                
        mtot.append(traj)
    #################################################################################################################################################

    #################################################################################################################################################    
    # Renumber trajectories -------------------------------------------------------------------------------------------------------------------------
    for itra, duma in enumerate(mtot):
        
        mtot[itra][:,3] = itra # Broadcast new trajectory number
    #################################################################################################################################################    
       
    #################################################################################################################################################
    # Store in array with single precision ----------------------------------------------------------------------------------------------------------
    mtot = np.vstack(mtot).astype('float32') 

    # Save trajectories -----------------------------------------------------------------------------------------------------------------------------
    ntmp = "".join((name, '_Connecting.txt'))   
    np.savetxt(os.path.join(oupt, ntmp), mtot, delimiter = '\t')
    #################################################################################################################################################
    
#####################################################################################################################################################
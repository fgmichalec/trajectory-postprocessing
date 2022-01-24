# trajectory-postprocessing

A set of MATLAB and Python scripts to build and process trajectories using *ptv* and *xuap* files from OpenPTV and its postprocessing extension, or using ASCII files from LaVision DaVis. 

**Building**

*Buildingptv* and *Buildingxuap* reconstruct trajectories from *ptv* and *xuap* files without linking broken segments.

*BuildingDaVis* was used to study the behavior of motile particles in turbulence. This script loads ASCII files from LaVision DaVis, reorganizes the data, filters trajectories, calculates the velocity and acceleration, and interpolates the flow velocity at the positions of the motile particles using scattered interpolation. Input files are typically very large, so the script relies on sparce matrices to move data around during processing.  

**SubVolume** 

*SubVolume* discards positions that are outside the investigation volume or that belong to very short trajectories. These short trajectories typically result from impurities in the water, bright elements in the image background, or reflections on the surface of the aquarium. The script flags blocks of coordinates and identifies the starting and ending indices of each trajectory within these blocks before renumbering trajectory segments. This approach is faster than looping over trajectories.

**Connecting**

*Connecting* connects broken segments to reconstruct long trajectories. It is helpful in the case of active particles that deviate strongly from the flow streamlines or whose intermittent motion is challenging for most tracking algorithms. 

The script predicts the positions of the particle in time steps *n* = 1, 2, …, *m* after the end of a broken segment based on the velocity and acceleration of the particle in the previous time steps, looks for nearby candidates (in space and time) that are close to the predicted positions, and attempts to connect the two segments by adding equally spaced coordinates between the end of the broken segment and the starting position of the candidate segment. The code looks for candidates in time step *n* = 1 then *n* = 2 and so on until a connection is made or until it reaches *m*. 

The connection is made only if the candidate segment is the best candidate among all the other segments in the dataset. In other words, it is not enough for a candidate to be close to a broken segment (specificaly, close to the predicted positions), it should also not be closer to any other segment from the dataset both in time and in space. 

**Gluing**

*Gluingptv* reconstructs trajectories from *ptv* files and connects broken segments on the fly. 

After it reaches the end of a segment (flag for missing link in the second column of *ptv* files), the script predicts the positions of the particle in time steps *n* = 1, 2, …, *m* based on the velocity and acceleration of the particle in the previous time steps, looks for coordinates that are close to the predicted positions (in space and time), and attempts to connect the segment that is already reconstructed to the candidate position by adding equally spaced coordinates. The code looks for candidates in time step *n* = 1 then *n* = 2 and so on until a connection is made or until it reaches *m*. 

The connection is made if a candidate is found: the script does not evaluate whether the candidate position would be a better fit for another segment in the dataset. 

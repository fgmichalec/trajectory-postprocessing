% This function reconstructs trajectories without linking broken segments, 
% using xuap files (txt) from the postprocessing code associated with 
% OpenPTV.    

% The full sequence of input files is divided into several subsequences
% that reflect changes in experimental conditions during the recording. 
% To process the entire sequence without creating subsequences, set 
% seqNumb to one. 

% Output ------------------------------------------------------------------
% [nb, ts, fg, cx, cy, cz] in file coor (txt)
% [nb, ts, fg, ux, uy, uz] in file velo (txt)
% [nb, ts, fg, ax, ay, az] in file acce (txt)

% where: nb is the trajectory number
%        ts is the time step, starting at one and without padding 
%        fg is the flag for filtered positions
%        cx, cy and cz are the coordinates (m)
%        u is the velocity (m/s)
%        a is the acceleration (m/s^2)

% Time indexing -----------------------------------------------------------
% Indicate the first and last frames of the sequence used to build 
% trajectories in the name of the output file. The code will retrieve 
% this information to determine which files to process.  

% This also allows to keep the correspondence between the time step in 
% column ts, which always starts at one, and the number contained in the 
% name of the xuap files, in case the first frame of the sequence used to 
% build the trajectories does not correspond to the first frame of the 
% sequence processed with OpenPTV.   

% The time step ts does not reset between subsequences, but the trajectory 
% number nb does reset between subsequences.   

function Buildingxuap

close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seqNumb = 003;     % Number of subsequences
colNumb = 015;     % Number of columns in xuap files
maxPart = 6e3;     % Max number of particles per frame
minSave = 002;     % Min trajectory duration for saving (in frames) 

% The minimal duration is two frames no matter the value of minSave, 
% because the initialization of a new trajectory requires an existing
% link to the next time step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inpt = {'D:\Path to res folder for Recording 01\';
        'D:\Path to res folder for Recording 02\';
        'D:\Path to res folder for Recording 03\'};
       
oupt = {'D:\Output folder for Recording 01\'; 
        'D:\Output folder for Recording 02\'; 
        'D:\Output folder for Recording 03\'};         
                         
nmat = {'Recording_01_100001_101000_Buildingxuap';   % Output file name
        'Recording_02_100001_101000_Buildingxuap';   % Output file name
        'Recording_03_100001_101000_Buildingxuap'};  % Output file name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               
global partMat
global trajInd
            
% Directory loop ----------------------------------------------------------         
for idir = 1:size(inpt,1)
    
    strp = strfind(nmat{idir}, '_') ;
    
    % Adjust here depending on the name of your output file ---------------
    nbeg = str2double(nmat{idir}(plus(strp(2),1):minus(strp(3),1)));
    nend = str2double(nmat{idir}(plus(strp(3),1):minus(strp(4),1))); 
    
% The first time step indicated in the column ts of the output file starts 
% at one and not at nbeg, however this first time step corresponds to nbeg 
% and also to svec(1). 

% Subsequence loop --------------------------------------------------------
svec = round(linspace(nbeg, nend, plus(seqNumb,1))); 

for iseq = 1:minus(size(svec,2),1)
    
    clearvars -except seqNumb colNumb maxPart minSave partMat trajInd ...
                      inpt oupt nmat idir svec iseq  
                                                               
fbeg = svec(iseq); % First frame of the subsequence (padded format)
fend = svec(plus(iseq,1)); % Last frame of the subsequence (padded format)

trajInd = ones; % Initialize new trajectory number for each subsequence

% Build partMat matrix ----------------------------------------------------    
partMat = zeros(fend - fbeg + 1, maxPart, plus(colNumb,1)); 

partMat(:,:,1) = (-1); % This works both for ptv and xuap files
partMat(:,:,2) = (-2); % This works both for ptv and xuap files

posiVec = zeros(size(partMat,1),1);

for indx = fbeg:fend
      
    fidx = fopen(fullfile(inpt{idir}, sprintf('xuap.%d', indx)), 'r');
    
    temp = cell2mat(textscan(fidx, repmat('%f', [1, colNumb])));
    
    fclose(fidx);

    posiVec(indx - fbeg + 1) = size(temp,1);    
    partMat(indx - fbeg + 1, 1:size(temp,1), 1:colNumb) = temp;
        
    % We can initialize all the trajectories of the current time 
    % step if it is the first time step of the subsequence

    if isequal(indx, fbeg) == true
        
       partMat(indx - fbeg + 1, 1:size(temp,1), 1) = (-1);
       % This works both for ptv and xuap files
    end    
       
end % End of partMat building loop

% Open output files -------------------------------------------------------
ouptCoor = fopen(fullfile(oupt{idir},...
           sprintf('%s_%02i_coor.txt', nmat{idir}, iseq)), 'a');
ouptVelo = fopen(fullfile(oupt{idir},...
           sprintf('%s_%02i_velo.txt', nmat{idir}, iseq)), 'a');             
ouptAcce = fopen(fullfile(oupt{idir},...
           sprintf('%s_%02i_acce.txt', nmat{idir}, iseq)), 'a');     

% Time step loop ----------------------------------------------------------
for ifra = 1:size(partMat,1) 

    fprintf('Processing frame number %0d out of %d\n', ...
             ifra, plus(minus(fend,fbeg),1))
    
    % Particle loop -------------------------------------------------------  
    for ipar = 1:posiVec(ifra) 

    % New trajectory ------------------------------------------------------            
    if partMat(ifra, ipar, plus(colNumb,1)) == 0 && ...
       partMat(ifra, ipar, 2) > 0 
   
% For xuap files it is partMat(ifra, ipar, 2) > 0 because the ranking  
% of the particles starts at one and the flag for missing link is (-1) 
% For ptv files it is partMat(ifra, ipar, 2) >= 0 because the ranking  
% of the particles starts at zero and the flag for missing link is (-2).   
                 
% If no initialization of all the trajectories at the first frame of each 
% subsequence, do not add partMat(ifra, ipar, 1) == (-1) as a condition
% for the if statement since it prevents the processing of the first time  
% step. This is valid both for ptv and xuap files because we have populated
% the first and second pages of partMat with (-1) and (-2), respectively. 
       
       % Start and build trajectory ---------------------------------------
       traj = [ifra, ipar];
       traj = buildTraj(traj);  
              
       % Update partMat ---------------------------------------------------
       lini = sub2ind(size(partMat), traj(:,1), traj(:,2), ...
                      repmat(plus(colNumb,1), [size(traj,1), 1]));
                  
       partMat(lini) = ones; % We flag particles already in a trajectory
              
       % Format and save trajectory ---------------------------------------
       if size(traj,1) >= minSave
           
          data = formatTraj(traj, fbeg, svec);
          
          % Save [nb, ts, fg, cx, cy, cz] ---------------------------------
          fprintf(ouptCoor, '%i %i %i %f %f %f \n', data{1}');        
            
          % Save [nb, ts, fg, ux, uy, uz] ---------------------------------
          fprintf(ouptVelo, '%i %i %i %f %f %f \n', data{2}');        
            
          % Save [nb, ts, fg, ax, ay, az] ---------------------------------
          fprintf(ouptAcce, '%i %i %i %f %f %f \n', data{3}');  
           
       end % End of format and save trajectory section 
        
    end % End of new trajectory section
    
    end % End of particle loop
    
end % End of time step loop 

% Close output files ------------------------------------------------------
fclose(ouptCoor);
fclose(ouptVelo);
fclose(ouptAcce);

end % End of subsequence loop

end % End of directory loop

end % End of main function
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function traj = buildTraj(traj) 

global partMat

ifra = traj(end,1);
ipar = traj(end,2);
     
while (partMat(ifra, ipar, 2) > 0 && ifra < size(partMat,1))
      % The ranking of the particles starts at one in xuap files      
        
      ipar = partMat(ifra, ipar, 2); 

      ifra = plus(ifra,1);      
        
      traj = [traj; ifra, ipar]; %#ok<AGROW>
           
end
end % End of buildTraj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = formatTraj(traj, fbeg, svec)

global partMat
global trajInd 

% Coordinates -------------------------------------------------------------
cx = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(06, size(traj,1),1)));                   
cy = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(07, size(traj,1),1)));
cz = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(08, size(traj,1),1)));
                                      
% Velocity ----------------------------------------------------------------
ux = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(09, size(traj,1),1)));                   
uy = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(10, size(traj,1),1)));
uz = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(11, size(traj,1),1)));

% Acceleration ------------------------------------------------------------
ax = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(12, size(traj,1),1)));                   
ay = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(13, size(traj,1),1)));
az = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(14, size(traj,1),1)));
                 
% Flag for filtered positions ---------------------------------------------
fg = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(15, size(traj,1),1)));                 

% Export data -------------------------------------------------------------
data{1} = cat(2, trajInd * ones(size(traj,1),1), ...
                 traj(:,1) + fbeg - svec(1), fg, cx, cy, cz);
data{2} = cat(2, trajInd * ones(size(traj,1),1), ...
                 traj(:,1) + fbeg - svec(1), fg, ux, uy, uz);
data{3} = cat(2, trajInd * ones(size(traj,1),1), ...
                 traj(:,1) + fbeg - svec(1), fg, ax, ay, az);
    
trajInd = plus(trajInd,1);

end % End of formatTraj
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
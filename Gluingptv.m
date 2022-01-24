% This function reconstructs trajectories using ptv files (txt) from 
% OpenPTV. It connects broken segments on the fly, meaning that it  
% may connect a broken fragment to a candidate fragment even if the 
% candidate fragment is a better fit for another broken fragment 
% in the dataset.

% The full sequence of input files is divided into several subsequences
% that reflect changes in experimental conditions during the recording. 
% To process the entire sequence without creating subsequences, set 
% seqNumb to one. 

% Output ------------------------------------------------------------------
% Trajectory file (txt) containing [cx, cy, cz, nb, fg, ts, rk] 

% where: cx, cy and cz are the coordinates (mm)
%        nb is the trajectory number
%        fg is the flag for added position (1 if added and 0 otherwise)
%        ts is the time step, starting at one and without padding
%        rk is the row index of the particle in the corresponding ptv file 
%                                   (only for particles that are not added)

% Time indexing -----------------------------------------------------------
% Indicate the first and last frames of the sequence used to build 
% trajectories in the name of the output file. The code will retrieve 
% this information to determine which files to process.  

% This also allows to keep the correspondence between the time step in 
% column ts, which always starts at one, and the number contained in the 
% name of the ptv files, in case the first frame of the sequence used to 
% build the trajectories does not correspond to the first frame of the 
% sequence processed with OpenPTV.   

% The time step ts does not reset between subsequences, but the trajectory 
% number nb does reset between subsequences.    

function Gluingptv

close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seqNumb = 003;       % Number of subsequences
colNumb = 005;       % Number of columns in ptv files
maxPart = 6e3;       % Max number of particles per frame
minSave = 002;       % Min trajectory duration for saving (in frames)

% The minimal duration is two frames no matter the value of minSave, 
% because the initialization of a new trajectory requires an existing
% link to the next time step.

maxiDisp = 05; % Max jump displacement (mm)               
maxiDura = 10; % Max jump duration (frames) 
polySpan = 20; % Max duration of the regression window (frames)
sepaRadi = 02; % Radius of search volume around predicted positions (mm)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inpt = {'D:\Path to res folder for Recording 01\';
        'D:\Path to res folder for Recording 02\';
        'D:\Path to res folder for Recording 03\'};
       
oupt = {'D:\Output folder for Recording 01\'; 
        'D:\Output folder for Recording 02\'; 
        'D:\Output folder for Recording 03\'};         
                         
nmat = {'Recording_01_100001_101000_Gluingptv';  % Output file name
        'Recording_02_100001_101000_Gluingptv';  % Output file name
        'Recording_03_100001_101000_Gluingptv'}; % Output file name 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                   
global partMat
global posiVec
global posiOri
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
    
    clearvars -except ...
               seqNumb colNumb maxPart minSave ....
               maxiDisp maxiDura polySpan sepaRadi ...
               partMat posiVec posiOri trajInd ...
               inpt oupt nmat idir svec iseq                                         
                                                               
fbeg = svec(iseq); % First frame of the subsequence (padded format)
fend = svec(plus(iseq,1)); % Last frame of the subsequence (padded format)

trajInd = ones; % Initialize new trajectory number for each subsequence

% Build partMat matrix ----------------------------------------------------
partMat = zeros(fend - fbeg + 1, maxPart, plus(colNumb,2)); 

partMat(:,:,1) = (-1); 
partMat(:,:,2) = (-2);

posiOri = zeros(size(partMat,1),1); 
posiVec = zeros(size(partMat,1),1); 

% posiOri stores the number of particles in the ptv files before adding
% positions, for each time step. posiVec stores the number of particles
% in partMat, i.e. after adding positions, for each time step. 

for indx = fbeg:fend
    
    fidx = fopen(fullfile(inpt{idir}, sprintf('ptv_is.%d', indx)), 'r');   
    pidx = fscanf(fidx, '%i', [1 1]);
    
    temp = cell2mat(textscan(fidx, repmat('%f', [1, colNumb])));
    % temp = transpose(fscanf(fidx, '%i %i %f %f %f', [colNumb, pidx]));
    % temp = transpose(fscanf(fidx, '%f', [colNumb, pidx]));
    
    fclose(fidx); 
    
    posiOri(indx - fbeg + 1) = pidx;
    posiVec(indx - fbeg + 1) = pidx;

    partMat(indx - fbeg + 1, 1:pidx, 1:colNumb) = temp;
    
    % We can initialize all the trajectories of the current time 
    % step if it is the first time step of the subsequence

    if isequal(indx, fbeg) == true
        
       partMat(indx - fbeg + 1, 1:size(temp,1), 1) = (-1);
       % partMat(indx - fbeg + 1, 1:maxPart, 1) = (-1);      
       % partMat(indx - fbeg + 1, 1:pidx, 1) = (-1); 
       
    end 
     
end % End of partMat building loop

% Open output file --------------------------------------------------------
ouptTota = fopen(fullfile(oupt{idir},...
           sprintf('%s_%02i.txt', nmat{idir}, iseq)), 'a');              

% Time step loop ----------------------------------------------------------
for ifra = 1:size(partMat,1) 
    
    fprintf('Processing frame number %i\n', ifra)
      
    % Particle loop -------------------------------------------------------  
    for ipar = 1:posiOri(ifra) 
    
    % New trajectory ------------------------------------------------------    
    if partMat(ifra, ipar, plus(colNumb,1)) == 0 && ...
       partMat(ifra, ipar, 2) >= 0 
   
    % If no initialization of all the trajectories at the first frame of 
    % each subsequence, do not add partMat(ifra, ipar, 1) == (-1) as a 
    % condition for the if statement since it prevents the processing of  
    % the first time step. 
       
    traj = [ifra, ipar]; % Start new trajectory
    
    ifraTemp = traj(end,1); % Define a temporary time index
    iparTemp = traj(end,2); % Define a temporary particle index        
      
    % Build trajectory (if no gluing, use this section only) --------------
%     while partMat(ifraTemp, iparTemp, 2) >= 0 && ...
%           ifraTemp < size(partMat,1) 
%               
%           iparTemp = plus(partMat(ifraTemp, iparTemp, 2),1);
%           ifraTemp = plus(ifraTemp,1);
%                       
%           traj = cat(1, traj, [ifraTemp, iparTemp]); 
%                                            
%     end % End of building section

    % Useless to add partMat(ifraTemp, iparTemp, plus(colNumb,1)) == 0 as 
    % a condition in the while statement because there is no risk to use
    % the same particle twice.
    
    % Build trajectory (with gluing) --------------------------------------
    while partMat(ifraTemp, iparTemp, 2) >= 0 && ifraTemp < size(partMat,1) 
              
          iparTemp = plus(partMat(ifraTemp, iparTemp, 2),1);
          ifraTemp = plus(ifraTemp,1);
                      
          traj = cat(1, traj, [ifraTemp, iparTemp]); 
          
          % Check whether we need to enter the gluing section ------------- 
          if partMat(ifraTemp, iparTemp, 2) < 0 
 
             traj = glueTraj(traj, maxiDura,...       
                                   sepaRadi, polySpan, maxiDisp);
                               
             % We need to retrieve the time and particle indices so that 
             % we can keep building the trajectory after a successful 
             % gluing event or exit the while loop if the gluing is
             % not successful.
           
             ifraTemp = traj(end,1);         
             iparTemp = traj(end,2);
                               
          end % End of check section
                                           
    end % End of building section   
       
    % Update partMat ------------------------------------------------------
    lini = sub2ind(size(partMat), traj(:,1), traj(:,2), ...
                   repmat(plus(colNumb,1), [size(traj,1), 1]));

    partMat(lini) = ones; % We flag particles already in a trajectory 
              
    % Format and save trajectory ------------------------------------------
    if size(traj,1) >= minSave
        
       data = formatTraj(traj, fbeg, svec);
                      
       % Save [cx, cy, cz, nb, fg, ts, rk] --------------------------------
       fprintf(ouptTota, '%f %f %f %i %i %i %i \n', transpose(data)); 
         
    end % End of format and save trajectory section  

    end % End of new trajectory section
    
    end % End of particle loop
    
end % End of time step loop 

% Close output file -------------------------------------------------------
fclose(ouptTota);

end % End of subsequence loop

end % End of directory loop

end % End of main function

function traj = glueTraj(traj, maxiDura, sepaRadi, polySpan, maxiDisp)

global partMat
global posiVec
global posiOri
 
ifraTemp = traj(end,1);
iparTemp = traj(end,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine order of polynomial and size of regression window -------------
order = min([minus(size(traj,1),1), 2]);
polySpan = min([polySpan, size(traj,1)]);

% The order is one for trajectories that are two time step long and 
% two for longer trajectories no matter their duration. There are no
% trajectories shorter than two time steps.

% Determine size of jump --------------------------------------------------
if ifraTemp > minus(size(partMat,1), maxiDura)
   maxiDura = minus(size(partMat,1), ifraTemp);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coordinates at the end of the fragment (jump onset) ---------------------
coor = reshape(partMat(ifraTemp, iparTemp, 3:5), 1, 3);

% Indices of training examples --------------------------------------------
pidx = ((size(traj,1) - polySpan + 1):size(traj,1))';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Polynomial fitting ------------------------------------------------------
be = pidx(1) + traj(1,1) - 1; 
en = pidx(end) + traj(1,1) - 1;

coma = [partMat(sub2ind(size(partMat), ...
               traj(pidx,1), traj(pidx,2), repmat(03, size(pidx,1),1))),...
        partMat(sub2ind(size(partMat), ...
               traj(pidx,1), traj(pidx,2), repmat(04, size(pidx,1),1))),...
        partMat(sub2ind(size(partMat), ...
               traj(pidx,1), traj(pidx,2), repmat(05, size(pidx,1),1)))];      

dema = (repmat((be:en)', [1 plus(order,1)]))... % Design matrix
                          .^(repmat((0:order), [(en - be + 1) 1])); 
                        
coef = dema\coma;

% Estimated positions -----------------------------------------------------
p{1} = fliplr(coef(:,1)');
p{2} = fliplr(coef(:,2)');
p{3} = fliplr(coef(:,3)');

c{1} = polyval(p{1}, ...
               plus(traj(pidx(end)),1):plus(traj(pidx(end)),maxiDura));
c{2} = polyval(p{2}, ...
               plus(traj(pidx(end)),1):plus(traj(pidx(end)),maxiDura));
c{3} = polyval(p{3},...
               plus(traj(pidx(end)),1):plus(traj(pidx(end)),maxiDura));          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for jind = 1:maxiDura
        
    % List all candidate points -------------------------------------------
    x = transpose(...
        partMat(plus(ifraTemp,jind), 1:posiOri(plus(ifraTemp,jind)), 03));
    y = transpose(...
        partMat(plus(ifraTemp,jind), 1:posiOri(plus(ifraTemp,jind)), 04));
    z = transpose(...
        partMat(plus(ifraTemp,jind), 1:posiOri(plus(ifraTemp,jind)), 05));
    
    % Calculate distances from coordinates at onset -----------------------
    dtmp = sqrt(sum((bsxfun(@minus, ...
                     [x, y, z], [coor(1), coor(2), coor(3)])).^2, 2));  
                                
    % Discard points ------------------------------------------------------
    disc = transpose(partMat(plus(ifraTemp,jind), ...
                     1:posiOri(plus(ifraTemp,jind)), 06)) == 1     | ... 
           transpose(partMat(plus(ifraTemp,jind), ...
                     1:posiOri(plus(ifraTemp,jind)), 01)) ~= (-1)  | ... 
           transpose(partMat(plus(ifraTemp,jind), ...
                     1:posiOri(plus(ifraTemp,jind)), 02)) == (-2)  | ... 
           dtmp > maxiDisp;
          
    % We discard points that do not start a fragment, that are already  
    % part of a previously reconstructed trajectory, that do not have 
    % a valid link to the next time step, and that are too far from 
    % the coordinates at onset.
    
    [x(disc), y(disc), z(disc)] = deal(NaN); 
    
    % Keep the NaNs to maintain the indexing
          
    % Calculate distance from predicted position --------------------------
    dtmp = sqrt(sum((bsxfun(@minus, ...
                     [x, y, z], ...
                     [c{1}(jind), c{2}(jind), c{3}(jind)])).^2, 2));
       
    [disp, indx] = nanmin(dtmp); % Minimal separation distance
    
    % Link to next fragment -----------------------------------------------
    if ~isnan(disp) && disp < sepaRadi 
       
        time = transpose((ifraTemp + 1):(ifraTemp + jind - 1));        
        posi = plus(posiVec(time),1);

        traj = cat(1, traj, [time, posi], ...
                            [plus(ifraTemp,jind), indx]);

        posiVec(time) = posi;
         
        % Add positions if the jump is over one time step -----------------
        if ~isempty(time) 

           incr = rdivide(bsxfun(@minus, ...
                          [x(indx), y(indx), z(indx)], coor), jind);
                            
           temp = repmat(coor, minus(jind,1), 1) + ...
                  times(transpose(1:minus(jind,1)), incr);

           % Update ptv array (coordinates) -------------------------------              
           partMat(sub2ind(size(partMat), ...
                   time, posi, repmat(03, size(time,1),1))) = temp(:,1);
           partMat(sub2ind(size(partMat), ...
                   time, posi, repmat(04, size(time,1),1))) = temp(:,2);
           partMat(sub2ind(size(partMat), ...
                   time, posi, repmat(05, size(time,1),1))) = temp(:,3);        

           % Update ptv array (flag for added positions) ------------------
           partMat(sub2ind(size(partMat), ...
                   time, posi, repmat(07, size(time,1),1))) = 1;                         

        end % End of add positions and update section 
        
        break % We have found a link and can exit the loop
 
    end % End of link to next fragment section   
  
end % End of jump loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
end % End of glueTraj
          
function data = formatTraj(traj, fbeg, svec)

global partMat
global trajInd 

% Coordinates -------------------------------------------------------------
cx = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(03, size(traj,1),1)));                   
cy = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(04, size(traj,1),1)));
cz = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(05, size(traj,1),1))); 
                    
% Flag for added positions ------------------------------------------------                    
fg = partMat(sub2ind(size(partMat), ...
                     traj(:,1), traj(:,2), repmat(07, size(traj,1),1)));
                    
% Row index in the corresponding ptv file ---------------------------------
rk = traj(:,2); % The row index starts at one instead of zero
rk(fg == 1) = NaN; % NaN for added positions 

% Export data -------------------------------------------------------------
data = [cx, cy, cz,...                          % Coordinates                
        trajInd * ones(size(traj,1),1),...      % nb
        fg,...                                  % fg 
        traj(:,1) + fbeg - svec(1), ...         % ts
        rk];                                    % rk

% Export data (if no gluing, use this section) ----------------------------
% data = [cx, cy, cz,...                          % Coordinates                
%         trajInd * ones(size(traj,1),1),...      % nb
%         zeros(size(traj,1),1),...               % fg 
%         traj(:,1) + fbeg - svec(1), ...         % ts
%         traj(:,2)];                             % rk
    
trajInd = plus(trajInd,1);

end % End of formatTraj

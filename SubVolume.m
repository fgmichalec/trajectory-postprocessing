% This script discards positions that are outside the investigation volume
% or that belong to short trajectories. These short trajectories typically
% result from impurities in the water, bright elements in the image 
% background, or reflections on the surface of the aquarium. 
   
% The code flags blocks of coordinates and identifies the starting and 
% ending indices of each trajectory within these blocks. This approach 
% is much faster than looping over trajectories.

% Input and output --------------------------------------------------------
% File (txt) containing [cx, cy, cz, nb, fg, ts, rk] 

% where: cx, cy and cz are the coordinates (mm)
%        nb is the trajectory number
%        fg is the flag for added positions (1 if added and 0 otherwise)
%        ts is the time step, starting at one and without padding
%        rk is the row index of the particle in the corresponding ptv file
%                                  (only for particles that are not added)

function SubVolume

close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inpt = {'D:\Input folder for Recording 01\';
        'D:\Input folder for Recording 02\';
        'D:\Input folder for Recording 03\'};
       
oupt = {'D:\Output folder for Recording 01\'; 
        'D:\Output folder for Recording 02\'; 
        'D:\Output folder for Recording 03\'};  

nmat = {'Recording_01_100001_101000_Buildingptv'; % Input file name
        'Recording_02_100001_101000_Buildingptv'; % Input file name
        'Recording_03_100001_101000_Buildingptv'; % Input file name
        };        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
leng = 05; % Minimal length (in mm) 
dura = 20; % Minimal duration (in frames) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Limits of the investigation volume (in mm) ------------------------------
volu{1} = [10, 90]; % X axis (OpenPTV reference frame)
volu{2} = [10, 95]; % Y axis (OpenPTV reference frame)
volu{3} = [15, 85]; % Z axis (OpenPTV reference frame)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% File loop ---------------------------------------------------------------
for fidx = 1:numel(nmat) 
    clearvars -except inpt oupt nmat fidx leng dura
    
    fprintf('Processing file <strong>%s</strong>\n', nmat{fidx})
  
    data = load(fullfile(inpt{fidx}, sprintf('%s.txt', nmat{fidx})));
                        
    % Flag coordinates outside the investigation volume -------------------
    disc = data(:,1) < volu{1}(1) | data(:,1) > volu{1}(2) |...
           data(:,2) < volu{2}(1) | data(:,2) > volu{2}(2) |... 
           data(:,3) < volu{3}(1) | data(:,3) > volu{3}(2);                   
       
    % Extract indices of coordinates inside the investigation volume ------    
    ivec = diff(cat(1, 0, ~disc, 0));
    
    stti = find(ivec == 1); 
    stpi = minus(find(ivec == (-1)),1); 
          
    keep = cell(0); 
    
    % Keep or reject based on length and duration -------------------------
    for segi = 1:size(stti,1) 
        
        tvec = find(diff(data(stti(segi):stpi(segi), 4))); 
        
        bidx = cat(1, stti(segi), (tvec + stti(segi)));               
        eidx = cat(1, (tvec + stti(segi) - 1), stpi(segi));
                 
        for posi = 1:size(bidx,1) 

        if numel(bidx(posi):eidx(posi)) >= dura && ... 
           ismember(0, (max(data(bidx(posi):eidx(posi), [1 2 3])) - ...
                        min(data(bidx(posi):eidx(posi), [1 2 3])) < ...
                        leng)) == true
                    
           % Fragment is longer than threshold in space in at least 
           % one dimension.            
                    
           keep{plus(numel(keep),1)} = bidx(posi):eidx(posi); 
           
        end % End of length and duration condition
        
        end % End of trajectory section
        
    end % End of segment section
         
% New trajectory numbering ------------------------------------------------
for coun = 1:numel(keep)
    data(keep{coun},4) = times(ones(numel(keep{coun}),1), coun);
end

% Write to disk -----------------------------------------------------------
temp = data(cat(2, keep{:}),:);  

save(fullfile(oupt{fidx}, sprintf('%s_SubVolume.txt', nmat{fidx})), ...
                          'temp', '-ascii', '-single', '-tabs')

end % End of file loop

end % End of function

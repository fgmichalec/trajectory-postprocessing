% This script filters copepod trajectories reconstructed with LaVision 
% DaVis, calculates the velocity and acceleration, and interpolates the 
% flow velocity at the position of the copepods using cua files.

% LaVision DaVis files can be very large depending on the recording 
% duration and particle density, and the information is organized along
% time steps, which is not convenient when studying particle motion in  
% a Lagrangian framework. This script relies on sparce matrices to 
% organize data along trajectories and to efficiently process these 
% large files.  

% Input -------------------------------------------------------------------
% (a) Files from LaVision DaVis in ASCII format (dat), containing 
%     the raw coordinates, velocity and acceleration of the copepods,  
%     and the trajectory number. 

% (b) cua files (txt) containing the filtered coordinates (columns  
%     1 to 3), the filtered velocity (columns 4 to 6) and the filtered 
%     acceleration (columns 7 to 9) of the tracer particles. 

% Output ------------------------------------------------------------------
% The structure array mtot(n).field (mat) with n the trajectory number. 

% (a) mtot(n).step is a vector containing the time step within the  
%     sequence, with origine at one.

% (b) mtot(n).coor is a matrix containing the filtered coordinates of 
%     the copepod for this trajectory.

% (c) mtot(n).velo and mtot(n).acce are matrices containing the three   
%     components of the filtered velocity and acceleration of the 
%     copepod, respectively.

% (d) mtot(n).flow is a matrix containing the three components of the 
%     flow velocity interpolated at the copepod position.

function BuildingDaVis

close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clus = parcluster('local'); % Adjust depending on local or cluster
% pool = parpool(clus, 48);   % Adjust depending on local or cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inpt.cope = '/cluster/scratch/Recording 01/'; % Path to trajectory file 
% inpt.cope = '/cluster/scratch/Recording 02/'; % Path to trajectory file
% inpt.cope = '/cluster/scratch/Recording 03/'; % Path to trajectory file

% inpt.trac = '/cluster/scratch/Recording 01/'; % Path to cua files
% inpt.trac = '/cluster/scratch/Recording 02/'; % Path to cua files
% inpt.trac = '/cluster/scratch/Recording 03/'; % Path to cua files

% oupt = '/cluster/scratch/Recording 01/'; % Path to output file
% oupt = '/cluster/scratch/Recording 02/'; % Path to output file
% oupt = '/cluster/scratch/Recording 03/'; % Path to output file

% name = 'Recording_01.dat';  % Name of input trajectory file
% name = 'Recording_02.dat';  % Name of input trajectory file
% name = 'Recording_03.dat';  % Name of input trajectory file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters --------------------------------------------------------------
ordr = 003;   % Order of Savitzky-Golay filter
span = 011;   % Span of Savitzky-Golay filter
freq = 300;   % Recording frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load input data ---------------------------------------------------------
flid.inpt = fopen(fullfile(inpt.cope, name), 'r');

% Retrieve headers TITLE and VARIABLES ------------------------------------
TITLE = textscan(flid.inpt, 'TITLE = %s', 1, 'Delimiter', '\n');
VARIABLES = textscan(flid.inpt, 'VARIABLES = %s', 1, 'Delimiter', '\n');

TITLE = cell2mat(TITLE{1});
VARIABLES = cell2mat(VARIABLES{1});

TITLE = cell2mat(regexp(TITLE, '"(.*?)"', 'tokens', 'once')); 

VARIABLES = cellfun(@(a) a{1}, regexp(VARIABLES, '"(.*?)"', 'tokens'),...
                    'UniformOutput', false); % Several instances
                
% TITLE is a string, VARIABLES is a cell array of strings
                
% Extract each block, one block corresponds to a time step ----------------
coun = ones; % Initialize counter

while (~feof(flid.inpt)) 
        
    % Retrieve other fields -----------------------------------------------
    a = textscan(flid.inpt, '%s', 1, 'Delimiter', '\n');
    b = textscan(flid.inpt, '%s', 1, 'Delimiter', '\n');
    c = textscan(flid.inpt, '%s', 1, 'Delimiter', '\n');   
    d = textscan(flid.inpt, '%s', 1, 'Delimiter', '\n'); 
    
    a = cell2mat(a{1});
    b = cell2mat(b{1});
    c = cell2mat(c{1});   
    d = cell2mat(d{1});  
    
    T = cell2mat(regexp(a, '"(.*?)"', 'tokens', 'once')); 
      
    temp = regexp(b, ...
           'STRANDID=(?<STRANDID>.*), SOLUTIONTIME=(?<SOLUTIONTIME>.*)', ...
           'names', 'once');

    STRANDID = str2double(temp.STRANDID);    
    SOLUTIONTIME = str2double(temp.SOLUTIONTIME);
      
    % T is a string, STRANDID and SOLUTIONTIME are scalars
    
    temp = regexp(c, ...
           'I=(?<I>.*), J=(?<J>.*), K=(?<K>.*), ZONETYPE = (?<ZONETYPE>.*)', ...
           'names', 'once'); % Careful as the spacing is not consistent
       
    I = str2double(temp.I);
    J = str2double(temp.J);
    K = str2double(temp.K);
    
    ZONETYPE = temp.ZONETYPE;
    
    DATAPACKING = char(regexp(d, 'DATAPACKING = (.*)', 'tokens', 'once'));
    
    % I, J and K are scalars, ZONETYPE and DATAPACKING are strings
        
    % Retrieve data points ------------------------------------------------
    form = repmat('%f', [1, 13]); % Number of columns in DaVis file

    DATA = cell2mat(...
           textscan(flid.inpt, form, plus(I,1), 'CollectOutput', true));
       
    % Store information ---------------------------------------------------   
    A(coun).TITLE = TITLE; %#ok<*AGROW>
    A(coun).VARIABLES = VARIABLES;  
    A(coun).T = T;        
    A(coun).STRANDID = STRANDID;
    A(coun).SOLUTIONTIME = SOLUTIONTIME;    
    A(coun).I = I;    
    A(coun).J = J;    
    A(coun).K = K; 
    A(coun).ZONETYPE = ZONETYPE;     
    A(coun).DATAPACKING = DATAPACKING; 
    A(coun).DATA = DATA;

    coun = plus(coun,1);

end

fclose(flid.inpt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Format input data into sparse matrices ----------------------------------
n = numel(A); % Number of time steps

% time = transpose({A.SOLUTIONTIME});

% time is a cell array that contains the time step (s), with the first
% frame starting at zero. We do not use this information. 

% fram = regexp(transpose({A.T}), 'Snapshot (.*)', 'tokens', 'once'); 
% fram = cellfun(@(var) var{1}, fram, 'UniformOutput', false);

% fram is a cell array that contains the frame number as a string, in 
% padded format and starting at zero. We do not use this information. 

[x, y, z, p] = deal(cell(n,1)); % Keep row vector format

parfor indx = 1:n
        
   temp = A(indx).DATA;
   
   % p{indx} = single(cat(2, repmat(indx, [size(temp,1), 1]), temp(:,9)));
   p{indx} = cat(2, repmat(indx, [size(temp,1), 1]), temp(:,9));  
           
   x{indx} = temp(:,1); % Keep double precision      
   y{indx} = temp(:,2); % Keep double precision                     
   z{indx} = temp(:,3); % Keep double precision 
   
   % x{indx} = single(temp(:,1)); % Single precision      
   % y{indx} = single(temp(:,2)); % Single precision                     
   % z{indx} = single(temp(:,3)); % Single precision                        
end

clearvars A

memo{1} = whos('p');
memo{2} = whos('x');
memo{3} = whos('y');
memo{4} = whos('z');

fprintf('Variable p uses %i bytes of memory\n', memo{1}.bytes)
fprintf('Variable x uses %i bytes of memory\n', memo{2}.bytes)
fprintf('Variable y uses %i bytes of memory\n', memo{3}.bytes)
fprintf('Variable z uses %i bytes of memory\n', memo{4}.bytes)

p = cell2mat(p); 
x = cell2mat(x); 
y = cell2mat(y); 
z = cell2mat(z); 

x(x == 0) = eps(0); % We will use sparse matrices
y(y == 0) = eps(0); % We will use sparse matrices
z(z == 0) = eps(0); % We will use sparse matrices

m = numel(unique(p(:,2))); % Total number of trajectories

fprintf('Variable m is equal to %i\n', m)
fprintf('Variable n is equal to %i\n', n)

x = sparse(p(:,2) + 1, p(:,1), x, m, n);
y = sparse(p(:,2) + 1, p(:,1), y, m, n);
z = sparse(p(:,2) + 1, p(:,1), z, m, n);

% x = sparse(double(p(:,2)) + 1, double(p(:,1)), double(x), m, n);
% y = sparse(double(p(:,2)) + 1, double(p(:,1)), double(y), m, n);
% z = sparse(double(p(:,2)) + 1, double(p(:,1)), double(z), m, n);

% sparse does not accept single as input variable

clearvars p

memo{1} = whos('x');
memo{2} = whos('y');
memo{3} = whos('z');

fprintf('Variable x (sparse) uses %i bytes of memory\n', memo{1}.bytes)
fprintf('Variable y (sparse) uses %i bytes of memory\n', memo{2}.bytes)
fprintf('Variable z (sparse) uses %i bytes of memory\n', memo{3}.bytes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert sparse matrices to cell arrays ----------------------------------
[ro, co, vx] = find(x); 
[~ , ~ , vy] = find(y);
[~ , ~ , vz] = find(z);

clearvars x y z

pnew = accumarray(ro, co, [], @(var) {var}); % Unsorted
xnew = accumarray(ro, vx, [], @(var) {var}); % Unsorted
ynew = accumarray(ro, vy, [], @(var) {var}); % Unsorted
znew = accumarray(ro, vz, [], @(var) {var}); % Unsorted

clearvars ro co vx vy vz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create output file ------------------------------------------------------
mtot = repmat(struct('step', [], 'coor', [], ...
                     'velo', [], 'acce', [], 'flow', []), [m, 1]);   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtering and calculation of derivatives --------------------------------
parfor indx = 1:m
    
[~, i] = sort(pnew{indx});

cx = xnew{indx}(i); % Retrieve in cell array
cy = ynew{indx}(i); % Retrieve in cell array
cz = znew{indx}(i); % Retrieve in cell array

a = pnew{indx}(i);  % Retrieve in cell array

% [~, a, cx] = find(x(indx,:)); % Directly in sparse matrix
% [~, ~, cy] = find(y(indx,:)); % Directly in sparse matrix
% [~, ~, cz] = find(z(indx,:)); % Directly in sparse matrix

[~, g] = sgolay(ordr, span); 
halfsp = rdivide(minus(span,1),2);

coor = [conv(cx, factorial(0) / (-1/freq)^0 * g(:, 0 + 1), 'same'),...
        conv(cy, factorial(0) / (-1/freq)^0 * g(:, 0 + 1), 'same'),...
        conv(cz, factorial(0) / (-1/freq)^0 * g(:, 0 + 1), 'same')];

velo = [conv(cx, factorial(1) / (-1/freq)^1 * g(:, 1 + 1), 'same'),...
        conv(cy, factorial(1) / (-1/freq)^1 * g(:, 1 + 1), 'same'),...
        conv(cz, factorial(1) / (-1/freq)^1 * g(:, 1 + 1), 'same')];        

acce = [conv(cx, factorial(2) / (-1/freq)^2 * g(:, 2 + 1), 'same'),...
        conv(cy, factorial(2) / (-1/freq)^2 * g(:, 2 + 1), 'same'),...
        conv(cz, factorial(2) / (-1/freq)^2 * g(:, 2 + 1), 'same')];

% coor = [conv(cx', factorial(0) / (-1/freq)^0 * g(:, 0 + 1), 'same'),...
%         conv(cy', factorial(0) / (-1/freq)^0 * g(:, 0 + 1), 'same'),...
%         conv(cz', factorial(0) / (-1/freq)^0 * g(:, 0 + 1), 'same')];
% 
% velo = [conv(cx', factorial(1) / (-1/freq)^1 * g(:, 1 + 1), 'same'),...
%         conv(cy', factorial(1) / (-1/freq)^1 * g(:, 1 + 1), 'same'),...
%         conv(cz', factorial(1) / (-1/freq)^1 * g(:, 1 + 1), 'same')];        
% 
% acce = [conv(cx', factorial(2) / (-1/freq)^2 * g(:, 2 + 1), 'same'),...
%         conv(cy', factorial(2) / (-1/freq)^2 * g(:, 2 + 1), 'same'),...
%         conv(cz', factorial(2) / (-1/freq)^2 * g(:, 2 + 1), 'same')];    

% If using conv it is necessary to remove positions at the beginning
% and end of the trajectory because they are not filtered.

% Indices of entries to keep ----------------------------------------------
keep = plus(halfsp,1):minus(numel(cx), halfsp);    

% Store data --------------------------------------------------------------
if isempty(keep) == false

    % mtot(indx).step = transpose(a(keep)); % Time step (in frame)
    mtot(indx).step = a(keep); % Time step (in frame)

    mtot(indx).coor = coor(keep,:); % Coordinates 
    mtot(indx).velo = velo(keep,:); % Velocity (components)                                            
    mtot(indx).acce = acce(keep,:); % Acceleration (components)   
end
       
end

mtot = mtot(~cellfun(@isempty, {mtot.step})); % Remove empty structures

clearvars pnew xnew ynew znew
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store filtered coordinates into sparse matrices -------------------------
a = num2cell((1:numel(mtot))');
b = arrayfun(@(var) numel(var.step), mtot, 'UniformOutput', false);

numb = cell2mat(...
       cellfun(@(i,j) repmat(i, [j 1]), a, b, 'UniformOutput', false));

smat = cat(1, mtot(:).step);
cmat = cat(1, mtot(:).coor);

cmat(cmat(:,1) == 0,1) = eps(0); % We will use sparse matrices
cmat(cmat(:,2) == 0,2) = eps(0); % We will use sparse matrices
cmat(cmat(:,3) == 0,3) = eps(0); % We will use sparse matrices

x = sparse(numb, smat, cmat(:,1), numel(mtot), n);
y = sparse(numb, smat, cmat(:,2), numel(mtot), n);
z = sparse(numb, smat, cmat(:,3), numel(mtot), n);

clearvars a b numb smat cmat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert sparse matrices to cell arrays with respect to time -------------
[ro, co, vx] = find(x); 
[~ , ~ , vy] = find(y);
[~ , ~ , vz] = find(z);

clearvars x y z

pnew = accumarray(co, ro, [n, 1], @(var) {var}); % Unsorted
xnew = accumarray(co, vx, [n, 1], @(var) {var}); % Unsorted
ynew = accumarray(co, vy, [n, 1], @(var) {var}); % Unsorted
znew = accumarray(co, vz, [n, 1], @(var) {var}); % Unsorted

clearvars ro co vx vy vz

% There are as many cells as there are time steps. Each cell of pnew
% contains the trajectory indices in row format. Each cell of xnew, ynew
% and znew contains the filtered coordinates along the corresponding 
% dimension, also in row format. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate flow velocity -----------------------------------------------
fmat = dir(fullfile(inpt.trac, 'cua*')); % Sorted because there is padding

[flox, floy, floz] = deal(cell(n,1));

parfor indx = 1:n 
    
    fprintf('Processing file %s\n', fmat(indx).name)

    fidx = fopen(fullfile(fmat(indx).folder, fmat(indx).name), 'r');
    data = cell2mat(textscan(fidx, repmat('%f', [1, 9])));
    
    fclose(fidx);
    
    [flox{indx}, floy{indx}, ... % In case there is no tracer data
                 floz{indx}] = deal(NaN(numel(pnew{indx}),1)); 
  
    if isempty(data) == false % There are tracer particles
    
       % Scattered interpolant --------------------------------------------                    
       fx = scatteredInterpolant(data(:,1), data(:,2), data(:,3), data(:,4),...
                                 'natural', 'none');
       fy = scatteredInterpolant(data(:,1), data(:,2), data(:,3), data(:,5),...
                                 'natural', 'none');                      
       fz = scatteredInterpolant(data(:,1), data(:,2), data(:,3), data(:,6),...
                                 'natural', 'none'); 
                                                                                      
       % Estimate velocity at multiple query points -----------------------
       flox{indx} = fx(xnew{indx}, ynew{indx}, znew{indx}); 
       floy{indx} = fy(xnew{indx}, ynew{indx}, znew{indx}); 
       floz{indx} = fz(xnew{indx}, ynew{indx}, znew{indx}); 

% flox{indx}, floy{indx} and floz{indx} contain the flow velocity 
% interpolated at the position of the multiple copepods in this time 
% step, in row  format and in the order indicated in pnew. 

% Cells may contain NaNs if the flow field cannot be interpolated because
% a copepod is outside the convex hull formed by the point cloud of tracer
% particles or because there is no tracer data at this given time step to
% construct a scattered interpolant (beginning and end of the sequence). 
% Cells may even be empty if there are no copepods in this time step.     

    end % End of isempty condition

end % End of parfor loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store flow velocity data into sparse matrices ---------------------------
a = num2cell((1:n)'); % num2cell((1:numel(pnew))');
b = cellfun(@(var) numel(var), pnew, 'UniformOutput', false);

smat = cell2mat(...
       cellfun(@(i,j) repmat(i, [j 1]), a, b, 'UniformOutput', false));

numb = cat(1, pnew{:});

flox = cat(1, flox{:});
floy = cat(1, floy{:});
floz = cat(1, floz{:});

flox(flox == 0) = eps(0); % We will use sparse matrices
floy(floy == 0) = eps(0); % We will use sparse matrices
floz(floz == 0) = eps(0); % We will use sparse matrices

flox = sparse(numb, smat, flox, numel(mtot), n);
floy = sparse(numb, smat, floy, numel(mtot), n);
floz = sparse(numb, smat, floz, numel(mtot), n);

clearvars a b numb smat
clearvars pnew xnew ynew znew
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store flow velocity data into cell arrays -------------------------------
[ro, co, vx] = find(flox); 
[~ , ~ , vy] = find(floy);
[~ , ~ , vz] = find(floz);

clearvars flox floy floz

pnew = accumarray(ro, co, [], @(var) {var}); % Unsorted
xnew = accumarray(ro, vx, [], @(var) {var}); % Unsorted
ynew = accumarray(ro, vy, [], @(var) {var}); % Unsorted
znew = accumarray(ro, vz, [], @(var) {var}); % Unsorted

clearvars ro co vx vy vz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store flow velocity data in output file ---------------------------------
parfor indx = 1:numel(mtot)
    
    [~, i] = sort(pnew{indx});

    mtot(indx).flow = cat(2, xnew{indx}(i),... % Retrieve in cell array
                             ynew{indx}(i),... % Retrieve in cell array
                             znew{indx}(i));   % Retrieve in cell array
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save output file --------------------------------------------------------
clearvars pnew xnew ynew znew
 
save(fullfile(oupt, ...
     sprintf('%s_BuildingDaVis.mat', erase(name, '.dat'))), ...
     'mtot', '-mat', '-v7.3')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% pool.delete()

end % End of function

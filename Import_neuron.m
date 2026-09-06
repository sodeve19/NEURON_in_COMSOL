
% File: Import_neuron.m
% Authors: Rana Jaylani, Sian McConchie

% Description:
% This script uses MATLAB LiveLink to create a time-dependent COMSOL model
% Imports neuron at specified location and angle into a COMSOL model

% Requirements:
%   - MATLAB with COMSOL LiveLink
%   - COMSOL Multiphysics (Tested on 6.2, 6.3 and 6.4)
%   - segment\_currents\_originalGold .csv
%   - segment\_morphology\_originalGold.csv

% Usage:
%   1. Update the file name and save path as suitable
%   2. Open the script in MATLAB LiveLink
%   3. Click Run

%% Import Single Neuron
ModelName = 'Model_Mircowire_Array'; % -> CHANGE NAME TO MODEL REQUIRED
% Neuron Configuration files
morph_csvPath = 'segment_morphology_originalGold.csv';
cur_csvPath = 'segment_currents_originalGold.csv';

%% Load Model
model = mphopen(ModelName);
compTag = 'comp1';
geomTag = 'geom1';
% Ensure component and geometry exist
try model.component(compTag); catch, model.component.create(compTag,true); end
try model.component(compTag).geom(geomTag); catch, model.component(compTag).geom.create(geomTag,3); end
geom = model.component(compTag).geom(geomTag);

%% Read Morphology CSV
scaleFactor = 1e-6;      % Coordinates are in µm
% create offsets
x_os = 0; % -> INSERT REQUIRED CELL POSITIONING IN um
y_os = 0;
z_os = 0;

% Rotations (degrees)
rotX_deg = 0;   % pitch
rotY_deg = 90;   % yaw
rotZ_deg = 0;   % roll

% Rotation order 
rotOrder = "ZYX";

% get segment start points
morphT   = readtable(morph_csvPath);
xs   = (morphT.xs(:) + x_os) * scaleFactor;
ys   = (morphT.ys(:) + y_os) * scaleFactor;
zs   = (morphT.zs(:) + z_os) * scaleFactor;

% get segment end points
xe   = (morphT.xe(:) + x_os) * scaleFactor;
ye   = (morphT.ye(:) + y_os ) * scaleFactor;
ze   = (morphT.ze(:) + z_os) * scaleFactor;

% Build point lists 
P1 = [xs, ys, zs];  
P2 = [xe, ye, ze];  

% rotate about global origin (0,0,0)
origin = [0 0 0];

% Apply rotation
P1r = rotateXYZ(P1, rotX_deg, rotY_deg, rotZ_deg, origin, rotOrder);
P2r = rotateXYZ(P2, rotX_deg, rotY_deg, rotZ_deg, origin, rotOrder);

% Unpack back into arrays
xs = P1r(:,1); ys = P1r(:,2); zs = P1r(:,3);
xe = P2r(:,1); ye = P2r(:,2); ze = P2r(:,3);

% get segment names and ids
segid = morphT.seg_id(:); 

%% Define interpolation function: I_Seg(t, seg)
tic
disp('Loading interpolation functions...\n');

n_segs = height(morphT);          % number of segments 
firstValueCol = 2;               
MaxSegPerFunc = 100;         % COMSOL limit: keep < 1000 columns total
nFuncs = ceil(n_segs / MaxSegPerFunc);
FuncStart = 2;
FuncEnd = MaxSegPerFunc+1;
fcount = 1;

curT = readtable(cur_csvPath);

% Multiply all columns except the time by 2.09
curT{:, 2:end} = curT{:, 2:end} * 2.09;

for k = 1:nFuncs
    CurFunc = curT(:, [1,FuncStart:FuncEnd]);
    FuncFileName = sprintf('cur_func_%d.csv',k);
    writetable(CurFunc,FuncFileName);
    WidthT= width(CurFunc);
    FuncLab = sprintf('int%d',k);

    try
        model.func.remove(FuncLab);  % remove if exists
    catch
    end

    model.func.create(FuncLab, 'Interpolation');
    model.func(FuncLab).set('source', 'file');
    model.func(FuncLab).set('filename', FuncFileName);
    n_segs = height(morphT);
    model.func(FuncLab).setIndex('argunit', 'ms', 0);
    indexcount = 0;

    for n = 2:WidthT
        col = sprintf('col%d',n);
        functag = sprintf('segi%d',fcount);
        model.func(FuncLab).setEntry('columnType', col, 'value');
        model.func(FuncLab).setEntry('funcnames', col, functag);
        model.func(FuncLab).setIndex('fununit', 'nA/um', indexcount);
        fcount = fcount + 1;
        indexcount = indexcount + 1;
        fprintf('(%d/%d) Function loaded for seg %d \n',fcount, n_segs, fcount);
    end
    model.func(FuncLab).set('interp', 'linear');
    model.func(FuncLab).set('extrap', 'const');  % hold outside range (or 'zero')

    FuncStart = FuncEnd+1;
    if FuncStart+MaxSegPerFunc-2 < n_segs
        FuncEnd = FuncStart+MaxSegPerFunc-1;
    else
        FuncEnd = n_segs+1;
    end

    fprintf('\n\n\n(%d/%d) FUNCTION %d LOADED \n\n\n',k);
end

fprintf('Step Complete:CSV Functions loaded in %.2f s\n\n', toc);
%% Loop to create lines and apply currents
nTotal = numel(xs);         % Total number of start points
S      = unique(segid(:))'; % List of seg_ids
nSeg   = numel(S);          % Total number of differnet segment ids
LinesMade = 0;
fprintf('Creating %d geometry lines across %d segments...\n', nTotal, nSeg);

tic
for j = 1:nSeg
    s   = S(j);                    % segment id
    idx = find(segid==s);          % rows for this segment

    for k = 1:numel(idx)
        i   = idx(k);

        % create geometry line
        gtag = sprintf('segid_%d', s);
        geom.create(gtag, 'LineSegment');
        geom.feature(gtag).set('specify1','coord');
        geom.feature(gtag).set('coord1',{num2str(xs(i),'%.16g'), ...
                                         num2str(ys(i),'%.16g'), ...
                                         num2str(zs(i),'%.16g')});
        geom.feature(gtag).set('specify2','coord');
        geom.feature(gtag).set('coord2',{num2str(xe(i),'%.16g'), ...
                                         num2str(ye(i),'%.16g'), ...
                                         num2str(ze(i),'%.16g')});

        % turn on selection for this feature
        SegLabel = sprintf('Segid_%d',s);
        geom.feature(gtag).label(SegLabel);
        geom.feature(gtag).set('selresult', true);
        % build feature so selection exists
        geom.run(gtag);

        seltag = sprintf('sel_%d',j);
        selname = sprintf('segid_%d',s);

        geom.create(seltag, 'ExplicitSelection');
        geom.feature(seltag).selection('selection').init(1);
        geom.feature(seltag).label(selname);
        geom.feature(seltag).selection('selection').set(gtag, 1);
        geom.run(seltag);

        curtag = sprintf('lsc_%d',j);
        autoSelTag = sprintf('%s_%s', geomTag, seltag);
        autoSelTagline = sprintf('%s_%s_edg', geomTag, gtag); 
        
        funclab = sprintf('segi%d',s);

        model.component(compTag).physics('ec').create(curtag, 'LineCurrentSource', 1);
        model.component(compTag).physics('ec').feature(curtag).set('Qjl', sprintf('%s(t)', funclab));
        model.component(compTag).physics('ec').feature(curtag).selection.named(autoSelTag);
        model.component(compTag).physics('ec').feature(curtag).label(sprintf('Cur_Segid_%d',s));
    end
    
    % Update process text
    LinesMade = LinesMade + numel(idx);
    fprintf('(%d/%d) seg %d loaded \n', ...
            j, nSeg, s);
end

geom.run('fin');
fprintf('Step Complete:All segments loaded in %.2f s\n', toc);
%% Save Model
disp('Saving updated model...');
mphsave(model, ModelName);
disp('Step Complete:Currents assigned per seg_id and model saved');

%% Functions
function P_rot = rotateXYZ(P, rotX_deg, rotY_deg, rotZ_deg, origin, order)
% Rotate Nx3 points about X/Y/Z around "origin"
% order: "ZYX" (default) means apply X then Y then Z in a stable convention
    if nargin < 5 || isempty(origin), origin = [0 0 0]; end
    if nargin < 6 || isempty(order),  order  = "ZYX";   end

    ax = deg2rad(rotX_deg);
    ay = deg2rad(rotY_deg);
    az = deg2rad(rotZ_deg);

    Rx = [ 1     0        0;
           0  cos(ax) -sin(ax);
           0  sin(ax)  cos(ax) ];

    Ry = [ cos(ay) 0 sin(ay);
              0    1    0;
          -sin(ay) 0 cos(ay) ];

    Rz = [ cos(az) -sin(az) 0;
           sin(az)  cos(az) 0;
              0        0    1 ];

    R = eye(3);
    for ch = char(order)
        switch ch
            case 'X', R = Rx * R;
            case 'Y', R = Ry * R;
            case 'Z', R = Rz * R;
            otherwise, error('Unknown rotation axis in order: %s', ch);
        end
    end

    P0 = P - origin;
    P_rot = (R * P0.').';
    P_rot = P_rot + origin;
end
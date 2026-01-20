
% File: Model_Code_SingleNeuron.m
% Authors: Rana Jaylani, Sian McConchie

% Description:
%   This script uses MATLAB LiveLink to create a time-dependent COMSOL 
%   Multiphysics model of a single neuron

% Requirements:
%   - MATLAB with COMSOL LiveLink
%   - COMSOL Multiphysics 6.3

% Usage:
%   1. Update the file name and save path as suitable
%   2. Open the script in MATLAB LiveLink
%   3. Click Run

%%
clear, clc
%% Starting COMSOL LiveLink
import com.comsol.model.*
import com.comsol.model.util.*

% Create new model
model = ModelUtil.create('Model');

% Name and path to save model
ModelName = 'Model_SingleNeuron';
path = 'D:\2026_CapstonePaper\2_MATLABCODE';    % path to folder to save model
fullpath = [path '\' ModelName '.mph'];         % path to folder with model name     

%% Create a component
compTag = 'comp1';
geomTag = 'geom1';
model.modelNode.create(compTag);
model.modelNode(compTag).label(ModelName);

% Study setup
start_time = 0;
step_time = 3.3e-5;     % ≈ 3.33e-5 s
stop_time = 0.01;       % Length of study in secs -> 10ms increasing this length may increase runtime and affect time-limited neuron currents.
timelist_str = sprintf('range(%g,%g,%g)', start_time, step_time, stop_time);
eventtolerance = 1e-5;  % This number changes accuracy:Smaller tolerance = higher accuracy & slower to solve

% Parameters
Electrode_r = 5e-6;     % Radius of the microwire
model.param.set('Electrode_r', Electrode_r);
Electrode_l = 100e-6;   % Height of active site
model.param.set('Electrode_l', Electrode_l);
model.param.set('elec_sa', '2*pi*Electrode_r*Electrode_l + pi*Electrode_r^2');
model.param.set('elec_sa_cm', 'elec_sa*10^4');
Tissue_r = 7e-4;        % Tissue Radius

% Electrode Position within the tissue
Elec_xpos = 50e-6;
Elec_ypos = -50e-6;
Elec_zpos = -75e-6;

%% Selections: 3 -domain, 2 -boundaries, 1 -edges, 0 -points
% ElectrodeDomain
% ElectrodeBoundries
% ElectrodeBoundries_NoTop
% TissueDomain
% TissueBoundries
% AllDomains
%% Geometry
model.geom.create(geomTag, 3);          % Creates a 3D geometry
model.geom(geomTag).lengthUnit('m');    % Units in this case m: meters
model.view('view1').set('transparency', true);

geom = model.component(compTag).geom(geomTag);

% Creates Electrode in the model using parameters defined previously
model.geom(geomTag).create('Electrode', 'Cylinder');
model.geom(geomTag).feature('Electrode').set('r', 'Electrode_r');
model.geom(geomTag).feature('Electrode').set('h', 'Electrode_l');
model.geom(geomTag).feature('Electrode').set('pos', [Elec_xpos Elec_ypos Elec_zpos]);
model.geom(geomTag).run('Electrode');
geom.run;

% Electrode Selection (Domain)
ElectrodeDomain = createDomainSelection(model, compTag, geomTag, 'Electrode', 'ElectrodeDomain');
% Electrode Selection (Boundries)
model.selection.create('ElectrodeBoundries', 'Explicit');
model.selection('ElectrodeBoundries').label('ElectrodeBoundries');
model.selection('ElectrodeBoundries').geom(geomTag, 2);     
model.selection('ElectrodeBoundries').set([1 2 3 4 5 6]);   % Selects all boudries on the electrode
% Electrode Selection for probe (Boundries - Top) 
model.selection.create('ElectrodeBoundries_NoTop', 'Explicit');
model.selection('ElectrodeBoundries_NoTop').label('ElectrodeBoundries_NoTop');
model.selection('ElectrodeBoundries_NoTop').geom(geomTag, 2);     
model.selection('ElectrodeBoundries_NoTop').set([1 2 3 5 6]);% Selects all boundaries on the electrode except the one that is not in contact with the tissue

% Add tissue
model.geom(geomTag).create('Tissue', 'Sphere');
model.geom(geomTag).feature('Tissue').set('r', Tissue_r);
model.geom(geomTag).run('Tissue');
geom.run;

% Tissue Selection (Domain)
TissueDomain = createDomainSelection(model, compTag, geomTag, 'Tissue', 'TissueDomain');
% Tissue Boundry
boundaryIDs = mphgetadj(model, geomTag, 'boundary', 'domain', TissueDomain);
allBnds = 1:length(boundaryIDs);
ElectrodeBoundries = model.selection('ElectrodeBoundries').entities;
% Boundaries not in the selection
TissueBoundries = setdiff(allBnds, ElectrodeBoundries);
% Tissue Selection - For grounding
model.selection.create('TissueBoundries', 'Explicit');
model.selection('TissueBoundries').label('TissueBoundries');
model.selection('TissueBoundries').geom(geomTag, 2);     
model.selection('TissueBoundries').set(TissueBoundries);

% Union
model.geom(geomTag).create('uni1', 'Union');
model.geom(geomTag).feature('uni1').selection('input').set({'Electrode' 'Tissue'});
model.geom(geomTag).feature('uni1').set('keep', true);
model.geom(geomTag).run('uni1');
geom.run;

model.geom('geom1').runPre('fin');
geom.run;
disp('Step Complete:Geometry Created')
%% Define properties to all Materials in the model
% Brain grey Matter
model.material.create('mat1', 'Common');
model.material('mat1').propertyGroup.create('MultipoleDebye', 'MultipoleDebye', 'Multipole Debye');
model.material('mat1').label('Brain Grey Matter');
model.material('mat1').propertyGroup('def').set('electricconductivity', {'0.627[S/m]' '0' '0' '0' '0.627[S/m]' '0' '0' '0' '0.627[S/m]'});
model.material('mat1').propertyGroup('def').set('relpermittivity', {'57.301' '0' '0' '0' '57.301' '0' '0' '0' '57.301'});
model.material('mat1').propertyGroup('MultipoleDebye').set('Tref', '37[degC]');
model.material('mat1').propertyGroup('MultipoleDebye').set('Gvm', {'34.036' '8.406' '7.798'});
model.material('mat1').propertyGroup('MultipoleDebye').set('tauvm', {'7.216e-12[s]' '2.267e-11[s]' '1.999e-10[s]'});

% Platinium
model.material.create('mat2', 'Common');
model.material('mat2').propertyGroup.create('Enu', 'Enu', 'Young''s modulus and Poisson''s ratio');
model.material('mat2').label('Pt - Platinum');
model.material('mat2').set('family', 'custom');
model.material('mat2').set('customspecular', [0.7843137254901961 1 1]);
model.material('mat2').set('diffuse', 'custom');
model.material('mat2').set('customdiffuse', [0.7843137254901961 0.7843137254901961 0.7843137254901961]);
model.material('mat2').set('ambient', 'custom');
model.material('mat2').set('customambient', [0.7843137254901961 0.7843137254901961 0.7843137254901961]);
model.material('mat2').set('noise', true);
model.material('mat2').set('fresnel', 0.9);
model.material('mat2').set('roughness', 0.1);
model.material('mat2').set('diffusewrap', 0);
model.material('mat2').set('reflectance', 0);
model.material('mat2').propertyGroup('def').set('electricconductivity', {'8.9e6[S/m]' '0' '0' '0' '8.9e6[S/m]' '0' '0' '0' '8.9e6[S/m]'});
model.material('mat2').propertyGroup('def').set('thermalexpansioncoefficient', {'8.80e-6[1/K]' '0' '0' '0' '8.80e-6[1/K]' '0' '0' '0' '8.80e-6[1/K]'});
model.material('mat2').propertyGroup('def').set('heatcapacity', '133[J/(kg*K)]');
model.material('mat2').propertyGroup('def').set('density', '21450[kg/m^3]');
model.material('mat2').propertyGroup('def').set('thermalconductivity', {'71.6[W/(m*K)]' '0' '0' '0' '71.6[W/(m*K)]' '0' '0' '0' '71.6[W/(m*K)]'});
model.material('mat2').propertyGroup('Enu').set('E', '168e9[Pa]');
model.material('mat2').propertyGroup('Enu').set('nu', '0.38');
model.material('mat2').propertyGroup('def').set('relpermittivity', {'1'});

% Assinging Materials
model.material('mat1').selection.all;
model.material('mat2').selection.named('ElectrodeDomain');

disp('Step Complete:Materials Assigned')
%% Electrode Boundry Probe
model.component(compTag).probe.create('bnd1', 'Boundary');
model.component(compTag).probe('bnd1').set('expr', 'V');
model.component(compTag).probe('bnd1').set('intsurface', true);
model.component(compTag).probe('bnd1').selection.named('ElectrodeBoundries_NoTop');
disp('Step Complete:Boundry Probe Assigned')
%% Extract all Domains
domains = geom.getNDomains;   % Total number of domains
domainIDs = 1:domains;

%% Physics selection
model.selection.create('Doms', 'Explicit');     % Create named selection
model.selection('Doms').label('All Domains');
model.selection('Doms').geom(geomTag, 3);      % 3 = domain dimension
model.selection('Doms').set(domainIDs);

% Check the selection
doms = mphgetselection(model, 'Doms');
disp('Physics domain IDs:');
disp(doms);

%% Electric Currents (EC)
model.physics.create('ec', 'ConductiveMedia', geomTag);
model.physics('ec').selection.named('Doms'); 

% Ground Entire Tissue Boundry
model.physics('ec').create('gnd1', 'Ground', 2);
model.physics('ec').feature('gnd1').selection.named('TissueBoundries');

disp('Step Complete:Electric Currents Physics added')
%% Electric Circuit (CIR)
model.physics.create('cir', 'Circuit', geomTag); 

% Voltage Source aka Electrode Probe
model.physics('cir').create('V1', 'VoltageSource', -1);
model.physics('cir').feature('V1').setIndex('Connections', 0, 1, 0);
model.physics('cir').feature('V1').set('value', 'bnd1');

% Surface Conditions
model.physics('cir').create('C1', 'Capacitor', -1);
model.physics('cir').feature('C1').setIndex('Connections', 1, 0, 0);
model.physics('cir').feature('C1').setIndex('Connections', 2, 1, 0);
model.physics('cir').feature('C1').set('C', '20*elec_sa_cm [uF]');
model.physics('cir').create('R1', 'Resistor', -1);
model.physics('cir').feature('R1').setIndex('Connections', 1, 0, 0);
model.physics('cir').feature('R1').setIndex('Connections', 2, 1, 0);
model.physics('cir').feature('R1').set('R', ['((1.33*10^4)/elec_sa_cm) [' 'ohm' ']']);
% Front End amplifier
model.physics('cir').create('R2', 'Resistor', -1);
model.physics('cir').feature('R2').setIndex('Connections', 2, 0, 0);
model.physics('cir').feature('R2').setIndex('Connections', 3, 1, 0);
model.physics('cir').feature('R2').set('R', ['450 [M' 'ohm' ']']);
model.physics('cir').create('C2', 'Capacitor', -1);
model.physics('cir').feature('C2').setIndex('Connections', 2, 0, 0);
model.physics('cir').feature('C2').setIndex('Connections', 3, 1, 0);
model.physics('cir').feature('C2').set('C', '15 [pF]');
% Volt meter
model.physics('cir').create('vm1', 'VoltMeter', -1);
model.physics('cir').feature('vm1').setIndex('Connections', 2, 0, 0);
model.physics('cir').feature('vm1').setIndex('Connections', 3, 1, 0);
model.physics('cir').create('gnd2', 'GroundNode', -1);
model.physics('cir').feature('gnd2').setIndex('Connections', 3, 0, 0);

disp('Step Complete:Electric Circuits Physics added')
%% Mesh
model.mesh.create('mesh1', geomTag);
model.mesh('mesh1').contribute('geom/detail', true);
model.mesh('mesh1').create('ftet1', 'FreeTet');

% hauto	refers to mesh size
model.mesh('mesh1').feature('size').set('hauto', 7);
model.mesh('mesh1').run;

disp('Step Complete:Mesh Added')
%% Study
model.study.create('std1');
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').setSolveFor('/physics/ec', true);
model.study('std1').feature('time').setSolveFor('/physics/cir', true);
model.study('std1').feature('time').set('tlist', timelist_str);

model.study('std1').createAutoSequences('all');

model.sol('sol1').feature('t1').set('tstepsbdf', 'manual');
model.sol('sol1').feature('t1').set('timestepbdf', num2str(step_time));
model.sol('sol1').feature('t1').set('eventtol', eventtolerance);
model.sol('sol1').feature('t1').set('stabcntrl', true);

model.study('std1').createAutoSequences('all');

% Set to automatic newton
model.sol('sol1').feature('t1').feature('se1').feature('ss2').set('subdtech', 'auto');

disp('Step Complete:Study Settings Added')

%% Save the model
disp('Saving updated model...');
mphsave(fullpath);
clear model

disp('Step Complete: Model Created and saved under:')
disp(fullpath)

%% Import Single Neuron
% Configuration
morph_csvPath = 'segment_morphologies.csv'; % -> CHANGE MORPHOLOGY FILE PATH
cur_csvPath = 'segment_currents_paramA_singlespike_10ms.csv'; % -> CHANGE CURRENT FILE PATH
cellName = 'pyr_paramA'; % Update with cell type, parameter type and recording no. 
scaleFactor = 1e-6;      % Coordinates are in µm

%% Load Model
model = mphopen(ModelName);
% Ensure component and geometry exist
try model.component(compTag); catch, model.component.create(compTag,true); end
try model.component(compTag).geom(geomTag); catch, model.component(compTag).geom.create(geomTag,3); end
geom = model.component(compTag).geom(geomTag);

%% Read Morphology CSV
% create offsets
x_os = 0; % -> INSERT REQUIRED CELL POSITIONING
y_os = 0;
z_os = 0;

% get segment start points
morphT   = readtable(morph_csvPath);
xs   = (morphT.xs(:) + x_os) * scaleFactor;
ys   = (morphT.ys(:) + y_os) * scaleFactor;
zs   = (morphT.zs(:) + z_os) * scaleFactor;

% get segment end points
xe   = (morphT.xe(:) + x_os) * scaleFactor;
ye   = (morphT.ye(:) + y_os ) * scaleFactor;
ze   = (morphT.ze(:) + z_os) * scaleFactor;

% get segment names and ids
segid = morphT.seg_id(:); 
%% Define interpolation function: I_Seg(t, seg)
tic
disp('Loading interpolation functions...\n');
functag = sprintf('%s_func',cellName);
model.func.create(cellName, 'Interpolation');
model.func(cellName).label(cellName);
model.func(cellName).set('source', 'file');
model.func(cellName).set('filename', cur_csvPath);
n_segs = height(morphT);
fcount = 1;
indexcount = 0;
model.func(cellName).setIndex('argunit', 'ms', 0);
for n = 2:n_segs + 1
col = sprintf('col%d',n);
functag = sprintf('%s_seg%d_I',cellName,fcount);
model.func(cellName).setEntry('columnType', col, 'value');
model.func(cellName).setEntry('funcnames', col, functag);
model.func(cellName).setIndex('fununit', 'pA/um', indexcount);
fcount = fcount + 1;
indexcount = indexcount + 1;
fprintf('(%d/%d) Function loaded for seg %d \n',indexcount, n_segs, indexcount);
end
model.func(cellName).set('interp', 'linear');
model.func(cellName).set('extrap', 'const');  % Hold outside range (or 'zero')

fprintf('CSV Functions loaded in %.2f s\n\n', toc);

%% Create Groups
linegroupname = cellName;
currentgroupname = cellName;

geom.nodeGroup.create(linegroupname);
geom.nodeGroup(linegroupname).placeAfter('uni1');
geom.nodeGroup(linegroupname).label(linegroupname)

model.nodeGroup.create(currentgroupname, 'Physics', 'ec');
model.nodeGroup(currentgroupname).placeAfter('gnd1');
model.nodeGroup(currentgroupname).label(cellName)

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
    %section_name = seg_label{j};  % section label from morphology

    for k = 1:numel(idx)
        i   = idx(k);

        % create geometry line
        gtag = sprintf('%s_ln_seg%d',cellName, s);
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
        geom.feature(gtag).label(gtag);
        geom.feature(gtag).set('selresult', true);
        % build feature so selection exists
        geom.run(gtag);

        seltag = sprintf('%s_sel_seg%d',cellName, s);
        selname = sprintf('%s_seg%d',cellName, s);

        geom.create(seltag, 'ExplicitSelection');
        geom.feature(seltag).selection('selection').init(1);
        geom.feature(seltag).label(selname);
        geom.feature(seltag).selection('selection').set(gtag, 1);
        geom.run(seltag);

        curtag = sprintf('%s_cur_seg%d',cellName,j);
        autoSelTag = sprintf('%s_%s', geomTag, seltag);
        autoSelTagline = sprintf('%s_%s_edg', geomTag, gtag);
        
        funclab = sprintf('%s_seg%d_I',cellName, s);

        model.component(compTag).physics('ec').create(curtag, 'LineCurrentSource', 1);
        model.component(compTag).physics('ec').feature(curtag).set('Qjl', sprintf('%s(t)', funclab));
        model.component(compTag).physics('ec').feature(curtag).selection.named(autoSelTag);
        model.component(compTag).physics('ec').feature(curtag).label(sprintf('%s_seg%d',cellName, s));
        
        geom.nodeGroup(linegroupname).add(seltag);
        geom.nodeGroup(linegroupname).add(gtag);
        model.nodeGroup(currentgroupname).add(curtag);
    end
    
    % Update process text
    LinesMade = LinesMade + numel(idx);
    fprintf('(%d/%d) seg %d loaded \n', ...
            j, nSeg, s);
end

geom.run('fin');
fprintf('All segments loaded in %.2f s\n', toc);

%% Save Model
disp('Saving updated model...');
mphsave(model, ModelName);
disp('Done. Currents assigned per seg_id and model saved.');

%% Functions
% Creates a new domain selection with all nonselected domains
function newDomain = createDomainSelection(model, compTag, geomTag, featureTag, domainSelName)
    % Run geometry feature to ensure it's added
    model.geom(geomTag).run(featureTag);

    % Get total number of domains
    geom = model.component(compTag).geom(geomTag);
    totalDomains = geom.getNDomains;
    allDomains = 1:totalDomains;

    % Collect already used domains (from existing selections)
    usedDomains = [];
    selTags = model.selection.tags;
    for i = 1:length(selTags)
        sel = model.selection(selTags(i));

        % Check if selection lives in 3D (domain)
        if sel.dim == 3
            entities = sel.entities;
            usedDomains = union(usedDomains, entities);
        end
    end

    % Find new domain(s)
    newDomain = setdiff(allDomains, usedDomains);

    % Create domain selection
    model.selection.create(domainSelName, 'Explicit');
    model.selection(domainSelName).label(domainSelName);
    model.selection(domainSelName).geom(geomTag, 3);
    model.selection(domainSelName).set(newDomain);
end
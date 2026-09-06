
% File: Model_HighdensityShank.m
% Authors: Rana Jaylani, Sian McConchie

% Description:
% This script uses MATLAB LiveLink to create a time-dependent COMSOL model
% Creates a COMSOL model with a high density shank of 188 electrodes

% Requirements:
%   - MATLAB with COMSOL LiveLink
%   - COMSOL Multiphysics (Tested on 6.2, 6.3 and 6.4)

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

% Name and path to save model - CHANGE AS REQUIRED
ModelName = 'Model_HighDensityShank1'; % -> CHANGE NAME TO MODEL REQUIRED
path = '';    % path to folder to save model
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
Electrode_w = 7e-6;     % Width pixel
model.param.set('Electrode_w', Electrode_w);
Electrode_l = 7.5e-6;   % Height of active site
model.param.set('Electrode_l', Electrode_l);
model.param.set('elec_sa', 'Electrode_w^2');
model.param.set('elec_sa_cm', 'elec_sa*10^4');
Tissue_r = 15e-4;        % Tissue Radius

% Electrodes Position within the tissue
Elec_xpos = 0;
model.param.set('Elec_xpos', Elec_xpos);
Elec_ypos = 0;
model.param.set('Elec_ypos', Elec_ypos);
Elec_zpos = 0;
model.param.set('Elec_zpos', Elec_zpos);

%% Geometry
model.geom.create(geomTag, 3);          % Creates a 3D geometry
model.geom(geomTag).lengthUnit('m');    % Units in this case m: meters
model.view('view1').set('transparency', true);

geom = model.component(compTag).geom(geomTag);

% Cumulative Selections Set up 
model.geom(geomTag).selection.create('cselElec','CumulativeSelection');
model.geom(geomTag).selection('cselElec').label('Electrode');
model.geom(geomTag).selection.create('cselIns','CumulativeSelection');
model.geom(geomTag).selection('cselIns').label('Insulated');
model.geom(geomTag).selection.create('cselTiss','CumulativeSelection');
model.geom(geomTag).selection('cselTiss').label('Tissue');

% Polymide Block of Neuropixel
model.geom(geomTag).create('InsulatedBlock', 'Block');
model.geom(geomTag).feature('InsulatedBlock').set('size', {'75 [um]' '15 [um]' '1600 [um]'});
model.geom(geomTag).feature('InsulatedBlock').set('pos', {'Elec_xpos' 'Elec_ypos' 'Elec_zpos + 100e-6'});
model.geom(geomTag).feature('InsulatedBlock').set('base', 'center');
model.geom(geomTag).feature('InsulatedBlock').set('contributeto','cselIns');
model.geom(geomTag).run('InsulatedBlock');
geom.run;

% Electrode Channel
model.geom(geomTag).create('Elec1', 'Block');
model.geom(geomTag).feature('Elec1').set('size', {'Electrode_w' 'Electrode_l' 'Electrode_w'});
model.geom(geomTag).feature('Elec1').set('base', 'center');
model.geom(geomTag).feature('Elec1').set('pos', {'Elec_xpos + 14 [um]' 'Elec_ypos + Electrode_l/2' 'Elec_zpos'});
model.geom(geomTag).feature('Elec1').set('contributeto','cselElec');
model.geom(geomTag).run('Elec1');
model.geom(geomTag).create('Elec2', 'Block');
model.geom(geomTag).feature('Elec2').set('size', {'Electrode_w' 'Electrode_l' 'Electrode_w'});
model.geom(geomTag).feature('Elec2').set('base', 'center');
model.geom(geomTag).feature('Elec2').set('pos', {'Elec_xpos - 14 [um]' 'Elec_ypos + Electrode_l/2' 'Elec_zpos'});
model.geom(geomTag).feature('Elec2').set('contributeto','cselElec');
model.geom(geomTag).run('Elec2');
geom.run;

% Duplicate electrode channesl into array
model.geom(geomTag).create('copy1', 'Copy');
model.geom(geomTag).feature('copy1').selection('input').set({'Elec1' 'Elec2'});
model.geom(geomTag).feature('copy1').set('displz', 'range(-600 [um],15 [um],800 [um])');
model.geom(geomTag).feature('copy1').set('contributeto','cselElec');
model.geom(geomTag).run('copy1');
geom.run;

model.geom(geomTag).create('dif1', 'Difference');
model.geom(geomTag).feature('dif1').selection('input').set({'InsulatedBlock'});
model.geom(geomTag).feature('dif1').selection('input2').set({'copy1'});
model.geom(geomTag).feature('dif1').set('keepsubtract', true);
model.geom(geomTag).run('dif1');
geom.run;

% Add tissue
model.geom(geomTag).create('Tissue', 'Sphere');
model.geom(geomTag).feature('Tissue').set('r', Tissue_r);
model.geom(geomTag).feature('Tissue').set('contributeto','cselTiss');
model.geom(geomTag).run('Tissue');
geom.run;

% Union
model.geom(geomTag).create('uni1', 'Union');
model.geom(geomTag).feature('uni1').selection('input').set({'copy1' 'Tissue' 'InsulatedBlock'});
model.geom(geomTag).feature('uni1').set('keep', true);
model.geom(geomTag).run('uni1');
geom.run;

model.geom('geom1').runPre('fin');
geom.run;
disp('Step Complete:Geometry Created')
%% Define properties to all Materials in the model
% Brain grey Matter
model.material.create('mat1', 'Common');
% Apply COMSOL generated material properties
model.material('mat1').propertyGroup.create('MultipoleDebye', 'MultipoleDebye', 'Multipole Debye');
model.material('mat1').label('Brain Grey Matter');
model.material('mat1').propertyGroup('def').set('electricconductivity', {'0.627[S/m]' '0' '0' '0' '0.627[S/m]' '0' '0' '0' '0.627[S/m]'});
model.material('mat1').propertyGroup('def').set('relpermittivity', {'57.301' '0' '0' '0' '57.301' '0' '0' '0' '57.301'});
model.material('mat1').propertyGroup('MultipoleDebye').set('Tref', '37[degC]');
model.material('mat1').propertyGroup('MultipoleDebye').set('Gvm', {'34.036' '8.406' '7.798'});
model.material('mat1').propertyGroup('MultipoleDebye').set('tauvm', {'7.216e-12[s]' '2.267e-11[s]' '1.999e-10[s]'});

% Platinium
model.material.create('mat2', 'Common');
% Apply COMSOL generated material properties
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
% Apply custom specified permitivity
model.material('mat2').propertyGroup('def').set('relpermittivity', {'1'});

% Insulating Material - Polymide (Kapton)
model.material.create('mat3', 'Common');
% Apply COMSOL generated material properties
model.material('mat3').label('Kapton MT [solid,without adhesive]');
model.material('mat3').info.create('Composition');
model.material('mat3').info('Composition').body('polyimide');
model.material('mat3').info('Composition').title('Composition');
model.material('mat3').info.create('Note');
model.material('mat3').info('Note').body('surface resistivity > 10y Mohm, at 1 kHz: Dk = 4.2');
model.material('mat3').info('Note').title('Note');
model.material('mat3').propertyGroup('def').set('thermalconductivity', 'k_solid_without_adhesive_2(T)');
model.material('mat3').propertyGroup('def').setPropertyInfo('thermalconductivity', ['Reference: D.L. Rule, D.R. Smith, and L.L. Sparks, "Thermal Conductivity of a polymide Film Between 4.2 and 300 K, with and without Alumina Particles as Fillers", NIST Internal Report 3948 (1990) https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir3948.pdf' newline 'Note: without adhesive']);
model.material('mat3').propertyGroup('def').set('resistivity', ['1.0E16[' 'ohm' '*m]']);
model.material('mat3').propertyGroup('def').setPropertyInfo('resistivity', ['Reference: DuPont Engineering Polymers, Product brochure https://www.dupont.com/brands.html' newline 'Note: minimum value, measured in accordance with IEC 60093, room temperature value']);
model.material('mat3').propertyGroup('def').set('heatcapacity', 'C(T)');
model.material('mat3').propertyGroup('def').setPropertyInfo('heatcapacity', ['Reference: NIST Material Measurement Laboratory, "Properties of solid materials from cryogenic- to room-temperatures", Applied Chemicals and Materials Division https://trc.nist.gov/cryogenics/materials/materialproperties.htm' newline 'Note: material is only described as Kapton, 3% error']);
model.material('mat3').propertyGroup('def').set('electricconductivity', '1.0E-16[S/m]');
model.material('mat3').propertyGroup('def').setPropertyInfo('electricconductivity', ['Reference: DuPont Engineering Polymers, Product brochure https://www.dupont.com/brands.html' newline 'Note: maximum value, measured in accordance with IEC 60093, calculated as the reciprocal of the resistivity, room temperature value']);
model.material('mat3').propertyGroup('def').set('density', '1420[kg/m^3]');
model.material('mat3').propertyGroup('def').setPropertyInfo('density', 'Note: room temperature value');
model.material('mat3').propertyGroup('def').func.create('k_solid_without_adhesive_2', 'Piecewise');
model.material('mat3').propertyGroup('def').func('k_solid_without_adhesive_2').set('funcname', 'k_solid_without_adhesive_2');
model.material('mat3').propertyGroup('def').func('k_solid_without_adhesive_2').set('arg', 'T');
model.material('mat3').propertyGroup('def').func('k_solid_without_adhesive_2').set('extrap', 'constant');
model.material('mat3').propertyGroup('def').func('k_solid_without_adhesive_2').set('pieces', {'5.0' '20.0' '(-1.30918997E-04*log(T)^4 + 1.24927049E-03*log(T)^3 - 4.02241446E-03*log(T)^2 + 5.43345876E-03*log(T) - 2.53007600E-03)*1.000000000e+02'; '20.0' '300.0' '(-1.57395530E-05*log(T)^4 + 1.80652338E-04*log(T)^3 - 3.83608296E-04*log(T)^2 + 6.52163650E-05*log(T) + 3.48765353E-04)*1.000000000e+02'});
model.material('mat3').propertyGroup('def').func('k_solid_without_adhesive_2').label('Piecewise');
model.material('mat3').propertyGroup('def').func('k_solid_without_adhesive_2').set('fununit', 'W/(m*K)');
model.material('mat3').propertyGroup('def').func('k_solid_without_adhesive_2').set('argunit', 'K');
model.material('mat3').propertyGroup('def').func.create('C', 'Piecewise');
model.material('mat3').propertyGroup('def').func('C').set('funcname', 'C');
model.material('mat3').propertyGroup('def').func('C').set('arg', 'T');
model.material('mat3').propertyGroup('def').func('C').set('extrap', 'constant');
model.material('mat3').propertyGroup('def').func('C').set('pieces', {'4.0' '30.0' '2.809665621-1.394348543*T^1+0.2106638561*T^2+0.004752015656*T^3-3.279997847E-4*T^4+4.282248688E-6*T^5'; '30.0' '300.0' '-86.85311476+7.816916992*T^1-0.03664788194*T^2+9.912799744E-5*T^3-1.084409538E-7*T^4'});
model.material('mat3').propertyGroup('def').func('C').label('Piecewise 1');
model.material('mat3').propertyGroup('def').func('C').set('fununit', 'J/(kg*K)');
model.material('mat3').propertyGroup('def').func('C').set('argunit', 'K');
model.material('mat3').propertyGroup('def').addInput('temperature');
model.material('mat3').set('family', 'custom');
model.material('mat3').set('lighting', 'cooktorrance');
model.material('mat3').set('fresnel', 0.12);
model.material('mat3').set('roughness', 0.45);
model.material('mat3').set('metallic', 0);
model.material('mat3').set('pearl', 0.015);
model.material('mat3').set('diffusewrap', 0.4);
model.material('mat3').set('clearcoat', 0.15);
model.material('mat3').set('reflectance', 1);
model.material('mat3').set('ambient', 'custom');
model.material('mat3').set('customambient', [0.6901960784313725 0.6901960784313725 0.6901960784313725]);
model.material('mat3').set('diffuse', 'custom');
model.material('mat3').set('customdiffuse', [0.7529411764705882 0.7529411764705882 0.7529411764705882]);
model.material('mat3').set('specular', 'custom');
model.material('mat3').set('customspecular', [1 1 1]);
model.material('mat3').set('noisecolor', 'custom');
model.material('mat3').set('customnoisecolor', [0 0 0]);
model.material('mat3').set('transparency', 0);
model.material('mat3').set('specular', 'custom');
model.material('mat3').set('customspecular', [0.7843137254901961 0.7843137254901961 0.7843137254901961]);
model.material('mat3').set('diffuse', 'custom');
model.material('mat3').set('customdiffuse', [0.39215686274509803 0.39215686274509803 0.9803921568627451]);
model.material('mat3').set('ambient', 'custom');
model.material('mat3').set('customambient', [0.39215686274509803 0.39215686274509803 0.7843137254901961]);
model.material('mat3').set('shininess', 500);
% Apply custom specified permitivity
model.material('mat3').propertyGroup('def').set('relpermittivity', {'3.3'});

% Assinging Materials
model.material('mat1').selection.all;
model.material('mat2').selection.named([geomTag '_cselElec_dom']);
model.material('mat3').selection.named([geomTag '_cselIns_dom']);

disp('Step Complete:Materials Assigned')
%% Electrodes Boundry Probe
elecBnds = model.selection([geomTag '_cselElec_bnd']).entities;
elecBnds = double(elecBnds(:));

Recordingboundries = [(386:479)';(954:1047)'];
for i = 1:length(Recordingboundries)
    probeName = ['Electrode' num2str(i)];
    
    model.component(compTag).probe.create(probeName, 'Boundary');
    model.component(compTag).probe(probeName).set('expr', 'V');
    model.component(compTag).probe(probeName).set('intsurface', true);
    model.component(compTag).probe(probeName).selection.set(Recordingboundries(i));
end

disp('Step Complete:Boundry Probes Assigned');
%% Extract all Domains
domains = geom.getNDomains;   % Total number of domains
domainIDs = 1:domains;
% All domains in the finalised geometry
allDoms  = 1:double(model.geom(geomTag).getNDomains());
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
elecDoms = model.selection([geomTag '_cselElec_dom']).entities;
elecDoms = double(elecDoms(:))';
insDoms = model.selection([geomTag '_cselIns_dom']).entities;
insDoms = double(insDoms(:))';
tissDoms = setdiff(allDoms, union(elecDoms, insDoms)); % tissue domain

% outer tissue surface
outerTiss = setdiff( adjBnds(model,geomTag,tissDoms), ...
                     union(adjBnds(model,geomTag,elecDoms), ...
                           adjBnds(model,geomTag,insDoms)) );

model.selection.create('TissueOuter','Explicit');
model.selection('TissueOuter').label('TissueOuter');
model.selection('TissueOuter').geom(geomTag, 2);
model.selection('TissueOuter').set(outerTiss);

% Ground Entire Tissue Boundry
model.physics('ec').create('gnd1', 'Ground', 2);
model.physics('ec').feature('gnd1').selection.named('TissueOuter');

disp('Step Complete:Electric Currents Physics added')
%% Electric Circuit (CIR)
model.physics.create('cir', 'Circuit', geomTag); 
% Create circuit with surface conditions and front end amplifier for each ch
for i = 1:length(Recordingboundries) % Number of Electrodes
    % Unique names for each circuit
    Vname  = ['V' num2str(i)];
    C1name = ['C1_' num2str(i)];
    R1name = ['R1_' num2str(i)];
    R2name = ['R2_' num2str(i)];
    C2name = ['C2_' num2str(i)];
    VMname = ['vm' num2str(i)];
    GNDname= ['grnd' num2str(i)];

    % Node numbering (3 nodes per electrode)
    n1  = 3*(i-1) + 1;
    n2 = 3*(i-1) + 2;
    n3 = 3*(i-1) + 3;

    bndVar = ['Electrode' num2str(i)];

    % Voltage Source
    model.physics('cir').create(Vname, 'VoltageSource', -1);
    model.physics('cir').feature(Vname).setIndex('Connections', 0, 1, 0);
    model.physics('cir').feature(Vname).setIndex('Connections', n1, 0, 0);
    model.physics('cir').feature(Vname).set('value', bndVar);

    % Surface conditions
    model.physics('cir').create(C1name, 'Capacitor', -1);
    model.physics('cir').feature(C1name).setIndex('Connections', n1, 0, 0);
    model.physics('cir').feature(C1name).setIndex('Connections', n2, 1, 0);
    model.physics('cir').feature(C1name).set('C', '20*elec_sa_cm [uF]');

    model.physics('cir').create(R1name, 'Resistor', -1);
    model.physics('cir').feature(R1name).setIndex('Connections', n1, 0, 0);
    model.physics('cir').feature(R1name).setIndex('Connections', n2, 1, 0);
    model.physics('cir').feature(R1name).set('R', '((1.33e4)/elec_sa_cm) [ohm]');

    % Front-end amplifier
    model.physics('cir').create(R2name, 'Resistor', -1);
    model.physics('cir').feature(R2name).setIndex('Connections', n2, 0, 0);
    model.physics('cir').feature(R2name).setIndex('Connections', n3, 1, 0);
    model.physics('cir').feature(R2name).set('R', '450[Mohm]');

    model.physics('cir').create(C2name, 'Capacitor', -1);
    model.physics('cir').feature(C2name).setIndex('Connections', n2, 0, 0);
    model.physics('cir').feature(C2name).setIndex('Connections', n3, 1, 0);
    model.physics('cir').feature(C2name).set('C', '15[pF]');

    % Voltmeter
    model.physics('cir').create(VMname, 'VoltMeter', -1);
    model.physics('cir').feature(VMname).setIndex('Connections', n2, 0, 0);
    model.physics('cir').feature(VMname).setIndex('Connections', n3, 1, 0);

    % Ground
    model.physics('cir').create(GNDname, 'GroundNode', -1);
    model.physics('cir').feature(GNDname).setIndex('Connections', n3, 0, 0);
end

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

disp('Step Complete:Model Created and saved under:')
disp(fullpath)

%% Fucntions 
function b = adjBnds(model, geomTag, doms)
    tmp = mphgetadj(model, geomTag, 'boundary', 'domain', doms);
    if iscell(tmp)
        tmp = cellfun(@(c) double(c(:))', tmp, 'UniformOutput', false);
        tmp = [tmp{:}];
    end
    b = unique(double(tmp(:))');
end
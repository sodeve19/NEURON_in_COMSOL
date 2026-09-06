
% File: Model_Mircowire_Array.m
% Authors: Rana Jaylani, Sian McConchie

% Description:
% This script uses MATLAB LiveLink to create a time-dependent COMSOL model
% Creates a 9x9 array of mircowire electrodes

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
ModelName = 'Model_Mircowire_Array'; % -> CHANGE NAME TO MODEL REQUIRED
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
stop_time = 0.01;       % Length of study in secs -> 10ms increasing this length may increase runtime and affect time-limited neuron currents
timelist_str = sprintf('range(%g,%g,%g)', start_time, step_time, stop_time);
eventtolerance = 1e-5;  % This number changes accuracy:Smaller tolerance = higher accuracy & slower to solve

% Parameters
Electrode_r = 5e-6;     % Radius of the microwire
model.param.set('Electrode_r', Electrode_r);
InsElectrode_r = Electrode_r + 2e-6; % Radius of insulated part of microwire
model.param.set('InsElectrode_r', InsElectrode_r);
Electrode_l = 100e-6;   % Height of active site
model.param.set('Electrode_l', Electrode_l);
model.param.set('elec_sa', '2*pi*Electrode_r*Electrode_l + pi*Electrode_r^2');
model.param.set('elec_sa_cm', 'elec_sa*10^4');
Tissue_r = 15e-4;        % Tissue Radius

% Array Settings and Position
pitch = 60e-6; % Pitch between microwires 
model.param.set('pitch', pitch);
ArrayNum = 2; %3 by 3 % Array size

Elec_xpos = 0;
Elec_xpos = Elec_xpos - (ArrayNum -1)*pitch/2;
model.param.set('Elec_xpos', Elec_xpos);
Elec_ypos = 0;
Elec_ypos = Elec_ypos - (ArrayNum -1)*pitch/2;
model.param.set('Elec_ypos', Elec_ypos);
Elec_zpos = 0;
model.param.set('Elec_zpos', '0 - Electrode_l/2');

%% Geometry
model.geom.create(geomTag, 3);          % Creates a 3D geometry
model.geom(geomTag).lengthUnit('m');    % Units in this case m: meters
model.view('view1').set('transparency', true);

geom = model.component(compTag).geom(geomTag);

% Cumulative Selections Set up
model.geom(geomTag).selection.create('cselElec','CumulativeSelection');
model.geom(geomTag).selection('cselElec').label('AllElectrodes');
model.geom(geomTag).selection.create('cselIns','CumulativeSelection');
model.geom(geomTag).selection('cselIns').label('AllInsulated');
model.geom(geomTag).selection.create('cselTiss','CumulativeSelection');
model.geom(geomTag).selection('cselTiss').label('Tissue');

% Creates Electrode in the model using parameters defined previously
model.geom(geomTag).create('Electrode', 'Cylinder');
model.geom(geomTag).feature('Electrode').set('r', 'Electrode_r');
model.geom(geomTag).feature('Electrode').set('h', 'Electrode_l');
model.geom(geomTag).feature('Electrode').set('pos', {'Elec_xpos' 'Elec_ypos' 'Elec_zpos'});
model.geom(geomTag).feature('Electrode').set('contributeto','cselElec');
model.geom(geomTag).run('Electrode');
geom.run;

% Insulated part of electrode
model.geom(geomTag).create('ElectrodeInsulated', 'Cylinder');
model.geom(geomTag).feature('ElectrodeInsulated').set('r', 'InsElectrode_r');
model.geom(geomTag).feature('ElectrodeInsulated').set('h', Tissue_r*2);
model.geom(geomTag).feature('ElectrodeInsulated').set('pos', {'Elec_xpos' 'Elec_ypos' ('Elec_zpos + Electrode_l')});
model.geom(geomTag).feature('ElectrodeInsulated').set('contributeto','cselIns'); 
model.geom(geomTag).run('ElectrodeInsulated');
geom.run;

% Make Array
model.geom('geom1').create('arr1', 'Array');
model.geom('geom1').feature('arr1').selection('input').set({'Electrode' 'ElectrodeInsulated'});
model.geom('geom1').feature('arr1').set('fullsize', [ArrayNum ArrayNum 1]);
model.geom('geom1').feature('arr1').set('displ', {'pitch' 'pitch' '0'});
model.geom('geom1').run('arr1');

% Add tissue
model.geom(geomTag).create('Tissue', 'Sphere');
model.geom(geomTag).feature('Tissue').set('r', Tissue_r);
model.geom(geomTag).feature('Tissue').set('contributeto','cselTiss');
model.geom(geomTag).run('Tissue');
geom.run;

% Merge Array
model.geom(geomTag).create('unicyl','Union');
model.geom(geomTag).feature('unicyl').selection('input').set({'arr1'});
model.geom(geomTag).feature('unicyl').set('intbnd', true);

% Copy the tissue so it isn't destroyed by the intersection
model.geom(geomTag).create('copyTissue','Copy');
model.geom(geomTag).feature('copyTissue').selection('input').set({'Tissue'});

% Remove excess electrodes
model.geom(geomTag).create('int1', 'Intersection');
model.geom(geomTag).feature('int1').selection('input').set({'unicyl','copyTissue'});

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

% Insulating Material - Parylene
model.material.create('mat3', 'Common');
% Apply COMSOL generated material properties
model.material('mat3').label('Parylene C [solid,50% RH]');
model.material('mat3').info.create('Composition');
model.material('mat3').info('Composition').body('poly(para-xylylene)');
model.material('mat3').info('Composition').title('Composition');
model.material('mat3').info.create('Note');
model.material('mat3').info('Note').body(['surface resistivity = 10' native2unicode(hex2dec({'00' 'b9'}), 'unicode') 't ohm, breakdown strength = 220.5 kV/mm (5600 V/mil)']);
model.material('mat3').info('Note').title('Note');
model.material('mat3').propertyGroup('def').set('thermalconductivity', '0.084[W/(m*K)]');
model.material('mat3').propertyGroup('def').setPropertyInfo('thermalconductivity', ['Reference: Specialty Coating Systems, SCS Parylene Properties brochure (2007) https://scscoatings.com/' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 290 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (563 K), measured in accordance with ASTM C177, room temperature value']);
model.material('mat3').propertyGroup('def').set('resistivity', ['8.8E14[' 'ohm' '*m]']);
model.material('mat3').propertyGroup('def').setPropertyInfo('resistivity', ['Reference: Specialty Coating Systems, SCS Parylene Properties brochure (2007) https://scscoatings.com/' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 290 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (563 K), 50% relative humidity, room temperature value']);
model.material('mat3').propertyGroup('def').set('thermalexpansioncoefficient', '(alpha(T)+(Tempref-293[K])*if(abs(T-Tempref)>1e-3,(alpha(T)-alpha(Tempref))/(T-Tempref),d(alpha(T),T)))/(1+alpha(Tempref)*(Tempref-293[K]))');
model.material('mat3').propertyGroup('def').setPropertyInfo('thermalexpansioncoefficient', ['Reference: Specialty Coating Systems, SCS Parylene Properties brochure (2007) https://scscoatings.com/' newline 'Note: the reference temperature is 20 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (293 K), T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 290 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (563 K), limited data' newline 'Reference temperature: 293.00[K]']);
model.material('mat3').propertyGroup('def').set('heatcapacity', '711.999999[J/(kg*K)]');
model.material('mat3').propertyGroup('def').setPropertyInfo('heatcapacity', ['Reference: Specialty Coating Systems, SCS Parylene Properties brochure (2007) https://scscoatings.com/' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 290 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (563 K), room temperature value']);
model.material('mat3').propertyGroup('def').set('electricconductivity', '1.13636364E-15[S/m]');
model.material('mat3').propertyGroup('def').setPropertyInfo('electricconductivity', ['Reference: Specialty Coating Systems, SCS Parylene Properties brochure (2007) https://scscoatings.com/' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 290 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (563 K), 50% relative humidity, calculated as the reciprocal of the resistivity, room temperature value']);
model.material('mat3').propertyGroup('def').set('density', 'rho(T)');
model.material('mat3').propertyGroup('def').setPropertyInfo('density', ['Reference: Specialty Coating Systems, SCS Parylene Properties brochure (2007) https://scscoatings.com/' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 290 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (563 K), calculated from the linear expansion and the room temperature density']);
model.material('mat3').propertyGroup('def').set('TD', '9.15264E-8');
model.material('mat3').propertyGroup('def').setPropertyInfo('TD', ['Reference: Specialty Coating Systems, SCS Parylene Properties brochure (2007) https://scscoatings.com/' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 290 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (563 K), calculated from the thermal conductivity, density, and specific heat, room temperature value']);
model.material('mat3').propertyGroup('def').func.create('alpha', 'Piecewise');
model.material('mat3').propertyGroup('def').func('alpha').set('funcname', 'alpha');
model.material('mat3').propertyGroup('def').func('alpha').set('arg', 'T');
model.material('mat3').propertyGroup('def').func('alpha').set('extrap', 'constant');
model.material('mat3').propertyGroup('def').func('alpha').set('pieces', {'293.0' '373.16' '3.5E-5'});
model.material('mat3').propertyGroup('def').func('alpha').label('Piecewise');
model.material('mat3').propertyGroup('def').func('alpha').set('fununit', '1/K');
model.material('mat3').propertyGroup('def').func('alpha').set('argunit', 'K');
model.material('mat3').propertyGroup('def').func.create('rho', 'Piecewise');
model.material('mat3').propertyGroup('def').func('rho').set('funcname', 'rho');
model.material('mat3').propertyGroup('def').func('rho').set('arg', 'T');
model.material('mat3').propertyGroup('def').func('rho').set('extrap', 'constant');
model.material('mat3').propertyGroup('def').func('rho').set('pieces', {'293.0' '373.16' '1329.46325-0.140856421*T^1+9.40815456E-6*T^2'});
model.material('mat3').propertyGroup('def').func('rho').label('Piecewise 1');
model.material('mat3').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
model.material('mat3').propertyGroup('def').func('rho').set('argunit', 'K');
model.material('mat3').propertyGroup('def').addInput('temperature');
model.material('mat3').propertyGroup('def').addInput('strainreferencetemperature');
model.material('mat3').propertyGroup.create('ThermalExpansion', 'ThermalExpansion', 'Thermal expansion');
model.material('mat3').propertyGroup('ThermalExpansion').set('dL', '(dL(T)-dL(Tempref))/(1+dL(Tempref))');
model.material('mat3').propertyGroup('ThermalExpansion').setPropertyInfo('dL', ['Reference: Specialty Coating Systems, SCS Parylene Properties brochure (2007) https://scscoatings.com/' newline 'Note: the reference temperature is 20 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (293 K), T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 290 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (563 K), limited data, calculated from the mean coefficient of thermal expansion' newline 'Reference temperature: 293.00[K]']);
model.material('mat3').propertyGroup('ThermalExpansion').set('alphatan', 'CTE(T)');
model.material('mat3').propertyGroup('ThermalExpansion').setPropertyInfo('alphatan', ['Reference: Specialty Coating Systems, SCS Parylene Properties brochure (2007) https://scscoatings.com/' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 290 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (563 K), limited data, calculated from the mean coefficient of thermal expansion']);
model.material('mat3').propertyGroup('ThermalExpansion').func.create('dL', 'Piecewise');
model.material('mat3').propertyGroup('ThermalExpansion').func('dL').set('funcname', 'dL');
model.material('mat3').propertyGroup('ThermalExpansion').func('dL').set('arg', 'T');
model.material('mat3').propertyGroup('ThermalExpansion').func('dL').set('extrap', 'constant');
model.material('mat3').propertyGroup('ThermalExpansion').func('dL').set('pieces', {'293.0' '373.16' '-0.010255+3.5E-5*T^1'});
model.material('mat3').propertyGroup('ThermalExpansion').func('dL').label('Piecewise');
model.material('mat3').propertyGroup('ThermalExpansion').func('dL').set('fununit', '');
model.material('mat3').propertyGroup('ThermalExpansion').func('dL').set('argunit', 'K');
model.material('mat3').propertyGroup('ThermalExpansion').func.create('CTE', 'Piecewise');
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE').set('funcname', 'CTE');
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE').set('arg', 'T');
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE').set('extrap', 'constant');
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE').set('pieces', {'293.0' '373.16' '3.5E-5'});
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE').label('Piecewise 1');
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE').set('fununit', '1/K');
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE').set('argunit', 'K');
model.material('mat3').propertyGroup('ThermalExpansion').addInput('temperature');
model.material('mat3').propertyGroup('ThermalExpansion').addInput('strainreferencetemperature');
% Apply custom specified permitivity
model.material('mat3').propertyGroup('def').set('relpermittivity', {'3.1'});

% Assinging Materials
model.material('mat1').selection.all;
model.material('mat2').selection.named([geomTag '_cselElec_dom']);   % all electrodes
model.material('mat3').selection.named([geomTag '_cselIns_dom']);    % all insulated parts

disp('Step Complete:Materials Assigned')
%% Electrode Boundry Probe
elecDoms = model.selection([geomTag '_cselElec_dom']).entities;   % all 9 electrode domains
elecDoms = double(elecDoms(:))';
insDoms = model.selection([geomTag '_cselIns_dom']).entities;
insDoms = double(insDoms(:))';

% all boundaries that border any insulator domain (these are the caps to drop)
tmp = mphgetadj(model, geomTag, 'boundary', 'domain', insDoms);
if iscell(tmp)
    tmp = cellfun(@(c) double(c(:))', tmp, 'UniformOutput', false);
    tmp = [tmp{:}];
end
insBnds = unique(double(tmp(:))');

for k = 1:numel(elecDoms)
    d = elecDoms(k);
    
    tmp = mphgetadj(model, geomTag, 'boundary', 'domain', d);
    if iscell(tmp)
        tmp = cellfun(@(c) double(c(:))', tmp, 'UniformOutput', false);
        tmp = [tmp{:}];
    end
    elecBnds = unique(double(tmp(:))'); % all faces of this one electrode
    activeBnds = setdiff(elecBnds, insBnds); % keep only the faces shared with tissue
    
    % explicit boundary selection for this electrode
    selTag = sprintf('ElecActive%d', k);
    model.selection.create(selTag, 'Explicit');
    model.selection(selTag).label(selTag);
    model.selection(selTag).geom(geomTag, 2);
    model.selection(selTag).set(activeBnds);
    
    % its own boundary probe
    pTag = sprintf('bnd%d', k);
    model.component(compTag).probe.create(pTag, 'Boundary');
    model.component(compTag).probe(pTag).set('expr', 'V');
    model.component(compTag).probe(pTag).set('intsurface', true);
    model.component(compTag).probe(pTag).selection.named(selTag);
end

disp('Step Complete:Boundry Probes Assigned');
%% Extract all Domains
domains = geom.getNDomains;   % Total number of domains
domainIDs = 1:domains;

% all domains in the finalized geometry
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
nElec = numel(elecDoms);   % one circuit per electrode
% Create circuit with surface conditions and front end amplifier for each ch
for i = 1:nElec
    % Unique names for each circuit
    Vname  = ['V'   num2str(i)];
    C1name = ['C1_' num2str(i)];
    R1name = ['R1_' num2str(i)];
    R2name = ['R2_' num2str(i)];
    C2name = ['C2_' num2str(i)];
    VMname = ['vm'  num2str(i)];
    GNDname= ['grnd' num2str(i)];

    % Node numbering (3 nodes per electrode)
    n1 = 3*(i-1) + 1;
    n2 = 3*(i-1) + 2;
    n3 = 3*(i-1) + 3;

    bndVar = ['bnd' num2str(i)];

    % Voltage source
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
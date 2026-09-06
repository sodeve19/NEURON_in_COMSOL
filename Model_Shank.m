
% File: Model_Shank.m
% Authors: Rana Jaylani, Sian McConchie

% Description:
% This script uses MATLAB LiveLink to create a time-dependent COMSOL model
% Creates a COMSOL model with a shank of 8 electrodes

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
ModelName = 'Model_Shank'; % -> CHANGE NAME TO MODEL REQUIRED
path = ''; % path to folder to save model
fullpath = [path '\' ModelName '.mph']; % path to folder with model name  

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
Electrode_r = 25e-6;     % Radius of the microwire
model.param.set('Electrode_r', Electrode_r);
Electrode_l = 7.5e-6;   % Height of active site
model.param.set('Electrode_l', Electrode_l);
model.param.set('elec_sa', 'pi*(Electrode_r^2)');
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

% Silicon Block of Neurnexus
model.geom(geomTag).create('SiliconBlock', 'Block');
model.geom(geomTag).feature('SiliconBlock').set('size', {'75 [um]' '15 [um]' '1600 [um]'});
model.geom(geomTag).feature('SiliconBlock').set('pos', {'Elec_xpos' 'Elec_ypos' 'Elec_zpos + 100e-6'});
model.geom(geomTag).feature('SiliconBlock').set('base', 'center');
model.geom(geomTag).feature('SiliconBlock').set('contributeto','cselIns');
model.geom(geomTag).run('SiliconBlock');
geom.run;

% Electrode Channel
model.geom(geomTag).create('Elec1', 'Cylinder');
model.geom(geomTag).feature('Elec1').set('r', 'Electrode_r');
model.geom(geomTag).feature('Elec1').set('h', 'Electrode_l');
model.geom(geomTag).feature('Elec1').set('rot', 45);
model.geom(geomTag).feature('Elec1').set('axistype', 'y');
model.geom(geomTag).feature('Elec1').set('pos', {'Elec_xpos' 'Elec_ypos' 'Elec_zpos'});
model.geom(geomTag).feature('Elec1').set('contributeto','cselElec');
model.geom(geomTag).run('Elec1');
geom.run;

% Duplicate electrode channesl into array
model.geom(geomTag).create('copy1', 'Copy');
model.geom(geomTag).feature('copy1').selection('input').set({'Elec1'});
model.geom(geomTag).feature('copy1').set('displz', 'range(-600 [um],200 [um],800 [um])');
model.geom(geomTag).feature('copy1').set('contributeto','cselElec');
model.geom(geomTag).run('copy1');
geom.run;

model.geom(geomTag).create('dif1', 'Difference');
model.geom(geomTag).feature('dif1').selection('input').set({'SiliconBlock'});
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
model.geom(geomTag).feature('uni1').selection('input').set({'copy1' 'Tissue' 'SiliconBlock'});
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

% Insulating Material - Silicon
model.material.create('mat3', 'Common');
% Apply COMSOL generated material properties
model.material('mat3').label('Silicon [solid,ideally pure]');
model.material('mat3').propertyGroup('def').set('resistivity', 'res_solid_ideally_pure_1(T)');
model.material('mat3').propertyGroup('def').setPropertyInfo('resistivity', ['Reference: A. Goldsmith, H.J. Hirschhorn, and T.E. Waterman, "Thermophysical properties of solid materials. volume II: Alloys (Melting Temperature above 1000 F)", WADC Technical Report 58-476, v2, AD253710 (1960) http://contrails.iit.edu/reports/6803; G.L. Pearson and J. Bardeen, "Electrical Properties of Pure Silicon and Silicon Alloys Containing Boron and Phosphorus", Physical Review, v75, No. 5, p865 (1949) https://doi.org/10.1103/PhysRev.75.865' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 1414 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (1687 K)']);
model.material('mat3').propertyGroup('def').set('thermalexpansioncoefficient', '(alpha_solid_1(T)+(Tempref-293[K])*if(abs(T-Tempref)>1e-3,(alpha_solid_1(T)-alpha_solid_1(Tempref))/(T-Tempref),d(alpha_solid_1(T),T)))/(1+alpha_solid_1(Tempref)*(Tempref-293[K]))');
model.material('mat3').propertyGroup('def').setPropertyInfo('thermalexpansioncoefficient', ['Reference: C.A. Swenson, "Recommended Values for the Thermal Expansivity of Silicon from 0 to 1000 K", Journal of Physical and Chemical Reference Data, v12, No. 2, p179 (1983) https://srd.nist.gov/JPCRD/jpcrd220.pdf' newline 'Note: the reference temperature is 20 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (293 K), T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 1414 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (1687 K), less than 1% error' newline 'Reference temperature: 293.00[K]']);
model.material('mat3').propertyGroup('def').set('heatcapacity', 'C_solid_1(T)');
model.material('mat3').propertyGroup('def').setPropertyInfo('heatcapacity', ['Reference: P.D. Desai, "Electrical Resistivity of Iron and Silicon", Journal of Physical and Chemical Reference Data, v15, No. 3, p967 (1986) https://srd.nist.gov/JPCRD/jpcrd298.pdf; K.K. Kelley, "Contributions to the Data on Theoretical Metallurgy, Pt. XIII", US Bureau of Mines, Bulletin No. 584 (1960) http://pbadupws.nrc.gov/docs/ML1212/ML12124A257.pdf' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 1414 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (1687 K), 1.5% to 5% error']);
model.material('mat3').propertyGroup('def').set('electricconductivity', '10E-14 [S/m]');
model.material('mat3').propertyGroup('def').setPropertyInfo('electricconductivity', ['Reference: A. Goldsmith, H.J. Hirschhorn, and T.E. Waterman, "Thermophysical properties of solid materials. volume II: Alloys (Melting Temperature above 1000 F)", WADC Technical Report 58-476, v2, AD253710 (1960) http://contrails.iit.edu/reports/6803; G.L. Pearson and J. Bardeen, "Electrical Properties of Pure Silicon and Silicon Alloys Containing Boron and Phosphorus", Physical Review, v75, No. 5, p865 (1949) https://doi.org/10.1103/PhysRev.75.865' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 1414 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (1687 K), calculated as the reciprocal of the resistivity']);
model.material('mat3').propertyGroup('def').set('HC', 'HC_solid_1(T)');
model.material('mat3').propertyGroup('def').setPropertyInfo('HC', ['Reference: P.D. Desai, "Electrical Resistivity of Iron and Silicon", Journal of Physical and Chemical Reference Data, v15, No. 3, p967 (1986) https://srd.nist.gov/JPCRD/jpcrd298.pdf; K.K. Kelley, "Contributions to the Data on Theoretical Metallurgy, Pt. XIII", US Bureau of Mines, Bulletin No. 584 (1960) http://pbadupws.nrc.gov/docs/ML1212/ML12124A257.pdf' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 1414 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (1687 K), 1.5% to 5% error']);
model.material('mat3').propertyGroup('def').set('VP', 'VP_solid_1(T)');
model.material('mat3').propertyGroup('def').setPropertyInfo('VP', ['Reference: P.D. Desai, "Electrical Resistivity of Iron and Silicon", Journal of Physical and Chemical Reference Data, v15, No. 3, p967 (1986) https://srd.nist.gov/JPCRD/jpcrd298.pdf' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 1414 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (1687 K), total of all complexes']);
model.material('mat3').propertyGroup('def').set('density', 'rho_solid_1(T)');
model.material('mat3').propertyGroup('def').setPropertyInfo('density', ['Reference: C.A. Swenson, "Recommended Values for the Thermal Expansivity of Silicon from 0 to 1000 K", Journal of Physical and Chemical Reference Data, v12, No. 2, p179 (1983) https://srd.nist.gov/JPCRD/jpcrd220.pdf' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 1414 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (1687 K), calculated from the linear expansion and the room temperature density']);
model.material('mat3').propertyGroup('def').func.create('res_solid_ideally_pure_1', 'Piecewise');
model.material('mat3').propertyGroup('def').func('res_solid_ideally_pure_1').set('funcname', 'res_solid_ideally_pure_1');
model.material('mat3').propertyGroup('def').func('res_solid_ideally_pure_1').set('arg', 'T');
model.material('mat3').propertyGroup('def').func('res_solid_ideally_pure_1').set('extrap', 'constant');
model.material('mat3').propertyGroup('def').func('res_solid_ideally_pure_1').set('pieces', {'79.0' '364.0' 'exp(2.41553243E-14*T^6 - 4.53060537E-11*T^5 + 3.33138089E-08*T^4 - 1.25961454E-05*T^3 + 2.65814575E-03*T^2 - 3.13800137E-01*T + 1.66255491E+01)';  ...
'364.0' '584.0' 'exp(-7.48706678E-12*T^5 + 1.66560544E-08*T^4 - 1.49026159E-05*T^3 + 6.73181905E-03*T^2 - 1.54364035E+00*T + 1.42728132E+02)';  ...
'584.0' '909.0' 'exp(2.60884402E-10*T^4 - 7.15630353E-07*T^3 + 7.39495950E-04*T^2 - 3.57498909E-01*T + 6.64957193E+01)';  ...
'909.0' '1073.16' 'exp(-4.13090442E-03*T - 3.07044447E+00)'});
model.material('mat3').propertyGroup('def').func('res_solid_ideally_pure_1').label('Piecewise');
model.material('mat3').propertyGroup('def').func('res_solid_ideally_pure_1').set('fununit', 'ohm*m');
model.material('mat3').propertyGroup('def').func('res_solid_ideally_pure_1').set('argunit', 'K');
model.material('mat3').propertyGroup('def').func.create('alpha_solid_1', 'Piecewise');
model.material('mat3').propertyGroup('def').func('alpha_solid_1').set('funcname', 'alpha_solid_1');
model.material('mat3').propertyGroup('def').func('alpha_solid_1').set('arg', 'T');
model.material('mat3').propertyGroup('def').func('alpha_solid_1').set('extrap', 'constant');
model.material('mat3').propertyGroup('def').func('alpha_solid_1').set('pieces', {'0.0' '30.0' '7.35637E-7+2.453566E-9*T^1+1.20482E-11*T^2';  ...
'30.0' '130.0' '7.707972079E-7+2.098318E-10*T^1+4.628581E-11*T^2+7.569451E-14*T^3-8.713366E-16*T^4';  ...
'130.0' '293.0' '-3.24590075E-7+2.257142E-8*T^1-9.684044E-11*T^2+2.835316E-13*T^3-3.440569E-16*T^4';  ...
'293.0' '1000.0' '6.794245226E-7+9.501405E-9*T^1-1.271286E-11*T^2+8.584038E-15*T^3-2.241706E-18*T^4'});
model.material('mat3').propertyGroup('def').func('alpha_solid_1').label('Piecewise 1');
model.material('mat3').propertyGroup('def').func('alpha_solid_1').set('fununit', '1/K');
model.material('mat3').propertyGroup('def').func('alpha_solid_1').set('argunit', 'K');
model.material('mat3').propertyGroup('def').func.create('C_solid_1', 'Piecewise');
model.material('mat3').propertyGroup('def').func('C_solid_1').set('funcname', 'C_solid_1');
model.material('mat3').propertyGroup('def').func('C_solid_1').set('arg', 'T');
model.material('mat3').propertyGroup('def').func('C_solid_1').set('extrap', 'constant');
model.material('mat3').propertyGroup('def').func('C_solid_1').set('pieces', {'1.0' '7.0' '-4.832181096E-5+7.68448084E-5*T^1-3.418133862E-5*T^2+2.808307076E-4*T^3-3.12897302E-7*T^4';  ...
'7.0' '20.0' '0.05251912742-0.03964814884*T^1+0.01004609362*T^2-7.812515424E-4*T^3+3.961556803E-5*T^4';  ...
'20.0' '50.0' '-1.805175785+0.7619034712*T^1-0.08653737912*T^2+0.003735361395*T^3-3.333975631E-5*T^4';  ...
'50.0' '293.0' '-82.940126+2.712235323*T^1+0.01404751222*T^2-7.977691376E-5*T^3+1.079905462E-7*T^4';  ...
'293.0' '900.0' '62.99772383+3.7706731*T^1-0.00694853616*T^2+5.9532044E-6*T^3-1.914384179E-9*T^4';  ...
'900.0' '1685.0' '769.4802262+0.1871751311*T^1-3.183959566E-5*T^2'});
model.material('mat3').propertyGroup('def').func('C_solid_1').label('Piecewise 2');
model.material('mat3').propertyGroup('def').func('C_solid_1').set('fununit', 'J/(kg*K)');
model.material('mat3').propertyGroup('def').func('C_solid_1').set('argunit', 'K');
model.material('mat3').propertyGroup('def').func.create('sigma_solid_ideally_pure_1', 'Piecewise');
model.material('mat3').propertyGroup('def').func('sigma_solid_ideally_pure_1').set('funcname', 'sigma_solid_ideally_pure_1');
model.material('mat3').propertyGroup('def').func('sigma_solid_ideally_pure_1').set('arg', 'T');
model.material('mat3').propertyGroup('def').func('sigma_solid_ideally_pure_1').set('extrap', 'constant');
model.material('mat3').propertyGroup('def').func('sigma_solid_ideally_pure_1').set('pieces', {'79.0' '364.0' '1/(exp(2.41553243E-14*T^6-4.53060537E-11*T^5+3.33138089E-08*T^4-1.25961454E-05*T^3+2.65814575E-03*T^2-3.13800137E-01*T+1.66255491E+01))';  ...
'364.0' '584.0' '1/(exp(-7.48706678E-12*T^5+1.66560544E-08*T^4-1.49026159E-05*T^3+6.73181905E-03*T^2-1.54364035E+00*T+1.42728132E+02))';  ...
'584.0' '909.0' '1/(exp(2.60884402E-10*T^4-7.15630353E-07*T^3+7.39495950E-04*T^2-3.57498909E-01*T+6.64957193E+01))';  ...
'909.0' '1073.16' '1/(exp(-4.13090442E-03*T-3.07044447E+00))'});
model.material('mat3').propertyGroup('def').func('sigma_solid_ideally_pure_1').label('Piecewise 3');
model.material('mat3').propertyGroup('def').func('sigma_solid_ideally_pure_1').set('fununit', 'S/m');
model.material('mat3').propertyGroup('def').func('sigma_solid_ideally_pure_1').set('argunit', 'K');
model.material('mat3').propertyGroup('def').func.create('HC_solid_1', 'Piecewise');
model.material('mat3').propertyGroup('def').func('HC_solid_1').set('funcname', 'HC_solid_1');
model.material('mat3').propertyGroup('def').func('HC_solid_1').set('arg', 'T');
model.material('mat3').propertyGroup('def').func('HC_solid_1').set('extrap', 'constant');
model.material('mat3').propertyGroup('def').func('HC_solid_1').set('pieces', {'1.0' '7.0' '-1.357142742E-6+2.158225189E-6*T^1-9.59999972E-7*T^2+7.887270952E-6*T^3-8.787876952E-9*T^4';  ...
'7.0' '20.0' '0.001475023955-0.001113538128*T^1+2.821495881E-4*T^2-2.19418374E-5*T^3+1.112623087E-6*T^4';  ...
'20.0' '50.0' '-0.05069927193+0.0213984404*T^1-0.002430445434*T^2+1.049094905E-4*T^3-9.363637192E-7*T^4';  ...
'50.0' '293.0' '-2.329414948+0.07617448976*T^1+3.945314097E-4*T^2-2.240575095E-6*T^3+3.032968211E-9*T^4';  ...
'293.0' '900.0' '1.391571826+0.1090406047*T^1-2.043950944E-4*T^2+1.786898954E-7*T^3-5.888574152E-11*T^4';  ...
'900.0' '1685.0' '21.99407257+0.004585166104*T^1-6.452267736E-7*T^2'});
model.material('mat3').propertyGroup('def').func('HC_solid_1').label('Piecewise 4');
model.material('mat3').propertyGroup('def').func('HC_solid_1').set('fununit', 'J/(mol*K)');
model.material('mat3').propertyGroup('def').func('HC_solid_1').set('argunit', 'K');
model.material('mat3').propertyGroup('def').func.create('VP_solid_1', 'Piecewise');
model.material('mat3').propertyGroup('def').func('VP_solid_1').set('funcname', 'VP_solid_1');
model.material('mat3').propertyGroup('def').func('VP_solid_1').set('arg', 'T');
model.material('mat3').propertyGroup('def').func('VP_solid_1').set('extrap', 'constant');
model.material('mat3').propertyGroup('def').func('VP_solid_1').set('pieces', {'293.0' '1687.0' 'exp(-28983.23815+18690.8144*log(T)^1-4934.14352*log(T)^2+662.645789*log(T)^3-45.0342817*log(T)^4+1.23403537*log(T)^5)'});
model.material('mat3').propertyGroup('def').func('VP_solid_1').label('Piecewise 5');
model.material('mat3').propertyGroup('def').func('VP_solid_1').set('fununit', 'Pa');
model.material('mat3').propertyGroup('def').func('VP_solid_1').set('argunit', 'K');
model.material('mat3').propertyGroup('def').func.create('rho_solid_1', 'Piecewise');
model.material('mat3').propertyGroup('def').func('rho_solid_1').set('funcname', 'rho_solid_1');
model.material('mat3').propertyGroup('def').func('rho_solid_1').set('arg', 'T');
model.material('mat3').propertyGroup('def').func('rho_solid_1').set('extrap', 'constant');
model.material('mat3').propertyGroup('def').func('rho_solid_1').set('pieces', {'0.0' '30.0' '2331.510015-7.113612E-5*T^1+3.674386E-6*T^2';  ...
'30.0' '130.0' '2331.593513-0.005873649*T^1+1.206114E-4*T^2-5.479876E-7*T^3+1.606517E-10*T^4';  ...
'130.0' '293.0' '2330.438322+0.02130626*T^1-9.544145E-5*T^2+4.607415E-8*T^3+4.840886E-11*T^4';  ...
'293.0' '1000.0' '2332.565+0.003839515*T^1-5.433308E-5*T^2+4.287211E-8*T^3-1.366545E-11*T^4'});
model.material('mat3').propertyGroup('def').func('rho_solid_1').label('Piecewise 6');
model.material('mat3').propertyGroup('def').func('rho_solid_1').set('fununit', 'kg/m^3');
model.material('mat3').propertyGroup('def').func('rho_solid_1').set('argunit', 'K');
model.material('mat3').propertyGroup('def').addInput('temperature');
model.material('mat3').propertyGroup('def').addInput('strainreferencetemperature');
model.material('mat3').propertyGroup.create('ThermalExpansion', 'ThermalExpansion', 'Thermal expansion');
model.material('mat3').propertyGroup('ThermalExpansion').set('dL', '(dL_solid_1(T)-dL_solid_1(Tempref))/(1+dL_solid_1(Tempref))');
model.material('mat3').propertyGroup('ThermalExpansion').setPropertyInfo('dL', ['Reference: C.A. Swenson, "Recommended Values for the Thermal Expansivity of Silicon from 0 to 1000 K", Journal of Physical and Chemical Reference Data, v12, No. 2, p179 (1983) https://srd.nist.gov/JPCRD/jpcrd220.pdf' newline 'Note: the reference temperature is 20 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (293 K), T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 1414 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (1687 K), less than 1% error' newline 'Reference temperature: 293.00[K]']);
model.material('mat3').propertyGroup('ThermalExpansion').set('alphatan', 'CTE_solid_1(T)');
model.material('mat3').propertyGroup('ThermalExpansion').setPropertyInfo('alphatan', ['Reference: C.A. Swenson, "Recommended Values for the Thermal Expansivity of Silicon from 0 to 1000 K", Journal of Physical and Chemical Reference Data, v12, No. 2, p179 (1983) https://srd.nist.gov/JPCRD/jpcrd220.pdf' newline 'Note: T' native2unicode(hex2dec({'00' '98'}), 'unicode')  native2unicode(hex2dec({'00' '9a'}), 'unicode') ' = 1414 ' native2unicode(hex2dec({'00' 'b0'}), 'unicode') 'C (1687 K), less than 1% error']);
model.material('mat3').propertyGroup('ThermalExpansion').func.create('dL_solid_1', 'Piecewise');
model.material('mat3').propertyGroup('ThermalExpansion').func('dL_solid_1').set('funcname', 'dL_solid_1');
model.material('mat3').propertyGroup('ThermalExpansion').func('dL_solid_1').set('arg', 'T');
model.material('mat3').propertyGroup('ThermalExpansion').func('dL_solid_1').set('extrap', 'constant');
model.material('mat3').propertyGroup('ThermalExpansion').func('dL_solid_1').set('pieces', {'0.0' '30.0' '-2.157945075E-4';  ...
'30.0' '130.0' '-2.275622662E-4+8.396104E-7*T^1-1.724143E-8*T^2+7.834799E-11*T^3-2.303956E-14*T^4';  ...
'130.0' '293.0' '-5.217231294E-5-3.263667E-6*T^1+1.532991E-8*T^2-1.223001E-11*T^3';  ...
'293.0' '1000.0' '-5.844526349E-4+1.124129E-6*T^1+3.311476E-9*T^2-1.161022E-12*T^3'});
model.material('mat3').propertyGroup('ThermalExpansion').func('dL_solid_1').label('Piecewise');
model.material('mat3').propertyGroup('ThermalExpansion').func('dL_solid_1').set('fununit', '');
model.material('mat3').propertyGroup('ThermalExpansion').func('dL_solid_1').set('argunit', 'K');
model.material('mat3').propertyGroup('ThermalExpansion').func.create('CTE_solid_1', 'Piecewise');
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE_solid_1').set('funcname', 'CTE_solid_1');
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE_solid_1').set('arg', 'T');
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE_solid_1').set('extrap', 'constant');
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE_solid_1').set('pieces', {'0.0' '35.0' '-1.128708E-10+4.013852E-10*T^1-1.52474E-10*T^2+1.935922E-11*T^3-8.512984E-13*T^4+9.848876E-15*T^5';  ...
'35.0' '121.0' '-4.540190868E-7+5.088778E-8*T^1-1.847617E-9*T^2+2.348401E-11*T^3-1.242463E-13*T^4+2.442997E-16*T^5';  ...
'121.0' '293.0' '8.511291633E-8-3.514147E-8*T^1+4.332654E-10*T^2-1.449023E-12*T^3+1.630929E-15*T^4';  ...
'293.0' '1000.0' '-2.623074067E-6+3.020773E-8*T^1-5.538848E-11*T^2+4.755059E-14*T^3-1.534596E-17*T^4'});
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE_solid_1').label('Piecewise 1');
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE_solid_1').set('fununit', '1/K');
model.material('mat3').propertyGroup('ThermalExpansion').func('CTE_solid_1').set('argunit', 'K');
model.material('mat3').propertyGroup('ThermalExpansion').addInput('temperature');
model.material('mat3').propertyGroup('ThermalExpansion').addInput('strainreferencetemperature');
model.material('mat3').set('family', 'custom');
model.material('mat3').set('lighting', 'cooktorrance');
model.material('mat3').set('specular', 'custom');
model.material('mat3').set('customspecular', [0.7843137254901961 1 1]);
model.material('mat3').set('fresnel', 0.9);
model.material('mat3').set('roughness', 0.1);
model.material('mat3').set('metallic', 0);
model.material('mat3').set('pearl', 0);
model.material('mat3').set('diffusewrap', 0);
model.material('mat3').set('clearcoat', 0);
model.material('mat3').set('reflectance', 0);
model.material('mat3').set('shininess', 130);
% Apply custom specified permitivity
model.material('mat3').propertyGroup('def').set('relpermittivity', {'11.7'});

% Assinging Materials
model.material('mat1').selection.all;
model.material('mat2').selection.named([geomTag '_cselElec_dom']);
model.material('mat3').selection.named([geomTag '_cselIns_dom']);

disp('Step Complete:Materials Assigned')
%% Electrodes Boundry Probe
elecBnds = model.selection([geomTag '_cselElec_bnd']).entities;
elecBnds = double(elecBnds(:));

Recordingboundries = elecBnds(end-15:end-8);
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
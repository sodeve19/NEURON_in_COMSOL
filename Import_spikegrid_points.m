
% File: Import_spikegrid_points.m
% Authors: Rana Jaylani, Sian McConchie

% Description:
% Code that imports a specified grid of point probes into a COMSOL model

% Requirements:
%   - MATLAB with COMSOL LiveLink
%   - COMSOL Multiphysics (Tested on 6.2, 6.3 and 6.4)

% Usage:
%   1. Update the file name and save path as suitable
%   2. Open the script in MATLAB LiveLink
%   3. Click Run

%% Configuration 
modelPath  = 'Model';  % -> CHANGE NAME TO MODEL REQUIRED
Zpoint = [340 320 300 280 260 240 220 200 180 160 140 120 100 80 60 40 20 0 -20 -40 -60 -80 -100 -120 -140 -160];
Ypoint = [60 40 20 0 -20 -40 -60];
compTag    = 'comp1';
geomTag    = 'geom1';
scaleFactor = 1e-6;   % because coordinates are in µm

%% Load Existing Model 
import com.comsol.model.*
import com.comsol.model.util.*
model = mphopen(modelPath);

% Ensure component and geometry exist
try model.component(compTag); catch, model.component.create(compTag,true); end
try model.component(compTag).geom(geomTag); catch, model.component(compTag).geom.create(geomTag,3); end
geom = model.component(compTag).geom(geomTag);

for m = 1:length(Ypoint)
    for n = 1:length(Zpoint)
        Rowlab = sprintf('Row_%d_%d',m,n);
        Xindex = 0*scaleFactor;
        Yindex = Ypoint(1,m)*scaleFactor;
        Zindex = Zpoint(1,n)*scaleFactor;

        model.component('comp1').geom('geom1').create(Rowlab, 'Point');
        model.component('comp1').geom('geom1').feature(Rowlab).label(Rowlab);
        model.component('comp1').geom('geom1').feature(Rowlab).setIndex('p', Xindex, 1);
        model.component('comp1').geom('geom1').feature(Rowlab).setIndex('p', Yindex, 1);
        model.component('comp1').geom('geom1').feature(Rowlab).setIndex('p', Zindex, 2);
        geom.feature(Rowlab).set('selresult', true);
        geom.run(Rowlab);

        seltag = sprintf('%s_sel',Rowlab);
        selname = sprintf('%s_sellab',Rowlab);

        model.component('comp1').geom('geom1').create(seltag, 'ExplicitSelection');
        model.component('comp1').geom('geom1').feature(seltag).selection('selection').init(0);
        model.component('comp1').geom('geom1').feature(seltag).selection('selection').set(Rowlab, 1);
        model.component('comp1').geom('geom1').run(seltag);

        prtag = sprintf('%s_pr',Rowlab);
        autoSelTag = sprintf('%s_%s', geomTag, seltag);

        model.component('comp1').probe.create(prtag, 'Point');
        model.component('comp1').probe(prtag).selection.named(autoSelTag);
        model.component('comp1').probe(prtag).label(Rowlab);
        model.component('comp1').probe(prtag).set('probename', Rowlab);
        fprintf('%s electrode added \n',Rowlab);
    end
end
geom.run('fin');

mphsave(model, modelPath);
disp('Step Complete:Electrde grid added and model saved');

% File: Import_BioNoise.m
% Authors: Rana Jaylani, Sian McConchie

% Description:
% This script uses MATLAB LiveLink to create a time-dependent COMSOL model
% Imports required density of neurons far from electrode site as points

% Requirements:
%   - MATLAB with COMSOL LiveLink
%   - COMSOL Multiphysics (Tested on 6.2, 6.3 and 6.4)

% Usage:
%   1. Update the file name and save path as suitable
%   2. Open the script in MATLAB LiveLink
%   3. Click Run

%% Starting COMSOL LiveLink
import com.comsol.model.*
import com.comsol.model.util.*

% Name to of model
modelpath = 'Model_HighDensityShank'; % -> CHANGE NAME TO MODEL REQUIRED
model = mphopen(modelpath);

% Ensure component and geometry exist
compTag   = 'comp1';
geomTag   = 'geom1';
try model.component(compTag); catch, model.component.create(compTag,true); end
try model.component(compTag).geom(geomTag); catch, model.component(compTag).geom.create(geomTag,3); end
geom = model.component('comp1').geom(geomTag);
%% Setup to model needs
start_time = 0;
step_time = 3.3e-5;  % ≈ 3.33e-5 s
stop_time = 0.01;

refractory_period = 0.003;
BionoiseCSV = 'neuron_currents_nA.csv';

% Center Location of Zone A
LocZoneA_xpos = 0;
LocZoneA_ypos = 0;
LocZoneA_zpos = 0;

Tissue_r = 0.0007 - 0.00015; % Radius minus5.5e-4 buffer for tissue
Tissue_Volume_m3 = (4/3)*pi*((Tissue_r)^3); %meter cubed
Tissue_Volume_mm3 = Tissue_Volume_m3 * 1e9; %mm cubed

ZoneA_r = 100e-6;
ZoneA_Volume_m3 = (4/3)*pi*((ZoneA_r)^3); %meter cubed
ZoneA_Volume_mm3 = ZoneA_Volume_m3 * 1e9; %mm cubed

buffer = 10e-6; %Avoid touching, currently at 10um

% Brain Congifuration
NDensity = 300000; %per cubic mm - from paper
NActive = 2/100; %percentage - from paper

P_FiringRate = 1.18; 
IN_FiringRate = 5.9;

NeuronNu_ZoneA = round(NDensity* ZoneA_Volume_mm3* NActive);
NeuronNu_Tissue = round((NDensity* Tissue_Volume_mm3* NActive) - NeuronNu_ZoneA);
%% BioNoise neurons
positions = nan(NeuronNu_Tissue, 3);
posCount = 0; %valid neurons placed
excluded_point = [LocZoneA_xpos, LocZoneA_ypos, LocZoneA_zpos];
for i = 1:NeuronNu_Tissue
    valid = false; %To check if placement is okay
    attempts = 0; %Number of tries per neuron
    while ~valid && attempts < 1000 %Either till postion found or after x attempts
        % Random point inside sphere
        r = Tissue_r * rand^(1/3);   % radius scaling for uniform dist
        theta = acos(2*rand - 1);
        phi = 2*pi*rand;

        %Convert spherical postion to cartesian
        x = r * sin(theta) * cos(phi);
        y = r * sin(theta) * sin(phi);
        z = r * cos(theta);
        candidate = [x y z];
        
        % Check if inside (not too close to boundary)
        if norm(candidate) + buffer > Tissue_r
            attempts = attempts + 1;
            continue
        end
        %Check if in zone A
        if norm(candidate - excluded_point) <= ZoneA_r
            attempts = attempts + 1;
            continue
        end
        % Check overlap with existing points
        if posCount == 0 || all(sqrt(sum((positions(1:posCount,:) - candidate).^2, 2)) > buffer)
            posCount = posCount + 1;
            positions(posCount, :) = candidate;
            valid = true;
        else
            attempts = attempts + 1;
        end
    end
    if ~valid
        break % no more space
    end
end

%Trim  the unused preallocated rows
positions = positions(1:posCount, :);

fprintf('Step Complete:Positions set for %d Points\n', size(positions,1));
%% Generate currents BEFORE geometry — filter out silent neurons
nPoints = size(positions, 1);
nIN = round(0.2 * nPoints);
nP  = nPoints - nIN;

% Build logic for true false: 1 = IN, 0 = P
choices = [ones(1, nIN), zeros(1, nP)];
choices = choices(randperm(nPoints)); % Shuffle so they arent all togther

% Preallocate currents
timeVec = (start_time:step_time:stop_time).';  % column vector
allCurrents_temp = zeros(length(timeVec), nPoints);
activeIdx = false(1, nPoints);

for i = 1:nPoints
    rate = IN_FiringRate * (choices(i) == 1) + ...
           P_FiringRate * (choices(i) ~= 1);
    [~, currentVec, spike_times] = generatePoissonSpikeTrain( ...
        rate, stop_time, step_time, refractory_period, BionoiseCSV);
    rows = min(length(timeVec), length(currentVec));
    allCurrents_temp(1:rows, i) = currentVec(1:rows);

    % Use spike_times directly — no need to scan the whole vector
    if ~isempty(spike_times)
        activeIdx(i) = true;
    end
end

% Keep only active neurons (positions and currents)
positions   = positions(activeIdx, :);
allCurrents = allCurrents_temp(:, activeIdx);
nActive     = size(positions, 1);

fprintf('Filtered: %d / %d neurons are active (have spikes)\n', nActive, nPoints);
%% Create geometry points only for active neurons
% Create temporary "all points" selection
model.component('comp1').selection.create('AllPoints_temp', 'Explicit'); 
model.selection('AllPoints_temp').label('AllPoints_temp'); 
model.component('comp1').selection('AllPoints_temp').geom(0); 
model.component('comp1').selection('AllPoints_temp').all; 

% Get the current entities (numeric vertex IDs)
AllPoints_ents = model.selection('AllPoints_temp').entities;

% Explicit slection 
model.component('comp1').selection.create('OtherPoints', 'Explicit');
model.selection('OtherPoints').label('OtherPoints');
model.component('comp1').selection('OtherPoints').geom(0);
model.component('comp1').selection('OtherPoints').set(AllPoints_ents);

for i = 1:nActive
    tag = sprintf('cell%d', i);
    model.geom(geomTag).create(tag, 'Point');
    model.geom(geomTag).feature(tag).set('selresult', true);
    model.geom(geomTag).feature(tag).set('p', positions(i,:));  % 3×1 column vector
    model.geom(geomTag).run(tag);
end

model.geom('geom1').run('fin');
OtherPoints_ents = model.selection('OtherPoints').entities;
PostBNPoints_ents = model.selection('AllPoints_temp').entities;

%Create selection for BioNoisePoints
BioNoisePoints_sel = setdiff(PostBNPoints_ents, OtherPoints_ents);
model.component('comp1').selection.create('BioNoisePoints', 'Explicit');
model.selection('BioNoisePoints').label('BioNoisePoints');
model.component('comp1').selection('BioNoisePoints').geom(geomTag, 0);  % 0 = points
model.component('comp1').selection('BioNoisePoints').set(BioNoisePoints_sel);

disp('Step Complete:Geometries of Bionoise created');
model.geom('geom1').runPre('fin');
geom.run;
disp('Step Complete:Geometry Created')
%% Electric Currents (EC)
BioNoisePoints = model.selection('BioNoisePoints').entities;
maxCols  = 900;
nBatches = ceil(nActive / maxCols);

for b = 1:nBatches
    % Indices for this batch
    idxStart = (b-1)*maxCols + 1;
    idxEnd   = min(b*maxCols, nActive);
    batchPoints   = BioNoisePoints(idxStart:idxEnd);
    batchCurrents = allCurrents(:, idxStart:idxEnd);

    % CSV for this batch
    batchData = [timeVec, batchCurrents];
    batchCSV  = fullfile(pwd, sprintf('allCurrents_batch%d.csv', b));
    writematrix(batchData, batchCSV);

    % Create COMSOL interpolation function
    funcTag = sprintf('BioNoiseInterp%d', b);
    model.func.create(funcTag, 'Interpolation');
    model.func(funcTag).set('source','file');
    model.func(funcTag).set('filename', batchCSV);
    model.func(funcTag).setIndex('argunit', 's', 0);
    model.func(funcTag).set('interp', 'linear');
    model.func(funcTag).set('extrap', 'const');

    % Assign function names for each point in this batch
    for i = 1:length(batchPoints)
        colTag   = sprintf('col%d', i+1);
        funcName = sprintf('pointID%d', batchPoints(i));
        model.func(funcTag).setEntry('columnType', colTag, 'value');
        model.func(funcTag).setEntry('funcnames', colTag, funcName);
        model.func(funcTag).setIndex('fununit', 'A', i-1);
    end

    % Assign each point to the physics using its batch function
    for i = 1:length(batchPoints)
        seg_id  = double(batchPoints(i));
        termTag = sprintf('term_%d', seg_id);
        model.physics('ec').create(termTag, 'PointCurrentSource', 0);
        model.physics('ec').feature(termTag).selection.set(seg_id);
        model.physics('ec').feature(termTag).set('Qjp', ...
            sprintf('pointID%d(t)', seg_id));
    end
end

disp('Step Complete:Electric Currents Physics added')
%% Save the model
disp('Saving updated model...');
mphsave(modelpath);
clear model

disp('Step Complete:Model Created and saved under: ')
disp(modelpath)

%% Fucntions 
function [time, signal, spike_times] = generatePoissonSpikeTrain(FiringRate, stop_time, step_time, refractory_period, spike_csv)
    % GENERATEPOISSONSPIKETRAIN Generate a spike train using a refractory-modified Poisson process
    % Inputs:
    %   FiringRate        - target firing rate in Hz
    %   stop_time         - total simulation time in seconds
    %   step_time         - time step of output signal in seconds
    %   refractory_period - refractory period in seconds
    %   spike_csv         - path to CSV containing spike waveforms
    %
    % Outputs:
    %   time       - time vector for the signal
    %   signal     - generated spike train signal (sum of spikes)
    %   spike_times- times of each spike in seconds
    
    % Load neuron spike library
    spike_library = readmatrix(spike_csv);  
    spike_waveforms = spike_library(:, 2:end);     % all spike signals
    spike_dt = spike_library(2,1) - spike_library(1,1); % actual spike sampling interval
    spike_length = size(spike_waveforms, 1);       % number of samples per spike
    
    % Pick one spike waveform randomly
    rand_col = randi(size(spike_waveforms,2));
    spike_template = spike_waveforms(:, rand_col);
    
    % Generate spike times with refractory-modified Poisson process
    spike_times = [];
    t = 0;
    while true
        isi = refractory_period + (-log(rand) / FiringRate);  % refractory + exponential
        t = t + isi;
        if t > stop_time
            break;
        end
        spike_times = [spike_times t];
    end
    
    % Create signal - initialized to zeros
    time = 0:step_time:(stop_time-step_time);
    signal = zeros(size(time));
    
    % Insert spike waveforms
    for i = 1:length(spike_times)
        start_idx = round(spike_times(i)/step_time) + 1;
        end_idx = min(start_idx + spike_length - 1, length(signal));
        
        spike_to_insert = spike_template(1:(end_idx-start_idx+1));
        signal(start_idx:end_idx) = signal(start_idx:end_idx) + spike_to_insert';
    end
    signal = signal *1e-9; %Make it in Current rather than nanoCurrent
    % % Plot results
    % figure()
    % 
    % % Continuous signal
    % subplot(2,1,1)
    % plot(time, signal, 'b', 'LineWidth', 1.2)
    % xlabel('Time (s)')
    % ylabel('Amplitude (A)')
    % title(sprintf('Poisson Spike Train with Firing Rate (%.1f Hz)', FiringRate))
    % grid on
    % 
    % % Raster plot
    % subplot(2,1,2)
    % hold on
    % for i = 1:length(spike_times)
    %     line([spike_times(i) spike_times(i)], [0 1], 'Color','r','LineWidth',1)
    % end
    % xlim([0 stop_time])
    % ylim([0 1])
    % xlabel('Time (s)')
    % ylabel('Spikes')
    % title('Spike Raster')
    % grid on
    % hold off
end


% File: COMSOLformatting.m
% Authors: Rana Jaylani, Sian McConchie

% Description:
% This script uses MATLAB script to convert the neruon currents from Gold et al format
% Outputs 2 CSVs - morphology and currents

% Requirements:
%   - MATLAB
%   - Gold et al Output for a cell

% Usage:
%   1. Update the file directory and cell prefix as needed
%   2. Click Run

%%
clear, clc

%% User settings - edit these only
dataDir     = "D:\2026_CapstonePaper\17_DetailedNeurons\d151_A3_P\nrn";
cellPrefix  = "d151_0022";   % used to match files like d151_0001_geom.dat, d151_0001_tXX.dat, d151_0001_times.dat

geomFile    = cellPrefix + "_geom.dat";
timesFile   = cellPrefix + "_times.dat";

outGeomCSV  = "segment_morphology_originalGold.csv";
outCurrCSV  = "segment_currents_originalGold.csv";

%% Read geometry file - Tmorph
geomPath = fullfile(dataDir, geomFile);

lines = readlines(geomPath);
lines = strtrim(lines);

% Keep only lines that contain numeric data
isDataLine = false(size(lines));
for i = 1:numel(lines)
    isDataLine(i) = ~isempty(sscanf(lines(i), '%f'));
end
lines = lines(isDataLine);

% Parse numeric and label after %
nLines = numel(lines);
labels = strings(nLines,1);
data   = zeros(nLines,10);   % expecting 10 numeric columns

for i = 1:nLines
    L = lines(i);

    % Label after %
    if contains(L, "%")
        labels(i) = strtrim(extractAfter(L, "%"));
    end

    % Numeric values
    nums = sscanf(L, '%f');
    k = min(numel(nums), 10);
    data(i,1:k) = nums(1:k);
end

% Build table
Tmorph = array2table(data, 'VariableNames', ...
    {'xs','ys','zs','ds','as','xe','ye','ze','de','ae'});
Tmorph = addvars(Tmorph, (1:height(Tmorph))', 'Before', 1, 'NewVariableNames', 'seg_id');
Tmorph.segment = categorical(labels);

% segment length (um)
dx = Tmorph.xe - Tmorph.xs;
dy = Tmorph.ye - Tmorph.ys;
dz = Tmorph.ze - Tmorph.zs;
Tmorph.length = sqrt(dx.^2 + dy.^2 + dz.^2);

% Save geometry CSV
writetable(Tmorph, fullfile(dataDir, outGeomCSV));

%% Read all current files - Tcurr
% Match any .dat for this cellPrefix, then filter ONLY those with _t<digit>
files = dir(fullfile(dataDir, cellPrefix + "_*.dat"));

isValid = false(numel(files),1);
for i = 1:numel(files)
    isValid(i) = ~isempty(regexp(files(i).name, '_t\d', 'once'));
end
files = files(isValid);

old_names = {files.name}';
extract_numbers = string(regexp(old_names, "\d{2,3}\.\d{3}", 'match'));
as_number = str2double(extract_numbers);

% Sort by filename
%[~, idx] = sort({files.name});
[~, idx] = sort(as_number);
files = files(idx);

if isempty(files)
    error("No current files found matching pattern %s in %s", cellPrefix + "_t*.dat", dataDir);
end

nFiles = numel(files);

% Load first file to get number of segments
firstVec = load(fullfile(dataDir, files(1).name));
nSeg = numel(firstVec);

% Preallocate: rows=time, cols=segment
M = zeros(nFiles, nSeg);

for i = 1:nFiles
    v = load(fullfile(dataDir, files(i).name));

    if numel(v) ~= nSeg
        error("File %s has inconsistent length (expected %d, got %d).", files(i).name, nSeg, numel(v));
    end

    M(i,:) = v(:).';   % force row
end

segNames = "seg_" + string(1:nSeg);
Tcurr = array2table(M, 'VariableNames', segNames);

%% Convert A to pA, dvide by length for pA/um
L = Tmorph.length(:);   % (N x 1), um

if nSeg ~= numel(L)
    error("Mismatch: current files have %d segments but geom file has %d segments.", nSeg, numel(L));
end
if any(L <= 0)
    error("Found zero/negative lengths in Tmorph.length. Can't divide by length.");
end

I = table2array(Tcurr);                         % nA
I_pA_per_um = (I) ./ (L.');              % nA/um (broadcast across columns)

Tcurr_pA_per_um = array2table(I_pA_per_um, 'VariableNames', Tcurr.Properties.VariableNames);

% Rename columns to requested format
Tcurr_pA_per_um.Properties.VariableNames = "I_pA_um_seg" + string(1:width(Tcurr_pA_per_um));

%% Add time_ms
timePath = fullfile(dataDir, timesFile);
time_ms = load(timePath);
time_ms = time_ms(:);

if height(Tcurr_pA_per_um) ~= numel(time_ms)
    error("Time vector length (%d) does not match number of time files (%d).", numel(time_ms), height(Tcurr_pA_per_um));
end

Tcurr_pA_per_um = addvars(Tcurr_pA_per_um, time_ms, 'Before', 1, 'NewVariableNames', 'time_ms');

% Save final CSV
writetable(Tcurr_pA_per_um, fullfile(dataDir, outCurrCSV));


% File: Plot_spikegrid_points.m
% Authors: Rana Jaylani, Sian McConchie

% Description:
% Takes output csv from COMSOL and plots the spike grid on the neuron

% Requirements:
%   - MATLAB
%   csv of ouputed voltage from COMSOL model using Import_spikegrid_points

% Usage:
%   1. Update the file name
%   2. Click Run

%% EAP grid from CSV with BUCKET-NORMALIsED scaling to a constant legend height
dataFile = 'BlankModel.csv'; % -> CHANGE NAME TO MODEL REQUIRED

%% Parameters 
nRows = 7; nCols = 26;
gridSize = 20;                 % µm between electrodes
somaRow = 4; somaCol = 18;      % skip this electrode
start_time_ms = 5.4; end_time_ms = 9.4;
xgrid_plot_pcnt = 0.8; eap_line_width = 0.6;

% Buckets (in V) and colors. Each bucket will be scaled to SAME display height.
bucketVals   = [10 20 40 80 160]*1e-6;   % V (µV -> V)
bucketColors = {'#656FE9','#6c05a2','#1aa7a1','#f5963e','#C2437A'};

% Display height for EVERY bucket (how tall the legend bar is, and the p2p height on the grid)
legend_cell_height = 3*gridSize/4; % set to gridSize/2 if you want half-height traces

%% Load CSV
T = readtable(dataFile, 'PreserveVariableNames', true, 'NumHeaderLines', 7);
names = T.Properties.VariableNames;
t_ms = T{:,1} * 1e3;
mask_t = (t_ms >= start_time_ms) & (t_ms <= end_time_ms);
plot_times = t_ms(mask_t);

%% Map Columns to (row,col), skip soma
RC = [];
RC_all = [];

for k = 2:numel(names)
    raw = names{k};
    tok = regexp(raw, 'Row[_\s]*(\d+)[_\s]*(\d+)', 'tokens', 'once');
    if isempty(tok)
        raw2 = regexprep(raw, '[^\w\s;_]', '');
        tok = regexp(raw2, 'Row[_\s]*(\d+)[_\s]*(\d+)', 'tokens', 'once');
    end
    if isempty(tok), continue; end
    r = str2double(tok{1}); 
    c = str2double(tok{2});

    RC_all = [RC_all; r, c, k];   % includes soma

    if r==somaRow && c==somaCol, continue; end
    RC = [RC; r, c, k];
end
if isempty(RC), error('No electrode columns matched "Row_r_c" (after skipping soma).'); end

%% canvas
x_min = -340;
x_max = 160;
y_min = -60;
y_max = 60;
xyMax = [x_min, x_max, y_min, y_max];

main_fig_h = figure('Color','w'); 
main_ax_h = axes('Parent', main_fig_h);
axis(main_ax_h, [xyMax(1)-gridSize, xyMax(2)+gridSize, xyMax(3)-gridSize, xyMax(4)+gridSize]);
hold(main_ax_h,'on');

baseFontSize = 8;
set(main_ax_h,'XTick',xyMax(1):gridSize:xyMax(2), ...
              'YTick',xyMax(3):gridSize:xyMax(4), ...
              'DataAspectRatio',[1 1 1], ...
              'Color','none','FontSize',baseFontSize, ...
              'Layer', 'Top');
xlabel(main_ax_h,'\mu m', 'FontSize',baseFontSize); 
ylabel(main_ax_h,'\mu m', 'FontSize',baseFontSize);
title(main_ax_h, 'Point Electrode Grid of Pyramidal Cell, Parameter A');
%% Plot Background neuron
% Read the CSV file
data = readtable('segment_morphology_originalGold.csv');
% Loop through each row (each segment)
for i = 1:height(data)
    % Extract start and end 3D coordinates
    xs = data.xs(i);
    ys = data.ys(i);
    zs = data.zs(i);
    
    xe = data.xe(i);
    ye = data.ye(i);
    ze = data.ze(i);

    x_proj = [xs, xe];
    y_proj = [ys, ye];
    
    % Plot each line separately
    plot(x_proj, y_proj, 'LineWidth', eap_line_width -0.05, 'Color', '#d8d8d8');
end
%% Plot each electrode (bucket-normalized scaling)
for i = 1:size(RC,1)
    r = RC(i,1); c = RC(i,2); colIdx = RC(i,3);

    v = T{mask_t, colIdx};      % volts
    if ~any(isfinite(v)), continue; end

    % Peak-to-peak and bucket selection
    pk = max(v) - min(v);
    idx = find(bucketVals >= pk, 1, 'first');
    if isempty(idx), idx = numel(bucketVals); end
    S_bucket = bucketVals(idx);          % this trace's bucket (V)
    clr      = bucketColors{idx};

    % Per-bucket scaling: make p2p == legend_cell_height
    y_rescale_trace = legend_cell_height / max(S_bucket, eps);   % µm per Volt

    % Position in grid (Row_1_* at bottom)
    x_um = (c-1)*gridSize + xyMax(1);
    y_um = (nRows - r)*gridSize + xyMax(3);

    % Time scaling within cell width
    x_rescale = (gridSize / max(end_time_ms - start_time_ms, eps)) * xgrid_plot_pcnt;

    new_xs = (plot_times * x_rescale) + x_um - start_time_ms*x_rescale - gridSize/2;
    new_ys = (v * y_rescale_trace)    + y_um;
    
    plot(main_ax_h, new_xs, new_ys, 'Color', clr, 'LineWidth', eap_line_width);
end

%% Time Scale
line(main_ax_h,[xyMax(1)-gridSize/2, xyMax(1)+gridSize*xgrid_plot_pcnt-gridSize/2], ...
                 [xyMax(3)-gridSize/1.8,  xyMax(3)-gridSize/1.8], 'Color','black');
text(main_ax_h, xyMax(1)-gridSize/2, xyMax(3)-gridSize/1.3, ...
     sprintf(' %.2f ms', (end_time_ms - start_time_ms)), ...
     'Color','black','FontUnits','normalized','FontSize',0.04);

%% Voltage legend (constant height that defines bucket-normalization)
xGrid    = xyMax(1):gridSize:xyMax(2);
textGrid = max(1, floor(numel(xGrid)/numel(bucketVals)));
y_top    = xyMax(3) - gridSize/2;
bar_h    = legend_cell_height;     % same height for every bucket

for k = 1:numel(bucketVals)
    spacingFactor = 0.5;   % try 1.2–1.8 for more or less space
    x_pos = xyMax(1) - gridSize/4 + gridSize*k*textGrid*spacingFactor;    line(main_ax_h, [x_pos x_pos], [y_top - bar_h, y_top], ...
         'Color', bucketColors{k}, 'LineWidth', 2);
    text(main_ax_h, x_pos + gridSize/4, y_top - bar_h/2, ...
         sprintf('=%.1f \\muV', 1e6*bucketVals(k)), ...
         'Color', bucketColors{k}, 'FontUnits','normalized','FontSize',0.04, ...
         'HorizontalAlignment','left','VerticalAlignment','middle','Interpreter','tex');
end

legend_offset = gridSize ; % move legend further down (1.0 = same, 1.5 = lower, 2 = even lower)
axis(main_ax_h, [x_min-20 x_max+20 y_min-30 y_max+20]);

top_pad    = gridSize*0.5;                    % small space above the grid
legend_off = gridSize*0.1;                    % how far below the grid the legend sits
bottom_pad = legend_off + legend_cell_height + gridSize*0.8;  % room for legend + time bar

% also remove default outer whitespace around axes
set(main_ax_h, 'LooseInset', get(main_ax_h,'TightInset'));
T = readtable(dataFile, 'PreserveVariableNames', true);

% Preallocate output
colNames = T.Properties.VariableNames;
p2pVals = nan(1, numel(colNames)-1);   % skip time column

for k = 2:numel(colNames)
    v = T{:, k};
    v = v(isfinite(v));    % remove NaNs if any
    if isempty(v)
        p2pVals(k-1) = NaN;
    else
        p2pVals(k-1) = max(v) - min(v);
    end
end

% Display results in a readable table
Result = table(colNames(2:end)', p2pVals', ...
    'VariableNames', {'ColumnName', 'PeakToPeak_V'});
% Optional: convert to µV for clarity
Result.PeakToPeak_uV = Result.PeakToPeak_V * 1e6;
disp(Result(:, {'ColumnName', 'PeakToPeak_uV'}));

set(main_fig_h, 'Units', 'centimeters', 'Position', [2, 2, 18, 14]); % width × height

% Saving png output of Figure
outFile = 'EAP_grid.png'; % -> CHANGE NAME TO MODEL REQUIRED
exportgraphics(main_fig_h, outFile,'Resolution', 600);
%exportgraphics(main_fig_h, 'EAP_grid_transparent.png', 'BackgroundColor', 'none', 'ContentType', 'auto');

% Shun_analyzeRTPP

% 2025/01/31

% Plots trajectory of RTPP centroid and performs other analysis

%% Load the CSV file

clear; close all;
% addpath(genpath(osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Analysis/NeuroDAP/Methods')));
addpath(genpath("C:\Users\sallyx\Documents\GitHub\NeuroDAP\"));
addpath(genpath("C:\Users\sallyx\HMS Dropbox\Jia Yin Xiao\ForSally"));
addpath("C:\Users\sallyx\Documents\MATLAB\slanCM")
[twoColors,~,~,~,~,~,bluePurpleRed] = loadColors;

% filename = uipickfiles('FilterSpec',osPathSwitch('/Volumes/Neurobio/MICROSCOPE/Shun/Project misc/Recordings'),...
%                         'Prompt','Select an date folder');
filename = uipickfiles('FilterSpec','C:\Users\sallyx\HMS Dropbox\Jia Yin Xiao\ForSally', ...
    'Prompt','Select an date folder');

sessionList = dir(filename{1});
sessionList = sessionList(~ismember({sessionList.name},{'.','..'}));
summaryMask = contains({sessionList.name}, {'Summary','3mW','8mW'});  %Change for each experiment!!!
sessionList = sessionList(~summaryMask);
nSessions = length(sessionList);
mid_point = 350; %find the mid point
extractedNames = regexp(sessionList(1).folder, '(?<=ForSally\\)[^\\]+', 'match', 'once'); %for figure title
disp(extractedNames)
%% Process data

session_cutoffs = [NaN, NaN; NaN, NaN; NaN, NaN; NaN, NaN];%[0, 0, 0, 0]; % in min; 10.56
Fs = 20; % frame rate

sessions = struct([]);
output_results = struct([]);
for s = 1:nSessions
    cur_session = dir(fullfile(sessionList(s).folder,sessionList(s).name));
    tableName = cur_session(contains({cur_session.name},'times-')).name;
    tablePath = cur_session(contains({cur_session.name},'times-')).folder;
    cur_data = readtable(fullfile(tablePath,tableName));

    dirsplit = split(tablePath,filesep);
    sessions(s).name = dirsplit{end};
    sessions(s).path = tablePath;
    sessions(s).data = cur_data;
    sessions(s).cutoff = session_cutoffs(s,:);
    % sessions(s).stimside = {};
    % sessions(s).totalT = {};
    % sessions(s).updata ={};
    output_results(s).name = dirsplit{end};

    % Add time column
    sessions(s).data.Item1 = datetime(sessions(s).data.Item1, ...
        'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSSSSSZ','TimeZone', 'UTC');
    startTime = sessions(s).data.Item1(1);
    minuteSinceStart = minutes(sessions(s).data.Item1 - startTime);
    sessions(s).data.time = minuteSinceStart;

    % Align to cutoff
    % if session_cutoffs(s) >= 0; offset = session_cutoffs(s);
    % else; offset = max(minuteSinceStart) + session_cutoffs(s);
    % end
    % sessions(s).data.time = minuteSinceStart - offset;

    %Detect stim side if controlled mid point
    idx_true = find(strcmp(cur_data.Item3, 'True'));
    if median(cur_data{idx_true,3}) <= mid_point
        sessions(s).stimside = 'left';
    else
        sessions(s).stimside = 'right'; 
    end

end

disp('Finished: animal data loaded');

%% Plot summary figure

initializeFig(0.5,1); tiledlayout(1,nSessions+1);
%stim_side = 'left';
noMovementThreshold = 20;

leftColor = [156, 219, 17]./255;
rightColor = [144, 126, 171]./255;
stimColor = [7, 162, 222]./255;
ctrlColor = [.3 .3 .3];



for s = 1:nSessions
    %cur_data = sessions(s).data(sessions(s).data.time >= 0,:);
    if ~any(isnan(sessions(s).cutoff))
        session_filter_start = session_cutoffs(s,1);
        session_filter_end = session_cutoffs(s,2);
        cur_data = sessions(s).data(sessions(s).data.time >= session_filter_start & sessions(s).data.time <= session_filter_end,:);
    else
        cur_data = sessions(s).data;
    end
    cur_name = sessions(s).name;

    if ~contains(cur_name,'RTPP')
        stim_side = sessions(s+1).stimside;
    else
    stim_side = sessions(s).stimside;
    end

    output_results(s).stimside = stim_side;

    X_raw = cur_data.Item2_X;
    Y_raw = cur_data.Item2_Y;
    Y_midpoint = mid_point; %(min(Y_raw)+max(Y_raw))/2;

    % Drop potential sleep time
    % staticWindow = getStaticPeriod(X_raw, Y_raw,noMovementThreshold=3,windowDuration=30);
    % cur_data_clean = cur_data(~staticWindow,:);
    % X = X_raw(~staticWindow);
    % Y = Y_raw(~staticWindow);

    %Check if the animal start at the stim side
    isTrue = strcmp(cur_data.Item3, 'True');
    streak = conv(double(isTrue), ones(5,1), 'valid') == 5;
    idx_r = find(streak, 1);
    if contains(sessions(s).name,'RTPP') & ~isempty(idx_r)
        X = X_raw(idx_r:end);
        Y = Y_raw(idx_r:end);
        updated_data = cur_data(idx_r:end,:);   
    else
        X = X_raw;
        Y = Y_raw;
        updated_data = cur_data;
    end

    % Check if there is background contaminating recognition
    % Create logical masks for both conditions; only affect X Y plotting
    if strcmp(sessions(s).stimside, 'right')
        cond1 = strcmp(updated_data.Item3, 'True')  & updated_data.Item2_Y < mid_point;
        cond2 = strcmp(updated_data.Item3, 'False') & updated_data.Item2_Y > mid_point;
    else
        cond1 = strcmp(updated_data.Item3, 'True')  & updated_data.Item2_Y > mid_point;
        cond2 = strcmp(updated_data.Item3, 'False') & updated_data.Item2_Y < mid_point;
    end
    rows_to_remove = cond1 | cond2;
    if contains(cur_name,'RTPP')
        X = X(~rows_to_remove);
        Y = Y(~rows_to_remove);
    end
    %

    sessions(s).totalT = length(Y);
    sessions(s).updata = updated_data(~rows_to_remove,:);

    % X = X_raw; Y = Y_raw;
    
    % Plot the trajectory
    nexttile;
    if ~contains(cur_name,'RTPP')
        plot(X(Y>=Y_midpoint), Y(Y>=Y_midpoint), Color=rightColor, LineWidth=2); hold on;
        plot(X(Y<Y_midpoint), Y(Y<Y_midpoint), Color=leftColor, LineWidth=2);
        legend({'Right', 'Left'}, 'Location', 'northeast');
    elseif strcmpi(stim_side,'right')
        plot(X(Y<=Y_midpoint), Y(Y<=Y_midpoint), Color=ctrlColor, LineWidth=2); hold on;
        plot(X(Y>Y_midpoint), Y(Y>Y_midpoint), Color=stimColor, LineWidth=2);
        legend({'Stim OFF', 'Stim ON'}, 'Location', 'northeast');
    elseif strcmpi(stim_side,'left')
        plot(X(Y<=Y_midpoint), Y(Y<=Y_midpoint), Color=stimColor, LineWidth=2); hold on;
        plot(X(Y>Y_midpoint), Y(Y>Y_midpoint), Color=ctrlColor, LineWidth=2);
        legend({'Stim ON', 'Stim OFF'}, 'Location', 'northeast');
    end
    
    xlim([300,650]); xlabel('X Position');
    ylim([0,700]);ylabel('Y Position');
    title(cur_name);

    % Plot time in each chamber
    if ~contains(cur_name,'RTPP')
        if strcmpi(stim_side,'right')
            stim = length(find(Y>=Y_midpoint));
        elseif strcmpi(stim_side,'left')
            stim = length(find(Y<=Y_midpoint));
        end
    else
        stim = sum(strcmp(updated_data.Item3,'True'));
    end
    
    stim_pct(s) = stim/length(updated_data.Item3) * 100;

    right_side = find(Y>=Y_midpoint);
    right_pct(s) = length(right_side)/length(Y)*100;


    output_results(s).stim_percentage = stim/length(updated_data.Item3) * 100;
    output_results(s).right_percentage = length(right_side)/length(Y)*100;
    
end

% Plot bar plot

%nexttile;
figure(2)
set(figure(2), 'Units', 'pixels', 'Position', [100, 100, 500, 800]);
mergedLabels = cellfun(@(n, s) [n, ' - ', s], {sessions.name}, {sessions.stimside}, 'UniformOutput', false);
for s = 1:nSessions
    if ~contains(sessions(s).name,'RTPP'); color = ctrlColor;
    else; color = stimColor; end
    plotScatterBar(s,stim_pct(s),Color=color,style='bar'); hold on
end
xticks(1:nSessions); xticklabels(mergedLabels);
ylim([0,100])
ylabel('Time spent in stimulated side (%)');


title(extractedNames)

%% Plot left and right for direct comparison (conditioned)
figure(3)
set(figure(3), 'Units', 'pixels', 'Position', [100, 100, 1000, 800]);
stackData = [100 - right_pct; right_pct]'; 
b = bar(stackData, 'stacked');
b(1).FaceColor =[0.9290 0.6940 0.1250];
b(2).FaceColor = [0 0.4470 0.7410];  
ylim([0 100]);
ylabel('Percentage (%)');
xticks(1:nSessions); xticklabels(mergedLabels);
legend({'Left side', 'Right side'}, 'Location', 'northeastoutside');
title(extractedNames)


%% Heatmap (Individual heatmap skip this one)
filenames = {sessions.name};
fileIndex = 1; % Change for each plot

X = sessions(fileIndex).updata.Item2_X;
Y = sessions(fileIndex).updata.Item2_Y;

% Define heatmap grid
squareSize = 15; % in pixels


edgesX = linspace(300, 640, round(range(X)/squareSize));
edgesY = linspace(0, 700, round(range(Y)/squareSize));

% Count occurences of points in each grid bin
[counts, edgesX, edgesY] = histcounts2(X, Y, edgesX, edgesY, 'Normalization','probability');
% timeInBins = counts * 1/Fs;
% log_counts = log(counts);
% log_times = log(timeInBins);
% gamma_counts = counts .^ 0.5;

% counts(counts >= 0.1) = 0.1;

%4. plot
figure;
% imagesc(edgesX(1:end-1), edgesY(1:end-1), counts');%gamma_counts');
cmap = turbo(); % your custom colormap function
colormap(cmap);
imagesc(edgesX(1:end-1), edgesY(1:end-1), counts');
set(gca, 'YDir', 'normal'); 
% ax = gca;
% ax.Colormap = cmap;
% ax.Color = 'w'; % background color
colorbar;

%cap the color bar
%clim([0 0.3]); 

% set(gca,'ColorScale','log');
%colormap(sky); colorbar;
hold on;
plot([edgesX(1) edgesX(end)], [350 350], 'w--', 'LineWidth', 2);  % draw the line
text(edgesX(1) - 10, 350, 'Separation', ...
    'Color', 'black', 'FontSize', 12, 'HorizontalAlignment', 'right');
%hold on;
%plot(X, Y, '-', 'LineWidth', 1.5, 'Color', [1 1 1 0.2]);

hold off;
xlim([edgesX(1) edgesX(end-1)]); xlabel('X Position');
ylim([edgesY(1) edgesY(end-1)]); ylabel('Y Position');
title (['Heatmap plot of ', filenames{fileIndex}]);
axis equal tight;
axis off;

filepath = fullfile(sessionList(fileIndex).folder,sessionList(fileIndex).name);
svg_name = 'Heatmap Plot.svg';
tif_name = 'Heatmap Plot.tif';
svg_path = fullfile(filepath,svg_name);
tif_path = fullfile(filepath,tif_name);

%exportgraphics(gcf, svg_path, 'ContentType', 'vector');
print(svg_path, '-dsvg');
exportgraphics(gcf, tif_path, 'Resolution', 300);
%% Heatmaps for all sessions
squareSize = 15;
edgesX = linspace(300, 640, round((640 - 300) / squareSize));
edgesY = linspace(0, 700, round(700 / squareSize));

numSessions = length(sessions);
counts_all = cell(numSessions, 1);

% First pass: compute and store heatmaps
for i = 1:numSessions
    X = sessions(i).updata.Item2_X;
    Y = sessions(i).updata.Item2_Y;

    counts = histcounts2(X, Y, edgesX, edgesY, 'Normalization', 'probability');
    counts_all{i} = counts;
end

% Find global max for capping
all_values = cell2mat(cellfun(@(x) x(:), counts_all, 'UniformOutput', false));
cmax = prctile(all_values, 99.5);  % Use 99.5th percentile as a stable cap (avoids outliers)
% Or set manually: cmax = 0.1;

% Create tiled layout
figure(4);
set(figure(4), 'Units', 'pixels', 'Position', [100, 100, 1200, 800]); 
tiledlayout('flow');
cmap = turbo();

for i = 1:numSessions
    nexttile;

    counts = counts_all{i};
    imagesc(edgesX(1:end-1), edgesY(1:end-1), counts');
    set(gca, 'YDir', 'normal');
    colormap(cmap);
    caxis([0 cmax]);%0.2]);  % Shared color scale
    t = title(sessions(i).name, 'Interpreter', 'none');
    t.Units = 'normalized';
    t.Position(2) = t.Position(2) + 0.1;
    axis equal tight;
    axis off;

    % Optional: add white dashed line
    hold on;
    plot([edgesX(1) edgesX(end)], [mid_point mid_point], 'w--', 'LineWidth', 2);
    hold off;
end

% Add shared colorbar
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'Probability';

%% Save images

figHandles = [figure(1), figure(2), figure(3),figure(4)];  % replace with your actual figure variables
figureNames = {
    'Track plot';
    'Percentage of time on the stimulated side bar plot'; 
    'Percentage of time left and right bar plot'; 
    'Heatmap plot'
};

filepath = fullfile(sessionList(1).folder,'/Summary');
for i = 1:length(figHandles)
    fig = figHandles(i);
    name = figureNames{i};

    % Define output folder
    % Create file names with appropriate extensions
    svg_path = fullfile(filepath, [name, '.svg']);
    tif_path = fullfile(filepath, [name, '.tif']);

    % Ensure folder exists
    if ~exist(filepath, 'dir')
        mkdir(filepath);
    end

    % Save files
    print(fig, svg_path, '-dsvg');
    exportgraphics(fig, tif_path, 'Resolution', 300);
end

filename_output = 'output_results.mat';
save(fullfile(filepath, filename_output), 'output_results');
%% Velocity

% velocity calculations
dx = diff(X);
dy = diff(Y);
dt = 1/Fs;
velocity = sqrt(dx.^2 + dy.^2) ./ dt;

% normalize velocity
velocity_norm = (velocity - min(velocity)) / (max(velocity) - min(velocity));

% colormap
cmap = sky;
colors = interp1(linspace(0, 1, size(cmap, 1)), cmap, velocity_norm);

% make plot
figure;
hold on;

% color gradient
for i = 1:length(velocity)
    plot(X(i:i+1), Y(i:i+1), 'Color', colors(i, :), 'LineWidth', 2)
end

% 7. axis label and titles
xlabel('X-Position')
ylabel('Y-Position')
title('RTPP Practice Analysis')
clim([min(velocity) max(velocity)]); % Set the color axis based on velocity range
colormap(sky);
h = colorbar;  % Create a colorbar
ylabel(h, 'Velocity (units per second)');

hold off;

%% (Test) Show removed static periods

windowDuration = 30;
noMovementThrehsold = 5;
staticWindow = getStaticPeriod(X_raw, Y_raw,...
    noMovementThreshold=noMovementThrehsold,windowDuration=windowDuration);

figure; tiledlayout(1,2);
nexttile;
% Plot the full trajectory in blue
plot(X_raw, Y_raw, 'b.-');  
hold on;
% Overlay the points marked for removal in red circles
plot(X_raw(staticWindow), Y_raw(staticWindow), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
xlabel('X Position');
ylabel('Y Position');
legend('Recorded Trajectory','Removed Data');
title('Animal Trajectory with Removed Data Highlighted');

nexttile;
movingWindowInSamples = windowDuration * 20;
diffXinWindow = movsum(diff(X_raw),movingWindowInSamples);
diffYinWindow = movsum(diff(Y_raw),movingWindowInSamples);
distanceInWindow = diffXinWindow + diffYinWindow;
time = (1:length(distanceInWindow))./20 / 60;
plot(time, distanceInWindow); hold on;
yline(50,'r',LineWidth=2); hold on 
yline(-50,'r',LineWidth=2); hold on; 
plot(time(staticWindow),distanceInWindow(staticWindow),'r');
xlabel('Time (min)');
ylabel('Distane travel in 30s moving window');


%% (Test) Show static period criteria on result

cur_data = sessions(2).data(sessions(2).data.time >= 0,:);
for s = 1:100
    % Drop potential sleep time
    staticWindow = getStaticPeriod(X_raw,Y_raw,noMovementThreshold=s,windowDuration=30);

    Y_raw = cur_data.Item2_Y;
    Y = Y_raw(~staticWindow);
    stim = find(Y<=Y_midpoint);
    stim_pct(s) = length(stim)/length(Y) * 100;
end

plot(stim_pct)
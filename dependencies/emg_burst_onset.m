%% EVA EMG Sliding Average – EMG Burst Onset
% University of Rochester School of Medicine & Dentistry
% Movement and Plasticity Lab (MAPL)
% Author: Aaron Huynh | PI: Dr. Ania Busza
% Collaborators:
%   Dr. Ania Busza
%   Dr. David Cunningham (Case Western University)
% Last commit: 02 December 2024 ANH

function emg_burst_onset(initials, subject, affected_arm, timepoint, recording_arm, PlotOnsets)
% v1 gets the onset based on the threshold of the minimum baseline average
%% CHANGE THESE
% Define baseline window size and step size
% sr * window (s) = # of samples
window_ms = 500; % ms
step_ms = 50; % ms
startSearchTime = -1.250;

aFq = 10; % lower bound of filter (Hz)
bFq = 500; % upper bound of filter (Hz)
envelope_lowpassHz = 20; % low-pass filter cut-off for EMG envelope

thresholdSize = 2; % how many SDs above the baseline mean to determine onset threshold
dur = 0.030; % s; duration signal needs to be held 2SD above baseline average to determine onset
dur2 = 0.070; % s; subsequent duration signal needs to be held 1SD above baseline average to confirm onset

%% ---------------------------------------------------------------------%%
%% Set Analysis Directory
%% ---------------------------------------------------------------------%%
root = pwd;
% parent_dir = fileparts(root);  % One level up from root
conditions_dir = sprintf('%s/Conditions (-3-3)', root);

% Get today's date
today = datetime('today');
formattedDate = datestr(today, 'yyyymmdd'); %#ok<DATST> % format date
formattedDate = char(formattedDate); % convert to string

% Determine affected and unaffected directories
% Read patient log to determine this
affected_arm_lower = lower(affected_arm);
if strcmp(affected_arm_lower, 'l')
    unaffected_arm_lower = 'r';
elseif strcmp(affected_arm_lower, 'r')
    unaffected_arm_lower = 'l';
else
    error("Incorrect affected arm input. Enter 'R' or 'L'.")
end
recording_arm_lower = lower(recording_arm);
analysis_dir = sprintf('%s/%s_%s_EMG-Burst-Onset_%s_%s (%dms Window, %dSD)', root, timepoint, subject, formattedDate, initials, window_ms, thresholdSize);

% if strcmp(affected_arm_lower, 'l')
%     affected_dir = sprintf('%s/Left (Affected)', analysis_dir);
%     unaffected_dir = sprintf('%s/Right (Unaffected)', analysis_dir);
% elseif strcmp(affected_arm_lower, 'r')
%     affected_dir = sprintf('%s/Right (Affected)', analysis_dir);
%     unaffected_dir = sprintf('%s/Left (Unaffected)', analysis_dir);
% else
%     error('\nInvalid input argument for ''affected_arm''.\nEnter ''L'' for Left arm or ''R'' for Right arm');
% end
% 
% if strcmp(recording_arm_lower, affected_arm_lower)
%     recording_dir = affected_dir;
%     recording_analysis_dir = affected_dir;
%     figures_dir = sprintf('%s/Figures', affected_dir);
%     % raw_dir = sprintf('%s/Raw', affected_dir);
%     % processed_dir = sprintf('%s/Processed', affected_dir);
% elseif strcmp(recording_arm_lower, unaffected_arm_lower)
%     recording_dir = unaffected_dir;
%     recording_analysis_dir = unaffected_dir;
%     figures_dir = sprintf('%s/Figures', unaffected_dir);
%     % raw_dir = sprintf('%s/Raw', unaffected_dir);
%     % processed_dir = sprintf('%s/Processed', unaffected_dir);
% end

figures_dir = sprintf('%s/Figures', analysis_dir);

% List of directories to create
% dir_list = {figures_dir, raw_dir, processed_dir, recording_dir, recording_analysis_dir, analysis_dir};
dir_list = {figures_dir, analysis_dir};

% Loop through and create directories if they do not exist
fprintf("\nCreating directories...")
for dir = 1:length(dir_list)
    % if the directory doesn't already exist, create it
    if ~exist(dir_list{dir}, 'dir')
        mkdir(dir_list{dir});
    end
    % add to path to access for analysis
    addpath(dir_list{dir})
end
fprintf("\nDirectories created!\n")

% save path so you don't have to add directories each time
addpath(conditions_dir)
savepath

%% ---------------------------------------------------------------------%%
%% Load Condition File
%% ---------------------------------------------------------------------%%
condition_file = sprintf('%s/%s_%s_%s', conditions_dir, subject, timepoint, recording_arm);
load(condition_file);

%% ---------------------------------------------------------------------%%
%% Format Data
% Channel 1 = Flexor
% Channel 2 = Extensor
%% ---------------------------------------------------------------------%%
fprintf("\nFormatting data for analysis...")

% ConditionInfo is stored in the condition files (.mat) acquired from DIG
events = ConditionInfo(:, 2);
RawEMG = mat2struct(events, ConditionData);

if strcmp(recording_arm_lower, 'l')
    eventRF = find(strcmp(events, 'Ch1LeRF'));
    eventRE = find(strcmp(events, 'Ch2LeRE'));
else
    eventRF = find(strcmp(events, 'Ch1RiRF'));
    eventRE = find(strcmp(events, 'Ch2RiRE'));
end

RF_IDX = eventRF; % default = 8
RE_IDX = eventRE; % default = 16
target_eventIDX = [RF_IDX RE_IDX];

fprintf("\nRaw data formatted!\n")

%% ---------------------------------------------------------------------%%
%% Pre-Processing
%% ---------------------------------------------------------------------%%
%% Bandpass Filter & Rectify
fprintf("\nApplying Bandpass Filter and Rectifying EMG signal...")
% [ProcessedEMG, events_processed] = preprocess_emg(ConditionData, events, [aFq bFq], sampleRate, 'butter');
[ProcessedEMG2, events_processed] = preprocess_emg(ConditionData, events, [aFq bFq], sampleRate, 'butter');
fprintf("\nEMG Bandpass filtered and Rectified!\n")

%% ---------------------------------------------------------------------%%
%% Notch Filter
% Account for power line interference at 60 Hz
% https://www.mathworks.com/help/signal/ug/remove-the-60-hz-hum-from-a-signal.html
%% ---------------------------------------------------------------------%%
fprintf("\nApplying 60Hz notch Filter...")

d = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1', 59, 'HalfPowerFrequency2', 61, ...
    'DesignMethod', 'butter', 'SampleRate', sampleRate);

ProcessedEMG = struct();
for e = 1:length(events_processed)
    [~, cols] = size(ProcessedEMG2.(events_processed{e}));
    for trial = 1:cols
        processed = ProcessedEMG2.(events_processed{e})(:, trial);
        notch = filtfilt(d, processed);
        ProcessedEMG.(events_processed{e})(:, trial) = notch;
    end
end
fprintf("\nEMG notched filtered!\n")

%% ---------------------------------------------------------------------%%
%% Smooth Signal
%% ---------------------------------------------------------------------%%
fprintf("\nSmoothing Processed EMG signal...")
SmoothedEMG = struct();
% EnvelopeEMG = struct();

events_smoothed = events + "_smoothed";
% events_envelope = events + "_envelope";
% Ch1EventsSmoothed = events_smoothed(1:length(events_smoothed) / 2);
% Ch2EventsSmoothed = events_smoothed((length(events_smoothed) / 2) + 1:end);

% for i = 1:length(ConditionData)
for e = 1:length(events_processed)
    [~, cols] = size(ProcessedEMG.(events_processed{e}));

    % [b, a] = butter(6, 5 / nyq_freq, 'low'); % butterworth filter
    [b, a] = butter(6, envelope_lowpassHz / (sampleRate / 2), 'low'); % butterworth filter
    for trial = 1:cols
        processed_data = ProcessedEMG.(events_processed{e})(:, trial); % EMG recordings for each condition in a given trial
        LowPassEMG = filtfilt(b, a, processed_data);
        SmoothedEMG.(events_smoothed{e})(:, trial) = LowPassEMG;
    end
end
fprintf("\nEMG Smoothed!\n\n")


%% ---------------------------------------------------------------------%%
%% Get Max EMG 
%% ---------------------------------------------------------------------%%
fprintf("Get max EMG (mV) for each trial...")
MaxEMG = struct();

for e = 1:length(events_smoothed)
    [~, cols] = size(SmoothedEMG.(events_smoothed{e}));

        for trial = 1:cols
            tempEMG = SmoothedEMG.(events_smoothed{e})(:, trial);
            maxEMG = max(tempEMG, [], 'all');
            maxIDX = find(tempEMG == maxEMG, 1);
            MaxEMG.(events{e})(1, trial) = maxIDX;
            MaxEMG.(events{e})(2, trial) = maxEMG;            
        end

end

fprintf("\nMax EMG (mV) stored in MaxEMG\n\n")

%% ---------------------------------------------------------------------%%
%% Sliding Average to detect Baseline
% Sanity check:
% Number of Samples = 2000Hz (sr) × 3s (pre-stim dur) = 6000 samples
%% ---------------------------------------------------------------------%%
BaselineEMG = struct();

events_baseline = events + "_baseline";

for e = 1:length(events_baseline)

    % Baseline [-3, 0)
    baselineIDX = time <= 0;
    BaselineEMG.(events_baseline{e}) = SmoothedEMG.(events_smoothed{e})(baselineIDX, :);
    % BaselineEMG.(events_baseline{e}) = ProcessedEMG.(events_processed{e})(baselineIDX, :);
    % BaselineEMG.(events_baseline{e}) = RawEMG.(events{e})(baselineIDX, :);

end

BaselineStats = struct();

% Define window size and step size
% sr * window (s) = # of samples
window = sampleRate * (window_ms/1000);
step = sampleRate * (step_ms/1000);

% total_trials = 0;
for e = 1:length(events_baseline)
    [~, cols] = size(BaselineEMG.(events_baseline{e}));

    % total_trials = total_trials + size(BaselineEMG.(events_baseline{e}), 2);

    for trial = 1:cols
        data = BaselineEMG.(events_baseline{e})(:, trial);

        % Initialize an array to store the averages
        num_windows = round((length(data) - window) / step);
        windowAverages = zeros(1, num_windows);
        % windowMedians = zeros(1, num_windows);
        windowSTD = zeros(1, num_windows);
        windowKeys = [];
        windowKeysStart = zeros(num_windows, 1);
        windowKeysStop = zeros(num_windows, 1);
        windowMin = zeros(1, num_windows);
        windowMax = zeros(1, num_windows);

        % windowStartValues = zeros(total_trials, 1);

        % Loop through the data using the sliding window
        for i = 1:num_windows
            start = ((i - 1) * step) + 1;
            stop = start + window;

            % Compute the mean for the current window
            windowAverages(i) = mean(data(start:stop));
            % windowMedians(i) = median(data(start:stop), 1);
            windowSTD(i) = std(data(start:stop));
            windowKeysStart(i) = time(start);
            windowKeysStop(i) = time(stop);
            windowMin(i) = min(data(start:stop), [], 'all');
            windowMax(i) = max(data(start:stop), [], 'all');

            key_format = sprintf("(%.4f, %.4f)", time(start), time(stop));
            windowKeys = [windowKeys; key_format];

        end

        % Find the minimum average
        min_average = min(windowAverages, [], 'all');
        baselineIDX = find(windowAverages == min_average);
        % min_median = min(windowMedians, [], 'all');
        % baselineIDX = find(windowMedians == min_median, 1); % need to specify take the first index where min_median is found

        BaselineStats.(events_baseline{e})(1, trial) = min_average;
        % BaselineStats.(events_baseline{e})(1, trial) = min_median;
        BaselineStats.(events_baseline{e})(2, trial) = windowSTD(baselineIDX);
        BaselineStats.(events_baseline{e})(3, trial) = windowKeysStart(baselineIDX);
        BaselineStats.(events_baseline{e})(4, trial) = windowKeysStop(baselineIDX);
        BaselineStats.(events_baseline{e})(5, trial) = windowMin(baselineIDX);
        BaselineStats.(events_baseline{e})(6, trial) = windowMax(baselineIDX);

    end
end

%% ---------------------------------------------------------------------%%
%% Get Onset
% Theshold = z-score; 2 SD
% dur = 30ms
%% ---------------------------------------------------------------------%%
EMGBurstOnsets = struct();

close all;

% for e = 1:length(events_smoothed)
for e = target_eventIDX
    [~, cols] = size(SmoothedEMG.(events_smoothed{e}));

    for trial = 1:cols
        clear tempEMG tempAvg tempSD
        tempEMG = SmoothedEMG.(events_smoothed{e})(:, trial);
        tempAvg = BaselineStats.(events_baseline{e})(1, trial);
        tempSD = BaselineStats.(events_baseline{e})(2, trial);

        clear onsetThreshold aboveThreshold
        onsetThreshold = tempAvg + (tempSD * thresholdSize);
        % tempMax = BaselineStats.(events_baseline{e})(6, trial);
        % tempMin = BaselineStats.(events_baseline{e})(5, trial);
        % onsetThreshold = percentBaselineRange * (tempMax - tempMin);
        % onsetThreshold = onsetThreshold * zSize;
        samplesToSearch = ceil(dur * sampleRate) + 1;
        samplesToSearch2 = ceil(dur2 * sampleRate) + samplesToSearch;
        
        aboveThreshold = tempEMG > onsetThreshold;
        aboveThreshold2 = tempEMG > (tempAvg + tempSD);        
        startSearchIDX = find(time == startSearchTime); % start search at -500ms

        EMGBurstOnsets.(events{e})(1, trial) = 0;
        EMGBurstOnsets.(events{e})(2, trial) = 0;
        
        firstThresholdIndex = [];
        for i = startSearchIDX:(length(tempEMG) - samplesToSearch2)
            if all(aboveThreshold(i:i + samplesToSearch))
                if all(aboveThreshold2(i + samplesToSearch + 1: i + samplesToSearch2))
                    firstThresholdIndex = i;
                    EMGBurstOnsets.(events{e})(1, trial) = firstThresholdIndex;
                    EMGBurstOnsets.(events{e})(2, trial) = time(firstThresholdIndex);
                    break;
                end

            end
        end
        %% ---------------------------------------------------------------------%%
        %% Plotting
        %% ---------------------------------------------------------------------%%
        if PlotOnsets
            OnsetPlotTitle = sprintf("%s %s Trial %d (%s) %s-Arm – EMG Burst onset", subject, events{e}, trial, timepoint, recording_arm);
            OnsetFigureTitle = sprintf("%s %s Trial %d – EMG Burst onset", subject, events{e}, trial);
            eventEMGBurstOnset = figure('Name', OnsetFigureTitle);

            clear ax1 ax2 y_lim_ax1 y_lim_ax2
            % First subplot for Smoothed EMG
            ax1 = subplot(9, 5, 1:20, 'Parent', eventEMGBurstOnset);
            hold on;
            plot(time, SmoothedEMG.(events_smoothed{e})(:, trial), 'Color', 1/255*[0,104,87], 'LineWidth', 1, 'DisplayName', 'Smoothed EMG', 'HandleVisibility', 'off'); % Plot test data
            % Plot baseline rectangle for onset threshold optimization
            baselineStart = BaselineStats.(events_baseline{e})(3, trial);
            maxIDX = MaxEMG.(events{e})(1, trial);
            maxEMG = MaxEMG.(events{e})(2, trial);
            y_lim_ax1 = get(ax1, 'YLim');
            ylim(ax1, [0, y_lim_ax1(2)])
            % Define the x and y coordinates for the filled region
            x_fill_ax1 = [baselineStart, baselineStart + (window_ms/1000), baselineStart + (window_ms/1000), baselineStart];
            y_fill_ax1 = [0, 0, y_lim_ax1(2), y_lim_ax1(2)];

            % Use fill to create the shaded area
            fill(x_fill_ax1, y_fill_ax1, [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'Parent', ax1, 'DisplayName', 'Baseline Window');
            plot(time, SmoothedEMG.(events_smoothed{e})(:, trial), 'Color', 1/255*[0,104,87], 'LineWidth', 1, 'HandleVisibility', 'off'); % Plot test data

            % Plot horizontal line for onsetThreshold
            yline(onsetThreshold, 'LineStyle', '-.', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1, 'DisplayName', 'Onset Threshold');
            % yline(tempAvg, 'LineStyle', '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1, 'DisplayName', 'Baseline Average');
            % yline(onsetThreshold, 'o--', 'LineWidth', 1, 'DisplayName', 'Low Threshold');
            % yline(tempAvg, '--', 'Color', 1/255*[1 1 1], 'LineWidth', 1, 'DisplayName', 'Baseline Average');

            maxText = sprintf('Max = %.3f', maxEMG);
            xline(time(maxIDX), 'LineStyle', '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1, 'DisplayName', maxText)

            % Plot vertical line at time(firstIndex), if it exists
            if ~isempty(firstThresholdIndex)
                xline(time(firstThresholdIndex), 'k--', 'LineWidth', 1, 'DisplayName', 'Onset Time');
                fprintf('\t%s | Trial %02d:\t%.1f ms\n', events{e}, trial, (time(firstThresholdIndex)*1000));
            else
                fprintf('\t%s | Trial %02d:\t0.0 ms – NO ONSET DETECTED\n', events{e}, trial);
            end

            % Add labels and legend
            % xlabel('Time (s)');
            % ylabel('EMG Signal (mV)');
            lgd = legend('show', 'Location', 'northoutside', 'Orientation', 'horizontal');
            lgd.FontSize = 10;
            title(ax1, "Smoothed EMG", 'FontSize', 12, 'FontWeight', 'bold');

            % Second subplot for Processed EMG
            ax2 = subplot(9, 5, 26:45, 'Parent', eventEMGBurstOnset);
            hold on;
            plot(time, ProcessedEMG.(events_processed{e})(:, trial), 'DisplayName', 'Processed Data', 'HandleVisibility', 'off')
            % Plot baseline rectangle for onset threshold optimization
            % Use fill to create the shaded area
            y_lim_ax2 = get(ax2, 'YLim');
            ylim(ax2, [0, y_lim_ax2(2)])
            % Define the x and y coordinates for the filled region
            x_fill_ax2 = [baselineStart, baselineStart + (window_ms/1000), baselineStart + (window_ms/1000), baselineStart];
            y_fill_ax2 = [0, 0, y_lim_ax2(2), y_lim_ax2(2)];
            fill(x_fill_ax2, y_fill_ax2, [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'Parent', ax2, 'DisplayName', 'Baseline Window');
            plot(time, ProcessedEMG.(events_processed{e})(:, trial), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1, 'HandleVisibility', 'off')

            % Add labels and legend
            % xlabel('Time (s)');
            % ylabel('EMG Signal (mV)');
            % title('Processed EMG');
            % lgd = legend('show');
            % lgd.FontSize = 8; % Smaller font size

            hold on;
            yline(onsetThreshold, 'LineStyle', '-.', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1, 'DisplayName', 'Onset Threshold');
            % yline(tempAvg, 'LineStyle', '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1, 'DisplayName', 'Baseline Average');
            % yline(onsetThreshold, 'r--', 'LineWidth', 1, 'DisplayName', 'Onset Threshold');
            % yline(tempAvg, '--', 'Color', 1/255*[1 1 1], 'LineWidth', 1, 'DisplayName', 'Baseline Average');

            xline(time(maxIDX), 'LineStyle', '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1, 'DisplayName', maxText)

            % Plot vertical line at time(firstIndex), if it exists
            if ~isempty(firstThresholdIndex)
                xline(time(firstThresholdIndex), 'k--', 'LineWidth', 1, 'DisplayName', 'Onset Time');
            end
            hold off;
            title(ax2, "Processed EMG", 'FontSize', 12, 'FontWeight', 'bold');

            % Title of whole figure
            sgtitle(OnsetPlotTitle, 'FontSize', 14, 'FontWeight', 'bold')

            % Shared labels
            l = axes(eventEMGBurstOnset, 'visible', 'off');
            l.XLabel.Visible = 'on';
            xlabel(l, 'Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
            l.YLabel.Visible = 'on';
            ylabel(l, 'EMG (mV)', 'FontSize', 14, 'FontWeight', 'bold',...
                'Position', [-0.085, 0.5, 0]);

            % Save the figure
            filename_png = sprintf("%s_%s_Trial-%02d_%s_%s-Arm_EMGBurstOnset.png", events{e}, subject, trial, timepoint, recording_arm);
            dir_name = events{e}(end - 1:end);
            FigDir = sprintf('%s/%s-%s', figures_dir, recording_arm, dir_name);
            if ~exist(FigDir, 'dir')
                mkdir(FigDir);
            end
            saveas(eventEMGBurstOnset, fullfile(FigDir, filename_png));
            % filename_fig = sprintf("%s_%s_Trial-%d_%s_%s-Arm_RawEMG.fig", eventName, subject, trial, timepoint, affected_arm);
            % saveas(eventRawEMG, fullfile(EE_fig_dir, filename_fig));
            close(eventEMGBurstOnset); % comment out for sanity check
        end

    end
    if PlotOnsets
        fprintf("\n")
    end
end

% %% ---------------------------------------------------------------------%%
% %% Plotting Histogram of Baseline Windows
% %% ---------------------------------------------------------------------%%
% % create a dictionary with each key as a window for the baseline search
% % the value of the key will just be frequency
% 
% windowStartValues = [];
% for i = 1:length(events_baseline)
% 
%     % clear tempArray tempData startTimes
%     startTimesArray = [];
%     tempData = BaselineStats.(events_baseline{i});  % Access the field
%     startTimes = tempData(3, :);  % Extract the 3rd row
%     startTimesArray = [startTimesArray; startTimes];
%     windowStartValues = [windowStartValues; startTimesArray'];
% 
% end
% 
% windowStartValues = sort(windowStartValues);
% 
% % Find unique values and their counts
% [windowStarts, ~, idx] = unique(windowStartValues);
% counts = histcounts(idx, numel(windowStarts));
% 
% % Combine unique values and counts into a nx2 array
% window_counts = [round(windowStarts, 4), counts'];
% 
% samples = sampleRate * 3;
% bins = round((samples - window) / step);
% 
% close all
% figure();
% windowDistribution = histogram(windowStartValues, bins);
% ylabel("Count")
% xlabel("Baseline Windows")
% title(["Histogram of Baseline Windows EMG Burst Onset Detection" recording_arm "-Arm"])
% 
% % Save the figure
% filename_png = sprintf("%s_%s_%s-Arm_EMGBurstOnset-Histogram.png", timepoint, subject, recording_arm);
% saveas(windowDistribution, fullfile(figures_dir, filename_png));
% 
% figure();
% windowDistribution2 = histogram(Ch1RF_Onsets(:, 2), bins);
% ylabel("Count")
% xlabel("Baseline Windows")
% title(["Histogram of Baseline Windows EMG Burst Onset Detection" recording_arm "-Arm for Ch1RF"])
% 
% % Save the figure
% filename_png2 = sprintf("%s_%s_%s-Arm_EMGBurstOnset-Ch1RF-Histogram.png", timepoint, subject, recording_arm);
% saveas(windowDistribution2, fullfile(figures_dir, filename_png2));
% 
% figure();
% windowDistribution3 = histogram(Ch2RE_Onsets(:, 2), bins);
% ylabel("Count")
% xlabel("Baseline Windows")
% title(["Histogram of Baseline Windows EMG Burst Onset Detection" recording_arm "-Arm for Ch2RE"])
% 
% % Save the figure
% filename_png3 = sprintf("%s_%s_%s-Arm_EMGBurstOnset-Ch2RE-Histogram.png", timepoint, subject, recording_arm);
% saveas(windowDistribution3, fullfile(figures_dir, filename_png3));

%% ---------------------------------------------------------------------%%
%% Calculate AUC
%% ---------------------------------------------------------------------%%
fprintf("Calculating AUC...")
AUC = struct(); % row 1 = 500ms after onset, row 2 = onset to end (3000ms)
window500ms = ceil(0.500 * sampleRate) + 1;
window3000ms = ceil(3.000 * sampleRate) + 1;

for e = target_eventIDX
    [~, cols] = size(EMGBurstOnsets.(events{e}));

    for trial = 1:cols
        onsetIDX = EMGBurstOnsets.(events{e})(1, trial);
        % onsetTime = EMGBurstOnsets.(events{e})(2, trial);
        tempEMG = SmoothedEMG.(events_smoothed{e})(:, trial);

        tempEMGMatrix = [time' tempEMG];

        stop500 = onsetIDX + window500ms;
        % stop3000 = length(time);
        stop3000 = onsetIDX + window3000ms;

        AUC.(events{e})(1, trial) = 0;
        AUC.(events{e})(2, trial) = 0;

        if onsetIDX == 0
            AUC.(events{e})(1, trial) = 0;
            AUC.(events{e})(2, trial) = 0;
        elseif stop500 <= length(time)
            auc500 = trapz(tempEMGMatrix(onsetIDX:stop500, 1), tempEMGMatrix(onsetIDX:stop500, 2));
            cp500 = ChoiceProbability(tempEMGMatrix(onsetIDX:stop500, 1), tempEMGMatrix(onsetIDX:stop500, 2));
            AUC.(events{e})(1, trial) = auc500 * 1000; % multiple by 1000 to get in ms
            AUC.(events{e})(3, trial) = cp500;
        end

        if onsetIDX == 0
            AUC.(events{e})(1, trial) = 0;
            AUC.(events{e})(2, trial) = 0;
        elseif stop3000 <= length(time)
            auc3000 = trapz(tempEMGMatrix(onsetIDX:stop3000, 1), tempEMGMatrix(onsetIDX:stop3000, 2));
            cp3000 = ChoiceProbability(tempEMGMatrix(onsetIDX:stop3000, 1), tempEMGMatrix(onsetIDX:stop3000, 2));            
            AUC.(events{e})(2, trial) = auc3000 * 1000; % multiple by 1000 to get in ms
            AUC.(events{e})(4, trial) = cp3000;

        end
    end
end

fprintf("\nAUC calculated\n\n")

%% ---------------------------------------------------------------------%%
%% Time from onset to peak
%% ---------------------------------------------------------------------%%
Time2Peak = struct();

for e = target_eventIDX
    [~, cols] = size(EMGBurstOnsets.(events{e}));

    for trial = 1:cols
        onset = EMGBurstOnsets.(events{e})(2, trial) * 1000;
        peak = time(MaxEMG.(events{e})(1, trial)) * 1000;
        Time2Peak.(events{e})(1, trial) = peak;
        Time2Peak.(events{e})(2, trial) = onset;
        Time2Peak.(events{e})(3, trial) = peak - onset;
    end

end

%% ---------------------------------------------------------------------%%
%% Save Everything
%% ---------------------------------------------------------------------%%
SummaryDataCh1RF(:, 1) = EMGBurstOnsets.(events{RF_IDX})(2, :)' * 1000; % Onset Times
SummaryDataCh1RF(:, 2) = Time2Peak.(events{RF_IDX})(1, :)'; % Peak Time
SummaryDataCh1RF(:, 3) = Time2Peak.(events{RF_IDX})(3, :)'; % Time to Peak (peak - onset)
SummaryDataCh1RF(:, 4) = AUC.(events{RF_IDX})(1, :)'; % AUC 500ms after onset
SummaryDataCh1RF(:, 5) = AUC.(events{RF_IDX})(2, :)'; % AUC from onset to end of stimulus block
SummaryDataCh1RF(:, 6) = BaselineStats.(events_baseline{RF_IDX})(1, :); % baseline average
SummaryDataCh1RF(:, 7) = MaxEMG.(events{RF_IDX})(2, :)'; % max sEMG signal
SummaryDataCh1RF(:, 8) = SummaryDataCh1RF(:, 7) - SummaryDataCh1RF(:, 6);% normalized EMG (Max-Baseline Average)


SummaryDataCh2RE(:, 1) = EMGBurstOnsets.(events{RE_IDX})(2, :)' * 1000; % Onset Times
SummaryDataCh2RE(:, 2) = Time2Peak.(events{RE_IDX})(1, :)'; % Peak Time
SummaryDataCh2RE(:, 3) = Time2Peak.(events{RE_IDX})(3, :)'; % Time to Peak (peak - onset)
SummaryDataCh2RE(:, 4) = AUC.(events{RE_IDX})(1, :)'; % AUC 500ms after onset
SummaryDataCh2RE(:, 5) = AUC.(events{RE_IDX})(2, :)'; % AUC from onset to end of stimulus block
SummaryDataCh2RE(:, 6) = BaselineStats.(events_baseline{RE_IDX})(1, :); % baseline average
SummaryDataCh2RE(:, 7) = MaxEMG.(events{RE_IDX})(2, :)'; % max sEMG signal
SummaryDataCh2RE(:, 8) = SummaryDataCh2RE(:, 7) - SummaryDataCh2RE(:, 6);% normalized EMG (Max-Baseline Average)

fprintf("\nSaving everything...")
data_filename = sprintf('%s_%s_%s-Arm_EMG-Burst-Onset', timepoint, subject, recording_arm);

% if strcmp(recording_arm_lower, affected_arm_lower)
%     output_file = sprintf('%s/%s.mat', affected_dir, data_filename);
% elseif strcmp(recording_arm_lower, unaffected_arm_lower)
%     output_file = sprintf('%s/%s.mat', unaffected_dir, data_filename);
% end

output_file = sprintf('%s/%s.mat', analysis_dir, data_filename);


%% ---------------------------------------------------------------------%%
%% Clean up things you don't need
%% ---------------------------------------------------------------------%%
vars = {
    'a', 'affected_arm_lower', 'affected_dir', 'analysis_dir', 'ax1', 'ax2',...
    'b', 'baselineIDX', 'baselineStart', 'cols', 'condition_file',...
    'condition_name_baseline', 'condition_name_processed', 'condition_name_smoothed',...
    'conditions_dir', 'counts', 'data', 'data_filename', 'dir', 'dir_list',...
    'e', 'eventEMGBurstOnset', 'events_baseline', 'events_processed', 'events_smoothed',...
    'EventsToConditions', 'FigDir', 'figures_dir', 'filename_png', 'filename_png2',...
    'filename_png3', 'firstOnset', 'i', 'idx', 'l', 'lgd', 'LowPassEMG', 'min_average',...
    'OnsetFigureTitle', 'OnsetPlotTitle', 'processed_data', 'recording_analysis_dir',...
    'recording_arm_lower', 'root', 'samplesToSearch', 'start', 'startSearchIDX', 'step',...
    'stop', 'target_eventIDX', 'tempAvg', 'tempData', 'tempEMG', 'tempSD', 'today',...
    'trial', 'unaffected_arm_lower', 'window', 'windowMax', 'windowMin', 'windowStarts',...
    'windowStartValues', 'windowSTD', 'x_fill_ax1', 'x_fill_ax2', 'y_fill_ax1',...
    'y_fill_ax2', 'y_lim_ax1', 'y_lim_ax2', 'bins', 'samples'
    };
clear(vars{:});
clear vars

save(output_file);

fprintf("\nScript complete!\n\n")
end



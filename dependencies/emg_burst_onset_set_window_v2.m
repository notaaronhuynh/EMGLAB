%% EVA EMG Burst Onset Automation – Set Window
% University of Rochester School of Medicine & Dentistry
% Movement and Plasticity Lab (MAPL)
% Author: Aaron Huynh | PI: Dr. Ania Busza
% Collaborators:
%   Dr. Ania Busza
%   Dr. David Cunningham (Case Western University)
% Last commit: 22 November 2024 ANH

function emg_burst_onset_set_window_v2(initials, subject, affected_arm, timepoint, recording_arm, PlotOnsets)
% Manually set baseline window period

%% CHANGE THESE
windowStart = -2350; % ms
windowStop = -1350; % ms
windowSize = windowStop - windowStart; % leave this

dur = 0.03; % s; 30ms
dur2 = 0.070;
zSize = 2;
startSearchTime = -1.250; %s

% filter information
aFq = 10;
bFq = 500;
envelope_lowpassHz = 20;

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
analysis_dir = sprintf('%s/%s_%s_EMG-Burst-Onset_%s_%s (%d – %d)', root, timepoint, subject, formattedDate, initials, windowStart, windowStop);

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
% elseif strcmp(recording_arm_lower, unaffected_arm_lower)
%     recording_dir = unaffected_dir;
%     recording_analysis_dir = unaffected_dir;
%     figures_dir = sprintf('%s/Figures', unaffected_dir);
% end
figures_dir = sprintf('%s/Figures', analysis_dir);

% List of directories to create
% dir_list = {figures_dir, recording_dir, recording_analysis_dir, analysis_dir};
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
[Processed, events_processed] = preprocess_emg(ConditionData, events, [aFq bFq], sampleRate, 'butter');
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

% This doesn't change the signal that much but do it anyways for sanity
for e = 1:length(events_processed)
    [~, cols] = size(Processed.(events_processed{e}));
        for trial = 1:cols
            processed = Processed.(events_processed{e})(:, trial);
            notch = filtfilt(d, processed);
            ProcessedEMG.(events_processed{e})(:, trial) = notch;
        end
end

fprintf("\nEMG notched filtered!\n")

%% ---------------------------------------------------------------------%%
%% Smooth Signal
%% ---------------------------------------------------------------------%%
fprintf("\nSmoothing Processed EMG signal...")

% SmoothedEMG = struct();
% % EnvelopeEMG = struct();
% 
% events_smoothed = events + "_smoothed";
% 
% % for i = 1:length(ConditionData)
% for e = 1:length(events_processed)
%     [~, cols] = size(ProcessedEMG.(events_processed{e}));
% 
%     % [b, a] = butter(6, 5 / nyq_freq, 'low'); % butterworth filter
%     [b, a] = butter(6, envelope_lowpassHz / (sampleRate / 2), 'low'); % butterworth filter
%         for trial = 1:cols
%             processed_data = ProcessedEMG.(events_processed{e})(:, trial); % EMG recordings for each condition in a given trial
%             LowPassEMG = filtfilt(b, a, processed_data);     
%             SmoothedEMG.(events_smoothed{e})(:, trial) = LowPassEMG;
%         end
% end

events_smoothed = events + "_smoothed";
padLength = (0.100 * sampleRate) + 1;
pad = zeros(1, padLength);

for e = 1:length(events_processed)
    [~, cols] = size(ProcessedEMG.(events_processed{e}));

    % [b, a] = butter(6, 5 / nyq_freq, 'low'); % butterworth filter
    [b, a] = butter(6, envelope_lowpassHz / (sampleRate / 2), 'low'); % butterworth filter
    for trial = 1:cols
        processed_data = ProcessedEMG.(events_processed{e})(:, trial); % EMG recordings for each condition in a given trial
        processed_data = [pad processed_data' pad];
        
        processed_data = processed_data';

        filteredEMG = filtfilt(b, a, processed_data);
        LowPassEMG = filteredEMG((padLength + 1):(length(filteredEMG) - padLength));
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
%% Get Baseline
% Sanity check:
% Number of Samples = 2000Hz (sr) × 3s (pre-stim dur) = 6000 samples
%% ---------------------------------------------------------------------%%
BaselineEMG = struct();

events_baseline = events + "_baseline";

for e = 1:length(events_baseline)

    % Baseline [-3, 0)
    baselineIDX = time <= 0;
    BaselineEMG.(events_baseline{e}) = SmoothedEMG.(events_smoothed{e})(baselineIDX, :);

end

BaselineStats = struct();

windowStartIDX = find(time == (windowStart/1000));
windowStopIDX = find(time == (windowStop/1000));

for e = 1:length(events_baseline)
    [~, cols] = size(BaselineEMG.(events_baseline{e}));

    % total_trials = total_trials + size(BaselineEMG.(events_baseline{e}), 2);

    for trial = 1:cols
        data = BaselineEMG.(events_baseline{e})(:, trial);

        % Find the Baseline EMG average
        baselineAVG =  mean(data(windowStartIDX:windowStopIDX));
        baselineSTD =  std(data(windowStartIDX:windowStopIDX));

        BaselineStats.(events_baseline{e})(1, trial) = baselineAVG;
        BaselineStats.(events_baseline{e})(2, trial) = baselineSTD;

    end
end

%% ---------------------------------------------------------------------%%
%% Get Onset
% Theshold = z-score; 2 SD
% dur = 30ms
%% ---------------------------------------------------------------------%%
EMGBurstOnsets = struct();

close all;

for e = target_eventIDX
    [~, cols] = size(SmoothedEMG.(events_smoothed{e}));

    for trial = 1:cols
        clear tempEMG tempAvg tempSD
        tempEMG = SmoothedEMG.(events_smoothed{e})(:, trial);
        tempAvg = BaselineStats.(events_baseline{e})(1, trial);
        tempSD = BaselineStats.(events_baseline{e})(2, trial);

        clear onsetThreshold aboveThreshold
        onsetThreshold = tempAvg + (tempSD*zSize);
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
            baselineStart = windowStart/1000;
            maxIDX = MaxEMG.(events{e})(1, trial);
            maxEMG = MaxEMG.(events{e})(2, trial);
            y_lim_ax1 = get(ax1, 'YLim');
            ylim(ax1, [0, y_lim_ax1(2)])
            % Define the x and y coordinates for the filled region
            x_fill_ax1 = [baselineStart, baselineStart + (windowSize/1000), baselineStart + (windowSize/1000), baselineStart];
            y_fill_ax1 = [0, 0, y_lim_ax1(2), y_lim_ax1(2)];

            % Use fill to create the shaded area
            fill(x_fill_ax1, y_fill_ax1, [0.85 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'Parent', ax1, 'DisplayName', 'Baseline Window');
            plot(time, SmoothedEMG.(events_smoothed{e})(:, trial), 'Color', 1/255*[0,104,87], 'LineWidth', 1, 'HandleVisibility', 'off'); % Plot test data

            % Plot horizontal line for onsetThreshold
            yline(onsetThreshold, 'LineStyle', '-.', 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1, 'DisplayName', 'Onset Threshold');
            % yline(tempAvg, 'LineStyle', '-.', 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 1, 'DisplayName', 'Baseline Average');
            % yline(onsetThreshold, 'o--', 'LineWidth', 1, 'DisplayName', 'Low Threshold');
            % yline(tempAvg, '--', 'Color', 1/255*[1 1 1], 'LineWidth', 1, 'DisplayName', 'Baseline Average');

            maxText = sprintf('Max = %.4f', maxEMG);
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
            x_fill_ax2 = [baselineStart, baselineStart + (windowSize/1000), baselineStart + (windowSize/1000), baselineStart];
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
            % filename_fig = sprintf("%s_MAPS%s_Trial-%d_%s_%s-Arm_RawEMG.fig", eventName, subject, trial, timepoint, affected_arm);
            % saveas(eventRawEMG, fullfile(EE_fig_dir, filename_fig));
            close(eventEMGBurstOnset); % comment out for sanity check
        end

    end
    if PlotOnsets
        fprintf("\n")
    end
end

%% ---------------------------------------------------------------------%%
%% Calculate AUC
%% ---------------------------------------------------------------------%%
fprintf("Calculating AUC...")
AUC = struct(); % row 1 = 500ms after onset, row 2 = onset to end (3000ms)
window500ms = ceil(0.500 * sampleRate) + 1;

for e = target_eventIDX
    [~, cols] = size(EMGBurstOnsets.(events{e}));

    for trial = 1:cols
        onsetIDX = EMGBurstOnsets.(events{e})(1, trial);
        % onsetTime = EMGBurstOnsets.(events{e})(2, trial);
        tempEMG = SmoothedEMG.(events_smoothed{e})(:, trial);

        tempEMGMatrix = [time' tempEMG];

        stop500 = onsetIDX + window500ms;
        stop3000 = length(time);

        AUC.(events{e})(1, trial) = 0;

        if onsetIDX == 0
            AUC.(events{e})(1, trial) = 0;
            AUC.(events{e})(2, trial) = 0;
        elseif (onsetIDX + window500ms) <= length(time)
            auc500 = trapz(tempEMGMatrix(onsetIDX:stop500, 1), tempEMGMatrix(onsetIDX:stop500, 2));
            AUC.(events{e})(1, trial) = auc500;

            auc3000 = trapz(tempEMGMatrix(onsetIDX:stop3000, 1), tempEMGMatrix(onsetIDX:stop3000, 2));
            AUC.(events{e})(2, trial) = auc3000;

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

if strcmp(recording_arm_lower, affected_arm_lower)
    output_file = sprintf('%s/%s.mat', affected_dir, data_filename);
elseif strcmp(recording_arm_lower, unaffected_arm_lower)
    output_file = sprintf('%s/%s.mat', unaffected_dir, data_filename);
end

%% ---------------------------------------------------------------------%%
%% Clean up things you don't need
%% ---------------------------------------------------------------------%%
vars = {
    'a', 'aboveThreshold', 'affected_arm_lower', 'affected_dir', 'analysis_dir',...
    'auc3000', 'auc500', 'ax1', 'ax2', 'b', 'baselineAVG', 'baselineStart',...
    'baselineSTD', 'cols', 'condition_file', 'conditions_dir', 'd', 'data',...
    'data_filename', 'dir', 'dir_list', 'dir_name', 'e', 'eventEMGBurstOnset',...
    'FigDir', 'figures_dir', 'filename_png', 'firstThresholdIndex', 'formattedDate',...
    'i', 'l', 'lgd', 'LowPassEMG', 'maxEMG', 'maxIDX', 'maxText', 'onset',...
    'OnsetFigureTitle', 'onsetIDX', 'OnsetPlotTitle', 'onsetThreshold', 'peak',...
    'processed', 'processed_data', 'recording_analysis_dir', 'recording_arm_lower',...
    'recording_dir', 'root', 'samplesToSearch', 'startSearchIDX', 'stop3000',...
    'stop500', 'tempAvg', 'tempEMG', 'tempEMGMatrix', 'tempSD', 'today', 'trial',...
    'unaffected_arm_lower', 'unaffected_dir', 'x_fill_ax1', 'x_fill_ax2',...
    'y_fill_ax1', 'y_fill_ax2', 'y_lim_ax1', 'y_lim_ax2'
    };
clear(vars{:});
clear vars

save(output_file);

fprintf("\nScript complete!\n\n")
end



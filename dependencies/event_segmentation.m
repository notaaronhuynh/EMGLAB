clear all
close all

%% Load file
[file, path] = uigetfile('*.xdf'); % File Selection
file = [path,file];
streams = load_xdf(file);

cd(path)

parts = strsplit(path, '/');
file_id = parts{end-2};
file_parts = strsplit(file_id, '_');
id = strjoin(file_parts, ' ');

%% Find Streams
for s = 1:length(streams)

    type = streams{1, s}.info.type;

    if isequal(type, 'Events')
        eventsIDX = s;
    elseif isequal(type, 'Level')
        levelsIDX = s;
    elseif isequal(type, 'Runner') % target or actual position? not where block is
        runnerIDX = s;
        % elseif isequal(type, 'Settings')
        %     settingsIDX = s;
        % elseif isequal(type, 'Data')
        %     dataIDX = s;
    elseif isequal(type, 'EMG')
        emgIDX = s;
        % elseif isequal(type, 'Smooth EMG')
        %     smoothIDX = s;
    end

end

%% Store streams data
events = streams{1, eventsIDX}.time_series; % event labels
events_time = streams{1, eventsIDX}.time_stamps; % event times
events_sr = str2double(streams{1, eventsIDX}.info.nominal_srate);

startIDX = find(strcmp(events, 'START'));
stopIDX = find(strcmp(events, 'STOP'));

start_time = events_time(startIDX);
stop_time = events_time(stopIDX);
game_start_idx = find(strcmp(events, 'game start'));

calibrationStop = startIDX - 1;
tutorialStop = game_start_idx - 1;

calibration = events_time(1:calibrationStop);
tutorial = events_time(startIDX + 1:tutorialStop);

blockStartIDX = find(strcmp(events,'expected zone start'));
% blockStartIDX = blockStartIDX(4:end); % ignore tutorial events
blockStartTimes = events_time(blockStartIDX);
blockStopIDX = find(strcmp(events,'expected zone stop'));
% blockStopIDX = blockStopIDX(4:end); % ignore tutorial events
blockStopTimes = events_time(blockStopIDX);

levels = streams{1, levelsIDX}.time_series(1,:); % level labels
levels_time = streams{1, levelsIDX}.time_stamps; % levels times

runner_position = streams{1, runnerIDX}.time_series(2,:); % block position
runner_time = streams{1, runnerIDX}.time_stamps;
runner_sr = str2double(streams{1, runnerIDX}.info.nominal_srate);
runner_dt = diff(runner_time);
runner_dt = [0 runner_dt];

emg = streams{1, emgIDX}.time_series;
emg_ch1 = emg(2,:); % channel 1
emg_ch2 = emg(4,:); % channel 2
emg_time = streams{1, emgIDX}.time_stamps;
emg_sr = str2double(streams{1, emgIDX}.info.nominal_srate);
emg_dt = diff(emg_time);
emg_dt = [0 emg_dt];

%% Manually adjust event time stamps
events_dt = diff(events_time);
events_dt = [0 events_dt];

events_time_adjusted = [];

for i = 1:length(events_dt)
    temp = round(events_dt(i), 2);

    if temp == 3
        events_time_adjusted(i) = temp;
    elseif temp == 0.5
        events_time_adjusted(i) = temp;
    else
        events_time_adjusted(i) = events_dt(i);
    end

end

events_fixed = [events_time(1)];
i = 2;
while i <= length(events_time)

    adj = events_fixed(i - 1) + events_time_adjusted(i);
    events_fixed(i) = adj;
    i = i + 1;
end

events_fixed(end) = emg_time(end);

level_first_block_diff = nan(1, numel(levels_time)); 
level_first_block_idx = nan(1, numel(levels_time));  

for i = 1:numel(levels_time)
    % Find first event after the level_time
    idx = find(events_time > levels_time(i), 1, 'first');
    level_first_block_idx(i) = idx;
    if ~isempty(idx)
        level_first_block_diff(i) = events_time(idx) - levels_time(i);
    end
end

% events_fixed(level_diffs) = events_time(level_diffs);


%% Extract tutorial and calibration
start_time = events_fixed(startIDX);
stop_time = events_fixed(stopIDX);
game_start_idx = find(strcmp(events, 'game start'));

calibrationStop = startIDX - 1;
tutorialStop = game_start_idx - 1;

calibration = events_fixed(1:calibrationStop);
calibration_blockStartIDX = find(strcmp(calibration,'expected zone start'));
calibration_blockStartTimes = events_fixed(calibration_blockStartIDX);
calibration_blockStopIDX = find(strcmp(calibration,'expected zone stop'));
calibration_blockStopTimes = events_fixed(calibration_blockStopIDX);

tutorial = events_fixed(startIDX + 1:tutorialStop);
tutorial_blockStartIDX = find(strcmp(tutorial,'expected zone start'));
tutorial_blockStartTimes = events_fixed(tutorial_blockStartIDX);
tutorial_blockStopIDX = find(strcmp(tutorial,'expected zone stop'));
tutorial_blockStopTimes = events_fixed(tutorial_blockStopIDX);

% Start with index 2, keep two, skip one: pattern = [2 3 skip 5 6 skip ...]
clear idx
idx = [];

i = 2;
while i <= length(tutorial)
    idx = [idx, i];
    if i+1 <= length(tutorial)
        idx = [idx, i+1];
    end
    i = i + 3;
end

tutorial_clean = tutorial(idx);

clear tutorialExpStop tutorialExpStart
tutorialExpStart = [];
tutorialExpStop = [];

for i = 1:length(tutorial_clean)
    if mod(i, 2) == 0
        tutorialExpStop = [tutorialExpStop, i];
    else
        tutorialExpStart = [tutorialExpStart, i];
    end
end

tutorialStartTimes = events_fixed(tutorialExpStart);

%% Adjust runner position values for visualization
if max(emg_ch1) >= max(emg_ch2)
    max_emg = max(emg_ch1);
else
    max_emg = max(emg_ch2);
end

if mean(emg_ch1) >= mean(emg_ch2)
    med = median(emg_ch1);
else
    med = median(emg_ch2);
end

runner_position_adjusted = zeros(size(runner_position));
for i = 1:length(runner_position)
    switch round(runner_position(i))
        case 2
            runner_position_adjusted(i) = med + (max_emg * 0.4);
        case 1
            runner_position_adjusted(i) = med;
        case 0
            runner_position_adjusted(i) = med - (max_emg * 0.4);
    end
end

%% Manually adjust runner time stamps
blockStartTimes_fixed = events_fixed(blockStartIDX);
blockStopTimes_fixed = events_fixed(blockStopIDX);

% Interpolate runner_position onto events_time:
start_time = events_fixed(startIDX);
offset = start_time - runner_time(1);
runner_time_offset = runner_time + offset;

target_times = [blockStartTimes_fixed blockStopTimes_fixed]; % add the tutorial times
% target_times = events_fixed;
% target_times(repeated) = [];
[target_times_sorted, ~] = sort(target_times);

% % Create interpolant
% F = griddedInterpolant(runner_time, runner_position, 'nearest');
% % Get runner position at event times
% runner_at_events = F(events_fixed);
% F_event = griddedInterpolant(target_times_sorted, runner_at_events, 'nearest');
% runner_position_stretched = F_event(runner_time_offset);

states = interp1(runner_time, runner_position, target_times_sorted, 'previous');
states(isnan(states)) = 1;
% states = states(2:end);
% states = [states 1];
% states(1) = 1;

% change first state following the start of the level to be 'rest' (1)
for i = 1:length(levels_time)
    % Find first index in blockStopTimes where time > levels_time(i)
    idx = find(target_times_sorted > levels_time(i), 1, 'first');

    if ~isempty(idx) && idx > 1
        % Set previous index to 1
        states(idx-1) = 1;
    elseif idx == 1
        % If first target_times is already > levels_time(i), do nothing
        % (cannot go to previous index)
        continue
    else
        % If levels_time is after all blockStopTimes, set the last one
        states(end) = 1;
    end
end

stretched_states = interp1(target_times_sorted, states, runner_time, 'previous', 'extrap');
% stretched_states = interp1(target_times_sorted, states, runner_time, 'previous');
% % adding NaNs for first state for some reason when stretched
stretched_states(isnan(stretched_states)) = 1;

%% ChatGPT
% aligned_runner_time = runner_time;                      % Copy original
% aligned_runner_position = nan(size(runner_position));  % Preallocate the same length
%
% % Get the stop times
% expected_stop_times = events_fixed(strcmp(events, 'expected zone stop'));
%
% % Initial segment (before first stop) — just copy the original positions
% initial_mask = runner_time < expected_stop_times(1);
% aligned_runner_position(initial_mask) = runner_position(initial_mask);
%
% % Segments between stops
% for seg = 1:length(expected_stop_times)-1
%     t_start = expected_stop_times(seg);
%     t_end = expected_stop_times(seg + 1);
%     runner_mask = runner_time > t_start & runner_time < t_end;
%
%     runner_segment = runner_position(runner_mask);
%     target_value = mode(runner_segment);
%
%     aligned_runner_position(runner_mask) = target_value;
% end
%
% % Final segment after the last stop
% last_mask = runner_time >= expected_stop_times(end);
% aligned_runner_position(last_mask) = runner_position(last_mask);

expected_stop_indices = find(strcmp(events, 'expected zone start'));
expected_stop_times = events_fixed(expected_stop_indices(:));
expected_stop_times = unique(expected_stop_times(:));
expected_stop_times = expected_stop_times(expected_stop_times <= runner_time(end));

%% Preallocate Aligned Output
aligned_runner_time = runner_time;                 % Same length
aligned_runner_position = nan(size(runner_position));

%% 1️⃣ Handle the START Segment (before first stop)
if ~isempty(expected_stop_times)
    first_boundary = expected_stop_times(1);
    start_mask = runner_time < first_boundary;

    if any(start_mask)
        target_value = mode(runner_position(start_mask));
        aligned_runner_position(start_mask) = target_value;
    end
end

%% 2️⃣ Handle the MAIN Segments
for seg = 1:length(expected_stop_times) - 1
    t_start = expected_stop_times(seg);
    t_end = expected_stop_times(seg + 1);

    runner_mask = runner_time >= t_start & runner_time < t_end;

    if any(runner_mask)
        target_value = mode(runner_position(runner_mask));
        aligned_runner_position(runner_mask) = target_value;
    end
end

%% 3️⃣ Handle the END Segment (after the final stop)
last_boundary = expected_stop_times(end);
tail_mask = runner_time >= last_boundary;

if any(tail_mask)
    target_value = mode(runner_position(tail_mask));
    aligned_runner_position(tail_mask) = target_value;
end

% adjusted_mask = runner_time <= levels_time(1);
% aligned_runner_position(adjusted_mask) = stretched_states(adjusted_mask);

% Loop from second level
for i = 1:length(levels_time)
    % Get the start and end times
    start_time = levels_time(i);
    end_time = events_fixed(level_first_block_idx(i));
    
    % Find indices in aligned_runner_time within the range
    idx_range = aligned_runner_time >= start_time & aligned_runner_time < end_time;

    % Set aligned_runner_position to 1
    aligned_runner_position(idx_range) = 1;
end


%% Adjust states for easier visualization
states_adjusted = zeros(size(aligned_runner_position));
for i = 1:length(aligned_runner_position)
    switch round(aligned_runner_position(i))
        case 2
            states_adjusted(i) = med + (max_emg * 0.4);
        case 1
            states_adjusted(i) = med;
        case 0
            states_adjusted(i) = med - (max_emg * 0.4);
    end
end

%% Plot the data
close all

% Plot raw flexor EMG
figure(1);
% subplot(2, 1, 1)
plot(emg_time, emg_ch1, 'DisplayName', 'Channel 1 – Flexor')
plot_title_ch1 = sprintf("%s — Ch1 (FLEXOR)", id);
title(plot_title_ch1)
hold on
plot(runner_time, runner_position_adjusted, 'LineWidth', 1, 'Color', 'r', 'DisplayName', 'Target Position RAW')
xline(events_time(blockStartIDX), 'Color', [0.1 0.1 0.1], 'HandleVisibility','off')
xline(events_time(blockStopIDX), 'LineStyle', '-.', 'Color', [0.1 0.1 0.1], 'HandleVisibility','off')
xline(levels_time, 'Color', [0.1 0.1 0.1], 'LineWidth', 3, 'HandleVisibility', 'off')
xline(events_time(game_start_idx), 'Color', [0.1 0.1 0.1], 'LineWidth', 3, 'HandleVisibility', 'off')
xline(calibration, 'Color', [0.7 0.4 0.5], 'HandleVisibility', 'off')
% xline(tutorial, 'Color', [0.7 0.5 0.2], 'HandleVisibility', 'off')
xline(tutorial_clean(tutorialExpStart), 'Color', [0.7 0.5 0.2], 'HandleVisibility', 'off')
xline(tutorial_clean(tutorialExpStop), 'LineStyle', '-.', 'Color', [0.7 0.5 0.2], 'HandleVisibility', 'off')
xline(start_time, 'Color', [0.2 0.7 0.5], 'DisplayName', 'Start (Tutorial + Game)')
xlabel("Time (s)")
ylabel("EMG (mV)")
legend()
hold off

% Plot adjusted flexor EMG
figure(2);
% subplot(2, 1, 1)
plot(emg_time, emg_ch1, 'DisplayName', 'Channel 1 – Flexor')
plot_title_ch1 = sprintf("%s — Ch1 (FLEXOR)", id);
title(plot_title_ch1)
hold on
% plot(emg_time, emg_ch2 / 10, 'DisplayName', 'Channel 2 – Extensor')
% plot(runner_time_offset_adjusted, runner_position_adjusted, 'LineWidth', 1, 'Color', 'r', 'DisplayName', 'Target Position OFFSET')
% adjusted events
xline(events_fixed(blockStartIDX), 'Color', [0.1 0.1 0.1], 'HandleVisibility','off')
xline(events_fixed(blockStopIDX), 'LineStyle', '-.', 'Color', [0.1 0.1 0.1], 'HandleVisibility','off')
% % raw events
% xline(events_time(blockStartIDX), 'Color', 'b', 'HandleVisibility','off')
% xline(events_time(blockStopIDX), 'LineStyle', '-.', 'Color', 'b', 'HandleVisibility','off')
xline(levels_time, 'Color', [0.1 0.1 0.1], 'LineWidth', 3, 'HandleVisibility', 'off')
xline(events_fixed(game_start_idx), 'Color', [0.1 0.1 0.1], 'LineWidth', 3, 'HandleVisibility', 'off')
% xline(all_changes_runner_time, 'Color', 'magenta', 'LineWidth', 1.5, 'HandleVisibility', 'off')
xline(calibration, 'Color', [0.7 0.4 0.5], 'HandleVisibility', 'off')
% xline(tutorial, 'Color', [0.7 0.5 0.2], 'HandleVisibility', 'off')
xline(tutorial_clean(tutorialExpStart), 'Color', [0.7 0.5 0.2], 'HandleVisibility', 'off')
xline(tutorial_clean(tutorialExpStop), 'LineStyle', '-.', 'Color', [0.7 0.5 0.2], 'HandleVisibility', 'off')
xline(start_time, 'Color', [0.2 0.7 0.5], 'DisplayName', 'Start (Tutorial + Game)')
plot(aligned_runner_time, states_adjusted, 'LineWidth', 1.25, 'Color', 'r', 'DisplayName', 'Target Position INTERPOLATED')
% plot(runner_time, runner_position_adjusted, 'LineStyle', '--', 'LineWidth', 1.25, 'Color', 'black', 'DisplayName', 'Target Position RAW')
xlabel("Time (s)")
ylabel("EMG (mV)")
legend()
hold off

figure(3)
subplot(2, 1, 1)
plot(aligned_runner_time, runner_position_adjusted, 'LineWidth', 1.25, 'DisplayName', 'Target Position RAW')
title("Target Position RAW", 'FontSize', 15)
xlabel("Time (s)")
ylabel("Position")
hold on
xline(levels_time, 'Color', [0.1 0.1 0.1], 'LineWidth', 3, 'HandleVisibility', 'off')
hold off
subplot(2, 1, 2)
plot(aligned_runner_time, states_adjusted, 'LineWidth', 1.25, 'Color', 'r', 'DisplayName', 'Target Position INTERPOLATED')
title("Target Position INTERPOLATED", 'FontSize', 15)
xlabel("Time (s)")
ylabel("Position")
hold on
xline(levels_time, 'Color', [0.1 0.1 0.1], 'LineWidth', 3, 'HandleVisibility', 'off')
hold off

rEMGb=round(emg_time(:),3);
rEMGa=round(events_fixed,3);
emg_events_idx = find(rEMGa == rEMGb);

positionActual = aligned_runner_position;
rPositionb=round(aligned_runner_time(:),1);
rPositiona=round(events_fixed,1);
position_events_idx = find(rPositionb == rPositiona);

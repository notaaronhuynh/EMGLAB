clear all;
close all;

% load('MAPS004_T0_L (2 runs).mat') % test case

[file, path] = uigetfile('*.mat'); % File Selection
file = [path,file];
load(file);

cd(path)

[~,name,ext] = fileparts(file);

%% Segment Events and EMG Data
first_stop = find(strcmp(Events, "STOP"), 1, "first");
last_start = find(strcmp(Events, "START"), 1, "last");
last_game_start = find(strcmp(Events, "game start"), 1, "last");

first_stop_time = EventsTime(first_stop);
last_start_time = EventsTime(last_start);
last_game_start_time = EventsTime(last_game_start);

newEvents = Events([1:first_stop - 1, last_game_start + 1:end]);
newEventsTime = EventsTime([1:first_stop - 1, last_game_start + 1:end]);

emg_mask = (EMG_Time < EventsTime(first_stop - 1)) | (EMG_Time > EventsTime(last_game_start + 1));
newEMG_Time = EMG_Time(emg_mask); % keep aligned times
newCh1EMG = Ch1EMG(emg_mask);
newCh2EMG = Ch2EMG(emg_mask);
newEMG_Data = EMG_Data(:, emg_mask);

%% Segment ConditionInfo
% Extract column 2 as a string array
col7 = string(ConditionInfo(:,7));
info_idx = find(startsWith(col7, '1'));

newConditionInfo = ConditionInfo(info_idx, :);

%% Segment ConditionData
% Extract column 2 as a string array
col2 = string(ConditionInfo(:,2));

pairs = struct();

for i = 1:numel(col2)
    key = col2(i);
    match_idx = find(col2 == key);

    if numel(match_idx) >= 1
        fieldName = char(col2(i)); % full string as field name
        pairs.(fieldName) = match_idx;
    end
end

% fn = fieldnames(pairs);
% 
% for i = 1:numel(fn)
%     val = pairs.(fn{i});
%     if ~isequal(size(val), [2 1])   % If not exactly 2x1
%         pairs = rmfield(pairs, fn{i});
%     end
% end

fn = fieldnames(pairs); % update fn for trimmed pairs struct
mergedData = cell(length(fn), 1);

for k = 1:numel(fn)
    idx = pairs.(fn{k});
    tmp = {}; % temp holder for merged data
    m = [];

    for ii = 1:numel(idx)
        % Append the contents of ConditionData{row,1}{1,1}
        tmp = ConditionData{idx(ii),1}{1,1};
        m = [m; tmp];
    end

    mergedData{k, 1} = m; % store merged result
end

newConditionData = cell(length(fn), 1);

for i = 1:numel(mergedData)
    temp = mergedData{i, 1};
    newConditionData{i, 1}{:} = temp;
end

%% Segment EventThresholds
events_fn = fieldnames(EventThresholds);

% Extract last 2 characters from *all* fieldnames
fn_str = string(events_fn);
last2 = extractBetween(fn_str, strlength(fn_str)-1, strlength(fn_str));

% Group by suffix
[G, suffixes] = findgroups(last2);

newEventThresholds = struct();

for g = 1:max(G)
    members_rel = find(G == g);  % all matching suffix
    names = events_fn(members_rel);

    % Only create a field if there's a Ch1... in this group
    ch1_idx = find(startsWith(names, 'Ch1'), 1);
    if isempty(ch1_idx)
        continue; % skip if no Ch1 version
    end

    merged = [];
    for k = 1:numel(names)
        val = EventThresholds.(names{k});

        % Ensure it's a double and 2D
        val = double(val);

        % First time: just store
        if isempty(merged)
            merged = val;
        else
            % Check row count
            if size(val,1) ~= size(merged,1)
                error('Row size mismatch for %s and %s', names{1}, names{k});
            end
            % Horizontal concat
            merged = [merged val];
        end
    end


    % Store under Ch1's fieldname
    newEventThresholds.(names{ch1_idx}) = merged;
end


ConditionData = newConditionData;
Events = newEvents;
EventsTime = newEventsTime;
EMG_Time = newEMG_Time;
Ch1EMG = newCh1EMG;
Ch2EMG = newCh2EMG;
EMG_Data = newEMG_Data;
ConditionInfo = newConditionInfo;
EventThresholds = newEventThresholds;

new_filename = sprintf('%s (cleaned)%s', name, ext);

save(new_filename, "ConditionData", "ConditionInfo", "time", "sampleRate", "EMG_Data", "EventThresholds", "streams",...
    "Events", "EventsTime", "Thresholds", "ThresholdsTime", "EMG_Time")


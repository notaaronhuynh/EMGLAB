%% Bandpass Filter and Rectify EMG Signal
% University of Rochester School of Medicine & Dentistry
% Movement and Plasticity Lab (MAPL)
% Author: Aaron Huynh
% Last commit: 26 October 2024

%% Usage Notes
%%-----------------------------------------------------------------------%%
%% Example calls to function
%%-----------------------------------------------------------------------%%
%% If using EVA data...
% fields = ConditionInfo(:, 2);
% rawData = ConditionData;
% aFQ = 10; % Hz
% bFQ = 500; % Hz
% Fq = [aFq bFq];
% sr = sampleRate;
% preprocess_emg(rawData, fields, Fq, sr, 'butter')

%% Output is a struct with the fieldnames stored in fields, each field has
%% data that has been filtered and rectified, and the output variable is 
%% stored as 'ProcessedEMG'

%% General usage..
% fields = {'Ch1LeEE', 'Ch1LeEF', 'Ch1LeER', 'Ch1LeFE', 'Ch1LeFF',...
%     'Ch1LeFR', 'Ch1LeRE', 'Ch1LeRF', 'Ch1LeRR', 'Ch2LeEE', 'Ch2LeEF',...
%     'Ch2LeER', 'Ch2LeFE', 'Ch2LeFF', 'Ch2LeFR', 'Ch2LeRE', 'Ch2LeRF',...
%     'Ch2LeRR'};
% rawData = ConditionData;
% aFQ = 10; Hz
% bFQ = 500; Hz
% Fq = [aFq bFq];
% sr = 2000; % Hz
% preprocess_emg(rawData, fields, Fq, sr, 'butter')

%% Output is a struct with the fieldnames stored in fields, each field has
%% data that has been filtered and rectified, and the output variable is 
%% stored as 'ProcessedEMG'

function [ProcessedEMG] = preprocess_emg_app(rawData, fields, Fq, sr, filter, varargin)
%%-----------------------------------------------------------------------%%
%% Handle varargin
%%-----------------------------------------------------------------------%%
writeStruct = false; % default is to not write data to csv
if nargin > 0 % if there are any arguments in the function
    for i = 1:length(varargin)
        % look for a specific varargin and the argument after it
        if strcmp(varargin{i}, 'WriteStruct') && i + 1 <= length(varargin)
            writeStruct = true;
            fullfileName = varargin{i + 1};

        end
    end
else
    error("varargin not recognized.")
end

%%-----------------------------------------------------------------------%%
%% Bandpass Filter & Rectify
%%-----------------------------------------------------------------------%%
nyq_freq = sr/2;
filter_lower = lower(filter);

ProcessedEMG = struct(); % empty struct to store bandpassed trial data

events = fields;
% events_processed = events + "_processed";

aFq = Fq(1);

if length(Fq) > 1
    bFq = Fq(2);
end

for i = 1:length(rawData)

    processed_names = events{i};
    data = rawData{i, 1}{1}'; % EMG recordings for each condition in a given trial

    if strcmp(filter_lower, 'butter')
        [b, a] = butter(4, [(aFq/nyq_freq) (bFq/nyq_freq)], 'bandpass'); % bandpass butterworth filter
        filtered = filtfilt(b, a, data);
    elseif strcmp(filter_lower, 'bandpass')
        filtered = bandpass(data, [aFq bFq], fs); % simple bandpass
    elseif strcmp(filter_lower, 'low')
        [b, a] = butter(4, (aFq/nyq_freq)); % low-pass butterworth filter
        filtered = filtfilt(b, a, data);
    elseif strcmp(filter_lower, 'high')
        [b, a] = butter(4, (aFq/nyq_freq), 'high'); % high-pass butterworth filter
        filtered = filtfilt(b, a, data);
    else
        error("Filter option unavailable. Choose butter or bandpass.")
    end

    rectified = abs(filtered);

    ProcessedEMG.(processed_names) = rectified;

    % assign values to variable in a specific workspace
    assignin('base', 'ProcessedEMG', ProcessedEMG);

    if writeStruct
        % outfile_name = sprintf('MAPS%s_%s_%s_%s.csv', subject, timepoint, recording_date, condition_name_processed);
        % output_processed = fullfile(processed_dir, outfile_name);
        output = fullfileName;

        % Create headers as a cell array
        num_trials = size(data, 2); % count number of trials for a given event
        headers = arrayfun(@(x) sprintf('Trial%d', x), 1:num_trials, 'UniformOutput', false); % print header for correlation table

        % Write the headers to the CSV file
        fid = fopen(output, 'w'); % Open the file for writing
        fprintf(fid, '%s,', headers{1:end-1}); % Write headers
        fprintf(fid, '%s\n', headers{end}); % Write the last header and move to a new line
        fclose(fid); % Close the file

        % Write the data to the CSV file
        writematrix(ProcessedEMG.(processed_names), output, 'WriteMode', 'append');
    end

end

end
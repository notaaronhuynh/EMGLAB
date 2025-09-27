%% Format data from mat file into struct
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
% mat2struct(fields, ConditionData, 'RawEMG')
%
%% Output is a struct with the fieldnames stored in fields, each field has
%% data from ConditionData, and the output variable is stored as 'RawEMG'
%
%% General usage..
% fields = {'Ch1LeEE', 'Ch1LeEF', 'Ch1LeER', 'Ch1LeFE', 'Ch1LeFF',...
%     'Ch1LeFR', 'Ch1LeRE', 'Ch1LeRF', 'Ch1LeRR', 'Ch2LeEE', 'Ch2LeEF',...
%     'Ch2LeER', 'Ch2LeFE', 'Ch2LeFF', 'Ch2LeFR', 'Ch2LeRE', 'Ch2LeRF',...
%     'Ch2LeRR'};
% rawData = ConditionData;
% mat2struct(fields, ConditionData, 'RawEMG')
%
%% Output is a struct with the fieldnames stored in fields, each field has
%% data from ConditionData, and the output variable is stored as 'RawEMG'

function dataStruct = mat2struct(fields, rawData, varargin)

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

dataStruct = struct(); % empty struct to store data

for i = 1:numel(rawData)
    events = fields{i};
    data = rawData{i, 1}{1}';

    % Store data in the structure
    dataStruct.(events) = data;

    % assign values to variable in a specific workspace
    % assignin('base', workspaceName, dataStruct);

    if writeStruct
        % outfile_name = sprintf('MAPS%s_%s_%s_%s_Raw.csv', subject, timepoint, recording_date, condition_name);
        % output_raw = fullfile(raw_dir, outfile_name);
        output = fullfileName;

        % Create headers as a cell array
        num_trials = size(rawData, 2); % Assuming each column is a trial
        headers = arrayfun(@(x) sprintf('Trial%d', x), 1:num_trials, 'UniformOutput', false);

        % Write the headers to the CSV file
        fid = fopen(output, 'w'); % Open the file for writing
        fprintf(fid, '%s,', headers{1:end-1}); % Write headers separated by commas
        fprintf(fid, '%s\n', headers{end}); % Write the last header and move to a new line
        fclose(fid); % Close the file

        % Write the data to the CSV file
        writematrix(dataStruct.(events), output, 'WriteMode', 'append');
    end

end
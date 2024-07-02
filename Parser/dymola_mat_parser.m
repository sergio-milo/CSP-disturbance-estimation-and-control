close all
clearvars
clc

% Open the text file for reading
fileID = fopen('ray_tracing.txt', 'r');

% Check if the file is open successfully
if fileID == -1
    error('Cannot open the file.');
end

% Initialize a cell array to store lines
lines = {};

% Read lines until the end of the file
while ~feof(fileID)
    line = string(fgetl(fileID));
    line = strrep(line, ' ', '');
    line = strrep(line, "'", '');
    line = strrep(line, '[', '(');
    line = strrep(line, ']', ')');
    line = strrep(line, '.', '_');
    lines = [lines; string({line})]; %#ok<AGROW>
end

% Close the file
fclose(fileID);

parse_strings = ["Time","flow2DFV_H_Maff_feedback_T(", "sourceMassFlow_w", "heatSource2DFV_flux(", "heatSource2DFV_Q_in_tot"];
rename_strings = ["time","temperature", "flow", "heatmap", "d_tot"];
mask = false(size(lines));

for i = 1:length(parse_strings)
    mask = mask | contains(lines, parse_strings(i));
end

for i = 1:length(parse_strings)
    lines = strrep(lines, strrep(parse_strings(i), '(', ''), rename_strings(i));
end

indexes = find(mask);

load('ray_tracing.mat');

for i = 1:length(indexes)
    var_name = lines(indexes(i));
    if contains(var_name, '(')
        tmp = strsplit(var_name, '(');
        var_name = tmp(1);
        tmp = strsplit(tmp(2), ')');
        num = round(str2double(tmp(1)));
        out_struct.(var_name)(num, :) = s(:, indexes(i));
    else       
        out_struct.(var_name) = s(:, indexes(i));
    end
end

save("m" + "ray_tracing.mat", '-struct', 'out_struct');

%function to import data from .dat files (marquardt)

function [f,re,Im,r,h] = read_data(path,filename)

%% Import data from text file.
% Script for importing data from the following text file:
%
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2016/03/20 01:31:58

%% Initialize variables.
delimiter = {',',' '};
startRow = 28;

%% Format string for each line of text:
%   column6: text (%s)
%	column13: text (%s)
%   column20: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%*s%*s%*s%s%*s%*s%*s%*s%*s%*s%s%*s%*s%*s%*s%*s%*s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');
%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
f = (dataArray{:, 1});
re = dataArray{:, 2};
Im = dataArray{:, 3};

f=cellfun(@str2num,f);
re=cellfun(@str2num,re);
Im=cellfun(@str2num,Im);

%% Clear temporary variables
clearvars delimiter startRow formatSpec fileID dataArray ans;

%importing r

%% Import data from text file.
% Script for importing data from the following text file:
%
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2016/03/20 01:38:12

%% Initialize variables.
delimiter = {',',' '};
startRow = 8;
endRow = 12;

%% Format string for each line of text:
%   column3: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%f%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'HeaderLines', startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
r = dataArray{:, 1};

%% Clear temporary variables
clearvars delimiter startRow endRow formatSpec fileID dataArray ans;

%importing h
%% Import data from text file.
% Script for importing data from the following text file:
%
%    E:\project\correct code\Hem data\hem_1_marq.dat
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2016/03/20 01:39:14

%% Initialize variables.
delimiter = {',',' '};
startRow = 4;
endRow = 4;

%% Format string for each line of text:
%   column4: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*s%*s%*s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'HeaderLines', startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
h = (dataArray{:, 1});
h=cellfun(@str2num,h);

%% Clear temporary variables
clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;




end
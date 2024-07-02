function tnlist(n, filename)
%tnlist   list trajectory names with preceding row numbers
%
%         tnlist(n) lists string matrix "n" with preceding row numbers
%
%See also: traj, nhead, tnindex.

%    Copyright (c) 1995,1996,1997 by DLR.
%    Copyright (C) 1997-2001 Dynasim AB.
%    All rights reserved.

%TNLIST   List trajectory names with preceding row numbers.
%
%   TNLIST(N) lists string matrix N with preceding row numbers.
%   TNLIST(N, FILENAME) writes the output to the file specified by FILENAME.

% Check the number of input arguments
if nargin < 1
    error('tnlist requires at least one argument');
end

% If filename is not provided, print to the command window
if nargin < 2
    fid = 1; % Standard output (command window)
else
    % Open the file for writing
    fid = fopen(filename, 'wt');
    if fid == -1
        error('Cannot open file for writing');
    end
end

sizen = size(n,1);
for i = 1:sizen
    fprintf(fid,'%s\n', n(i,:));
end
fprintf(fid, '\n');

% Close the file if filename is provided
if nargin == 2
    fclose(fid);
end
end
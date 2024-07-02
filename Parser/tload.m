function [s,n] = tload(matfile,quiet)
%tload    load trajectory into Matlab work space
%
%         [s,n]=tload;       loads trajectory data 'dsres.mat'.
%         [s,n]=tload(file); loads trajectory data '<file>.mat'.
%         where,
%            s: trajectory matrix (column i is data of trajectory i) 
%            n: trajectory names (row i is from column i of "s").
%
%Example:  [s,n] = tload('result');
%          [s,n] = tload('result.mat');
%          [s,n] = tload('result.mat',1); % Same but no messages.
%
%See also: traj, tsave.

%    Copyright (c) 1995,1996,1997 by DLR.
%    Copyright (C) 1997-2001 Dynasim AB.
%    All rights reserved.

% For version 1.0, it is expected that three matrices are on file:
%    Aclass (contains class-name), names and data.
% For version 1.1 the following matrices are expected:
%    Aclass, name, description, dataInfo, data_1, data_2, ...

% determine file name
if nargin<2
  quiet=0;
end
  if nargin < 1
     file = 'dsres.mat';
  else
     matfile = lower(matfile);
     ii = findstr(matfile,'.');
     if isempty(ii)
        file = [matfile,'.mat'];
     elseif strcmp(matfile(1,ii:size(matfile,2)),'.mat')
        file = matfile;
     else
       error( ['filename (= "',matfile,'") has not extension ".mat"'])
     end
  end

% check file existence
   if exist( file ) ~= 2
     error( sprintf( ['"', file, '" does not exist in \n   %s'], pwd) )
   end

% load data
  eval( [ 'load ' file ] );

% read Aclass variable
  matlabVersion = version;
  if exist('Aclass') ~= 1
     if matlabVersion(1,1)=='4'
        if exist('class') ~= 1
          error( ['no traj. on file "' file '" ("Aclass" is missing).'])
        else
           Aclass = class;
        end
     else
       error( ['no trajectory on file "' file '" ("Aclass" is missing).'])
     end
  end

% check whether file has correct class name
  classReq = 'Atrajectory';
  ncol1 = size(classReq,2);
  [nrow2,ncol2] = size(Aclass);
  if ncol1 < ncol2 
     classReq = [ classReq, blanks(ncol2-ncol1) ];
  elseif ncol1 > ncol2
     Aclass = [ Aclass, blanks(ncol1-ncol2) ];
     ncol2  = size(Aclass, 2);
  end
  if nrow2 < 2 then
     error( [ 'file "' file '" is not of class ' classReq ] )
  elseif Aclass(1,:) ~= classReq(1,:)
     error( [ 'file "' file '" is not of class ' classReq ] )
  end

% Check version number
  if ['1.0'] == Aclass(2,1:3)
     vers = 0;
  elseif ['1.1'] == Aclass(2,1:3)
     vers = 1;
  else
    error( [ 'file "' file '" has wrong version number ' Aclass(2,:) ] )
  end

% Determine whether matrices have to be transposed
  if nrow2 < 4
     trans = 0;
  elseif Aclass(4,1:8) == ['binTrans']
     trans = 1;
  else
     trans = 0;
  end

% Action according to version number
  if vers == 0 
     % Check existance of data 
       if exist('data') ~= 1
         error( ['no traj. on file "' file '" ("data" is missing).'] )
       elseif exist('names') ~= 1
         error( ['no traj. on file "' file '" ("names" is missing).'] )
       end

     % copy data according to storage type
       if trans == 0
          s = data;
          n = names;
       else
          s = data';
          n = names';
       end
  else
     % Check existance of name and dataInfo matrix
       if exist('name') ~= 1
         error( ['no traj. on file "' file '" (matrix "name" missing).'])
       elseif exist('dataInfo') ~= 1
         error( ['no traj. on file "' file '" ("dataInfo" missing).'] )
       end

     % Copy name
       if trans == 0
          n = name;
       else
          n = name';
          dataInfo = dataInfo';
       end

     % Determine number of data matrices 
       ndata = max(dataInfo(:,1));
       ncol  = size(n,1);
       if ndata > 1 
          if trans == 0
             [nrow1,ncol1] = size(data_2);
             t1 = data_2(:,1);
          else
             [ncol1,nrow1] = size(data_2);
             t1 = data_2(1,:)';
          end
          s = [t1, zeros(nrow1,ncol-1)];
       else
         if trans == 0
             [nrow1,ncol1] = size(data_1);
             t1 = data_1(:,1);
          else
             [ncol1,nrow1] = size(data_1);
             t1 = data_1(1,:)';
          end
          s = [t1, zeros(nrow1,ncol-1)];
       end

     % Store matrices data_i
       for i=1:ndata
          % Determine matrix name
            if trans == 0
               eval( ['data = data_' int2str(i) ';'] );
            else
               eval( ['data = data_' int2str(i) ''';'] );
            end

          if size(data,1) > 0 
             % Make equal time axis
               [s,data] = tsame(s,data);

              % Store data in matrix s
              I = find( dataInfo(:,1)==i & (1:ncol)'>1);
              for j=I.'
                 s(:,j) = sign(dataInfo(j,2))*data(:,abs(dataInfo(j,2)));
              end
           end
       end
  end
if ~quiet
% print info message
  if file(1,1) == '/' | file(1,1) == '\' | file(1,2) == ':'
     disp( [ '> ', file, ' loaded.' ] )
  else
     machine = computer;
     if machine(1,1:2) == 'PC'
        disp( [ '> ', lower(pwd), '\', lower(file), ' loaded.'] )
     else
        disp( [ '> ', pwd, '/', file, ' loaded.' ] )
     end
  end
end

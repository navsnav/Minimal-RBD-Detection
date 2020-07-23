% RBD Sleep Detection Toolbox, version 1.0, November 2020
% Released under the GNU General Public License
%
% Copyright (C) 2020  Navin Cooray
% Institute of Biomedical Engineering
% Department of Engineering Science
% University of Oxford
% navin.cooray@eng.ox.ac.uk
%
%
% Referencing this work
% Navin Cooray, Fernando Andreotti, Christine Lo, Mkael Symmonds, Michele T.M. Hu, & Maarten De % Vos (in review). Detection of REM Sleep Behaviour Disorder by Automated Polysomnography Analysis. Clinical Neurophysiology.
%
% Last updated : 23-7-2020
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
function [result,descriptor]=histogram2(x,y,descriptor)
%HISTOGRAM2 Computes the two dimensional frequency histogram of two
%           row vectors x and y.
%   [RESULT,DESCRIPTOR] = HISTOGRAM2(X,Y) or
%   [RESULT,DESCRIPTOR] = HISTOGRAM2(X,Y,DESCRIPTOR) or
%where
%   DESCRIPTOR = [LOWERX,UPPERX,NCELLX;
%                 LOWERY,UPPERY,NCELLY]
%
%   RESULT       : A matrix vector containing the histogram
%   DESCRIPTOR   : The used descriptor
%
%   X,Y          : The row vectors to be analyzed
%   DESCRIPTOR   : The descriptor of the histogram
%     LOWER?     : The lowerbound of the ? dimension of the histogram
%     UPPER?     : The upperbound of the ? dimension of the histogram
%     NCELL?     : The number of cells of the ? dimension of the histogram
%
%   See also: http://www.cs.rug.nl/~rudy/matlab/

%   R. Moddemeijer 
%   Copyright (c) by R. Moddemeijer
%   $Revision: 1.2 $  $Date: 2001/02/05 09:54:29 $

if nargin <1
   disp('Usage: RESULT = HISTOGRAM2(X,Y)')
   disp('       RESULT = HISTOGRAM2(X,Y,DESCRIPTOR)')
   disp('Where: DESCRIPTOR = [LOWERX,UPPERX,NCELLX;')
   disp('                     LOWERY,UPPERY,NCELLY]')
   return
end

% Some initial tests on the input arguments

[NRowX,NColX]=size(x);

if NRowX~=1
  error('Invalid dimension of X');
end;

[NRowY,NColY]=size(y);

if NRowY~=1
  error('Invalid dimension of Y');
end;

if NColX~=NColY
  error('Unequal length of X and Y');
end;

if nargin>3
  error('Too many arguments');
end;

if nargin==2
  minx=min(x);
  maxx=max(x);
  deltax=(maxx-minx)/(length(x)-1);
  ncellx=ceil(length(x)^(1/3));
  miny=min(y);
  maxy=max(y);
  deltay=(maxy-miny)/(length(y)-1);
  ncelly=ncellx;
  descriptor=[minx-deltax/2,maxx+deltax/2,ncellx;miny-deltay/2,maxy+deltay/2,ncelly];
end;

lowerx=descriptor(1,1);
upperx=descriptor(1,2);
ncellx=descriptor(1,3);
lowery=descriptor(2,1);
uppery=descriptor(2,2);
ncelly=descriptor(2,3);

if ncellx<1 
  error('Invalid number of cells in X dimension')
end;

if ncelly<1 
  error('Invalid number of cells in Y dimension')
end;

if upperx<=lowerx
  error('Invalid bounds in X dimension')
end;

if uppery<=lowery
  error('Invalid bounds in Y dimension')
end;

result(1:ncellx,1:ncelly)=0;

xx=round( (x-lowerx)/(upperx-lowerx)*ncellx + 1/2 );
yy=round( (y-lowery)/(uppery-lowery)*ncelly + 1/2 );
for n=1:NColX
  indexx=xx(n);
  indexy=yy(n);
  if indexx >= 1 & indexx <= ncellx & indexy >= 1 & indexy <= ncelly
    result(indexx,indexy)=result(indexx,indexy)+1;
  end;
end;

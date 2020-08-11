% Copyright (C) 2020  Navin Cooray
% Institute of Biomedical Engineering
% Department of Engineering Science
% University of Oxford
% navin.cooray@eng.ox.ac.uk
%
%
% Referencing this work
% Navin Cooray, Fernando Andreotti, Christine Lo, Mkael Symmonds, Michele T.M. Hu, & Maarten De % Vos (in review). Proof of Concept: Screening for REM Sleep Behaviour Disorder with a Minimal Set of Sensors. Clinical Neurophysiology.
%
% Last updated : 11-8-2020
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

function     cost_opt = find_optimise_cost(rf,posterior,Ytrn)

    
firstparam = linspace(1,0,6);     %list of places to search for first parameter
secondparam = linspace(1,0,6);         %list of places to search for second parameter
thirdparam = linspace(1,2,6);     %list of places to search for third parameter
fourthparam = linspace(1,2,6);         %list of places to search fourth second parameter


[F,S,Q,R] = ndgrid(firstparam, secondparam,thirdparam,fourthparam);
fitresult = arrayfun(@(p1,p2,p3,p4) fittingfunction(p1,p2,p3,p4,rf,posterior,Ytrn), F,S,Q,R); %run a fitting on every pair fittingfunction(F(J,K), S(J,K))



[maxval, maxidx] = max(fitresult(:));
bestFirst = F(maxidx);
bestSecond = S(maxidx);
bestThird = Q(maxidx);
bestFourth = R(maxidx);


cost_opt = [0,1,bestFirst;1,0,bestSecond;bestThird,bestFourth,0];
% cost_opt = [0,1,0;1,0,0;1,1,0];

end



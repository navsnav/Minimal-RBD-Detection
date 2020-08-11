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

function fitresult = fittingfunction(p1,p2,p3,p4,rf,posterior,Ytrn)

    cost = [0,1,p1;1,0,p2;p3,p4,0];

    expCost = posterior*cost;
    [~,classIndex] = min(expCost,[],2);
    Yhat2 = rf.ClassNames(classIndex);
    Yhat2 = str2num(cell2mat(Yhat2));
    
    metrics  = process_classification_results2(Yhat2==5, Ytrn==5);
    
    ConfMat_Class_Summary = confusionmat(Yhat2, Ytrn, 'order', [0 2 5]);
    kappa = kappa_result(ConfMat_Class_Summary);    
%     fitresult = metrics(6); %F1
%     fitresult = metrics(7); %Kappa
    fitresult = metrics(4); %Precision
%     fitresult = kappa; %Precision
    
    

end
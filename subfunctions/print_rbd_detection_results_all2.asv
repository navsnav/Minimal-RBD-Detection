function    print_rbd_detection_results_all2(results_f_all,names_all,tablename,print_figures,print_folder)
% This function compares rbd detection perforance between 2 sets of results
% and produces a table.
%
% Inputs:
%  results_f_est - 1st set of results with mean in first column and
%                  standard  deviation in the 2nd column. In the following
%                  row order: 'Accuracy','Sensitivity','Specificity','Precision','Recall','F1'
%  results_f_new - 2nd set of results, same as above. 
%  name1 - text name of 1st set of results
%  name2 - text name of 2nd set of results
%  tablename - name of table
%  print_figures - flag to print/save table
%  print_folder -  folder where table is saved
%
% --
% RBD Sleep Detection Toolbox, version 1.0, November 2018
% Released under the GNU General Public License
%
% Copyright (C) 2018  Navin Cooray
% Institute of Biomedical Engineering
% Department of Engineering Science
% University of Oxford
% navin.cooray@eng.ox.ac.uk
%
%
% Referencing this work
% Navin Cooray, Fernando Andreotti, Christine Lo, Mkael Symmonds, Michele T.M. Hu, & Maarten De % Vos (in review). Detection of REM Sleep Behaviour Disorder by Automated Polysomnography Analysis. Clinical Neurophysiology.
%
% Last updated : 15-10-2018
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
%% RBD Detection table (annotated)




for t=1:length(names_all)
overall_performance_rbd_all = [mean(results_f_all{t},1);std(results_f_all{t},1)];

 op_t_rbd_all_cell{t} = {[num2str(overall_performance_rbd_all(1,1),2),repmat('±',1,1),num2str(overall_performance_rbd_all(2,1),2)],...
[num2str(overall_performance_rbd_all(1,2),2),repmat('±',1,1),num2str(overall_performance_rbd_all(2,2),2)],...
[num2str(overall_performance_rbd_all(1,3),2),repmat('±',1,1),num2str(overall_performance_rbd_all(2,3),2)],...
[num2str(overall_performance_rbd_all(1,4),2),repmat('±',1,1),num2str(overall_performance_rbd_all(2,4),2)],...
[num2str(overall_performance_rbd_all(1,5),2),repmat('±',1,1),num2str(overall_performance_rbd_all(2,5),2)],...
[num2str(overall_performance_rbd_all(1,6),2),repmat('±',1,1),num2str(overall_performance_rbd_all(2,6),2)],...
[num2str(overall_performance_rbd_all(1,7),2),repmat('±',1,1),num2str(overall_performance_rbd_all(2,7),2)]};
end


op_t_rbd_det = cell2table([op_t_rbd_stream_cell;op_t_rbd_ai_cell;op_t_rbd_qma_cell;op_t_rbd__est_cell;op_t_rbd_new_cell;op_t_rbd_ecg_cell;op_t_rbd_all_cell],'VariableNames',{'Accuracy','Sensitivity','Specificity','Precision','Recall','F1','Kappa'},...
    'RowNames',{name5,name0,name6,name1,name2,name3,name4});

fig_rbd_det = figure('units','normalized','outerposition',[0 0 1 1]);

uitable(fig_rbd_det,'Data', [op_t_rbd_stream_cell;op_t_rbd_ai_cell;op_t_rbd_qma_cell;op_t_rbd__est_cell;op_t_rbd_new_cell;op_t_rbd_ecg_cell;op_t_rbd_all_cell],'ColumnName',op_t_rbd_det.Properties.VariableNames,...
'RowName',op_t_rbd_det.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
if (print_figures), saveas(fig_rbd_det,strcat(print_folder,'\',tablename),'png'), end

end
function [Yhat_Results,RBD_Yhat_Results,RBD_Auto_Yhat_Results,All_Confusion] = RBD_Detection(Sleep_table_Pre,Sleep_Struct,rbd_group,indices,folds,SS_Features,RBD_Detection_Feats,n_trees,view_results,print_figures,print_folder,save_data,outfilename,display_flag,prior_mat)


% Copyright (c) 2018, Navin Cooray (University of Oxford)
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%
% 3. Neither the name of the University of Oxford nor the names of its
%    contributors may be used to endorse or promote products derived
%    from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
% LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
% DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
% THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%	Contact: navsnav@gmail.com
%	Originally written by Navin Cooray 19-Sept-2018
% Input:
%       Sleep_table_Pre:    Sleep table with all features and subjects (preprocessed to have no nans/infs).
%       Sleep_Struct: 	Structure containing all features and feature names
%       rbd_group:      Participant condition/diagnosis (0: Healthy control, 1: RBD participant). 
%       indices:        Indicised for cross-fold validation for sleep staging and RBD detection.
%       folds:          Number of folds for cross-fold validation
%       SS_Features:    Features indicies to be used for automated sleep staging
%       EMG_est_feats:  Indicies of Established EMG features in RBD detection 
%       EMG_feats:      Indicies of features to compare for RBD detection
%       n_trees:        Number of trees for Random Forest training
%       view_results:   Flag for displaying results (0: no, 1: yes).
%       print_figures:  Flag for saving figures (0: no, 1: yes).
%       save_data:      Save data and results in mat format
%       outfilename:    Filename for saved data  
% Output:
%       Yhat_Results:           Automated sleep staging results. 
%       EMG_Yhat_Results:       RBD Detection results using new features.  
%       EMG_Yhat_Results:       RBD Detection results using new
%                               features and manually annoated sleep
%                               stages.
%       EMG_est_Yhat_Results:   RBD Detection results using established
%                               features and manually annoated sleep
%                               stages.
%       EMG_Auto_Yhat_Results:  RBD Detection results using new
%                               features and automated sleep staging.
%       EMG_Auto_est_Yhat_Results:   RBD Detection results using established
%                               features and automated sleep staging.
%       All_Confusion:          Confusion matrices of automated sleep staging.  

Sleep_table = Sleep_table_Pre;
Save_Data_Name = outfilename;
numStates = size(unique(Sleep_table_Pre.AnnotatedSleepStage),1);
if print_figures, mkdir(print_folder), end

Sleep = table2array(Sleep_table);
[patients,ia,ic] = unique(Sleep_table.SubjectIndex);

num_rbd_feats = size(RBD_Detection_Feats.feats,2);
rbd_rf = {}; 
rf_importance_all = {};
rbd_rf_importance_all = {};
%%

%% Initialise

Yhat_Results =[];
Yhat_REM_Results = [];
votes_Results = [];
votes_REM_Results = [];
importance_Results = [];
Yhat_Results = zeros(size(Sleep,1),1);
votes_Results = zeros(size(Sleep,1),numStates);
EMG_Metric = table;

RBD_Yhat_CV = table;
RBD_Auto_Yhat_CV = table;
    
RBD_Yhat_Results = ones(num_rbd_feats,size(rbd_group,1),1)*-1;
RBD_votes_Results = zeros(num_rbd_feats,size(rbd_group,1),2);


EMG_Auto_Metric = table;

RBD_Auto_Yhat_Results = ones(num_rbd_feats,size(rbd_group,1))*-1;
RBD_Auto_votes_Results = zeros(num_rbd_feats,size(rbd_group,1),2);

RBD_Auto_Yhat = table;

results_f_ai = zeros(folds,7);
results_f_stream = zeros(folds,7);
results_f_qma = zeros(folds,7);


results_f_ai_auto = zeros(folds,7);
results_f_stream_auto = zeros(folds,7);
results_f_qma_auto = zeros(folds,7);

results_f_rem_auto = zeros(folds,7);

results_f = zeros(num_rbd_feats,folds,7);
results_f_auto = zeros(num_rbd_feats,folds,7);

posterior_struct = {};
if isempty(prior_mat)
    all_cost_opt = [];    
else 
    all_cost_opt = zeros(3,3,folds);   
end

for out=1:folds
    disp(['Fold: ',num2str(out)]);
    PatientTest = (indices==out); %patient id for testing
    PatientTrain = (indices~=out);%patient id for training
    
    PatientTest_idx = ismember(Sleep(:,1),patients(PatientTest)); %patient index for testing
    PatientTrain_idx = ismember(Sleep(:,1),patients(PatientTrain)); %patient index for training
    
    %% Train set (Sleep Staging % RBD Detection)    
    Xtrn = Sleep_table_Pre(PatientTrain_idx,:);
    Ytrn = table2array(Sleep_table_Pre(PatientTrain_idx,7));
    
    %% Testing set (Sleep Staging & RBD Detection)    
    Xtst = Sleep(PatientTest_idx,SS_Features);
    Ytst = Sleep(PatientTest_idx,7);
    tst_condition = Sleep(PatientTest_idx,6);       

    %% Train Sleep Stage RF

% Configure  parameters
%     predict_all=true;
%     extra_options.predict_all = predict_all;
%     extra_options.importance = 1; %(0 = (Default) Don't, 1=calculate)    
%     mtry = floor(sqrt(length(SS_Features))); %number of features used to creates trees
%     rf = classRF_train(Xtrn, Ytrn,n_trees,mtry,extra_options);  

    %Matlab Trees
%     prior_mat = 1;

    [rf,rf_importance,posterior] = Train_SleepStaging_RF(n_trees,Xtrn,SS_Features,Ytrn,prior_mat);    
    posterior_struct{out} = posterior;

    % Optimise cost matrix
    if isempty(prior_mat)
        cost_opt = [];          
    else
        cost_opt = find_optimise_cost(rf,posterior,Ytrn);
        all_cost_opt(:,:,out) = cost_opt;

    end
    
    %% Train RBD Detection RF from annotated sleep stages (training set)
       
    EMG_Table = Calculate_RBD_Values_table(Xtrn);
    EMG_Table_Pre = EMG_Table;
    % Training Set for RBD Detection
    EMG_Ytrn = table2array(EMG_Table_Pre(:,2));
    
 
    for t = 1:num_rbd_feats
       [rbd_rf{t} rbd_rf_importance{t}] = Train_RBDDetection_RF(n_trees,EMG_Table,RBD_Detection_Feats.feats{t},EMG_Ytrn);
    end
 %% Test Sleep Staging 

    % Matlab Trees
    [Yhat,votes] = Predict_SleepStaging_RF(rf,Xtst);

    if ~isempty(prior_mat)
        expCost = votes*cost_opt;
        [~,classIndex] = min(expCost,[],2);
        Yhat2 = rf.ClassNames(classIndex);
        Yhat = str2num(cell2mat(Yhat2));    
    end
    
 %% Test RBD Detection using Annotated Sleep Staging
    EMG_Annotated_Test_Table = Calculate_RBD_Values_table(Sleep_table_Pre(PatientTest_idx,:));    
    EMG_Ytst = table2array(EMG_Annotated_Test_Table(:,2));
    RBD_Xtst = {};
    RBD_Yhat = {};
    RBD_votes = {};
    for t = 1:num_rbd_feats
        RBD_Xtst{t} = table2array(EMG_Annotated_Test_Table(:,RBD_Detection_Feats.feats{t})); 

        [RBD_Yhat{t},RBD_votes{t}] = Predict_RBDDetection_RF(rbd_rf{t},RBD_Xtst{t},strrep(RBD_Detection_Feats.labels{t},' ','_'));    
    end
    
      
 %% Test RBD Detection using Automatic Sleep Staging
    
    Sleep_table_automatic = Sleep_table_Pre(PatientTest_idx,:);
    Sleep_table_automatic.AnnotatedSleepStage = Yhat; %Automatic sleep staging
    
    % Generate Test values based on automatic classified Sleep Staging
    EMG_Auto_Test_Table = Calculate_RBD_Values_table(Sleep_table_automatic);

    
    RBD_Auto_Yhat = {};
    RBD_Auto_votes = {};
    RBD_Auto_Xtst = {};
    for t = 1:num_rbd_feats
        RBD_Auto_Xtst{t} = table2array(EMG_Auto_Test_Table(:,RBD_Detection_Feats.feats{t})); 

        [RBD_Auto_Yhat{t},RBD_Auto_votes{t}] = Predict_RBDDetection_RF(rbd_rf{t},RBD_Auto_Xtst{t},strrep(RBD_Detection_Feats.labels{t},' ','_'));    
    end    
    
    
    %% Store Results  
    % Automated Sleep Staging
    Yhat_Results(PatientTest_idx) =  Yhat;
    votes_Results(PatientTest_idx,:) = votes;
    importance_Results(:,:,out) = [rf_importance];            
    % RBD Detection using Annoated Sleep Staging
    
    for t=1:num_rbd_feats
        RBD_Yhat_Results(t,PatientTest) = table2array(RBD_Yhat{t});    
        RBD_votes_Results(t,PatientTest,:) = RBD_votes{t} ;      
    end
    EMG_Metric = [EMG_Metric;EMG_Annotated_Test_Table];

    for t=1:num_rbd_feats
        rbd_rf_importance_all{t,out} = rbd_rf_importance{t};    
    end    
    
    % RBD Detection using Automatic Sleep Staging
  
    for t=1:num_rbd_feats
        RBD_Auto_Yhat_Results(t,PatientTest) = table2array(RBD_Auto_Yhat{t});    
        RBD_Auto_votes_Results(t,PatientTest,:) = RBD_Auto_votes{t} ;      
    end    
    EMG_Auto_Metric = [EMG_Auto_Metric;EMG_Auto_Test_Table];
    new_row = [];
    new_row2 = [];
    for t=1:num_rbd_feats
        new_row = [new_row,RBD_Yhat{t}];
        new_row2 = [new_row2,RBD_Auto_Yhat{t}];
    end
    RBD_Yhat_CV = [RBD_Yhat_CV;new_row];
    RBD_Auto_Yhat_CV = [RBD_Auto_Yhat_CV;new_row2];
    
 
%% RBD Detection Results

    results_f_ai(out,:)  = process_classification_results(EMG_Annotated_Test_Table.AI_REM<0.9, rbd_group(PatientTest)==1);
    results_f_stream(out,:)  = process_classification_results(EMG_Annotated_Test_Table.Stream>25, rbd_group(PatientTest)==1);
    results_f_qma(out,:)  = process_classification_results(max([EMG_Annotated_Test_Table.MAD_Dur,EMG_Annotated_Test_Table.MAD_Per],[],2)>0.07, rbd_group(PatientTest)==1);
    
    results_f_ai_auto(out,:)  = process_classification_results(EMG_Auto_Test_Table.AI_REM<0.9, rbd_group(PatientTest)==1);
    results_f_stream_auto(out,:)  = process_classification_results(EMG_Auto_Test_Table.Stream>25, rbd_group(PatientTest)==1);
    results_f_qma_auto(out,:)  = process_classification_results(max([EMG_Annotated_Test_Table.MAD_Dur,EMG_Annotated_Test_Table.MAD_Per],[],2)>0.07, rbd_group(PatientTest)==1);  
    results_f_rem_auto(out,:) = process_classification_results(Yhat==5, Ytst==5);     
    
    for t=1:num_rbd_feats
        results_f(t,out,:) = process_classification_results(table2array(RBD_Yhat{t})==1, rbd_group(PatientTest)==1);
        results_f_auto(t,out,:) = process_classification_results(table2array(RBD_Auto_Yhat{t})==1, rbd_group(PatientTest)==1);        
    end
    
end

%% Save Data
Sleep_names = Sleep_table.Properties.VariableNames;
EMG_Table_Names = EMG_Table.Properties.VariableNames;


%% Print Sleep Stage Results
states = unique(Sleep_table_Pre.AnnotatedSleepStage);
if (view_results)
   %Print Sleep Staging Results 
    print_results(Sleep,Yhat_Results,states,print_figures,print_folder,display_flag);
end

%% Print RBD Detection Results
if (view_results)
   %Print Comparison of RBD Detection (Annotated)
   rbd_detect_name0 = 'Atonia Index (Annotated)';
   rbd_detect_name5 = 'STREAM (Annotated)';
   rbd_detect_name6 = 'QMA (Annotated)';

   tablename = 'Summary_RBD_Detection_Annotated';  
   print_rbd_detection_results_all(results_f_ai,results_f_stream,results_f_qma,squeeze(results_f(1,:,:)),squeeze(results_f(2,:,:)),squeeze(results_f(3,:,:)),squeeze(results_f(4,:,:)),rbd_detect_name0,RBD_Detection_Feats.labels{1},RBD_Detection_Feats.labels{2},RBD_Detection_Feats.labels{3},RBD_Detection_Feats.labels{4},rbd_detect_name5,rbd_detect_name6,tablename,print_figures,print_folder); 

   %Print Comparison of RBD Detection (Automated)
   rbd_detect_name0 = 'Atonia Index (Automated)';
   rbd_detect_name5 = 'STREAM (Automated)';
   rbd_detect_name6 = 'QMA (Automated)';   
   tablename = 'Summary_RBD_Detection_Automated';  
   print_rbd_detection_results_all(results_f_ai_auto,results_f_stream_auto,results_f_qma_auto,squeeze(results_f_auto(1,:,:)),squeeze(results_f_auto(2,:,:)),squeeze(results_f_auto(3,:,:)),squeeze(results_f_auto(4,:,:)),rbd_detect_name0,RBD_Detection_Feats.labels{1},RBD_Detection_Feats.labels{2},RBD_Detection_Feats.labels{3},RBD_Detection_Feats.labels{4},rbd_detect_name5,rbd_detect_name6,tablename,print_figures,print_folder); 

   %Compare RBD Detection (annotated)
    label_name = 'Annotated';
    compare_rbd_detection_results(EMG_Metric,RBD_Yhat_CV,label_name,print_figures,print_folder,display_flag);
    label_name = 'Automated';
    compare_rbd_detection_results(EMG_Auto_Metric,RBD_Auto_Yhat_CV,label_name,print_figures,print_folder,display_flag);
end


%% Print Feature Importance Results
if (view_results)
   %RBD Importance (Gini)
   important_mat = [];
    for t = 1:folds
        important_mat(:,:,t) = [rbd_rf_importance_all{3,t}];
    end
    order_idx = size(important_mat ,2); %Mean Decrease in Gini 
    titlename = 'Ranked Feature Importance (ECG)';
    xname = 'Mean Permuted Prediction Error (Importance)';
    print_feature_importance(important_mat,order_idx,EMG_Table_Names,RBD_Detection_Feats.feats{3},titlename,xname,print_figures,print_folder);       

    important_mat = [];
    for t = 1:folds
        important_mat(:,:,t) = [rbd_rf_importance_all{2,t}];
    end    
    titlename = 'Ranked Feature Importance (EMG)';
    order_idx = size(important_mat,2); %Mean Decrease in Gini 
    print_feature_importance(important_mat,order_idx,EMG_Table_Names,RBD_Detection_Feats.feats{2},titlename,xname,print_figures,print_folder);       

    important_mat = [];
    for t = 1:folds
        important_mat(:,:,t) = [rbd_rf_importance_all{4,t}];
    end    
    titlename = 'Ranked Feature Importance (All)';
    order_idx = size(important_mat,2); %Mean Decrease in Gini 
    print_feature_importance(important_mat,order_idx,EMG_Table_Names,RBD_Detection_Feats.feats{4},titlename,xname,print_figures,print_folder);       
    

    titlename = 'Ranked Feature Importance (Sleep Staging)';
    order_idx = size(importance_Results,2); %Mean Decrease in Gini 
    print_feature_importance(importance_Results,order_idx,Sleep_table_Pre.Properties.VariableNames,SS_Features,titlename,xname,print_figures,print_folder);       
    
end

%% Print Annotated Vs Automatic RBD Metrics
if (view_results)
    print_annotated_vs_auto(EMG_Table_Names,RBD_Detection_Feats.feats{3},EMG_Metric,EMG_Auto_Metric,print_figures,print_folder);
    print_annotated_vs_auto(EMG_Table_Names,RBD_Detection_Feats.feats{2},EMG_Metric,EMG_Auto_Metric,print_figures,print_folder);

end

%% Print Confusion Matrices/Hypnograms
if (view_results)
    All_Confusion = print_confusion_mats(Sleep,Sleep_Struct,Yhat_Results,print_figures,print_folder);
end

%%
if (save_data),save(strcat(print_folder,'\',Save_Data_Name,'.mat'),'Sleep','Sleep_table','Sleep_Struct',...
        'Sleep_names','Yhat_Results','indices','folds',...
        'votes_Results','indices','folds',...
        'importance_Results','SS_Features','rbd_rf_importance_all',...
        'EMG_Auto_Metric','EMG_Metric','EMG_Table_Names',...
        'RBD_Yhat_CV','RBD_Auto_Yhat_CV',...
        'results_f_ai','results_f_ai_auto','results_f_rem_auto','All_Confusion','posterior_struct','all_cost_opt');
end


end
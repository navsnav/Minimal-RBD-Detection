% Main code to run algorithms to perform feature extraction ->
% train/test cross-fold evaluation for automated sleep staging while also assessing RBD
% Detection

%% Add paths
slashchar = char('/'*isunix + '\'*(~isunix));
main_dir = strrep(which(mfilename),[mfilename '.m'],'');

addpath(genpath([main_dir, 'libs', slashchar])) % add external libraries folder to path
addpath(genpath([main_dir, 'subfunctions', slashchar])) % add subfunctions folder to path
addpath(genpath([main_dir, 'dataprep', slashchar])) % add data preparation folder to path
addpath(genpath([main_dir, 'models', slashchar])) % add classifiers folder to path


%% Attain PSG Signals
% There are several options to get PSG signals
% (A) A folder containing all edf files and annotations
% (B) Download files eg using CAP sleep database
% (C) A folder containing all 'prepared' mat files of all PSG signals
% (D) Load Features matrix saved from ExtractFeatures

%% (A) Extract PSG Signals  - Use this section if you have a folder of edf files with annotations
% current_dir = pwd;
% data_folder = 'I:\Data\CAP Sleep Database'; %Choose file location
% outputfolder = [current_dir,'\data'];
% prepare_capslpdb(data_folder,outputfolder);
% 
% data_folder = 'I:\Data\Navin PSG recordings\Both Nights'; %Choose file location
% outputfolder = [current_dir,'\data\'];
% prepare_JR_data(data_folder,outputfolder);
% 
% 
% data_folder = 'I:\Data\MASS\SS1\'; %Choose file location
% outputfolder = [current_dir,'\data'];
% prepare_massdb(data_folder,outputfolder,0);

%% (B) Download PSG Signals  - Use this section to download example edf files and annotations as a test
outputfolder = [pwd,'\data\'];
% The following data will be downloaded from the CAPS database from
% physionet
list_of_files = {
    'n1';
    'n2';
    'n3';
    'n10';
    'n5'
    'rbd1';
    'rbd2';
    'rbd3';
    'rbd4';
    'rbd5'};

download_CAP_EDF_Annotations(outputfolder,list_of_files);
%Prepare mat files with PSG signals and annotations
std_flg = 0;
prepare_capslpdb(outputfolder,outputfolder,std_flg);

%% (C) Extract PSG Signals - Use this section if you have a dataset of mat files with hypnogram datasets
current_dir = pwd;
data_folder = [main_dir,'data\'];
output_filename = 'Features_CAP_Test.mat';

signals_for_processing = {'EEG','EOG','EMG','EEG-EOG','ECG'};

disp(['Extracting Features:',signals_for_processing]);
% % % Generate Features

std_flag = 0;
filt_flag = 1;
[Sleep, Sleep_Struct, Sleep_table] = ExtractFeatures_mat(data_folder,signals_for_processing,std_flag,filt_flag);

output_folder = [pwd,'\features'];
% Create a destination directory if it doesn't exist
if exist(output_folder, 'dir') ~= 7
    fprintf('WARNING: Output directory does not exist. Creating new directory ...\n\n');
    mkdir(output_folder);
end
cd(output_folder);

save(output_filename,'Sleep','Sleep_Struct','Sleep_table', '-v7.3');

disp('Feature Extraction Complete and Saved');
cd(current_dir);

%% (D) Load Features matrix saved from ExtractFeatures
current_dir = pwd;
data_folder = [main_dir, 'data', slashchar, 'features', slashchar];
cd(data_folder);
filename = 'Features_CAP_Test.mat';
load(filename);
cd(current_dir);
%% Balance RBD & HC Cohort if needed

[patients,ia,ic] = unique(Sleep_table.SubjectIndex);

% Remove dupicates, and participants with PD/MSA
rmv_id = [];

rmv_idx = ismember(Sleep_table.SubjectIndex,patients(rmv_id));

Sleep(rmv_idx,:) = [];
Sleep_table(rmv_idx,:) = [];
subject_names = fieldnames(Sleep_Struct);
Sleep_Struct_Old = Sleep_Struct;
Sleep_Struct = rmfield(Sleep_Struct,subject_names(rmv_id));

%% Parameters for generating results
outfilename = 'RBD_Detection_Test_Results'; %Filename/Folder to be created
view_results = 1; %Produce Graphs/figures (set to 1 to observe figures)
print_figures= 1; %Save Graphs/figures (set to 1 to save figures)
print_folder = [data_folder,'Graphs_', outfilename, ];
display_flag = 1; %Diplay results on command window
save_data = 1; %Save Data
cd(main_dir);

%% Preprocess Data

SS_Features = [11:122,123:151,152:167,325,326]; % All new EEG/EOG/EMG Features used for sleep staging (HR only 150s)
% SS_Features = [123:145,146:160]; % All new EOG/EMG Features used for sleep staging 


EMG_est_feats = [3,7,8,9]; %AI, MAD_Dur, MAD_Per, Stream
EMG_feats = [EMG_est_feats,10,11,14,15,32,33]; % For 3 Stage: AI ratios + N2% + N3% + Fractal Exponenet Ratios (REM:N2/N3)
ECG_feats = [14,131,127,153,156,160,182,186,188,136,138]; %For 3 stage: ECG 150s features without EMG
All_feats = [EMG_feats,127,131,127,153,156,160,182,136,138,186,188]; %ECG features without EMG

RBD_Detection_Feats = struct;
RBD_Detection_Feats.feats = {EMG_est_feats,EMG_feats,ECG_feats,All_feats};
RBD_Detection_Feats.labels = {'Established Metrics','New Features','ECG Features','All Features'}; 

numStates = 3;
% Preprocess Sleep Staging Features
disp('Precprocessing Features...');
[Sleep_table_Pre] = RBD_RF_Preprocess(Sleep_table,[-1,0,5],SS_Features,numStates);

disp('Precprocessing Complete.');

% Random Forest paramters
n_trees = 50;  %Paper used 500 trees (done to save time/space)

%% Cross Fold Indexing
Sleep = table2array(Sleep_table_Pre);
[patients,ia,ic] = unique(Sleep_table_Pre.SubjectIndex);
rbd_group = Sleep(ia,6)==5;

folds = ceil(log2(length(rbd_group))/5)*5; %find appropriate number of folds
indices = crossvalind('Kfold',rbd_group, folds);

 %% SS & RBD Detection
% Apply cross fold validation for automated sleep staging followed by RBD
% detection using established metrics and new metrics
prior_mat = [];
% prior_mat = 1;

disp(['Initiating ',num2str(folds),' fold cross validation with ',num2str(n_trees),' trees.']);
[RBD_Yhat_Results] = SS_RBD_Detection(Sleep_table_Pre,Sleep_Struct,rbd_group,indices,folds,SS_Features,RBD_Detection_Feats,n_trees,view_results,print_figures,print_folder,save_data,outfilename,display_flag,prior_mat);

%% SS Detection Only
% Apply cross fold validation for automated sleep staging followed by RBD
% detection using established metrics and new metrics
prior_mat = 0;
if isempty(prior_mat)
    disp(['Initiating ',num2str(folds),' fold cross validation with ',num2str(n_trees),' trees.']);
else    
    disp(['Initiating ',num2str(folds),' fold cross validation with ',num2str(n_trees),' trees and priors.']);
end

[Auto_SS_Results,RBD_New_Results,EMG_Est_Results,EMG_Auto_New_Results,EMG_Auto_Est_Results,All_Confusion] = SS_Detection(Sleep_table_Pre,Sleep_Struct,rbd_group,indices,folds,SS_Features,EMG_est_feats,EMG_feats,ECG_feats,n_trees,view_results,print_figures,print_folder,save_data,outfilename,display_flag,prior_mat);

%% RBD Detection Only

% prior_mat = [];
% [RBD_Yhat_Results] = RBD_Detection(Sleep_table_Pre,Sleep_Struct,rbd_group,indices,folds,SS_Features,RBD_Detection_Feats,n_trees,view_results,print_figures,print_folder,save_data,outfilename,display_flag,prior_mat);



%% Train Model Only

% [ss_rf,ss_rf_importance,posterior] = Train_SleepStaging_RF3(n_trees,Sleep_table_Pre,SS_Features,Sleep_table_Pre.AnnotatedSleepStage,[]);
% model_folder = 'C:\Users\scro2778\Documents\GitHub\RBD-Sleep-Detection\models'
% model_name = 'ThesisModel_27_5_2020';
% save(strcat(model_folder,'\',model_name,'.mat'),'ss_rf','SS_Features','-v7.3'); 
%  
 
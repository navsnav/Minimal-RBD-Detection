function prepare_new_JR_data_raw(filename,outputfolder,std_flag,raw_flag)
% This function prepares data from raw edf/txt files of PSG recordings.
% Specifically preparing data for each subject that provides neccessary
% signals (eg. EEG, EOG, EMG etc). Data is saved in mat format within
% folder
%
% Inputs:
%  filename    - location folder of PSG recordings
%  outputfolder   - folder to save prepared PSG signals in mat format for
%                   each subject
%
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
% Specify details of data to be saved after processing and extracting
% features.

%%
addpath(pwd)
cd(filename)
% cd('/data/datasets/Navin_Data/JR Data/')
% cd('I:\Data\Navin PSG recordings\P518')
% genderage = readtable('gender-age-8-1-2018.csv');
% genderage = readtable('gender-age.csv');
% genderage = readtable('gender-age-12-4-2018.csv');
cd([pwd,'\..\..\..\']);
genderage = readtable('SleepProjectRecordings.xlsx');
cd(filename);
CLASSES = {'W','S1','S2','S3','R'};

%% Parameters

%% Get the raw "xls" and "edf" data -- if they are not identical we have a problem!
flsedf = dir('*.edf');
flsedf = arrayfun(@(x) x.name,flsedf,'UniformOutput',false);
hypfldcsv = dir('*.csv');
hypfldcsv = arrayfun(@(x) x.name,hypfldcsv,'UniformOutput',false);

%% Main loop

for f = 1:length(flsedf)
    disp(['Loading ' flsedf{f}])
    
    [hdr, record] = edfread(flsedf{f});
    total_time = length(record)/max(hdr.frequency);
    if ~range(hdr.frequency(~isnan(hdr.frequency))) == 0
        zero_flag = 1;        
    else 
        zero_flag = 0;
    end    
    
    subject = flsedf{f}(1:end-4);
    
    patid = upper(flsedf{f}(1:7));
    idxrow = ~cellfun(@isempty, strfind(upper(genderage.ID),patid));
    patinfo.gender = genderage.Gender{find(idxrow,1,'first')}(1);
    patinfo.age = genderage.Age(idxrow);
    patinfo.fs_orig = nan(2,9);
    patinfo.fs = -1;
    patinfo.classes = CLASSES;
    patinfo.chlabels = {'EEG','EOG','EMG','ECG','EMG_l','EMG_r','MIC','CHEST','FLOW'};
    
    % Figuring out start point
    startdate = eval(['genderage.Start_N{idxrow}']);
    psgstarttime = eval(['genderage.PSGStart_N{idxrow}']);
    recstart = datetime([hdr.startdate,'.',hdr.starttime],'InputFormat','dd.MM.yy.HH.mm.ss');
    
    annotstart = datetime(strrep(startdate,',','.'),'InputFormat','yyyy.MM.dd.HH.mm.ss');
    try 
        recstart2 =  datetime(strrep(psgstarttime,',','.'),'InputFormat','yyyy.MM.dd.HH.mm.ss');
    catch
        recstart2 = datetime([hdr.startdate,'.',hdr.starttime],'InputFormat','dd.MM.yy.HH.mm.ss');
    end    
    %     if (str2double(hdr.starttime(1:2)) < 12)&&(str2double(hdr.starttime(1:2)) > 3) % recoding 514 messes up am and pm times
    %         recstart = recstart + hours(12);
    %     end
    skip_s = seconds(annotstart-recstart);
%     skip = skip_s*patinfo.fs;
    if skip_s < 0
        skip_s = seconds((annotstart+1)-recstart);
%         skip = seconds((annotstart+1)-recstart)*patinfo.fs;
    end
    %% EEG
    idxeeg = strmatch('C4',hdr.label);
%     if isempty(idxeeg)
%         idxeeg = strmatch('C3',hdr.label);        
%     end
    idxeeg = [idxeeg,strmatch('A1',hdr.label)];
    if isempty(idxeeg)
        idxeeg = strmatch('Fp1',hdr.label);
    end      
    fseeg = hdr.frequency(idxeeg);
    if length(idxeeg) > 1
        eeg1 = record(idxeeg(1),:);
        eeg2 = record(idxeeg(2),:);
        eeg = (eeg1 - eeg2);
        if fseeg(1) ~= fseeg(2)
            error('Different fs for EEG!?')
        end
    else
        eeg = record(idxeeg,:);
    end
    if any(isnan(eeg))
        error('Whats up here!?')
    end
    
    if zero_flag
       eeg = eeg(1:total_time*fseeg); 
    end
    if ~raw_flag
       eeg(1:skip_s*fseeg(1)) = [];
    end    
    patinfo.fs_orig(:,1) = fseeg;
    
    
    %% EOG
    idxeog = strmatch('ROC',hdr.label);
    idxeog = [idxeog, strmatch('LOC',hdr.label)];
    if isempty(idxeog)
        idxeog = strmatch('EOGl',hdr.label);
        idxeog = [idxeog, strmatch('EOGr',hdr.label)];             
    end
    if isempty(idxeog)    
        idxeog = strmatch('EOG',hdr.label);
        idxeog = [idxeog, strmatch('eog',hdr.label)];                
    end    
    fseog = hdr.frequency(idxeog);
    eog1 = record(idxeog(1),:);
    if length(idxeog) > 1
        eog2 = record(idxeog(2),:);
        eog = (eog1 - eog2);
        if fseog(1) ~= fseog(2)
            error('Different fs for EOG!?')
        end
    else 
        eog = eog1;
    end
    if zero_flag
       eog = eog(1:total_time*fseog); 
    end    
    if ~raw_flag
       eog(1:skip_s*fseog(1)) = [];
    end
    patinfo.fs_orig(:,2) = fseog;
    clear idxeeg idxeog
    
    %% EMG
    idxemg = strmatch('CHIN1',hdr.label);
    idxemg = [idxemg, strmatch('CHIN2',hdr.label)];
    if isempty(idxemg)
        idxemg = strmatch('CHIN',hdr.label);
        idxemg = [idxemg, strmatch('chin',hdr.label)];        
    end
    if isempty(idxemg)
        idxemg = strmatch('SEMG',hdr.label);
    end 
    if isempty(idxemg)
        idxemg = strmatch('PLM1',hdr.label);
    end     
    fsemg = hdr.frequency(idxemg);
    emg1 = record(idxemg(1),:);    
    if length(idxemg) > 1
        emg2 = record(idxemg(2),:);
        emg = (emg1 - emg2);
        if fsemg(1) ~= fsemg(2)
            error('Different fs for emg!?')
        end
    else 
       emg =  emg1;
    end
    if any(isnan(emg))
        error('Whats up here!?')
    end
    patinfo.fs_orig(:,3) = fsemg;
    if zero_flag
       emg = emg(1:total_time*fsemg); 
    end
    if ~raw_flag
       emg(1:skip_s*fsemg(1)) = [];
    end        
    clear  idxeeg idxemg
    %% ECG
    idxecg = strmatch('ECGL',hdr.label);
    idxecg = [idxecg, strmatch('ECGR',hdr.label)];
    if isempty(idxecg)
        idxecg = strmatch('ECG1',hdr.label);
        idxecg = [idxecg, strmatch('ECG2',hdr.label)];        
    end
    if isempty(idxecg)
        idxecg = strmatch('ECGl',hdr.label);
        idxecg = [idxecg, strmatch('ECGr',hdr.label)];        
    end
    if isempty(idxecg)
        idxecg = strmatch('X1',hdr.label);
        idxecg = [idxecg, strmatch('X2',hdr.label)];        
    end    
    if isempty(idxecg)
        disp('No ECG daata!?');
        ecg1 = [];
        ecg2 = [];
        ecg = [];
        fsecg = [NaN;NaN];
    else        
        fsecg = hdr.frequency(idxecg);
        if length(idxecg) > 1
            ecg1 = record(idxecg(1),:);
            ecg2 = record(idxecg(2),:);
            ecg = (ecg1 - ecg2);
            if fsecg(1) ~= fsecg(2)
                error('Different fs for ecg!?')
            end
        else 
            ecg = record(idxecg(1),:);
        end
        if any(isnan(ecg))
            error('Whats up here!?')
        end
        if zero_flag
           ecg = ecg(1:total_time*fsecg); 
        end
        if ~raw_flag
           ecg(1:skip_s*fsecg(1)) = [];
        end          
    end

    patinfo.fs_orig(:,4) = fsecg;    
    %% Leg EMGs
    idxemg = strmatch('LAT1',hdr.label);
    idxemg = [idxemg, strmatch('LAT2',hdr.label)];
    if isempty(idxemg)
        idxemg = strmatch('LLEMG',hdr.label);
        idxemg = [idxemg, strmatch('EMGref',hdr.label)];        
    end
    if isempty(idxemg)
        idxemg = strmatch('LTA',hdr.label);
        idxemg = [idxemg, strmatch('lta',hdr.label)];        
    end   
    if isempty(idxemg)
        idxemg = strmatch('PLM2',hdr.label);
    end       
    fsemg_l = hdr.frequency(idxemg);
    emg1 = record(idxemg(1),:);
    if length(idxemg) > 1
        emg2 = record(idxemg(2),:);
        emg_l = (emg1 - emg2);
        if fsemg_l(1) ~= fsemg_l(2)
            error('Different fs for emg!?')
        end
    else
       emg_l = emg1;
    end
    if any(isnan(emg))
        error('Whats up here!?')
    end
    patinfo.fs_orig(:,5) = fsemg_l;
    if zero_flag
       emg_l = emg_l(1:total_time*fsemg_l); 
    end 
    if ~raw_flag
       emg_l(1:skip_s*fsemg_l(1)) = [];
    end     
    %%
    idxemg = strmatch('RAT1',hdr.label);
    idxemg = [idxemg, strmatch('RAT2',hdr.label)];
    if isempty(idxemg)
        idxemg = strmatch('RLEMG',hdr.label);
        idxemg = [idxemg, strmatch('EMGREF',hdr.label)];        
    end
    if isempty(idxemg)
        idxemg = strmatch('RTA',hdr.label);
        idxemg = [idxemg, strmatch('rta',hdr.label)];        
    end       
    if isempty(idxemg)
        disp('No EMGr daata!?');
        emg1 = [];
        emg2 = [];
        emg_r = [];
        fsemg_r = [NaN;NaN];
    else     
        fsemg_r = hdr.frequency(idxemg);
        emg1 = record(idxemg(1),:);
        if length(idxemg) > 1
            emg2 = record(idxemg(2),:);
            emg_r = (emg1 - emg2);
            if fsemg_r(1) ~= fsemg_r(2)
                error('Different fs for emg!?')
            end
        else 
            emg_r = emg1;
        end
        if any(isnan(emg))
            error('Whats up here!?')
        end
        if zero_flag
           emg_r = emg_r(1:total_time*fsemg_r); 
        end
        if ~raw_flag
           emg_r(1:skip_s*fsemg_r(1)) = [];
        end           
    end
    patinfo.fs_orig(:,6) = fsemg_r;

    %% MIC
    idxmic = strmatch('SNORE',hdr.label);
    idxmic2 = strmatch('snore',hdr.label);


    if isempty(idxmic)
        disp('No MIC daata!?');
        mic = [];
        fsmic = [];
    else
        fsmic = hdr.frequency(idxmic);
        mic = record(idxmic(1),:);
        if isempty(idxmic2)
            patinfo.fs_orig(:,7) = fsmic;                
        else
            mic = mic -  record(idxmic2(1),:);
            fsmic(2) = hdr.frequency(idxmic2);
            if fsmic(1) ~= fsmic(2)
                error('Different fs for mic!?')
            end            
            patinfo.fs_orig(:,7) = fsmic;               
        end
        if zero_flag
            mic = mic(1:total_time*fsmic); 
        end    
        if ~raw_flag
           mic(1:skip_s*fsmic(1)) = [];
        end            
    end

    
    %% CHEST
    idxchest = strmatch('CHEST',hdr.label);
    idxchest2 = strmatch('Chest',hdr.label);
    if isempty(idxchest)
        idxchest = strmatch('Abdomen',hdr.label);        
    end
    if isempty(idxchest)
        disp('No Chest Data!?');
        chest = [];      
        fschest = [];
    else
        fschest = hdr.frequency(idxchest);
        chest = record(idxchest(1),:);       
        if isempty(idxchest2)
            patinfo.fs_orig(:,8) = fschest;    
        else
            chest = chest -  record(idxchest2(1),:);
            fschest(2) = hdr.frequency(idxchest2);
            if fschest(1) ~= fschest(2)
                error('Different fs for chest!?')
            end            
            patinfo.fs_orig(:,8) = fschest;             
        end
        if zero_flag
           chest = chest(1:total_time*fschest); 
        end      
        if ~raw_flag
           chest(1:skip_s*fschest(1)) = [];
        end          
    end

    
    %% FLOW
    idxflow = strmatch('FLOW',hdr.label);
    idxflow2 = strmatch('flow',hdr.label);
    
    if isempty(idxflow)
        disp('No Flow Data!?');
        flow = [];   
        fsflow = [];
    else
        fsflow = hdr.frequency(idxflow);
        flow = record(idxflow(1),:);         
        if isempty(idxflow2)
            patinfo.fs_orig(:,9) = fsflow;                
        else
            flow = flow -  record(idxflow2(1),:);
            fsflow(2) = hdr.frequency(idxflow2);
            if fsflow(1) ~= fsflow(2)
                error('Different fs for flow!?')
            end            
            patinfo.fs_orig(:,9) = fsflow;               
        end
        if zero_flag
            flow = flow(1:total_time*fsflow); 
        end    
        if ~raw_flag
           flow(1:skip_s*fsflow(1)) = [];
        end           
    end

    %% Filtering
    fss = round(patinfo.fs_orig(1,:));
    record(:,1:skip_s*max(hdr.frequency)) = [];
    

 
%%  Load annotations

    [epoch, labels]  = importcsv_new(sprintf('%s epochs.csv',string(patid)));
    skip_end_s = epoch(end)*30 - length(record)/max(hdr.frequency);
    
    clear record idxeeg idxemg    
    
    if skip_end_s > 0 %annotated hyp is larger than psg recording (cut-off signal and hyp)
       skip_end = skip_end_s*patinfo.fs; 
%        signal(:,end-skip_end:end) = [];
       epoch(end-floor(skip_end_s/30):end) = [];
       labels(end-floor(skip_end_s/30):end,:) = [];          

    end
    
    % Standardizing signals (does not help and not neccessary)
    if std_flag == 1
        for s=1:size(signal,1)
            signal(s,:) = signal(s,:) - nanmean(signal(s,:));
            signal(s,:) = signal(s,:)./nanstd(signal(s,:));
        end
    end
    

%     data = zeros(length(epoch),30*patinfo.fs,size(signal,1));
    data = {};
    tstamp = zeros(length(epoch),9);
    for seg = 1:length(epoch)
        if ~isempty(eeg)
            tstamp(seg,1) = (epoch(seg)-1)*30*fseeg(1) + 1; %not sure about +1            
            data_eeg(seg,:,:) = eeg(:,tstamp(seg,1):tstamp(seg,1)+30*fseeg(1)-1)';            
        end
        if ~isempty(eog)
            tstamp(seg,2) = (epoch(seg)-1)*30*fseog(1) + 1; %not sure about +1            
            data_eog(seg,:,:) = eog(:,tstamp(seg,2):tstamp(seg,2)+30*fseog(1)-1)';            
        end 
        if ~isempty(emg)
            tstamp(seg,3) = (epoch(seg)-1)*30*fsemg(1) + 1; %not sure about +1            
            data_emg(seg,:,:) = emg(:,tstamp(seg,3):tstamp(seg,3)+30*fsemg(1)-1)';            
        end  
        if ~isempty(ecg)
            tstamp(seg,4) = (epoch(seg)-1)*30*fsecg(1) + 1; %not sure about +1                        
            data_ecg(seg,:,:) = ecg(:,tstamp(seg,4):tstamp(seg,4)+30*fsecg(1)-1)';            
        end   
        if ~isempty(emg_r)
            tstamp(seg,5) = (epoch(seg)-1)*30*fsemg_r(1) + 1; %not sure about +1                        
            data_emg_r(seg,:,:) = emg_r(:,tstamp(seg,5):tstamp(seg,5)+30*fsemg_r(1)-1)';            
        end
        if ~isempty(emg_l)
            tstamp(seg,6) = (epoch(seg)-1)*30*fsemg_l(1) + 1; %not sure about +1                        
            data_emg_l(seg,:,:) = emg_l(:,tstamp(seg,6):tstamp(seg,6)+30*fsemg_l(1)-1)';            
        end    
        if ~isempty(mic)
            tstamp(seg,7) = (epoch(seg)-1)*30*fsmic(1) + 1; %not sure about +1                        
            data_mic(seg,:,:) = mic(:,tstamp(seg,7):tstamp(seg,7)+30*fsmic(1)-1)';            
        end    
        if ~isempty(chest)
            tstamp(seg,8) = (epoch(seg)-1)*30*fschest(1) + 1; %not sure about +1                        
            data_chest(seg,:,:) = chest(:,tstamp(seg,8):tstamp(seg,8)+30*fschest(1)-1)';            
        end 
        if ~isempty(flow)
            tstamp(seg,9) = (epoch(seg)-1)*30*fsflow(1) + 1; %not sure about +1                        
            data_flow(seg,:,:) = flow(:,tstamp(seg,9):tstamp(seg,9)+30*fsflow(1)-1)';            
        end         
%         data(seg,:,:) = signal(:,tstamp(seg):tstamp(seg)+30*patinfo.fs-1)';
    end
       
    if ~isempty(eeg)
        data.eeg = data_eeg;            
    end
    if ~isempty(eog)
        data.eog = data_eog;            
    end 
    if ~isempty(emg)
        data.emg = data_emg;            
    end  
    if ~isempty(ecg)
        data.ecg = data_ecg;            
    end   
    if ~isempty(emg_r)
        data.emg_r = data_emg_r;            
    end
    if ~isempty(emg_l)
        data.emg_l = data_emg_l;            
    end    
    if ~isempty(mic)
        data.mic = data_mic;            
    end    
    if ~isempty(chest)
        data.chest = data_chest;            
    end 
    if ~isempty(flow)
        data.flow = data_flow;            
    end         
    
    clear fseeg fseog fsemg fss fsflow fschest    
    clear eeg eog emg skip fseog fsemg fseeg fss idxrow eeg1 eeg2
    clear emg1 emg2 eog1 eog2       
    clear data_eeg data_eog data_emg data_emg_r data_emg_l data_mic data_flow data_chest data_ecg
    
    if any(isempty(data))
        error('NaNs')
    end
    idx = ~any(labels,2);  %rows
    labels(idx,:) = [];
    data(idx,:,:) = [];
    tstamp(idx,:) = [];
    
%     if size(data,1) ~= size(labels,1)
%         error('Different length for recording..')
%     end
    epoch = round(tstamp/patinfo.fs);
    

    %     savefig(['plot_' flsedf{f} '.fig'])
    %     drawnow
    %     close
    if raw_flag
        save([outputfolder,subject,'_raw.mat'],'data','labels','epoch','patinfo')
    else
            save([outputfolder,subject,'.mat'],'data','labels','epoch','patinfo')    
    end
    clear data labels patinfo epoch signal
end
end

% clear all
%
% fls = dir('*.mat');
% fls = arrayfun(@(x) x.name,fls,'UniformOutput',false);
%
% for i = 1:2:length(fls)
%     m1 = load(fls{i});
%     m2 = load(fls{i+1});
%     sprintf('%s with %s', fls{i}, fls{i+1})
%     data1 = m1.data;
%     data2 = m2.data;
%     labels1 = m1.labels;
%     labels2 = m2.labels;
%     save([fls{i}(1:end-11) '.mat'],'data1','data2','labels1','labels2');
%     clear data1 data2 m1 m2 labels1 labels2
% end



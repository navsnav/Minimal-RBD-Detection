function prepare_new_JR_data(filename,outputfolder,std_flag)
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
    patinfo.fs = 200;
    patinfo.classes = CLASSES;
    patinfo.chlabels = {'EEG','EOG','EMG','ECG','EMG_l','EMG_r','MIC','CHEST','FLOW'};
    
    % Figuring out start point
    startdate = eval(['genderage.Start_N{idxrow}']);
    psgstarttime = eval(['genderage.PSGStart_N{idxrow}']);
    recstart = datetime([hdr.startdate,'.',hdr.starttime],'InputFormat','dd.MM.yy.HH.mm.ss');
    
    annotstart = datetime(strrep(startdate,',','.'),'InputFormat','yyyy.MM.dd.HH.mm.ss');
%     try 
%         recstart =  datetime(strrep(psgstarttime,',','.'),'InputFormat','yyyy.MM.dd.HH.mm.ss');
%     catch
%         recstart = datetime([hdr.startdate,'.',hdr.starttime],'InputFormat','dd.MM.yy.HH.mm.ss');
%     end    
    %     if (str2double(hdr.starttime(1:2)) < 12)&&(str2double(hdr.starttime(1:2)) > 3) % recoding 514 messes up am and pm times
    %         recstart = recstart + hours(12);
    %     end
    skip = seconds(annotstart-recstart)*patinfo.fs;
    if skip < 0
        skip = seconds((annotstart+1)-recstart)*patinfo.fs;
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
    patinfo.fs_orig(:,1) = fseeg;
    
    
    %% EOG
    idxeog = strmatch('ROC',hdr.label);
    idxeog = [idxeog, strmatch('LOC',hdr.label)];
    if isempty(idxeog)
        idxeog = strmatch('EOGl',hdr.label);
        idxeog = [idxeog, strmatch('EOGr',hdr.label)];             
    else
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
    end
    if zero_flag
       ecg = ecg(1:total_time*fsecg); 
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
    fsemg = hdr.frequency(idxemg);
    emg1 = record(idxemg(1),:);
    if length(idxemg) > 1
        emg2 = record(idxemg(2),:);
        emg_l = (emg1 - emg2);
        if fsemg(1) ~= fsemg(2)
            error('Different fs for emg!?')
        end
    else
       emg_l = emg1;
    end
    if any(isnan(emg))
        error('Whats up here!?')
    end
    patinfo.fs_orig(:,5) = fsemg;
    if zero_flag
       emg_l = emg_l(1:total_time*fsemg); 
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
        fsemg = [NaN;NaN];
    else     
        fsemg = hdr.frequency(idxemg);
        emg1 = record(idxemg(1),:);
        if length(idxemg) > 1
            emg2 = record(idxemg(2),:);
            emg_r = (emg1 - emg2);
            if fsemg(1) ~= fsemg(2)
                error('Different fs for emg!?')
            end
        else 
            emg_r = emg1;
        end
        if any(isnan(emg))
            error('Whats up here!?')
        end
        if zero_flag
           emg_r = emg_r(1:total_time*fsemg); 
        end
    end
    patinfo.fs_orig(:,6) = fsemg;

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
    end

    %% Filtering
    fss = round(patinfo.fs_orig(1,:));
    

    
    eeg = resample(eeg,patinfo.fs,fss(1));
    eog = resample(eog,patinfo.fs,fss(2));
    emg = resample(emg,patinfo.fs,fss(3));
    emg_l = resample(emg_l,patinfo.fs,fss(5));
    
    
    if isempty(emg_r)
        emg_r = zeros(size(eeg));
    else
        emg_r = resample(emg_r,patinfo.fs,fss(6));
    end
    if isempty(mic)
        mic = zeros(size(eeg));
    else
        mic = resample(mic,patinfo.fs,fss(7));        
    end
    if isempty(chest)
        chest = zeros(size(eeg));        
    else
        chest = resample(chest,patinfo.fs,fss(8));
    end
    if isempty(flow)
        flow = zeros(size(eeg));            
    else
        flow = resample(flow,patinfo.fs,fss(9));
    end
    if isempty(ecg)
        ecg = zeros(size(eeg));            
    else
        ecg = resample(ecg,patinfo.fs,fss(4));
    end    
    clear fseeg fseog fsemg fss fsflow fschest
    clear record idxeeg idxemg    

    
    Nfir = 500;
    b_band = fir1(Nfir,[0.3 40].*2/patinfo.fs,'bandpass'); % bandpass
    eeg = filtfilt(b_band,1,eeg);
    b_band = fir1(Nfir,[0.3 40].*2/patinfo.fs,'bandpass'); % bandpass
    eog = filtfilt(b_band,1,eog);    
    
%     ecg = filtfilt(b_band,1,ecg);      
    
    pwrline1 = 50; %Hz
    pwrline2 = 60; %Hz
    b_notch1 = fir1(Nfir,[(pwrline1-1) (pwrline1+1)].*2/patinfo.fs,'stop');
    b_notch2 = fir1(Nfir,[(pwrline2-1) (pwrline2+1)].*2/patinfo.fs,'stop');
%     b_band = fir1(Nfir,[10 100].*2/patinfo.fs,'bandpass'); % bandpass  
  b_band = fir1(Nfir,10.*2/patinfo.fs,'high'); % bandpass
    
%      subplot(2,1,1)
%     N = length(emg);
%     freq = linspace(-patinfo.fs/2, patinfo.fs/2, N+1); freq(end) = [];
%     ydft = fft(emg);
%     shiftSpectrum = abs(fftshift(ydft));
%     plot(freq,shiftSpectrum,'r');
%     xlabel ('Frequency (Hz)');
%      ylabel('Magnitude');
%     xlim([0,100])
    
    emg = filtfilt(b_notch1,1,emg);
    emg = filtfilt(b_notch2,1,emg);
    emg = filtfilt(b_band,1,emg);
    
    emg_l = filtfilt(b_notch1,1,emg_l);
    emg_l = filtfilt(b_notch2,1,emg_l);
    emg_l = filtfilt(b_band,1,emg_l);
    
    emg_r = filtfilt(b_notch1,1,emg_r);
    emg_r = filtfilt(b_notch2,1,emg_r);
    emg_r = filtfilt(b_band,1,emg_r);    
    
    b_band = fir1(Nfir,3.*2/patinfo.fs,'high'); % bandpass
    mic = filtfilt(b_notch1,1,mic);
    mic = filtfilt(b_notch2,1,mic);
    mic = filtfilt(b_band,1,mic);       

    b_band = fir1(Nfir,[0.1 0.5].*2/patinfo.fs,'bandpass'); % bandpass  
    flow = filtfilt(b_band,1,flow);       
    
    
%     subplot(2,1,2)
%      N = length(emg);
%     freq = linspace(-patinfo.fs/2, patinfo.fs/2, N+1); freq(end) = [];
%     ydft = fft(emg);
%     shiftSpectrum = abs(fftshift(ydft));
%     plot(freq,shiftSpectrum,'r');
%     xlabel ('Frequency (Hz)');
%     ylabel('Magnitude');
%     xlim([0,100])
%     title(flsedf{f})
%%  Resample and merge into matrix
    signal = [eeg;eog;emg;ecg;emg_r;emg_l;mic;chest;flow];
    signal(:,1:skip) = [];

    [epoch, labels]  = importcsv_new(sprintf('%s epochs.csv',string(patid)));
    skip_end_s = epoch(end)*30 - length(signal)/patinfo.fs;
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
    
    clear eeg eog emg skip fseog fsemg fseeg fss idxrow eeg1 eeg2
    clear emg1 emg2 eog1 eog2
    % Loading annotations
    % sanity plot
    %     tm = [0:30*patinfo.fs:length(signal)];
    %     plot(signal')
    %     hold on
    %     for i = 1:length(tm)
    %         line([tm(i) tm(i)],[-100,100],'Color','k')
    %     end
    %     classes = {'W' 'R' 'N1', 'N2', 'N3'};
    %     txt = cell(length(epoch),1);
    %     for l = 1:length(epoch)
    %         txt(l) = classes(labels(l,:));
    %     end
    %     text(tm(epoch),100*ones(size(epoch)),txt)
    % Putting all together
    data = zeros(length(epoch),30*patinfo.fs,size(signal,1));
    tstamp = zeros(length(epoch),1);
    for seg = 1:length(epoch)
        tstamp(seg) = (epoch(seg)-1)*30*patinfo.fs + 1; %not sure about +1
        data(seg,:,:) = signal(:,tstamp(seg):tstamp(seg)+30*patinfo.fs-1)';
    end
    
    if any(isnan(data))
        error('NaNs')
    end
    idx = ~any(labels,2);  %rows
    labels(idx,:) = [];
    data(idx,:,:) = [];
    tstamp(idx,:) = [];
    
    if size(data,1) ~= size(labels,1)
        error('Different length for recording..')
    end
    epoch = round(tstamp/patinfo.fs);
    
    
    %     savefig(['plot_' flsedf{f} '.fig'])
    %     drawnow
    %     close
    save([outputfolder,subject,'.mat'],'data','labels','epoch','patinfo')
    
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



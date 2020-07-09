function prepare_JR_data(filename,outputfolder,std_flg)
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
genderage = readtable('gender-age-8-1-2018.csv');
% genderage = readtable('gender-age.csv');
% genderage = readtable('gender-age-12-4-2018.csv');

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
    subject = flsedf{f}(1:end-4);
    
    patid = flsedf{f}(10:12);
    idxrow = ~cellfun(@isempty, strfind(genderage.ID,patid));
    patinfo.gender = genderage.Gender{idxrow}(1);
    patinfo.age = genderage.Age(idxrow);
    patinfo.fs_orig = nan(2,9);
    patinfo.fs = 200;
    patinfo.classes = CLASSES;
    patinfo.chlabels = {'EEG','EOG','EMG','ECG','EMG_l','EMG_r','MIC','CHEST','FLOW'};
    
    % Figuring out start point
    night=regexp(flsedf{f},'_Night(\d).edf','tokens');
    startdate = eval(['genderage.Start_N' char(night{1}) '{idxrow}']);
    try 
        psgstarttime = eval(['genderage.PSGStart_N' char(night{1}) '{idxrow}']);
        recstart =  datetime(strrep(psgstarttime,',','.'),'InputFormat','yyyy.MM.dd.HH.mm.ss');        
    catch 
        psgstarttime = [hdr.startdate,',',hdr.starttime];
        recstart =  datetime(strrep(psgstarttime,',','.'),'InputFormat','dd.MM.yy.HH.mm.ss');                
    end
    
    annotstart = datetime(strrep(startdate,',','.'),'InputFormat','yyyy.MM.dd.HH.mm.ss');

    
    
    %     if (str2double(hdr.starttime(1:2)) < 12)&&(str2double(hdr.starttime(1:2)) > 3) % recoding 514 messes up am and pm times
    %         recstart = recstart + hours(12);
    %     end
    skip = seconds(annotstart-recstart)*patinfo.fs;
    if skip < 0
        skip = seconds((annotstart+1)-recstart)*patinfo.fs;
    end
    %% EEG
    idxeeg = strmatch('C4',hdr.label);
    idxeeg = [idxeeg,strmatch('A1',hdr.label)];
    fseeg = hdr.frequency(idxeeg);
    eeg1 = record(idxeeg(1),:);
    eeg2 = record(idxeeg(2),:);
    eeg = (eeg1 - eeg2);
    if fseeg(1) ~= fseeg(2)
        error('Different fs for EEG!?')
    end
    if any(isnan(eeg))
        error('Whats up here!?')
    end
    patinfo.fs_orig(:,1) = fseeg;
    
    
    %% EOG
    idxeog = strmatch('ROC',hdr.label);
    idxeog = [idxeog, strmatch('LOC',hdr.label)];
    fseog = hdr.frequency(idxeog);
    eog1 = record(idxeog(1),:);
    eog2 = record(idxeog(2),:);
    eog = (eog1 - eog2);
    if fseog(1) ~= fseog(2)
        error('Different fs for EOG!?')
    end
    patinfo.fs_orig(:,2) = fseog;
    clear idxeeg idxeog
    
    %% EMG
    idxemg = strmatch('CHIN1',hdr.label);
    idxemg = [idxemg, strmatch('CHIN2',hdr.label)];
    fsemg = hdr.frequency(idxemg);
    emg1 = record(idxemg(1),:);
    emg2 = record(idxemg(2),:);
    emg = (emg1 - emg2);
    if fsemg(1) ~= fsemg(2)
        error('Different fs for emg!?')
    end
    if any(isnan(emg))
        error('Whats up here!?')
    end
    patinfo.fs_orig(:,3) = fsemg;
    clear  idxeeg idxemg
    %% ECG
    idxecg = strmatch('ECGL',hdr.label);
    idxecg = [idxecg, strmatch('ECGR',hdr.label)];
    if isempty(idxecg)
        idxecg = strmatch('ECG1',hdr.label);
        idxecg = [idxecg, strmatch('ECG2',hdr.label)];        
    end
    
    fsecg = hdr.frequency(idxecg);
    ecg1 = record(idxecg(1),:);
    ecg2 = record(idxecg(2),:);
    ecg = (ecg1 - ecg2);
    if fsecg(1) ~= fsecg(2)
        error('Different fs for ecg!?')
    end
    if any(isnan(ecg))
        error('Whats up here!?')
    end
    patinfo.fs_orig(:,4) = fsecg;    
    %% Leg EMGs
    idxemg = strmatch('LAT1',hdr.label);
    idxemg = [idxemg, strmatch('LAT2',hdr.label)];
    if isempty(idxemg)
        idxemg = strmatch('LLEMG',hdr.label);
        idxemg = [idxemg, strmatch('EMGref',hdr.label)];        
    end
    fsemg = hdr.frequency(idxemg);
    emg1 = record(idxemg(1),:);
    emg2 = record(idxemg(2),:);
    emg_l = (emg1 - emg2);
    if fsemg(1) ~= fsemg(2)
        error('Different fs for emg!?')
    end
    if any(isnan(emg))
        error('Whats up here!?')
    end
    patinfo.fs_orig(:,5) = fsemg;
    
    idxemg = strmatch('RAT1',hdr.label);
    idxemg = [idxemg, strmatch('RAT2',hdr.label)];
    if isempty(idxemg)
        idxemg = strmatch('RLEMG',hdr.label);
        idxemg = [idxemg, strmatch('EMGREF',hdr.label)];        
    end    
    fsemg = hdr.frequency(idxemg);
    emg1 = record(idxemg(1),:);
    emg2 = record(idxemg(2),:);
    emg_r = (emg1 - emg2);
    if fsemg(1) ~= fsemg(2)
        error('Different fs for emg!?')
    end
    if any(isnan(emg))
        error('Whats up here!?')
    end
    patinfo.fs_orig(:,6) = fsemg;
    
    %% MIC
    idxmic = strmatch('SNORE',hdr.label);


    if isempty(idxmic)
        disp('No MIC daata!?');
        mic = [];
        fsmic = [];
    else
        fsmic = hdr.frequency(idxmic);
        mic = record(idxmic(1),:);        
        patinfo.fs_orig(:,7) = fsmic;                
    end
    
    %% CHEST
    idxchest = strmatch('CHEST',hdr.label);
    if isempty(idxchest)
        disp('No Chest Data!?');
        chest = [];      
        fschest = [];
    else
        fschest = hdr.frequency(idxchest);
        chest = record(idxchest(1),:);
        patinfo.fs_orig(:,8) = fschest;           
    end

    
    %% FLOW
    idxflow = strmatch('FLOW',hdr.label);
    if isempty(idxflow)
        disp('No Flow Data!?');
        flow = [];   
        fsflow = [];
    else
        fsflow = hdr.frequency(idxflow);
        flow = record(idxflow(1),:); 
        patinfo.fs_orig(:,9) = fsflow;                 
    end

    %% Filtering
    fss = patinfo.fs_orig(1,:);
    eeg = resample(eeg,patinfo.fs,fss(1));
    eog = resample(eog,patinfo.fs,fss(2));
    emg = resample(emg,patinfo.fs,fss(3));
    emg_l = resample(emg_l,patinfo.fs,fss(5));
    emg_r = resample(emg_r,patinfo.fs,fss(6));
    ecg = resample(ecg,patinfo.fs,fss(4));
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
    clear fseeg fseog fsemg fss fsflow fschest
    clear record idxeeg idxemg    

%% Ensure signals is calibrated to microvolts
    Factor = 1;
    if max(eeg) > 10   
        Factor = 1;
    else
        %MilliVolts, change to micro
        Factor = 1000;
        eeg = eeg*Factor;
    end    
    if max(eog) > 10   
        Factor = 1;
    else
        %MilliVolts, change to micro
        Factor = 1000;
        eog = eog*Factor;
    end      
    if max(emg) > 10   
        Factor = 1;
    else
        %MilliVolts, change to micro
        Factor = 1000;
        emg = emg*Factor;
    end       
    if max(emg_l) > 10   
        Factor = 1;
    else
        %MilliVolts, change to micro
        Factor = 1000;
        emg_l = emg_l*Factor;
    end       
    if max(emg_r) > 10   
        Factor = 1;
    else
        %MilliVolts, change to micro
        Factor = 1000;
        emg_r = emg_r*Factor;
    end     
    
    Nfir = 500;
    b_band = fir1(Nfir,[0.3 40].*2/patinfo.fs,'bandpass'); % bandpass
    eeg = filtfilt(b_band,1,eeg);
    b_band = fir1(Nfir,[0.3 40].*2/patinfo.fs,'bandpass'); % bandpass
    eog = filtfilt(b_band,1,eog);    
    
%   ecg = filtfilt(b_band,1,ecg);      
    
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
    %% Resample and merge into matrix
    signal = [eeg;eog;emg;ecg;emg_r;emg_l;mic;chest;flow];
    signal(:,1:skip) = [];

    [epoch, labels]  = importcsv(sprintf('P%s_Night%s.csv',string(patid),char(night{1})));
    skip_end_s = epoch(end)*30 - length(signal)/patinfo.fs; %remove +30?
    if skip_end_s > 0 %annotated hyp is larger than psg recording (cut-off signal and hyp)
       skip_end = skip_end_s*patinfo.fs; 
       signal(:,end-skip_end:end) = [];
       epoch(end-floor(skip_end_s/30):end) = [];
       labels(end-floor(skip_end_s/30):end,:) = [];          
    end
    
    % Standardizing signals (does not help and not neccessary)
    if std_flg == 1
        for s=1:size(signal,1)
            signal(:,s) = signal(:,s) - nanmean(signal(:,s));
            signal(:,s) = signal(:,s)./nanstd(signal(:,s));
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
        tstamp(seg) = (epoch(seg))*30*patinfo.fs + 1; %not sure about +1
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



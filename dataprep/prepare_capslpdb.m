function prepare_capslpdb(filename,outputfolder,std_flg)
% This function prepares data from raw edf/txt files of PSG recordings.
% Specifically preparing data for each subject that provides neccessary
% signals (eg. EEG, EOG, EMG etc). Data is saved in mat format within
% folder
%
% Inputs:
%  filename    - location folder of PSG recordings
%  outputfolder   - folder to save prepared PSG signals in mat format for
%                   each subject
%  std_flg     - Standardise Flag: standardise each signal 
%                (mean = 0 std_dec = 1)
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

cd(filename)
fls = dir('*.edf');
fls = arrayfun(@(x) x.name,fls,'UniformOutput',false);
genderage = readtable('gender-age.xlsx');
CLASSES = {'W','S1','S2','S3','R'};

if isempty(std_flg)
   std_flg = 0; 
end

for f = 1:length(fls)
    disp(fls{f})
    
    % Ignoring following recordings:
    % brux1 EMG fs (100Hz) < 200Hz
    % 8: n16 has no EOG or EMG
    % 11: n4 has no EOG or EMG (EOG sin, tib sin)
    % 4: n12 EMG fs (100Hz) < 200Hz
    % 15: n8 has no EMG
    if regexp(fls{f},'(brux1|n4|n8|n12|n16)')
        sprintf('Skipping %s due to inconsistencies', fls{f})
        continue
    end
    
    % Loading data
    [hdr, record] = edfread(fls{f});
    idxrow = strcmpi(genderage.Pathology,fls{f}(1:end-4));
    patinfo.gender = genderage.Gender{idxrow};
    patinfo.age = genderage.Age(idxrow);
    patinfo.ch_orig = cell(2,6);
    patinfo.fs_orig = nan(2,6);
    patinfo.fs = 200;
    patinfo.classes = CLASSES;
%     patinfo.chlabels = {'EEG','EOG','EMG'};
    patinfo.chlabels = {'EEG','EOG','EMG','ECG','MIC','CHEST'};
    
    %% Getting EEG
    idxeeg = find(~cellfun(@isempty, regexp(hdr.label,'C4-A1|C4A1')));
    if isempty(idxeeg)
        idxeeg = find(~cellfun(@isempty, regexp(hdr.label,'C3-A2|C3A2')));
    end
    if isempty(idxeeg)
        idxeeg = find(~cellfun(@isempty, regexp(hdr.label,'(C4)')));
        idxeeg = [idxeeg,find(~cellfun(@isempty, regexp(hdr.label,'(A1)')))];
    end
    
    % sanity check for multiple channels
    if length(idxeeg)>1
        eeg = diff(record(idxeeg,:))';
    elseif length(idxeeg)== 1
        eeg = record(idxeeg,:)';
    else
        warning('Skipping record, no EEG')
        continue
    end
    % Getting fs
    fseeg = hdr.frequency(idxeeg);
    if length(fseeg)>1
        if fseeg(1) ~= fseeg(2)
            error('Different fs for EEG!?')
        end
    end
    if any(fseeg < 100), error('Low sampling frequency?'), end
    for l = 1:length(idxeeg)
        patinfo.ch_orig(l,1) = cellstr(hdr.label(idxeeg(l)));
        patinfo.fs_orig(l,1) = fseeg(l);
    end
    clear fseeg idxeeg l
    %% Getting EOG
    idxeog = find(~cellfun(@isempty, regexp(hdr.label,'(ROC-LOC|EOG dx|ROCLOC)')));
    if isempty(idxeog)
        idxeog = find(~cellfun(@isempty, regexp(hdr.label,'(ROC|EOG-R|ROC-A2)')));
        idxeog = [idxeog,find(~cellfun(@isempty, regexp(hdr.label,'(LOC|EOG-L|LOC-A1)')))];
    end
    
    % sanity check for multiple channels
    if length(idxeog)>1
        eog = diff(record(idxeog,:))';
    elseif length(idxeog)== 1
        eog = record(idxeog,:)';
    else
        warning('Skipping record, no EOG')
        continue
    end
    
    % Getting fs
    fseog = hdr.frequency(idxeog);
    if length(fseog)>1
        if fseog(1) ~= fseog(2)
            error('Different fs for EOG!?')
        end
    end
    if any(fseog < 100), error('Low sampling frequency?'), end
    for l = 1:length(idxeog)
        patinfo.ch_orig(l,2) = cellstr(hdr.label(idxeog(l)));
        patinfo.fs_orig(l,2) = fseog(l);
    end
    clear fseog idxeog l
    %% Getting EMG
    idxemg = find(~cellfun(@isempty, regexp(hdr.label,'(EMG1-EMG2|EMG-EMG|CHIN-CHIN)')));
    if isempty(idxemg)
        idxemg = find(~cellfun(@isempty, regexp(hdr.label,'CHIN|EMG','match')));
    end
    
    % sanity check for multiple channels
    if length(idxemg)>1
        emg = diff(record(idxemg,:))';
    elseif length(idxemg)== 1
        emg = record(idxemg,:)';
    else
        warning('Skipping record, no EMG')
        continue
    end
    
    % Getting fs
    fsemg = hdr.frequency(idxemg);
    if length(fsemg)>1
        if fsemg(1) ~= fsemg(2)
            error('Different fs for EMG!?')
        end
    end
    if any(fsemg < 100), error('Low sampling frequency?'), end
    for l = 1:length(idxemg)
        patinfo.ch_orig(l,3) = cellstr(hdr.label(idxemg(l)));
        patinfo.fs_orig(l,3) = fsemg(l);
    end
    clear fsemg idxemg l
    %% Getting ECG
    
    idxecg = find(~cellfun(@isempty, regexp(hdr.label,'(ECG1-ECG2|ECG|EKG)')));
   
    % sanity check for multiple channels
    if length(idxecg)>1
        ecg = diff(record(idxecg,:))';
    elseif length(idxecg)== 1
        ecg = record(idxecg,:)';
    else
        warning('Skipping record, no ECG')
        continue
    end    
    
    % Getting fs
    fsecg = hdr.frequency(idxecg);
    if length(fsecg)>1
        if fsecg(1) ~= fsecg(2)
            error('Different fs for ECG!?')
        end
    end
    if any(fsecg < 100), error('Low sampling frequency?'), end
    for l = 1:length(idxecg)
        patinfo.ch_orig(l,4) = cellstr(hdr.label(idxecg(l)));
        patinfo.fs_orig(l,4) = fsecg(l);
    end
    clear fsecg idxecg l     
    
    %% MIC
    idxmic = find(~cellfun(@isempty, regexp(hdr.label,'(MIC)')));
   
    % sanity check for multiple channels
    if length(idxmic)>1
        mic = diff(record(idxmic,:))';
    elseif length(idxmic)== 1
        mic = record(idxmic,:)';
    else
        warning('Skipping record, no MIC')
        mic = [];
%         continue
    end    
    
    % Getting fs
    fsmic = hdr.frequency(idxmic);
    if length(fsmic)>1
        if fsmic(1) ~= fsmic(2)
            error('Different fs for MIC!?')
        end
    end
    if any(fsmic < 100), error('Low sampling frequency?'), end
    for l = 1:length(idxmic)
        patinfo.ch_orig(l,5) = cellstr(hdr.label(idxmic(l)));
        patinfo.fs_orig(l,5) = fsmic(l);
    end
    clear fsmic idxmic l    

    %% CHEST
    idxchest = find(~cellfun(@isempty, regexp(hdr.label,'(TORACE)')));
   
    % sanity check for multiple channels
    if length(idxchest)>1
        chest = diff(record(idxchest,:))';
    elseif length(idxchest)== 1
        chest = record(idxchest,:)';
    else
        warning('Skipping record, no CHEST')
        chest = [];
%         continue
    end    
    
    % Getting fs
    fschest = hdr.frequency(idxchest);
    if length(fschest)>1
        if fschest(1) ~= fschest(2)
            error('Different fs for CHEST!?')
        end
    end
%     if any(fschest < 100), error('Low sampling frequency?'), end
    for l = 1:length(idxchest)
        patinfo.ch_orig(l,6) = cellstr(hdr.label(idxchest(l)));
        patinfo.fs_orig(l,6) = fschest(l);
    end
    clear fschest idxchest l        
    %%
    
    clear record info idxrow
    
    %% Resampling signals
    fss = patinfo.fs_orig(1,:);    
    
    %% Preprocessing Filter coefficiens
    % Resampling to 100 Hz
    eeg = resample(eeg,patinfo.fs,fss(1));
    eog = resample(eog,patinfo.fs,fss(2));
    emg = resample(emg,patinfo.fs,fss(3));
    ecg = resample(ecg,patinfo.fs,fss(4));
    
    if (~isempty(mic))
        mic = resample(mic,patinfo.fs,fss(5));
    end
    
    if (~isempty(chest))
        chest = resample(chest,patinfo.fs,fss(6));
    end    
    clear fss
    
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
    
    Nfir = 500;
    % Preprocessing filters
    b_band = fir1(Nfir,[0.3 40].*2/patinfo.fs,'bandpass'); % bandpass
    eeg = filtfilt(b_band,1,eeg);
    b_band = fir1(Nfir,[0.3 40].*2/patinfo.fs,'bandpass'); % bandpass
    eog = filtfilt(b_band,1,eog);
    clear b_band
    
    % Preprocessing filters
    pwrline1 = 50; %Hz
    b_notch1 = fir1(Nfir,[(pwrline1-1) (pwrline1+1)].*2/patinfo.fs,'stop');
    pwrline2 = 60; %Hz
    b_notch2 = fir1(Nfir,[(pwrline2-1) (pwrline2+1)].*2/patinfo.fs,'stop');
    b_band = fir1(Nfir,10.*2/patinfo.fs,'high'); % bandpass
    emg = filtfilt(b_notch1,1,emg);
    emg = filtfilt(b_notch2,1,emg);
    emg = filtfilt(b_band,1,emg);
    
    if (~isempty(mic))
        b_band = fir1(Nfir,3.*2/patinfo.fs,'high'); % bandpass
        mic = filtfilt(b_notch1,1,mic);
        mic = filtfilt(b_notch2,1,mic);
        mic = filtfilt(b_band,1,mic);        
    end    
    % cut to shortest signal (ie stops flat-lining)
    rem = find(eeg(end:-1:1) ~= 0,1,'first');
    if (length(eeg) - rem)/patinfo.fs/60/60 < 5 % less than five hours available
        warning('This looks fishy, wrong fs? Skipping..')
        continue
    end
    eeg(end-rem:end) = [];
    rem = find(eog(end:-1:1) ~= 0,1,'first');
    if (length(eog) - rem)/patinfo.fs/60/60 < 5 % less than five hours available
        warning('This looks fishy, wrong fs? Skipping..')
        continue
    end
    eog(end-rem:end) = [];
    rem = find(emg(end:-1:1) ~= 0,1,'first');
    if (length(emg) - rem)/patinfo.fs/60/60 < 5 % less than five hours available
        warning('This looks fishy, wrong fs? Skipping..')
        continue
    end
    emg(end-rem:end) = [];
    rem = find(ecg(end:-1:1) ~= 0,1,'first');
    if (length(ecg) - rem)/patinfo.fs/60/60 < 5 % less than five hours available
        warning('This looks fishy, wrong fs? Skipping..')
        continue
    end    
    ecg(end-rem:end) = [];
    
    if ~isempty(mic)
        rem = find(mic(end:-1:1) ~= 0,1,'first');
        if (length(mic) - rem)/patinfo.fs/60/60 < 5 % less than five hours available
            warning('This looks fishy, wrong fs? Skipping..')
            continue
        end            
        mic(end-rem:end) = [];        
    end
    if ~isempty(chest)
        rem = find(chest(end:-1:1) ~= 0,1,'first');
        if (length(chest) - rem)/patinfo.fs/60/60 < 5 % less than five hours available
            warning('This looks fishy, wrong fs? Skipping..')
            continue
        end            
        chest(end-rem:end) = [];        
    end
    % merging elements into matrix    
    len = [length(eeg),length(eog),length(emg),length(ecg),length(mic),length(chest)];
    len(len<=0) = inf;
    len = min(len); %find min len that's not zero
    signals = [eeg(1:len)'];
    signals = [signals;eog(1:len)'];
    signals = [signals;emg(1:len)'];
    signals = [signals;ecg(1:len)'];
    if ~isempty(mic)
        signals = [signals;mic(1:len)'];
    end
    if ~isempty(chest)
        signals = [signals;chest(1:len)'];
    end    
    % Standardizing signals
    if std_flg == 1
        for s=1:size(signals,2)
            signals(s,:) = signals(s,:) - nanmean(signals(s,:));
            signals(s,:) = signals(s,:)./nanstd(signals(s,:));
        end
    end
            
    clear eeg eog emg ecg mic chest len rem
    %% Figuring out annotations
    infotxt = loadtxt([fls{f}(1:end-4) '.txt']);
    try
        anntm = datetime(infotxt.Timehhmmss,'InputFormat','HH:mm:ss');
    catch
        try
            anntm = datetime(cellstr(infotxt.Timehhmmss),'InputFormat','HH.mm.ss');
        catch
            try   % imnsonia ones dont have second column
                anntm = datetime(cellstr(infotxt.Position),'InputFormat','HH:mm:ss');
            catch
                anntm = datetime(cellstr(infotxt.Position),'InputFormat','HH.mm.ss');
            end
            
        end
    end
    
    recstart = datetime(hdr.starttime,'InputFormat','HH.mm.ss');
    
    % convert labels for neural networks
    stage = infotxt.SleepStage;
    stage = cellstr(stage);
    stage(cellfun(@(x) strcmp(x,'S4'), stage)) = {'S3'}; % converting to AASM annotation
    stage(cellfun(@(x) strcmp(x,'REM'), stage)) = {'R'}; % converting to AASM annotation
    rem = ~cellfun(@(x) ismember(x,CLASSES),stage); % remove S4 and MT
    stage(rem) = [];
    anntm(rem) = [];
    
    % in case annotation started after midnight
    for s = 1:length(stage)
        if seconds(anntm(s)-recstart) < 0
            anntm(s) = anntm(s)+1; %plus a day
        end
    end
    
    lastsamp = (seconds(anntm(end)-recstart)+30)*patinfo.fs+1; %number of total samples according to hyp (time starts from 0)
    signals(:,lastsamp:end) = [];
    data = zeros(length(stage)-1,30*patinfo.fs,size(signals,1));
    labels = zeros(length(stage)-1,5);
    epoch = zeros(length(labels),1);
    for s = 1:(length(labels)-1)
        try
            startsamp = seconds(anntm(s)-recstart)*patinfo.fs;
            if startsamp == 0
                startsamp = 1;
            end
            epoch(s) = round(startsamp/patinfo.fs);
            data(s,:,:) = signals(:,startsamp:startsamp+30*patinfo.fs-1)';
            labels(s,:) = ismember(CLASSES,stage(s));
        catch
            if (startsamp+30*patinfo.fs-1) > size(signals,2)
                continue
            end
        end
    end
    
    if any(isnan(data))
        error('NaNs')
    end
    idx = ~any(labels,2);  %rows
    labels(idx,:) = [];
    data(idx,:,:) = [];
    epoch(idx,:) = [];
    if size(data,1) ~= size(labels,1)
        error('Different length for recording..')
    end
    
    if ~exist(outputfolder, 'dir')
       mkdir(outputfolder);
    end     
    %% Update patinfo
    
    empty_idx = isnan(patinfo.fs_orig(1,:));
    patinfo.fs_orig(:,empty_idx) = [];
    patinfo.ch_orig(:,empty_idx) = [];
    patinfo.chlabels(:,empty_idx) = [];

    
    %Saving results
    save([outputfolder,fls{f}(1:end-3),'mat'],'data','labels','patinfo','epoch')
    clear infotxt startsamp lastsamp recstart stage anntm rem data labels patinfo signals r
end
end
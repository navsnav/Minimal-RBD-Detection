function prepare_massdb(filename,outputfolder,std_flg)
parentfold = filename;
cd(parentfold);

genderage = readtable([parentfold 'gender-age.csv'],'Format','%s%d%s');

fls = dir('*.edf');
fls = arrayfun(@(x) x.name,fls,'UniformOutput',false);
%%%genderage = readtable('gender-age.csv');
CLASSES = {'W','S1','S2','S3','R'};
    
    
    for f = 1:length(fls)
        info = wfdbdesc(fls{f});
        disp(fls{f})
        
        % Loading data
        [hdr, record] = edfread(fls{f});
        
        idxrow = strcmpi(genderage.Name,fls{f}(1:end-8));
        patinfo.gender = genderage.Sex{idxrow};
        patinfo.age = genderage.Age(idxrow);
        %%idxrow = strcmpi(genderage.Var1,fls{f}(1:end-4));
        %%%patinfo.gender = genderage.Var2{idxrow};
        %%%patinfo.age = genderage.Var3(idxrow);
        patinfo.ch_orig = cell(2,5);
        patinfo.fs_orig = nan(2,5);
        patinfo.fs = 200;
        patinfo.classes = CLASSES;
%         patinfo.chlabels = {'EEG','EOG','EMG'};
        patinfo.chlabels = {'EEG','EOG','EMG','ECG','CHEST'}; %No MIC
        
        %% Getting EEG
        idxeeg = zeros(2,1);
        if any(~cellfun(@isempty, regexp({info.Description},'(EEG A1)|(EEG A2)')))
            try
                idxeeg(1) = find(~cellfun(@isempty, regexp({info.Description},'(EEG C4-CLE)|(EEG C4-LER)')));
                idxeeg(2) = find(~cellfun(@isempty, regexp({info.Description},'EEG A1-CLE')));
            catch
                idxeeg(1) = find(~cellfun(@isempty, regexp({info.Description},'EEG C3-CLE')));
                idxeeg(2) = find(~cellfun(@isempty, regexp({info.Description},'EEG A2-CLE')));
            end
            eeg = diff(record(idxeeg,:))';
            
        else
            try
                idxeeg = find(~cellfun(@isempty, regexp({info.Description},'(EEG C4-LER)')));
            catch
                idxeeg = find(~cellfun(@isempty, regexp({info.Description},'EEG C3-LER')));
            end
            eeg = record(idxeeg,:)';
        end
        % sanity check for multiple channels
        if any(idxeeg==0)
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
            patinfo.ch_orig(l,1) = cellstr(info(idxeeg(l)).Description);
            patinfo.fs_orig(l,1) = fseeg(l);
        end
        clear fseeg idxeeg l
        %% Getting EOG
        idxeog = zeros(2,1);
        try
            idxeog(1) = find(~cellfun(@isempty, regexp({info.Description},'EOG Right Horiz')));
            idxeog(2) = find(~cellfun(@isempty, regexp({info.Description},'EOG Left Horiz')));
        catch
            %idxeog(1) = find(~cellfun(@isempty, regexp({info.Description},'EEG C3-CLE')));
            %idxeog(2) = find(~cellfun(@isempty, regexp({info.Description},'EEG A2-CLE')));
            error('Alternative name')
        end
        
        % sanity check for multiple channels
        if any(idxeog==0)
            warning('Skipping record, no EOG')
            continue
        end
        eog = diff(record(idxeog,:))';
        
        % Getting fspatinfo.fs
        fseog = hdr.frequency(idxeog);
        if length(fseog)>1
            if fseog(1) ~= fseog(2)
                error('Different fs for EOG!?')
            end
        end
        if any(fseog < 100), error('Low sampling frequency?'), end
        for l = 1:length(idxeog)
            patinfo.ch_orig(l,2) = cellstr(info(idxeog(l)).Description);
            patinfo.fs_orig(l,2) = fseog(l);
        end
        clear fseog idxeog l
        %% Getting EMG
        if (sum(~cellfun(@isempty, regexp({info.Description},'EMG Chin'))) < 2)
            
            idxemg = find(~cellfun(@isempty, regexp({info.Description},'EMG Chin')));
            emg = record(idxemg,:)';
        else
            try
                idxemg(1) = find(~cellfun(@isempty, regexp({info.Description},'EMG Chin1')));
                idxemg(2) = find(~cellfun(@isempty, regexp({info.Description},'EMG Chin2')));
            catch
                %idxeog(1) = find(~cellfun(@isempty, regexp({info.Description},'EEG C3-CLE')));
                %idxeog(2) = find(~cellfun(@isempty, regexp({info.Description},'EEG A2-CLE')));
                error('Alternative name')
            end
            emg = diff(record(idxemg,:))';
            
        end
        
        % sanity check for multiple channels
        if any(idxemg==0)
            warning('Skipping record, no EOG')
            continue
        end
        
        % Getting fs
        fsemg = hdr.frequency(idxemg);
        if length(fsemg)>1
            if fsemg(1) ~= fsemg(2)
                error('Different fs for EOG!?')
            end
        end
        if any(fsemg < 100), error('Low sampling frequency?'), end
        for l = 1:length(idxemg)
            patinfo.ch_orig(l,3) = cellstr(info(idxemg(l)).Description);
            patinfo.fs_orig(l,3) = fsemg(l);
        end
        clear fsemg idxemg l
        
        %% Getting ECG
        idxecg = zeros(2,1);
        try
            idxecg = find(~cellfun(@isempty, regexp({info.Description},'ECG I')));
        catch

        end
        ecg = diff(record(idxecg,:))';

        % sanity check for multiple channels
        if any(idxecg==0)
            warning('Skipping record, no ECG')
            continue
        end
        % Getting fs
        fsecg = hdr.frequency(idxecg);
        if length(fsecg)>1
            if fsecg(1) ~= fsecg(2)
                error('Different fs for EEG!?')
            end
        end
        if any(fsecg < 100), error('Low sampling frequency?'), end
        for l = 1:length(idxecg)
            patinfo.ch_orig(l,4) = cellstr(info(idxecg(l)).Description);
            patinfo.fs_orig(l,4) = fsecg(l);
        end
        clear fsecg idxecg l        
        %% Getting CHEST
        idxchest = zeros(2,1);
        try
            idxchest = find(~cellfun(@isempty, regexp({info.Description},'Resp Belt Thor')));
        catch

        end
        chest = diff(record(idxchest,:))';

        % sanity check for multiple channels
        if any(idxchest==0)
            warning('Skipping record, no CHEST')
            continue
        end
        % Getting fs
        fschest = hdr.frequency(idxchest);
        if length(fschest)>1
            if fschest(1) ~= fschest(2)
                error('Different fs for CHEST!?')
            end
        end
%         if any(fschest < 100), error('Low sampling frequency?'), end
        for l = 1:length(idxchest)
            patinfo.ch_orig(l,5) = cellstr(info(idxchest(l)).Description);
            patinfo.fs_orig(l,5) = fschest(l);
        end
        clear fschest idxchest l        
        
        clear record info idxrow
        
        %% Resampling signals
        fss = round(patinfo.fs_orig(1,:));
         % Resampling to 100 Hz
        eeg = resample(eeg,patinfo.fs,fss(1));
        eog = resample(eog,patinfo.fs,fss(2));
        emg = resample(emg,patinfo.fs,fss(3));
        ecg = resample(ecg,patinfo.fs,fss(4));
        chest = resample(chest,patinfo.fs,fss(5));

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
        
        clear fss
        %% Preprocessing Filter coefficiens
        Nfir = 500;
        
        % Preprocessing filters
        b_band = fir1(Nfir,[0.3 40].*2/patinfo.fs,'bandpass'); % bandpass
        eeg = filtfilt(b_band,1,eeg);
        
        clear b_notch1 b_notch2 b_band pwrline
        % Preprocessing filters
        b_band = fir1(Nfir,[0.3 40].*2/patinfo.fs,'bandpass'); % bandpass
        eog = filtfilt(b_band,1,eog);
        
        % Preprocessing filters
        pwrline = 50; %Hz
        b_notch1 = fir1(Nfir,[(pwrline-1) (pwrline+1)].*2/patinfo.fs,'stop');
        pwrline = 60; %Hz
        b_notch2 = fir1(Nfir,[(pwrline-1) (pwrline+1)].*2/patinfo.fs,'stop');
        b_band = fir1(Nfir,10.*2/patinfo.fs,'high'); % bandpass
        emg = filtfilt(b_notch1,1,emg);
        emg = filtfilt(b_notch2,1,emg);
        emg = filtfilt(b_band,1,emg);
        
     
        % cut to shortest signal
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
        % merging elements into matrix
        len = [length(eeg),length(eog),length(emg),length(ecg),length(chest)];
        len(len<=0) = inf;
        len = min(len); %find min len that's not zero
        signals = [eeg(1:len)'];
        signals = [signals;eog(1:len)'];
        signals = [signals;emg(1:len)'];
        signals = [signals;ecg(1:len)'];
        signals = [signals;chest(1:len)'];
        
        if std_flg == 1
            % Standardizing signals
            for s=1:size(signal,1)
                signals(s,:) = signals(s,:) - nanmean(signals(s,:));
                signals(s,:) = signals(s,:)./nanstd(signals(s,:));
            end
        end
        clear eeg eog emg ech chest len rem
        %% Figuring out annotations
        annfile = [fls{f}(1:end-8)];
        annfilesaf = [annfile '.saf'];
        annfileedf = dir([annfile '*.edf']);
        annfileedf = [annfileedf.name];
        if isunix()
            % Some shell magic
            derror = system(['sed "s/+/\n/g" ' annfilesaf ' > cleanedsaf.txt']);
            derror2 = system('sed -i "/Sleep stage/!d" cleanedsaf.txt');
            if (derror || derror2)
                error('Something went wront reading hypnogram')
            end
        else
            fid  = fopen(annfilesaf,'r');
            saf=fread(fid,'*char')';
            fclose(fid);
            saf = strrep(saf,'+',',');                    
            fid2  = fopen('cleanedsaf.txt','w');
            fwrite(fid2,saf,'*char');
            fclose(fid2);
   
        end
        
        if isunix()
            [fid,msg]=fopen('cleanedsaf.txt','r','n','UTF-8');
            hypno=textscan(fid,'%[^\n]','delimiter','\n');
            hypno = hypno{1};
            fclose(fid);
            delete('cleanedsaf.txt')
        else
            fid4 = fopen('cleanedsaf.txt');
            hypno=fread(fid4,'*char')';
            hypno = split(hypno,',');                     
            fclose(fid4);
            Index = find(contains(hypno,'Sleep stage '));
            hypno = hypno(Index); 
        end
        
        epoch = round(cellfun(@str2double,regexp(hypno,'^\d+\.?\d*','match')));
        nclasses = {'Sleep stage W' 'Sleep stage R' 'Sleep stage 1' 'Sleep stage 2' 'Sleep stage 3' 'Sleep stage 4'};
        stages = cellfun(@(x) regexp(x,nclasses),hypno,'UniformOutput',0);
        stages = cell2mat(cellfun(@(x) ~cellfun(@isempty, x), stages,'UniformOutput',0));
        stages = [stages(:,1:4), sum(stages(:,end-1:end),2)]; % converting R&K -> AASM
        stages = stages(:,[1, 3, 4, 5, 2]);  % reordering to match W, N1, N2, N3, R
        rem = ~any(stages,2);  %rows
        epoch(rem,:) = [];
        stages(rem,:) = [];
        clear nclasses hypno fid msg rem
        
        epoch = epoch*patinfo.fs; % including annotation delay on selected epochs
        if epoch(end) > length(signals)
            disp('More annotations than signal. Chop chop!')
            rem = epoch>length(signals);
            epoch(rem) = [];
            stages(rem,:) = [];
        end
        
        labels = logical(stages(1:length(epoch),:));
        clear stages delay ann hdrann hdr i l
        
        data = zeros(length(labels),30*patinfo.fs,3);
        annotinterv = roundn(mean(diff(epoch)/patinfo.fs),1); % annotation interval can be 20s or 30s
        
        if (annotinterv == 20)
            for s = 2:(length(epoch)-1)
                data(s,:,:) = signals(:,epoch(s)-5*patinfo.fs:epoch(s)+(annotinterv+5)*patinfo.fs-1)';
            end
        elseif annotinterv == 30
            for s = 1:(length(epoch))
                data(s,:,:) = signals(:,epoch(s):epoch(s)+annotinterv*patinfo.fs-1)';
            end
        end
        
        
        if any(isnan(data))
            error('NaNs')
        end
        idx = ~any(labels,2)|~any(all(data,2),3);  %rows
        labels(idx,:) = [];
        data(idx,:,:) = [];
        epoch(idx,:) = [];
        if size(data,1) ~= size(labels,1)
            error('Different length for recording..')
        end
        
%                 % Plotting for sanity check
%                 plot(signals')
%                 hold on
%                 for i = 1:length(epoch)
%                     line([epoch(i) epoch(i)],[-100,100],'Color','k')
%                 end
%                 txt = cell(length(epoch),1);
%                 for l = 1:length(epoch)
%                     txt(l) = CLASSES(labels(l,:));
%                 end
%                 text(epoch,100*ones(size(epoch)),txt)
%                 close
        %Saving results
        save([fls{f}(1:end-8)],'data','labels','patinfo','epoch')
        clear infotxt startsamp lastsamp recstart stage anntm rem data labels patinfo signals r
    end
end
function   [feats, features_struct, hr_data] = FeatExtract_HR_mat(hr_signal, fs,epoch_time,std_flag,buff_flag,subject,dbpath,hyp) %external function
% addpath('C:\Users\scro2778\Documents\CDT\Imaging\Practical 2 - Registration\Registration code 2015');
if buff_flag
    pad_zero_len = (epoch_time-30)/2; 
    signal = buffer([zeros(pad_zero_len*fs,1);hr_signal;zeros(pad_zero_len*fs,1)],epoch_time*fs,(epoch_time/30-1)*30*fs,'nodelay');    
else
    signal = reshape(hr_signal',fs*epoch_time,numel(hr_signal)/(fs*epoch_time));
end
NFEAT = 53; % number of features used
NFEAT_hrv = 47+7;
NFEAT_hrv = 73-5;
NFEAT_hrv = 73-5+3;

% Narrow BP
Fhigh = 5;  % highpass frequency [Hz]
Flow = 45;   % low pass frequency [Hz]
Nbut = 10;     % order of Butterworth filter
d_bp= design(fdesign.bandpass('N,F3dB1,F3dB2',Nbut,Fhigh,Flow,fs),'butter');
[b_bp,a_bp] = tf(d_bp);

if size(signal,1)<size(signal,2), signal = signal'; end % make sure it's column vector
signalraw =  signal;

%% Preprocessing
signal = filtfilt(b_bp,a_bp,signal);             % filtering narrow
signal = detrend(signal);    
% detrending (optional)
if std_flag==1
    signal = signal - mean(signal(:));
    signal = signal./std(signal(:));                     % standardizing
end
%% get mean for entire signal 

% reshape into 5min epochs
signal5 = reshape(signal,numel(signal),1);
if mod(numel(signal5),fs*5*60) > 0
   pad =  fs*5*60-mod(numel(signal5),fs*5*60);
   signal5 = [signal5;zeros(pad,1)];
end
signal5 = reshape(signal5',fs*5*60,numel(signal5)/(fs*5*60));
% Get qrs for each segment
all_qrs = [];
all_qrs_idx = [];

for n=1:size(signal5,2)
    % Get signal of interest
    sig_seg = signal5(:,n);

    % QRS detect
    [qrsseg,featqrs] = multi_qrsdetect(sig_seg,fs,['s_' num2str(n)]);
    
    qrs_t = qrsseg{2}./fs+(n-1)*5*60; %P%T: 2
    qrs_idx = qrsseg{2}+(n-1)*5*60*fs;
    qrs = [qrs_t];
    all_qrs_idx = [all_qrs_idx;qrs_idx];
    all_qrs = [all_qrs;qrs];
    % HRV features
    
end

% get rr features
if size(all_qrs,1) > size(all_qrs,2); all_qrs = all_qrs';end % check column or row vector

% HRV features
hrv_now=[all_qrs(2:end);diff(all_qrs)];
hrv_now_all = hrv_now;
%% Averaging Window

% raw_rr = hrv_now(2,:);
%  
% raw_rr_buff = buffer([zeros(1,10)+raw_rr(1),raw_rr,zeros(1,10)+raw_rr(end)],21,20,'nodelay');
% raw_rr_mean = mean(raw_rr_buff,1);
% % 
% % small_idx = raw_rr < 0.5*raw_rr_mean;
% % big_idx = raw_rr > 1.5*raw_rr_mean;
% % 
% rr_win = raw_rr;
% rr_win(small_idx) = raw_rr_mean(small_idx);
% rr_win(big_idx) = raw_rr_mean(big_idx);
% 
% all_HRV=get_hrv2(hrv_now');
all_HRV=get_hrv_all(hrv_now'); %use raw qrs rr's
% all_HRV=get_hrv_all([hrv_now(1,:);rr_win]'); %use windowed average qrs rr's
% % hrv_now_all_orig = hrv_now_all;
% hrv_now_all(2,:) = rr_win; % Replace with avg mRR

clear signal5;
% Plot to see detected QRS peaks
fig1 = figure;
plot(hr_signal);
hold on
plot(all_qrs_idx,hr_signal(all_qrs_idx),'r*');
xlabel('Sample');
ylabel('HR Signal');
title(subject,'interpreter','none');
% saveas(fig1,[dbpath,subject],'fig');
%% averaging window for RR

%% Create sliding 5min signal buffer
signal5min = reshape(signal,numel(signal),1);

% t = 0:1/fs:length(signal5min)/fs-1/fs;
% fig_hr = figure;
% plot(t,hr_signal)
% hold on
% plot(t(all_qrs_idx),signal5min(all_qrs_idx),'r*');
% title([subject,'_qrs'],'Interpreter','none');
% xlabel('Time (s)');
% ylabel('ECG');
% saveas(fig_hr,[dbpath,'\',subject,'_qrs'],'fig');
% close(fig_hr);

% rem_idx = ismember(hyp(:,1),5);
% rem_time = hyp(rem_idx,2);
% 
% rem_qrs_idx = (rem_time <= all_qrs )&(all_qrs <= rem_time+30);
% fig_l = figure;
% for i=1:size(rem_qrs_idx,1)
%     rem_all_qrs = all_qrs(rem_qrs_idx(i,:));
%     deltaRRi = diff(diff(rem_all_qrs));
%     deltaRRi2 = deltaRRi(2:end);
%     deltaRRi1 = deltaRRi(1:end-1);
%     plot(deltaRRi2,deltaRRi1,'r.');    
%     hold on;    
% end
% title([subject,'_REM_lorenz'],'Interpreter','none');
% xlabel('delta RR(i)');
% ylabel('delta RR(i-1)');
% xlim([-1,1]);
% ylim([-1,1]);
% grid on;
% saveas(fig_l,[dbpath,subject,'REM_lorenz'],'fig');
% close(fig_l);
% 
% nrem_idx = ismember(hyp(:,1),[1,2,3,4]);
% nrem_time = hyp(nrem_idx,2);
% nrem_qrs_idx = (nrem_time <= all_qrs )&(all_qrs <= nrem_time+30);
% fig_l = figure;
% for i=1:size(nrem_qrs_idx,1)
%     nrem_all_qrs = all_qrs(nrem_qrs_idx(i,:));
%     deltaRRi = diff(diff(nrem_all_qrs));
%     deltaRRi2 = deltaRRi(2:end);
%     deltaRRi1 = deltaRRi(1:end-1);
%     plot(deltaRRi2,deltaRRi1,'b.');    
%     hold on;    
% end
% title([subject,'_NREM_lorenz'],'Interpreter','none');
% xlabel('delta RR(i)');
% ylabel('delta RR(i-1)');
% xlim([-1,1]);
% ylim([-1,1]);
% grid on;
% saveas(fig_l,[dbpath,subject,'NREM_lorenz'],'fig');
% close(fig_l);
% 
% 
% wake_idx = ismember(hyp(:,1),0);
% wake_time = hyp(wake_idx,2);
% wake_qrs_idx = (wake_time <= all_qrs )&(all_qrs <= wake_time+30);
% fig_l = figure;
% for i=1:size(wake_qrs_idx,1)
%     wake_all_qrs = all_qrs(wake_qrs_idx(i,:));
%     deltaRRi = diff(diff(wake_all_qrs));
%     deltaRRi2 = deltaRRi(2:end);
%     deltaRRi1 = deltaRRi(1:end-1);
%     plot(deltaRRi2,deltaRRi1,'k.');    
%     hold on;    
% end
% title([subject,'_WAKE_lorenz'],'Interpreter','none');
% xlabel('delta RR(i)');
% ylabel('delta RR(i-1)');
% xlim([-1,1]);
% ylim([-1,1]);
% grid on;
% saveas(fig_l,[dbpath,subject,'WAKE_lorenz'],'fig');
% close(fig_l);
%%

% [I,rem_qrs_idx_col] = find(rem_qrs_idx);
% rem_all_qrs = all_qrs(rem_qrs_idx_col);
% fig_l = figure;
% deltaRRi = diff(diff(rem_all_qrs));
% deltaRRi2 = deltaRRi(2:end);
% deltaRRi1 = deltaRRi(1:end-1);
% plot(deltaRRi2,deltaRRi1,'r.');
% title([subject,'_REM_lorenz'],'Interpreter','none');
% xlabel('delta RR(i)');
% ylabel('delta RR(i-1)');
% xlim([-1,1]);
% ylim([-1,1]);
% grid on;

% fig_l = figure;
% deltaRRi = diff(diff(all_qrs));
% deltaRRi2 = deltaRRi(2:end);
% deltaRRi1 = deltaRRi(1:end-1);
% plot(deltaRRi2,deltaRRi1,'.');
% title([subject,'_lorenz'],'Interpreter','none');
% xlabel('delta RR(i)');
% ylabel('delta RR(i-1)');
% xlim([-1,1]);
% ylim([-1,1]);
% grid on;
% saveas(fig_l,[dbpath,subject,'_lorenz'],'fig');
% close(fig_l);


signal = buffer([zeros(135*fs,1);signal5min;zeros(135*fs,1)],5*60*fs,5*60*fs-30*fs,'nodelay' );

time_pts = 0:1/fs:length(signal5min)/fs-1/fs+135+135;

time_pts_buff = buffer(time_pts,5*60*fs,5*60*fs-30*fs,'nodelay' );

%% Get Features


fetbag = {};
feat_hrv = [];
feats = [];
%     parfor n = 1:nseg
for n=1:size(signal,2)
    % Get signal of interest
    sig_seg = signal(:,n);
   
    time_range = time_pts_buff(:,n);
    time_range_idx = (hrv_now_all(1,:)>time_range(1)) & (hrv_now_all(1,:)<time_range(end));
    hrv_now_orig = hrv_now_all(:,time_range_idx);
    hrv_now_orig(1,:) = hrv_now_orig(1,:)-time_range(1);    
    

    % HRV features (using P&T - 2)
    if length(hrv_now_orig(1,:))>5 % if too few detections, returns zeros
        try
            feat_basic=HRV_features4(sig_seg,hrv_now_orig(1,:),fs,all_HRV.mRR);
            feat_hrv = [feat_basic];

%                 feats_poincare = get_poincare(qrsseg{end}./fs,fs);
%                 feat_hrv = [feat_basic, feats_poincare];
            feat_hrv(~isreal(feat_hrv)|isnan(feat_hrv)|isinf(feat_hrv)) = 0; % removing not numbers
        catch
            warning('Some HRV code failed.')
            feat_hrv = zeros(1,NFEAT_hrv);
        end
    else
        disp('Skipping HRV analysis due to shortage of peaks..')
        feat_hrv = zeros(1,NFEAT_hrv);
    end

    % Heart Rate features
    HRbpm = median(60./(diff(qrsseg{2})));
    %obvious cases: tachycardia ( > 100 beats per minute (bpm) in adults)
    feat_tachy = normcdf(HRbpm,120,20); % sampling from normal CDF
    %See e.g.   x = 10:10:200; p = normcdf(x,120,20); plot(x,p)

    %obvious cases: bradycardia ( < 60 bpm in adults)
    feat_brady = 1-normcdf(HRbpm,60,20);

    % SQI metrics
%         feats_sqi = ecgsqi(sig_seg,qrsseg,fs);
    feats_sqi=[];
    % Features on residual
%         featsres = residualfeats(sig_segraw,fs,qrsseg{end});
    featsres=[];
    % Morphological features
%         feats_morph = morphofeatures(sig_segraw,fs,qrsseg,[fname '_s' num2str(n)]);
    feats_morph=[];

    feat_fer=[featqrs,feat_tachy,feat_brady,double(feats_sqi),featsres,feats_morph];
    feat_fer(~isreal(feat_fer)|isnan(feat_fer)|isinf(feat_fer)) = 0; % removing not numbers

    % Save features to table for training
    feats = [feats;feat_hrv,feat_fer];
%     fetbag{n} = [hyp(n),feats];
end
names = {'ECG_IrrIndex' 'ECG_IrrEv' 'ECG_OriginCount' 'ECG_PACEv' 'ECG_meanRR' 'ECG_medianRR' 'ECG_SDNN' 'ECG_RMSSD' 'ECG_SDSD' 'ECG_NN50' 'ECG_pNN50' 'ECG_LFpeak' 'ECG_HFpeak' 'ECG_totalpower' 'ECG_LFpower' ...
    'ECG_HFpower' 'ECG_nLF' 'ECG_nHF' 'ECG_LFHF' 'ECG_PoincareSD1' 'ECG_PoincareSD2' 'ECG_SampEn' 'ECG_ApEn' ...
    'ECG_RR' 'ECG_DET' 'ECG_ENTR' 'ECG_L' 'ECG_TKEO1'  'ECG_DAFa2' 'ECG_LZ' 'ECG_BD' 'ECG_PD' 'ECG_BDa' 'ECG_PDa' 'ECG_ZCI' 'ECG_mZCI' 'ECG_nsZCI'...
    'ECG_meanRR_norm' 'ECG_medianRR_norm' 'ECG_SDNN_norm' 'ECG_LFpeak_norm' 'ECG_HFpeak_norm' 'ECG_totalpower_norm' 'ECG_LFpower_norm' 'ECG_HFpower_norm'...
    'ECG_LFHF_norm' 'ECG_PoincareSD1_norm' 'ECG_PoincareSD2_norm' 'ECG_SampEn_norm' 'ECG_ApEn_norm' 'ECG_TKEO1_norm'...
    'ECG_Clvl1' 'ECG_Clvl2' 'ECG_Clvl3' 'ECG_Clvl4' 'ECG_Clvl5' 'ECG_Clvl6' 'ECG_Clvl7' 'ECG_Clvl8' 'ECG_Clvl9' ...
    'ECG_Clvl10' 'ECG_Dlvl1' 'ECG_Dlvl2' 'ECG_Dlvl3' 'ECG_Dlvl4' ...
    'ECG_Dlvl5' 'ECG_Dlvl6' 'ECG_Dlvl7' 'ECG_Dlvl8' 'ECG_Dlvl9' 'ECG_Dlvl10'};
names = [names 'ECG_amp_varsqi' 'ECG_amp_stdsqi' 'ECG_amp_mean'];
names = [names 'ECG_tachy' 'ECG_brady'];


allfeats = array2table(feats, 'VariableNames',names);
features_struct = table2struct(allfeats);
hr_data = reshape(signal,numel(signal),1);
% rmpath('C:\Users\scro2778\Documents\CDT\Imaging\Practical 2 - Registration\Registration code 2015');

    
end
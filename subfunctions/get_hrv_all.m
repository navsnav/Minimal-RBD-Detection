function HRV=get_hrv_all(hrv_now)
% Loads QRS peak detector location file and calculates all heart rate
% variability features (HRV).
%
% Time domain features:
%               SDNN: Standard deviation of all NN intervals
%              SDANN: Standard deviation of mean of NN intervals in 5 min
%                     windows
%              RMSSD: Square root of mean of squares of differences between
%                     adjacent NN intervals
%          SDNNindex: mean of standard deviation of all NN intervals in all
%                     5 mins windows
%               SDSD: Standard deviation of differences between adjacent NN
%                     intervals
%               NN50: Number of pairs of adjacent NN intervals differing by
%                     more than 50ms
%              pNN50: Percentage NN50 (NN50/totan number of NN intervals)
%
% Frequency domain features
%           VLF_peak: Frequency of maximum peak in very low frequency range
%            LF_peak: Frequency of maximum peak in low frequency range
%            HF_peak: Frequency of maximum peak in high frequency range
%     VLF_power_perc: Percentage of total power in VLF range
%      LF_power_perc: Percentage of total power in LF range
%      HF_power_perc: Percentage of total power in HF range
%      LF_power_norm: Percentage of low and high frequency power in the low
%                     frequency range
%      HF_power_norm: Percentage of low and high frequency power in the
%                     high frequency range
%        LF_HF_ratio: Ratio of low to high frequency power
%
% Non-linear features
%                SD1: Standard deviation in y=-x direction of Poincare plot
%                SD2: Standard deviation in y=x direction of Poincare plot
%               saen: Sample Entropy
%               apen: Approximate Entropy
%
%            Recurrence Quantification Analysis 
%                 RR: Recurrence Rate
%                Det: Determinism
%               ENTR: Shannon Entropy
%                  L: Average diagonal line length
%
%               TKEO: Mean Teager-Kaiser Energy Operator%   
%             DFA_a2: Detrended Fluctuation Analysis Exponent%                 
%                 LZ: Lempel Ziv Complexity
%
%
% --
% ECG classification from single-lead segments using Deep Convolutional Neural 
% Networks and Feature-Based Approaches - December 2017
% 
% Released under the GNU General Public License
%
% Copyright (C) 2017  Fernando Andreotti, Oliver Carr
% University of Oxford, Insitute of Biomedical Engineering, CIBIM Lab - Oxford 2017
% fernando.andreotti@eng.ox.ac.uk
%
% 
% For more information visit: https://github.com/fernandoandreotti/cinc-challenge2017
% 
% Referencing this work
%
% Andreotti, F., Carr, O., Pimentel, M.A.F., Mahdi, A., & De Vos, M. (2017). 
% Comparing Feature Based Classifiers and Convolutional Neural Networks to Detect 
% Arrhythmia from Short Segments of ECG. In Computing in Cardiology. Rennes (France).
%
% Last updated : December 2017
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


% load QRS time data

vlow_thresh=0.003;
low_thresh=0.04; % set frequency ranges for very low, low and high ranges (Hz)
mid_thresh=0.15;
high_thresh=0.4;

HRV=struct; % create HRV structure

%%

NN=1000*(hrv_now(:,1)-hrv_now(1,1));  % convert data points to time (ms)

NN_int_sec=1000*hrv_now(:,2); % find NN intervals (ms)

NN_mat=[NN,NN_int_sec]; %(ms)

HRV.mRR=mean(NN_int_sec);

HRV.medRR=median(NN_int_sec);

HRV.SDNN=std(NN_int_sec); 

D=diff(NN_int_sec);

D_sq=D.^2;

N=length(NN_int_sec);

mean_diff_int=sum((D_sq))/(N-1);

HRV.RMSSD=sqrt(mean_diff_int);

HRV.SDSD=std(D);
  
HRV.NN50=sum(D>0.05);   

HRV.pNN50=HRV.NN50/length(NN_int_sec);



end







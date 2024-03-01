function [TRsec, omega, vmax, vmin, vsig_wide] = FR_Preprocess(tseries)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: timecourse
% Output: 
% 1. power omega
% 2. gaussian optimization phase: max and min and its orginal timeseries -- vmax, vmin
% 3. gaussian optimization: its orginal timeseries -- vsig_wide

%% Frequency analysis
%--------------------------------------------------------------------------
%COMPUTE POWER SPECTRA FOR
%NARROWLY FILTERED DATA WITH LOW BANDPASS (0.04 to 0.07 Hz)
%WIDELY FILTERED DATA (0.04 Hz to justBelowNyquistFrequency)
%[justBelowNyquistFrequency depends on TR,
%for a TR of 2s this is 0.249 Hz]
%--------------------------------------------------------------------------
nNodes = size(tseries{1},1);%Number of parcels
TT=size(tseries{1},2);%The length of the timecourse
Subs = length(tseries);

band_pass = [0.04 0.07];
TRsec = 2;%for a TR of 2s this is 0.249 Hz]
Ts = TT*TRsec;
freq = (0:TT/2-1)/Ts;  %dTR
[~, idxMinFreq] = min(abs(freq-band_pass(1)));
[~, idxMaxFreq] = min(abs(freq-band_pass(2)));
nFreqs = length(freq);
delt = TRsec;                               % sampling interval: same as TR
fnq = 1/(2*delt);                           % Nyquist frequency

%filter construction: WIDE BANDPASS
k_order_filter = 2;                         % 2nd order butterworth filter
flp = band_pass(1);                         % lowpass frequency of filter,orignal 0.04Hz
fhi = fnq-0.001;%.249;                      % highpass needs to be limited by Nyquist frequency, which in turn depends on TR
Wn = [flp/fnq fhi/fnq];                     % butterworth bandpass non-dimensional frequency
[bfilt_wide, afilt_wide] = butter(k_order_filter,Wn);    % construct the filter
clear k_order_filter Wn fhi flp

%filter construction: NARROW LOW BANDPASS
k_order_filter = 2;                         % 2nd order butterworth filter
flp = band_pass(1);                         % lowpass frequency of filter
fhi = band_pass(2);                         % highpass
Wn=[flp/fnq fhi/fnq];                       % butterworth bandpass non-dimensional frequency
[bfilt_narrow,afilt_narrow] = butter(k_order_filter,Wn); % construct the filter
clear k_order_filter Wn fhi flp

% Frequency analysis
PowSpect_filt_narrow = zeros(nFreqs, nNodes, Subs);
PowSpect_filt_wide = zeros(nFreqs, nNodes, Subs);
for seed=1:nNodes
    for idxSub=1:Subs
        signaldata = tseries{idxSub};
        x=detrend(demean(signaldata(seed,:)));
        
        %Revise here since they are not filter signal
        ts_filt_narrow =zscore(filtfilt(bfilt_narrow,afilt_narrow,x));
        pw_filt_narrow = abs(fft(ts_filt_narrow));
        PowSpect_filt_narrow(:,seed,idxSub) = pw_filt_narrow(1:floor(TT/2)).^2/(TT/2);%dTR
        
        ts_filt_wide =zscore(filtfilt(bfilt_wide,afilt_wide,x));
        pw_filt_wide = abs(fft(ts_filt_wide));
        PowSpect_filt_wide(:,seed,idxSub) = pw_filt_wide(1:floor(TT/2)).^2/(TT/2);%dTR
        
        clear signaldata x ts_filt_narrow ts_filt_wide pw_filt_narrow pw_filt_wide
        
    end
    clear idxSub
end
clear seed
clear bfilt_narrow afilt_narrow bfilt_wide afilt_wide

% FC_emp_update_narrow = FC_matrix_update(tseries_update_filt_narrow);%Making functional matrix
% FC_emp_update_wide = FC_matrix_update(tseries_update_filt_wide);%Making functional matrix

Power_Areas_filt_narrow_unsmoothed = mean(PowSpect_filt_narrow,3);clear PowSpect_filt_narrow
Power_Areas_filt_wide_unsmoothed = mean(PowSpect_filt_wide,3);clear PowSpect_filt_wide

%Frequency optimization:gaussian optimization
Power_Areas_filt_narrow_smoothed = zeros(nFreqs, nNodes);
Power_Areas_filt_wide_smoothed = zeros(nFreqs, nNodes);
vsig = zeros(1, nNodes);
for seed=1:nNodes
    Power_Areas_filt_narrow_smoothed(:,seed)=gaussfilt(freq,Power_Areas_filt_narrow_unsmoothed(:,seed)',0.01);
    Power_Areas_filt_wide_smoothed(:,seed)=gaussfilt(freq,Power_Areas_filt_wide_unsmoothed(:,seed)',0.01);
    %relative power in frequencies of interest (band_pass(1) - band_pass(2) Hz) with respect
    %to entire power of bandpass-filtered data (band_pass(1) - just_below_nyquist)
    vsig_wide(seed) =...
        sum(Power_Areas_filt_wide_smoothed(idxMinFreq:idxMaxFreq,seed))/sum(Power_Areas_filt_wide_smoothed(:,seed));
    vsig_narrow(seed) =...
        sum(Power_Areas_filt_narrow_smoothed(idxMinFreq:idxMaxFreq,seed))/sum(Power_Areas_filt_narrow_smoothed(:,seed));
    
end
clear seed Power_Areas_filt_narrow_unsmoothed Power_Areas_filt_wide_unsmoothed
vmax=max(vsig_wide); %consider computing this later where needed
vmin=min(vsig_wide);%consider computing this later where needed

%a-minimization seems to only work if we use the indices for frequency of
%maximal power from the narrowband-smoothed data
[~, idxFreqOfMaxPwr]=max(Power_Areas_filt_narrow_smoothed);
f_diff = freq(idxFreqOfMaxPwr);

%f_diff  previously computed frequency with maximal power (of narrowly filtered data) by area
omega = repmat(2*pi*f_diff',1,2); %angular velocity
omega(:,1) = -omega(:,1);


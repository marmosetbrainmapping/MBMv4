% data_making

clear all
close all

%% fMRI data timecourse

path = 'D:\modeling\database\All_No_smooth_MBMV4_merge\';
dirOutput_file=dir(fullfile([path,'*.mat']));

% path = 'D:\modeling\database\All_Smooth_Paxinos\';
% dirOutput_file=dir(fullfile(['D:\modeling\database\All_Smooth_Paxinos\','*.mat']));
num_sub = 1;
for k = 1:size(dirOutput_file,1)
    load([path,dirOutput_file(k).name]);data = data_pool_merge;clear data_pool_merge
    tmp = strfind(dirOutput_file(k).name,'_');
    name{k,:} = dirOutput_file(k).name(1:tmp(1)-1);
    if k > 1
        if strcmp(name{k,:},name{k-1,:}) == 1
            data_pool = [data_pool;data];
            name_session = [name_session;[num_sub,k,size(data,1)]];
        else
            tseries_subject{num_sub} = data_pool;
            name_subject{num_sub} = name_session;
            num_sub = num_sub+1;
            data_pool = [];data_pool = [data_pool;data];
            name_session = [];name_session = [num_sub,k,size(data,1)];
            name_update{num_sub} = name{k,:}; 
        end
    else
        data_pool = data;
        name_session = [num_sub,k,size(data,1)];
        name_update{num_sub} = name{k,:}; 
    end
    clear data tmp
end
clear k data_pool name_session

%Number of subject
nSubs = num_sub-1;clear num_sub

%For every subject: logan
for k = 16:nSubs
    % for k = 9:nSubs
    
    if k ~=5
        name_temp = name{name_subject{k}(1,2)};
        
        %% Construction
        if k == 23 
            SC_conn_v4 = table2array(readtable(['D:\modeling\database\Tracking_All\',name_temp,'_T1T2DTI_postsurgery_invivoDTI_v4_0622.csv']));%High resolution connectivity
        else
            SC_conn_v4 = table2array(readtable(['D:\modeling\database\Tracking_All\',name_temp,'_T1T2DTI_invivoDTI_v4_0622.csv']));%High resolution connectivity
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Left = SC_conn_v4(:,1:126);Right = SC_conn_v4(:,127:end);
        Left_update = area_merge(Left);clear Left
        Right_update = area_merge(Right);clear Right
        temp = [Left_update,Right_update];clear Left_update Right_update
        temp = temp';
        Left = temp(:,1:126);Right = temp(:,127:end);
        Left_update = area_merge(Left);clear Left
        Right_update = area_merge(Right);clear Right
        SC_conn_v4_merge = [Left_update,Right_update];clear Left_update Right_update
        SC_conn_v4_merge = SC_conn_v4_merge';clear temp
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        SC_conn_v4 = SC_matrix(SC_conn_v4_merge);
        SC_conn_scale = SC_conn_v4;
        
        nNodes = size(SC_conn_scale,1);
        
        %     [tseries_update_paxinos,FC_emp,Tmax_paxinos,Subs] = FC_matrix(tseries_subject{k},Valid_paxinos,nNodes,name_subject{k});
        [tseries_update_paxinos,FC_emp,Tmax,Subs] = FC_matrix(tseries_subject{k},[],nNodes,name_subject{k});
        
        %% Frequency analysis
        %--------------------------------------------------------------------------
        %COMPUTE POWER SPECTRA FOR
        %NARROWLY FILTERED DATA WITH LOW BANDPASS (0.04 to 0.07 Hz)
        %WIDELY FILTERED DATA (0.04 Hz to justBelowNyquistFrequency)
        %[justBelowNyquistFrequency depends on TR,
        %for a TR of 2s this is 0.249 Hz]
        %--------------------------------------------------------------------------
        TT=Tmax;%The length of the timecourse
        band_pass = [0.04 0.07];
        TRsec = 2;%for a TR of 2s this is 0.249 Hz]
        Ts = TT*TRsec;
        freq = (0:TT/2-1)/Ts;  %dTR
        [~, idxMinFreq] = min(abs(freq-band_pass(1)));
        [~, idxMaxFreq] = min(abs(freq-band_pass(2)));
        nFreqs = length(freq);
        delt = TRsec;                               % sampling interval: same as TR
        fnq = 1/(2*delt);                           % Nyquist frequency
        
        %WIDE BANDPASS
        k_order_filter = 2;                         % 2nd order butterworth filter
        flp = band_pass(1);                         % lowpass frequency of filter,orignal 0.04Hz
        fhi = fnq-0.001;%.249;                      % highpass needs to be limited by Nyquist frequency, which in turn depends on TR
        Wn = [flp/fnq fhi/fnq];                     % butterworth bandpass non-dimensional frequency
        [bfilt_wide, afilt_wide] = butter(k_order_filter,Wn);    % construct the filter
        clear k_order_filter Wn fhi flp
        
        %NARROW LOW BANDPASS
        k_order_filter = 2;                         % 2nd order butterworth filter
        flp = band_pass(1);                         % lowpass frequency of filter
        fhi = band_pass(2);                         % highpass
        Wn=[flp/fnq fhi/fnq];                       % butterworth bandpass non-dimensional frequency
        [bfilt_narrow,afilt_narrow] = butter(k_order_filter,Wn); % construct the filter
        clear k_order_filter Wn fhi flp
        
        
        PowSpect_filt_narrow = zeros(nFreqs, nNodes, Subs);
        PowSpect_filt_wide = zeros(nFreqs, nNodes, Subs);
        for seed=1:nNodes
            for idxSub=1:Subs
                signaldata = tseries_update_paxinos{idxSub};
                x=detrend(demean(signaldata(seed,:)));
                
                %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             figure,subplot(2,3,1),
                %             plot(signaldata(seed,:));
                %             axis square
                %             subplot(2,3,4),
                %             plot(x);
                %             axis square
                %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Revise here since they are not filter signal
                ts_filt_narrow =zscore(filtfilt(bfilt_narrow,afilt_narrow,x));
                pw_filt_narrow = abs(fft(ts_filt_narrow));
                
                tseries_update_filt_narrow{idxSub}(seed,:) = ts_filt_narrow;
                
                %FOR EACH AREA AND TIMEPOINT COMPUTE THE INSTANTANEOUS PHASE IN THE NARROW BAND
                xFilt = filtfilt(bfilt_narrow,afilt_narrow,x);    % zero phase filter the data
                Xanalytic = hilbert(demean(xFilt));
                PhasesD(seed,:,idxSub) = angle(Xanalytic);
                
                
                clear Xanalytic xFilt
                
                %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             subplot(2,3,2),
                %             plot(ts_filt_narrow);
                %             axis square
                %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                PowSpect_filt_narrow(:,seed,idxSub) = pw_filt_narrow(1:floor(TT/2)).^2/(TT/2);  %dTR
                
                %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             subplot(2,3,5),
                %             plot(1/TRsec*(0:(TT/2))/TT,[pw_filt_narrow(1:floor(TT/2)).^2/(TT/2),0]);
                %             xlabel('Hz');
                %             axis square
                %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                ts_filt_wide =zscore(filtfilt(bfilt_wide,afilt_wide,x));
                tseries_update_filt_wide{idxSub}(seed,:) = ts_filt_wide;
                
                %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             subplot(2,3,3),
                %             plot(ts_filt_wide);
                %             axis square
                %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                pw_filt_wide = abs(fft(ts_filt_wide));
                PowSpect_filt_wide(:,seed,idxSub) = pw_filt_wide(1:floor(TT/2)).^2/(TT/2);   %dTR
                %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             subplot(2,3,6),
                %             plot(1/TRsec*(0:(TT/2))/TT,[pw_filt_wide(1:floor(TT/2)).^2/(TT/2),0]);
                %             xlabel('Hz');
                %             axis square
                %             close all
                %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                clear signaldata x ts_filt_narrow ts_filt_wide pw_filt_narrow pw_filt_wide
                
            end
            clear idxSub
        end
        clear seed
        clear bfilt_narrow afilt_narrow bfilt_wide afilt_wide
        
        FC_emp_update_narrow = FC_matrix_update(tseries_update_filt_narrow);%Making functional matrix
        FC_emp_update_wide = FC_matrix_update(tseries_update_filt_wide);%Making functional matrix
        
        Power_Areas_filt_narrow_unsmoothed = mean(PowSpect_filt_narrow,3);clear PowSpect_filt_narrow
        Power_Areas_filt_wide_unsmoothed = mean(PowSpect_filt_wide,3);clear PowSpect_filt_wide
        Power_Areas_filt_narrow_smoothed = zeros(nFreqs, nNodes);
        Power_Areas_filt_wide_smoothed = zeros(nFreqs, nNodes);
        vsig = zeros(1, nNodes);
        %Frequency optimization:gaussian optimization
        for seed=1:nNodes
            Power_Areas_filt_narrow_smoothed(:,seed)=gaussfilt(freq,Power_Areas_filt_narrow_unsmoothed(:,seed)',0.01);
            Power_Areas_filt_wide_smoothed(:,seed)=gaussfilt(freq,Power_Areas_filt_wide_unsmoothed(:,seed)',0.01);
            %relative power in frequencies of interest (band_pass(1) - band_pass(2) Hz) with respect
            %to entire power of bandpass-filtered data (band_pass(1) - just_below_nyquist)
%             vsig(seed) =...
%                 sum(Power_Areas_filt_wide_smoothed(idxMinFreq:idxMaxFreq,seed))/sum(Power_Areas_filt_wide_smoothed(:,seed));
            vsig(seed) =...
                sum(Power_Areas_filt_narrow_smoothed(idxMinFreq:idxMaxFreq,seed))/sum(Power_Areas_filt_narrow_smoothed(:,seed));
            
        end
        clear seed Power_Areas_filt_narrow_unsmoothed Power_Areas_filt_wide_unsmoothed
        vmax=max(vsig); %consider computing this later where needed
        vmin=min(vsig);%consider computing this later where needed
        
        %a-minimization seems to only work if we use the indices for frequency of
        %maximal power from the narrowband-smoothed data
        [~, idxFreqOfMaxPwr]=max(Power_Areas_filt_narrow_smoothed);
        f_diff = freq(idxFreqOfMaxPwr);
        
        %f_diff  previously computed frequency with maximal power (of narrowly filtered data) by area
        omega = repmat(2*pi*f_diff',1,2); %angular velocity
        omega(:,1) = -omega(:,1);
        
        
        %% Modeling
        %FROM HERE ON SIMULATIONS AND FITTING
        av = -0.5:0.01:0.5; %The abifurcation parameter (a)
        wG = 0:0.1:10;  %G value or vector: The coupling parameter vector. Change the range and step length as needed.
        nWeights = numel(wG);
        % for i = 1:length(av)
        parfor j = 1:length(wG)
            %1: Optimization by the gene algorithm
            [trackminm_op(:,j),best_vsigs_op(j,:),bestIter_op(j,:),a_op_bifurcation{j}] = BIFURCATION_opt_params(SC_conn_scale, Tmax, wG(j), Subs, nNodes, TRsec, omega, vmax, vmin, vsig);
            %         a_no_op_bifurcation(j,:) = av(i);%No Optimization
            a_no_op_bifurcation(j,:) = -0.5;%No Optimization
            bifpar(j,:) = a_op_bifurcation{j}(:,1)';%store them in an output variable
            %2: Simulation and fitting examine
            [simulated_timecourse{j},index_fitting(j,:)] = hbif_model_fitting(FC_emp, SC_conn_scale, Tmax, wG(j), Subs, nNodes, a_op_bifurcation{j}, a_no_op_bifurcation(j,:), TRsec, omega);
            %3: metastability
            %         [index_ksP(j,:),pcS{j},index_meta(j,:),index_cohe(j,:)] = hbif_model_metastability(simulated_timecourse{j}, Tmax, nSubs, nNodes, TRsec, PhasesD, band_pass);
        end
        clear j
        % end
        % clear i
        k
        save([name_temp,'_simulation_V4_merge.mat'],'-v7.3');
        
        clearvars -except k nSubs dirOutput_file name name_subject tseries_subject Valid_paxinos
    end
end
clear k

%% Functions
function SC_update = SC_matrix(SC)%Making structural matrix
% SC_update = SC+SC';SC_update=SC_update/2;
SC_update=SC/2;
%--------------------------------------------------------------------------
%SCALE STRUCTURAL CONNECTIVITY MATRIX
%--------------------------------------------------------------------------
SC_update = SC_update/max(max(SC_update))*0.2;%
end

function FC_update = FC_matrix_update(tseries)%Making structural matrix
for k = 1:length(tseries)
    %--------------------------------------------------------------------------
    %CALCULATE FUNCTIONAL CONNECTIVITY MATRIX
    %--------------------------------------------------------------------------
    r(:,:,k) = corrcoef(tseries{k}');
end
FC_update=mean(r,3);
clear k r
end

function [tseries_update,FC_emp,Tmax,Subs] = FC_matrix(tseries,Valid_paxinos,nNodes,Subj_list)%Making functional matrix

Tmax = 200;
% Tmax = 500;
% Tmax = min(Subj_list(:,3));

% Subs = floor(size(tseries,1)/Tmax);
Subs = size(Subj_list,1);

% targ = (1:Subs)*Tmax;targ = targ';
% targ = [[1;targ(1:end-1)+1],targ];

targ = cumsum(Subj_list(:,3));
targ = [[1;targ(1:end-1)+1],targ];

for k = 1:Subs
    tseries_update{k} = tseries(targ(k,1):targ(k,2),:);
    tseries_update{k} = tseries_update{k}(1:Tmax,:);
    tseries_update{k} = tseries_update{k}';
    if isempty(Valid_paxinos) ~= 1
        tseries_update{k} = tseries_update{k}(Valid_paxinos,:);
    end
    %--------------------------------------------------------------------------
    %CALCULATE FUNCTIONAL CONNECTIVITY MATRIX
    %--------------------------------------------------------------------------
    r(:,:,k) = corrcoef(tseries_update{k}');
end
FC_emp=mean(r,3);
clear k r
end

function data_update = area_merge(data)
combine = [3,123,120,0,0;
    4,114,0,0,0;
    32,6,110,0,0;
    30,11,0,0,0;
    12,125,0,0,0;
    16,56,60,0,0;
    17,28,0,0,0;
    18,19,0,0,0;
    23,44,119,0,0;
    24,84,0,0,0;
    26,33,0,0,0;
    27,74,0,0,0;
    36,38,0,0,0;
    37,98,121,0,0;
    47,50,95,0,0;
    49,77,97,63,80;
    51,102,0,0,0;
    52,82,0,0,0;
    53,69,0,0,0;
    54,93,0,0,0;
    55,78,0,0,0];

for k = 1:size(combine,1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp = nonzeros( combine(k,:) );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    data_update(:,k) = mean(data(:,temp),2);
    
    clear temp
end
Num = k;
clear k

merge = unique(combine(:));merge(1) = [];
difference = setdiff(1:126,merge);

for k = 1:length(difference)
    data_update(:,Num+k) = data(:,difference(k));
end
clear k
end
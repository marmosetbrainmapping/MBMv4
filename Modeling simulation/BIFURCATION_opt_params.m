% OPTIMIZATION OF BIFURCATION PARAMETERS

function [trackminm1,best_vsigs,bestIter,a1] = BIFURCATION_opt_params(C, Tmax, we, nSubs, nNodes, TR, omega, vmax, vmin, vsig)
%Cfg for optimization (change nIters if we don´t have the minimum error)
Cfg.simulID = 'MortenGroup_v3wide';
Cfg.opt_a.nIters = 100;  %n iterations 300
Cfg.opt_a.updateStrength = 0.1;
Cfg.opt_a.abortCrit = 0.1; % abortion criteria = maximally 2.5% error
Cfg.opt_a.gref = 3; %no less than 2
Cfg.plots.showOptimization = 1;

%filter
band_pass = [0.04 0.07];
delt = TR;                               % sampling interval: same as TR
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

%Timecourse
dt=0.1.*(TR/2);  %BEFORE dt = 0.1;
sig = 0.04; %was 0.04
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step

a = repmat(-0.05*ones(nNodes,1),1,2);
a1 = a; %JVS we start with the best a (a1) being a
trackminm1 = zeros(Cfg.opt_a.nIters, 1); %for tracking the minimization (good for debugging)


xs = zeros(Tmax*nSubs,nNodes);  % BEFORE xs = zeros(3000/2,nNodes);
wC = we*C; %structural connectivity matrix weighted with current global coupling parameter
sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj %HOW SHOULD WE TREAT ASYMMETRIC CONNECTIVITY MATRICES


fprintf(1, 'STARTING OPTIMIZATION OF BIFURCATION PARAMETERS\n');
minm=100;
bestIter = 1; %JVS for tracking purposes
for iter = 1:Cfg.opt_a.nIters
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    z = 0.1*ones(nNodes,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    for t=1:dt:1000
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
    end
    
    for t=1:dt:Tmax*nSubs*TR %BEFORE 3000 seconds worth of simulated data; why dont we limit it to Tmax*TRsec?
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
        
        if (abs(mod(t,TR))<0.01) %BEFORE mod(t,2)==0
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    vsigs = zeros(1, nNodes);
    for seed=1:nNodes
        
        x=detrend(demean(xs(1:nn,seed)'));%was x
        ts_filt_wide =zscore(filtfilt(bfilt_wide,afilt_wide,x));%zscore before!
        
        TT=length(x);
        Ts = TT*TR;
        freq = (0:TT/2-1)/Ts;
        [~, idxMinFreqS]=min(abs(freq-band_pass(1)));
        [~, idxMaxFreqS]=min(abs(freq-band_pass(2)));
        
        pw_filt_wide = abs(fft(ts_filt_wide));
        Pow1 = pw_filt_wide(1:floor(TT/2)).^2/(TT/2);
        Pow=gaussfilt(freq,Pow1,0.01);
        
        vsigs(seed)=sum(Pow(idxMinFreqS:idxMaxFreqS))/sum(Pow);
        
    end
    
    vsmin=min(vsigs);%vsig y vsigs max min points not equal
    vsmax=max(vsigs);
    bb=(vmax-vmin)/(vsmax-vsmin);
    aa=vmin-bb*vsmin;  %% adaptation of local bif parameters
    vsigs=aa+bb*vsigs;
    minm1=max(abs(vsig-vsigs)./vsig);
    trackminm1(iter,:) = minm1;%JVS tracking of minm1

    if minm1<minm
        minm=minm1;
        a1=a;
        bestIter = iter; %JVS: tracking
        best_vsigs = vsigs; %JVS: tracking
    end
     
    %CRITERION REACHED?
    if minm<Cfg.opt_a.abortCrit %default is 0.1
        break;
    end
    
    %UPDATE a VALUES FOR NEXT ITER
    if ~any(isnan(vsigs))
        a(:,1)=a(:,1)+Cfg.opt_a.updateStrength*(1-vsigs./vsig)';
        a(:,2)=a(:,2)+Cfg.opt_a.updateStrength*(1-vsigs./vsig)';
    else
        %JVS: this should not happen according to Gustavo, but it
        %can when a values run crazy
        warning('There are NaNs in the power spectra. Probably a-values too strong.');
        a = a1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------
%FEEDBACK
%--------------------------------------------------------------
if Cfg.plots.showOptimization
    idx_g = 1;
    showOptimPlot(idx_g, we, iter, a, a1, vsig, vsigs, bestIter, best_vsigs, trackminm1, Cfg)
end
fprintf(1, 'iter: %03d, minm1: %5.3f\n', iter, minm1);
%--------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function showOptimPlot(idx_g, we, iter, a, a1, vsig, vsigs, bestIter, best_vsigs, trackminm1, Cfg)
figure,
subplot(4,1,1)
ph_a = plot(a(:, 1), 'k', 'LineWidth', 2);
hold on
ph_a_best = plot(a1(:, 1), 'r', 'LineWidth', 2);

plot([0, length(vsigs)], [0 0], ':k')
hold off
set(gca, 'xlim', [0, length(vsigs)+1], 'ylim', [-3, 1])
ylabel('a')
xlabel('ROI')
box off
title('a')
legend([ph_a, ph_a_best], {'current a',  'best a'}, 'Location', 'NorthEast')

subplot(4,1,2)
ph = plot([vsig(:), vsigs(:), best_vsigs(:)], 'LineWidth', 1);
set(ph([1,3]), 'LineWidth', 3);
legend('vsig (data)', 'current vsigs (simulated)', 'best fitting vsigs (simulated)', 'Location', 'SouthEast')
set(gca, 'xlim', [0, length(vsigs)])
xlabel('ROI')
ylabel('relative power')
title('vsig and vsigs')
box off

subplot(4,1,3)
bar(abs(vsig-vsigs)./vsig)
hold on
plot([0, length(vsigs)-1], [Cfg.opt_a.abortCrit Cfg.opt_a.abortCrit], 'r')
hold off
set(gca, 'xlim', [0,  length(vsigs)], 'ylim', [0, Cfg.opt_a.abortCrit]*2)
xlabel('ROI')
ylabel('\delta power')
title('scaled difference of powers: abs(vsig-vsigs)./vsig')
box off

subplot(4,1,4)
ph = plot(trackminm1(:, idx_g));
xlabel('iter')
set(gca, 'xlim', [0, Cfg.opt_a.nIters+1], 'ylim', [0, Cfg.opt_a.abortCrit]*2);
hold on
plot([1, Cfg.opt_a.nIters], [Cfg.opt_a.abortCrit Cfg.opt_a.abortCrit], 'r')
plot(bestIter, trackminm1(bestIter, idx_g), 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'r')
hold off
legend(ph, 'max(abs(vsig-vsigs)./vsig)', 'location', 'NorthEast')
ylabel('max \delta of power')
box off

set(gcf, 'name', sprintf('G = %5.3f. Optimization of bifurcation parameters (iter %3d)', we, iter))
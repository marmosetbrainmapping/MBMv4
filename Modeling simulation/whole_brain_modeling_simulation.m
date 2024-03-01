%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The modeling is mainly referred Deco, G. et al. (2017). The dynamics of resting fluctuations in the brain: metastability and its dynamical cortical core. Sci Rep 7, 3095.
% If there is any question, just contact me without hesitation txgxp88@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [wG,bifpar,simulated_timecourse,FC_set1,FC_set2,index_fitting] = whole_brain_modeling_simulation(file_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load example dataset by the file name:
%%% Functional timecourse
% FC_set1: example dataset 1
% 1) Resting-state timecourse in each cell can be represented by a run/session, 
% 2) the matrix of Resting-state timecourse in each cell represents regions x time ,here we select 10min record. 

% FC_set2: example dataset 2
% same as above but for verification

%%% Structural connectivity
% SC_conn: the matrix of whole brain connectivity (Left+Right hemisphere)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
currentFolder = pwd;
load([currentFolder,'/Example_data/',file_name,'.mat'],'FC_set1','FC_set2','SC_conn');

%% structural data
SC_conn_scale = SC_matrix(SC_conn);%Normalized structural connectivity and scaled 20%

%% functional data
%Preprocess: calculate the statistical FC
[FC_emp_test,FC_emp_verify] = FC_matrix(FC_set1,FC_set2);

%% Necessary information
nNodes = size(SC_conn_scale,1);%Number of parcels
Sessions = length(FC_set1);%Number of runs/sessions
Tmax = size(FC_set1{1},1);%Tmax: the length of timecourse

%% Step 2: Frequency power preprocess
%Input: timecourse
[TRsec, omega, vmax, vmin, vsig] = FR_Preprocess(FC_set1);

%% Step 3: Modeling
%FROM HERE ON SIMULATIONS AND FITTING
wG = 0:0.1:8;  %G value or vector: The coupling parameter vector. ### You can change the range and step length.
nWeights = numel(wG);
parfor j = 1:length(wG)
    %1: Optimization of abifurcation param (a) by the gene algorithm
    [trackminm_op(:,j),best_vsigs_op(j,:),bestIter_op(j,:),a_op_bifurcation{j}] = BIFURCATION_opt_params(SC_conn_scale, Tmax, wG(j), Sessions, nNodes, TRsec, omega, vmax, vmin, vsig);
    bifpar(j,:) = a_op_bifurcation{j}(:,1)';%store them in an output variable
    %2: Simulation and fitting examine
    [simulated_timecourse{j},index_fitting(j,:)] = hbif_model_fitting(FC_emp_verify, SC_conn_scale, Tmax, wG(j), Sessions, nNodes, a_op_bifurcation{j}, TRsec, omega);
end
clear j

%% Functions
function SC_update = SC_matrix(SC)%Making structural matrix
SC_update = SC+SC';SC_update=SC_update/2;
% SC_update=SC/2;
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

function [FC_emp1,FC_emp2] = FC_matrix(tseries1,tseries2)%Making functional matrix

for k = 1:length(tseries1)
    %--------------------------------------------------------------------------
    %CALCULATE FUNCTIONAL CONNECTIVITY MATRIX
    %--------------------------------------------------------------------------
    r1(:,:,k) = corrcoef(tseries1{k}');
    %--------------------------------------------------------------------------
    %CALCULATE FUNCTIONAL CONNECTIVITY MATRIX
    %--------------------------------------------------------------------------
    r2(:,:,k) = corrcoef(tseries2{k}');
 
end

FC_emp1=mean(r1,3);
FC_emp2=mean(r2,3);

clear k r
end

end
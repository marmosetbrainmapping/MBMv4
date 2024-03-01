%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The modeling is mainly referred Deco, G. et al. (2017). The dynamics of resting fluctuations in the brain: metastability and its dynamical cortical core. Sci Rep 7, 3095.
% If there is any question, just contact me without hesitation txgxp88@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Names of Example files in the folder
dataset_atlas_V4 = 'Example_data_v4';
dataset_atlas_Paxinos = 'Example_data_Paxinos';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sending to the modeling functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input data includes:
%FC_set1 for training: example data is the timecourses with 22 runs (cell) in different atlas, and 10min (row 300 with TR 2s) for each run;
%FC_set2 for verification: example data is the timecourses with 22 runs (cell) in different atlas, and 10min (row 300 with TR 2s) for each run;
%SC_conn: structural connectivity can be calculated by the invivo or ex-vivo DTI or tracing data

% returning results
%wG: manual defined global counpling parameter
%bifpar_: optimized alpha parameter for every parcel
%simulated_timecourse_: simulated timecourse, which is totally same as the input data
%training_timecourse_: empirical timecourse for optimized alpha parameter, which is the same as the input dataset1
%verification_timecourse_: empirical timecourse for verification, which is the same as the input dataset2
%index_fitting_: the pearson correlation coefficient between simulated and verification datasets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[wG,bifpar_v4,simulated_timecourse_V4,training_timecourse_V4,verification_timecourse_V4,index_fitting_V4] = whole_brain_modeling_simulation(dataset_atlas_V4);
[~,bifpar_paxinos,simulated_timecourse_paxinos,training_timecourse_paxinos,verification_timecourse_paxinos,index_fitting_paxinos] = whole_brain_modeling_simulation(dataset_atlas_Paxinos);

%Draw the figures and compare the performances
figure,
subplot(2,2,1),
hold on
plot(wG,index_fitting_V4,'k.-');
plot(wG,index_fitting_paxinos,'g.-');
legend('V4','Paxinos')
hold off
xlim([0 10]);
xlabel('wG (global coupling)');
ylim([0 1]);
box off
axis square

subplot(2,2,2),
hold on
histogram(bifpar_v4(:),'Normalization','probability');
histogram(bifpar_paxinos(:),'Normalization','probability');
legend('V4','Paxinos');
hold off
xlabel('optimized alpha');
ylabel('normalized occurrence');
box off
axis square
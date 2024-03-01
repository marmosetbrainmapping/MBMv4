clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: all analysis will be implemented in seperate hemisphere
% If there is any question, just contact me without hesitation email: txgxp88@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load surface information: medial wall and space coordinates
currentFolder = pwd;
% We prepare all necessary information from the surface template, such as medial wall etc.
load([currentFolder,'\Surface_template\surface_data.mat']);

%% Boundary map

% load files: Gradient
%We load gradient map which could be either group or individual gradient map
file_lh = ft_read_cifti([currentFolder,'/Preprocessed_data/Example_Group_grad_lh.dconn.nii']);%here we use group gradient as the example
gradient_lh = file_lh.dconn;clear file_lh

% file_rh = ft_read_cifti([currentFolder,'Preprocess_data/Example_Group_grad_rh.dconn.nii']);
% gradient_rh = file_rh.dconn;clear file_rh

%a large number for the medial wall
gradient_lh(median_lh,:) = 1000;
% gradient_rh(median_rh,:) = 1000;

%finding the minimum points from gradient, here we define 3 neighbors to
%look for the minumum values
minimametrics_left = metric_minima_all(gradient_lh,3,vertexNbors_lh);
% minimametrics_right = metric_minima_all(gradient_rh,3,vertexNbors_rh);

%Watershed algorithm to find the boundaries
%here we define a searching step with every 300 points in the step of 1 point according to their neighbor information
labels_left = watershed_algorithm_all_par(gradient_lh,minimametrics_left,300,1,vertexNbors_lh);
% labels_right = watershed_algorithm_all_par(gradient_rh,minimametrics_right,300,1,vertexNbors_rh);

% Raw boundary map for view:
labels_left_binary = binary_label(labels_left);
% labels_right_binary = binary_label(labels_right);
labels_left_avg_raw = mean(labels_left_binary,2);
% labels_right_avg_raw = mean(labels_right_binary,2);

% Save it for workbench view
save(gifti(single(labels_left_avg_raw)),'Boundary_map_Left_raw.shape.gii','Base64Binary');
% save(gifti(single(labels_right_avg_raw)),'Boundary_map_Right_raw.shape.gii','Base64Binary');
% 

% Smooth boundary map for view:
%Here we define 4 times for the optimization iterations and you can change depending on your data quality 
[labels_left_avg_raw_smooth] = spatial_smooth(labels_left_avg_raw,vertexNbors_lh);
% % [labels_right_avg_raw_smooth] = spatial_smooth(labels_right_avg_raw,vertexNbors_rh);

% Save it for workbench view
save(gifti(single(labels_left_avg_raw_smooth(:,4))),'Boundary_map_Left_smooth.shape.gii','Base64Binary');
% % save(gifti(single(labels_right_avg_raw_smooth(:,4))),'Boundary_map_Right_smooth.shape.gii','Base64Binary');

% save('boundary_map_lh.mat','labels_left_avg_raw','labels_left_avg_raw_smooth','-v7.3');
% save('boundary_map_rh.mat','labels_right_avg_raw','labels_right_avg_raw_smooth','-v7.3);


%% Parcel generation
% Left
[minimametric_left_parcels,label_left_parcels,label_left_parcels_no_edge] = parcel_generation(labels_left_avg_raw_smooth(:,4),labels_left,left_cortex_length,[(1:size(vertexNbors_lh,1))',vertexNbors_lh],lh_vertices);
% Minimum point
save(gifti(single(minimametric_left_parcels)),'left_parcels.shape.gii','Base64Binary');
% Parcel map with edge show
save(gifti(single(label_left_parcels)),'location_min_area_left_parcels.shape.gii','Base64Binary');
% Parcel map without edge show
save(gifti(single(label_left_parcels_no_edge)),'location_min_area_no_edge_left_parcels.shape.gii','Base64Binary');
% Save results
% save('boundary_map_lh.mat','minimametric_left_parcels','label_left_parcels','label_left_parcels_no_edge','-append');

% Right
% [minimametric_right_parcels,label_right_parcels,label_right_parcels_no_edge] = parcel_generation(labels_right_avg_raw_smooth(:,4),labels_right,right_cortex_length,[(1:size(vertexNbors_rh,1))',vertexNbors_rh],rh_vertices);
% % Minimum point
% save(gifti(single(minimametric_right_parcels)),'location_minimum_right_parcels_NIH.shape.gii','Base64Binary');
% % Parcel map with edge show
% save(gifti(single(label_right_parcels)),'location_min_area_right_parcels_NIH.shape.gii','Base64Binary');
% % Parcel map without edge show
% save(gifti(single(label_right_parcels_no_edge)),'location_min_area_no_edge_right_parcels_NIH.shape.gii','Base64Binary');
% % Save results
% save('boundary_map_rh.mat','minimametric_right_parcels','label_right_parcels','label_right_parcels_no_edge','-append');


%% Binary boundary map
function update = binary_label(label)
update = label;
update(update~=0) = -1;
update(update==0) = 1;
update(update==-1) = 0;
end




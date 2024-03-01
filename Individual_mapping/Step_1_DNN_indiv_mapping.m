clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: all analysis will be implemented in seperate hemisphere
% If there is any question, just contact me without hesitation email: txgxp88@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Current step will finish the group correlation training
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load surface information: medial wall and spatial distances of vertex
% We prepare all necessary information from the surface template, such as medial wall etc.
currentFolder = pwd;
load([currentFolder,'/Surface_template/surface_data.mat']);

%% Load atlas
%Left
surf_lh = gifti([currentFolder,'/Surface_template/Atlas_example/MBMv4_lh.shape.gii']);
surf_lh = surf_lh.cdata;
surf_label_lh = surf_lh;surf_label_lh(median_lh,:) = [];
% %Right
% surf_rh = gifti([currentFolder,'\Surface_template\Atlas_example\MBMv4_rh.shape.gii']);
% surf_rh = surf_rh.cdata;
% surf_label_rh = surf_rh;surf_label_rh(median_rh,:) = [];

%% Load Group whole brain functional connectivity

%Left
corr_lh_group = ft_read_cifti([currentFolder,'/Preprocessed_data/Example_Group_lh.dconn.nii']);corr_lh_group = corr_lh_group.dconn;
corr_lh_group(median_lh,:) = [];
corr_lh_group(:,median_lh) = [];

% %Right
% % Group whole brain functional connectivity
% corr_rh_group = ft_read_cifti([currentFolder,'/Preprocessed_data/Example_Group_rh.dconn.nii']);corr_rh_group = corr_rh_group.dconn;
% corr_rh_group(median_rh,:) = [];
% corr_rh(:,median_rh) = [];

%If prefer, we can make dimension reduction, as below PCA METHOD
%Left
[~,score_lh_group,latent_lh_group] = pca(corr_lh_group);
contribute_lh_group = cumsum(latent_lh_group)./sum(latent_lh_group);clear latent_lh_group

% %Right
% [~,score_rh_group,latent_rh_group] = pca(corr_rh_group);
% contribute_rh_group = cumsum(latent_rh_group)./sum(latent_rh_group);clear latent_rh_group


%% Deep learning network, Load Group correlation and training
score_lh_pool = score_lh_group(:,1:100);%Covering 99% information
% score_rh_pool = score_rh_group(:,1:100);%Covering 99% information

%Left
Num_parcels_lh = unique(surf_label_lh);
[net_lh,Search_voxel_lh,Search_neighbor_lh] = deep_learning_network(Num_parcels_lh,score_lh_pool,surf_lh,surf_label_lh,vertexNbors_lh,median_lh);

% %Right 
% Num_parcels_rh = unique(surf_label_rh);
% [net_rh,Search_voxel_rh,Search_neighbor_rh] = deep_learning_network(Num_parcels_rh,score_rh_pool,surf_rh,surf_label_rh,vertexNbors_rh,median_rh);


save([currentFolder,'\Results_view\Group_training.mat'],'-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For the convinence we have a example result in the folder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% Functions
function B = neighbor_area(label,neighbor,area,medial_wall)
tmp_neigh = neighbor;tmp_neigh(medial_wall,:) = [];
tmp_neigh = tmp_neigh(area,:);
B = [];
for k = 1:size(tmp_neigh,1)
    A = nonzeros(tmp_neigh(k,:));
    B = [B;label(A,:)];
    clear A
end
clear k
B = unique(B);
end

function label_list = voxel_search(Search_neighbor,label)
label_list = [];
for k = 1:length(Search_neighbor)
    A = find(label==Search_neighbor(k));A(:,2) = Search_neighbor(k);
    label_list = [label_list;A];
end
clear k
end

%%% Deep network training
function [net,Search_voxel,Search_neighbor] = deep_learning_network(Num_parcels,score_lh_pool,surf_lh,surf_label_lh,vertexNbors_lh,median_lh)
for k = 1:length(Num_parcels)
    %%%%% Prepare for training
    training_goal = surf_label_lh;%Target labeling
    %%%%% Prepare for testing
    parcel_index = find(surf_label_lh == Num_parcels(k));
    Search_neighbor{k} = neighbor_area(surf_lh,vertexNbors_lh,parcel_index,median_lh);
    Search_neighbor{k} = setdiff(Search_neighbor{k},Num_parcels(k));
    Search_voxel{k} = voxel_search(Search_neighbor{k},surf_label_lh);
    parcel_index(:,2) = Num_parcels(k);parcel_index(:,3) = 1;
    Search_voxel{k}(:,3) = 0;Search_voxel{k} = [Search_voxel{k};parcel_index];
    %%%%%%%%%%%%%
    
    %Network
    input_feature = size(score_lh_pool,2);
    for i = 1:size(score_lh_pool,1)
        data_input_4D(:,:,:,i) = score_lh_pool(i,:)';
    end
    clear i
    num_output = length(unique(training_goal));
    
    layers = [
        imageInputLayer([size(score_lh_pool,2) 1 1]) % 22X1X1 refers to number of features per sample
        convolution2dLayer(3,16,'Padding','same')
        reluLayer
        fullyConnectedLayer(384) % 384 refers to number of neurons in next FC hidden layer
        fullyConnectedLayer(384) % 384 refers to number of neurons in next FC hidden layer
        fullyConnectedLayer(num_output) % refers to number of neurons in next output layer (0 or 1)
        softmaxLayer
        classificationLayer];
    
    options = trainingOptions('sgdm', ...
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropFactor',0.2, ...
        'LearnRateDropPeriod',5, ...
        'MaxEpochs',20, ...
        'MiniBatchSize',64);
    
    net{k} = trainNetwork(data_input_4D,categorical(training_goal)',layers,options);
    clear training_goal data_input_4D layers options
end
end


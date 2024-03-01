clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: all analysis will be implemented in seperate hemisphere
% If there is any question, just contact me without hesitation email: txgxp88@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Current step will finish the individual mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load surface information: medial wall and spatial distances of vertex
% We prepare all necessary information from the surface template, such as medial wall etc.
currentFolder = pwd;
load([currentFolder,'/Surface_template/surface_data.mat']);

%% Individual mapping
%Left
% individual whole brain functional connectivity
corr_lh_individual = read_cifti([currentFolder,'/Preprocessed_data/Example_Individual_lh.dconn.nii']);
corr_lh_individual = corr_lh_individual.cdata;%Whole brain connections, different from group level
corr_lh_individual = corr(corr_lh_individual');

% %Right
% % individual whole brain functional connectivity
% corr_rh_individual = ft_read_cifti([currentFolder,'/Preprocessed_data/Example_Individual_rh.dconn.nii']);corr_rh_individual = corr_rh_individual.dconn;
% corr_rh_individual(median_rh,:) = [];
% corr_rh_individual(:,median_rh) = [];

%If prefer, we can make dimension reduction
%Left
[~,score_lh_individual,latent_lh_individual] = pca(corr_lh_individual);
contribute_lh_individual = cumsum(latent_lh_individual)./sum(latent_lh_individual);clear latent_lh_individual
score_lh_indiv = score_lh_individual(:,1:100);%Covering 99% information

% %Right
% [~,score_rh_individual,latent_rh_individual] = pca(corr_rh_individual);
% contribute_rh_individual = cumsum(latent_rh_individual)./sum(latent_rh_individual);clear latent_rh_individual
% score_rh_indiv = score_rh_individual(:,1:100);%Covering 99% information

cortices_lh = 1:size(lh_vertices,1);cortices_lh = cortices_lh';cortices_lh(median_lh,:) = [];
% cortices_rh = 1:size(rh_vertices,1);cortices_rh = cortices_rh';cortices_rh(median_rh,:) = [];

%% Mapping by the trained Deep learning
%LOAD GROUP TRAINING RESULTS from previous training
load([currentFolder,'/Preprocessed_data/Example_Training_result.mat']);

%Left hemisphere
[brain_label_predict_lh,brain_label_predict_lh_max,brain_label_lh] = deep_learning_network_classify(Num_parcels_lh,net_lh,score_lh_indiv,Search_voxel_lh,cortices_lh);

% %Right hemisphere
% [brain_label_predict_rh,brain_label_predict_rh_max,brain_label_rh] = deep_learning_network_classify(Num_parcels_rh,net_rh,score_rh_indiv,Search_voxel_rh,cortices_rh);

%% Postwork after classification
%Load necessary information from such individual (boundary information and so on)
%The results come from the process of boundary map
load([currentFolder,'/Preprocessed_data/Example_Individual_boundary_map_lh.mat'],'labels_left_avg_raw_smooth');

%Left hemisphere
labels_pool_left_avg_raw_smooth = mean(labels_left_avg_raw_smooth,2);
% Return mapping results and its optimized version
[label_individual_lh,label_individual_lh_opt] = post_deep_learning_network(labels_pool_left_avg_raw_smooth,brain_label_predict_lh_max,brain_label_predict_lh,brain_label_lh,left_cortex_length,surf_lh,vertexNbors_lh,median_lh);
save(gifti(single(label_individual_lh)),'Indiv_Label_lh.shape.gii','Base64Binary');
save(gifti(single(label_individual_lh_opt)),'Indiv_Label_lh_opt.shape.gii','Base64Binary');

% %Right hemisphere
% labels_pool_right_avg_raw_smooth = mean(labels_pool_right_avg_raw_smooth,2);
% [label_individual_rh,label_individual_rh_opt] = post_deep_learning_network(labels_pool_right_avg_raw_smooth,brain_label_predict_rh_max,brain_label_predict_rh,brain_label_rh,right_cortex_length,surf_rh,vertexNbors_rh,median_rh);
% save(gifti(single(label_individual_rh)),'Indiv_Label_rh.shape.gii','Base64Binary');
% save(gifti(single(label_individual_rh_opt)),'Indiv_Label_rh_opt.shape.gii','Base64Binary');


%% Functions
%%% Classification by deep network
function [brain_label_predict,brain_label_predict_max,brain_label] = deep_learning_network_classify(Num_parcels,net,score_lh_indiv,Search_voxel,cortices_lh)

% for every parcel
for k = 1:length(Num_parcels) 
    data_test = score_lh_indiv(Search_voxel{k}(:,1),:);
    for i = 1:size(data_test,1)
        data_test_4D(:,:,:,i) = data_test(i,:)';
    end
    clear i
    predictedLabels = classify(net{k},data_test_4D);
    predictedLabels = [double(predictedLabels),Search_voxel{k}];
    clear data_test_4D

    brain_label_predict(:,k) = cortices_lh;brain_label_predict(:,k) = 0;
    brain_label_predict(predictedLabels(:,2),k) = predictedLabels(:,1);

    brain_label(:,k) = cortices_lh;brain_label(:,k) = 0;
    brain_label(predictedLabels(:,2),k) = predictedLabels(:,4);

    clear predictedLabels

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A = find(brain_label(:,k) == 1);B = brain_label_predict(A,k);
    C = tabulate(B);C(C(:,2)==0,:)=[];C = sortrows(C,2);
    brain_label_predict_max{k,:} = C;
    clear A B C
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %     %Write and check
    %     brain = lh_vertices(:,1);brain(:) = 0;
    %     brain(cortices_lh,:) = brain_label(:,k);
    %     save(gifti(single(brain)),'Deep_origin.shape.gii','Base64Binary');
    %
    %     brain_p = lh_vertices(:,1);brain_p(:) = 0;
    %     brain_p(cortices_lh,:) = brain_label_predict(:,k);
    %     save(gifti(single(brain_p)),'Deep_mapping.shape.gii','Base64Binary');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
clear k


end

%%% Post work after deep network mapping
function [label_individual_lh,label_individual_opt_lh] = post_deep_learning_network(labels_pool_left_avg_raw,brain_label_predict_max,brain_label_predict,brain_label,left_cortex_length,surf_lh,vertexNbors_lh,median_lh)
cortices_lh = 1:size(vertexNbors_lh,1);
cortices_lh = cortices_lh';
cortices_lh(median_lh,:) = [];

for k = 1:length(brain_label_predict_max)
    tmp = brain_label_predict_max{k};
    area_targ = tmp(end,1);
    tmp = find( brain_label_predict(:,k) == area_targ );
    tmp(:,2) = brain_label(tmp,k);
    tmp(:,3) = k;
    tmp(:,4) = cortices_lh(tmp(:,1),:);
    
    tmp( tmp(:,2)==0,: ) = [];
    
    Individual_area_match(:,k) = surf_lh;Individual_area_match(:,k) = 0;
    Individual_area_match(tmp(:,4),k) = k;
    Group_template(:,k) = surf_lh;Group_template(:,k) = 0;
    A = find(brain_label(:,k)~=0);Group_template(cortices_lh(A,:),k) = k;
    clear A tmp
end
clear k

Individual_area_match(median_lh,:) = NaN;

%dilate
Individual_area_match_smooth = Individual_area_match;
for i = 1:size(Individual_area_match_smooth,2)
    if sum(Individual_area_match(:,i)) ~= 0
        Loop = 0;
        while 1 ~= 0
            Loop = Loop+1;
            for j = 1:size(Individual_area_match_smooth,1)
                if isnan( Individual_area_match_smooth(j,i) ) ~= 1
                    neigh_tmp = nonzeros(vertexNbors_lh(j,:));
                    neigh_tmp_label = Individual_area_match_smooth(neigh_tmp,i);
                    if length(find(isnan(neigh_tmp_label)==1))~=length(neigh_tmp_label)
                        A = tabulate(neigh_tmp_label);A = sortrows(A,3);
                        if Individual_area_match_smooth(j,i) ~= A(end,1)
                            Individual_area_match_smooth(j,i) = A(end,1);
                        end
                        clear neigh_tmp_label neigh_tmp A
                    end
                end
            end
            clear j
            tmp = find(Individual_area_match_smooth(:,i)~=0 & isnan(Individual_area_match_smooth(:,i))~=1);
            totals(Loop) = length(tmp);
            clear tmp
            if Loop>2
                if totals(Loop)==totals(Loop-1)
                    break
                end
            end
        end
    end
end
clear i
Individual_area_match_smooth(median_lh,:) = 0;
Individual_area_match_smooth_temp = Individual_area_match_smooth;
Individual_area_match_smooth_temp(Individual_area_match_smooth_temp~=0) = 1;
A = sum(Individual_area_match_smooth_temp,2);B = find(A>1);
Individual_area_match_smooth(B,:) = 0;clear A B Individual_area_match_smooth_temp
% save(gifti(single(Individual_area_match_smooth)),'Individual_deep_learning.shape.gii','Base64Binary');

label_template = Group_template(:,1);label_template(:) = 0;
for k = 1:size(Group_template,2)
    A = find(Group_template(:,k)~=0);
    label_template(A,:) = k;
    clear A
end
clear k

Individual_area_match_smooth_update = opt_training(Individual_area_match_smooth,label_template);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Template opt
% label_template_opt = label_template;
% for k = 1:4
%     label_template_opt(median_lh,:) = NaN;
%     for k = 1:size(label_template_opt,1)
%         if isnan(label_template_opt(k,:)) ~= 1
%             A = vertexNbors_lh(k,:);
%             B = nonzeros(A);
%             C = label_template_opt(B,:);
%             D = tabulate(C);
%             D(D(:,2)==0,:)=[];
%             D = sortrows(D,3);                
%             if label_template_opt(k,:) ~= D(end,1) && D(end,1)~= 0
%                 label_template_opt(k,:) = D(end,1);
%             end
%             clear A B C D
%         end
%     end
%     clear k
% end
% save(gifti(single(label_template_opt)),'Group_template_opt.shape.gii','Base64Binary');
% save(gifti(single(label_template)),'Group_template.shape.gii','Base64Binary');
% save(gifti(single(Group_template)),'Group_template_matrix.shape.gii','Base64Binary');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[label_individual_lh,label_individual_opt_lh] = parcel_generation_individual(labels_pool_left_avg_raw,Individual_area_match_smooth,left_cortex_length,[(1:size(vertexNbors_lh,1))',vertexNbors_lh],label_template);



% save(gifti(single(label_individual_lh)),'Individual_Label_lh.shape.gii','Base64Binary');
% [label_individual_opt_lh] = parcel_generation_individual_opt(label_individual_lh,label_template,vertexNbors_lh,median_lh);
% save(gifti(single(label_individual_opt_lh)),'Individual_Label_optimal_lh.shape.gii','Base64Binary');

end

%%% Post optimization
function Individual_area_match_smooth_update = opt_training(Individual_area_match_smooth,label_template)
Individual_area_match_smooth_list = sum(Individual_area_match_smooth,2);
Individual_area_match_smooth_update = Individual_area_match_smooth_list;
temp = find(Individual_area_match_smooth_list~=0);
A = [Individual_area_match_smooth_list(temp,:),label_template(temp,:)];
Individual_area_match_smooth_update(temp,:) = label_template(temp,:);

% save(gifti(single(Individual_area_match_smooth_update)),'Individual_Label_lh_b.shape.gii','Base64Binary');
end


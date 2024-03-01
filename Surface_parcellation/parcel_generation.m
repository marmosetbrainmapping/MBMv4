function [minimametric,label,label_no_edge] = parcel_generation(edgemetric,edgemetric_matrix,medialmaskdata,neighbors,temp_location)

%%
row = [];
for k = 1:size(edgemetric_matrix,1)
    A = find(edgemetric_matrix(k,:)==0);
    row = [row,length(A)];
    clear A
end
clear k
erodedmedialmaskdata1 = find(row>median(row));clear row

%-----------------------------------------------------------------
%Get cortical indices
corticalindices = medialmaskdata{1,2}+1;

%Get medial wall indices
medialindices = setdiff(medialmaskdata{1,3},medialmaskdata{1,2}+1);

% %Load an eroded medial mask
% erodedmedialmaskdata = gifti(eroded_medial_maskfile);
% erodedmedialmaskdata = erodedmedialmaskdata.cdata;

%% Preparation for local minima: 1st filter
disp('Finding minima')

%make a copy of the edgemap to use
metric = edgemetric;

%nodes in eroded medial mask can't be minima
metric(medialindices) = 100000;
metric(erodedmedialmaskdata1) = 100000;

%filter independent nodes
Scatter = 0.3;
erodedmedialmaskdata2 = [];
for i = 1:size(metric,1)
    if metric(i) ~= 100000
        %get this point's neighbors
        nodeneigh = neighbors(i,:);
        nodeneigh(nodeneigh == 0) = [];
        nodeneigh_value = metric(nodeneigh);
        A = find(nodeneigh_value==100000);
        if length(A) > round( length(nodeneigh_value)*Scatter )
            erodedmedialmaskdata2 = [erodedmedialmaskdata2;i];
        end
        clear A nodeneigh nodeneigh_value
    end
end
clear i
metric(erodedmedialmaskdata2) = 100000;
% save(gifti(single(metric)),'location_min_test.shape.gii','Base64Binary');

%% Iteration for the minimum points and dialation areas: : 2nd filter
% Get minima
%Define various thresholds
area_thr = 1.3;
metric_thr = 0.1;

B = [find(metric~=100000),metric( find(metric~=100000),: )];
N = 0;
while 1 ~= 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = N+1;
    B_temp(N) = size(B,1);
    
    if  N>2
        if B_temp(N) == B_temp(N-1)
            break;
        end
    end
    
    minimametric = metric;
    minimametric(:)=0;
    minimametric(B(:,1),:)=1;
    %     save(gifti(single(minimametric)),['location_min_',num2str(N),'.shape.gii'],'Base64Binary');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    neighbors_TEMP = squareform(pdist(temp_location(B(:,1),:)));
    metric_TEMP = squareform(pdist(B(:,2)));
    minimametric_label = B;
    clear B
    
    %Update I: based on distance and metric difference
    remove =[];num = 0;
    for k = 1:size(neighbors_TEMP,1)
        if isempty(remove) == 1
            A = find(neighbors_TEMP(k,:)<area_thr & neighbors_TEMP(k,:)~=0);%
            B = metric_TEMP(k,A);
            C = A(find(B<metric_thr));
            remove =[remove,C];
            num = num+1;
            cluster{num} = [[k,minimametric_label(k,:)];[C',minimametric_label(C,:)]];
            clear A B C
        else
            mark = find(remove==k);
            if isempty(mark) ~= 1
                continue;
            else
                A = find(neighbors_TEMP(k,:)<area_thr & neighbors_TEMP(k,:)~=0);%
                B = metric_TEMP(k,A);
                C = A(find(B<1));
                remove =[remove,C];
                num = num+1;
                cluster{num} = [[k,minimametric_label(k,:)];[C',minimametric_label(C,:)]];
                clear A B C
            end
        end
    end
    clear k num remove
    
    A = [];
    for k = 1:length(cluster)
        B = find(cluster{k}(:,3)==min(cluster{k}(:,3)));
        A = [A;cluster{k}(B(1),2:3)];
        clear B
    end
    clear k
    B = unique(A,'rows');clear A
    
end

% erodedmedialmaskdata3 = [];Scatter=1;
% for i = 1:size(B,1)
%     %get this point's neighbors
%     nodeneigh = neighbors(B(i),:);
%     nodeneigh(nodeneigh == 0) = [];
%     nodeneigh_value = metric(nodeneigh);
%     A = find( nodeneigh_value>metric(B(i)) );
%     if length(A) < round( length(nodeneigh_value)*Scatter )
%         erodedmedialmaskdata3 = [erodedmedialmaskdata3;i];
%     end
%     clear A nodeneigh nodeneigh_value
% end
% clear i
% 
% B(erodedmedialmaskdata3,:) = [];
minimametric = metric;
minimametric(:)=0;
minimametric(B(:,1),:)=1;
% save(gifti(single(minimametric)),'location_min_final.shape.gii','Base64Binary');
clear B

%Grow parcels from minima using watershed technique
while 1~=0
    %Remove NaNs from edgemap (probably not any)
    edgemetric(isnan(edgemetric)) = 0;
    
    % Wateshild dialation:getting the area
    label = watershield_Dilation(edgemetric,minimametric,corticalindices,neighbors);

    %Remove values in medial wall
    label(medialindices,:) = 0;
    
    % Examine the size
    temp1 = find(minimametric == 1);
    temp2 = label(temp1,:);
    for k = 1:size(temp1,1)
        label_area(k,:) = length(find(label==temp2(k)));
    end
    clear k
    label_area = [temp1,label_area];
    A = find(label_area(:,2)<50);
    if isempty(A) ~= 1
        minimametric(label_area(A,1),:) = 0;
        clear label_area
    else
        break
    end
    clear temp1 temp2 A label
end

%Preparation:remove the edges
A = squareform( pdist( temp_location(corticalindices,:) ) );
temp_lh_update = label(corticalindices,:);
iter = 0;distance = 1;
while 1~=6
    iter = iter+1;
    B = find(temp_lh_update==0);count(iter) = length(B);
    for k = 1:length(B)
        C = A(B(k),:);C(B(k))=10000;
        C=C';C(:,2)=1:length(C);
        C = sortrows(C,1);
        index_edge(k,:) = [B(k),C(distance,2),temp_lh_update(C(distance,2))];
        clear C
    end
    clear k
    temp_lh_update(B,:) = index_edge(:,end);
    if length(count) > 1
        if isempty( find(temp_lh_update==0) ) == 1
            break;
        else
            if count(iter) == count(iter-1)
                if count(iter-1) ~= count(iter-2)
                    distance = distance+1;
                end
            end
        end
    end
    clear B index_edge
end
clear A B count iter distance
label_no_edge = label;
label_no_edge(corticalindices,:) = temp_lh_update;
clear A

% %%
% save(gifti(single(minimametric)),'location_minimum.shape.gii','Base64Binary');
% save(gifti(single(label)),'location_min_area.shape.gii','Base64Binary');
% save(gifti(single(label_no_edge)),'location_min_area_no_edge.shape.gii','Base64Binary');


function label = watershield_Dilation(edgemetric,minimametric,corticalindices,neighbors)

%Remove NaNs from edgemap (probably not any)
edgemetric(isnan(edgemetric)) = 0;

%Initialize the parcels and put unique values at all local minima
label = zeros(size(minimametric));
labelnum = nnz(minimametric);
labelpos = find(minimametric==1);
[ign sortorder] = sort(edgemetric(labelpos));
for j = 1:labelnum
    label(labelpos(j)) = sortorder(j);
end
clear j
%Initialize a variable keeping track of final border positions
watershedzone = zeros(size(label));

hiter = unique(edgemetric(corticalindices));

%Iterate through the edgemap values
for i = 1:length(hiter)
    
    string{i} = ['Growing parcels through ' num2str(i) ' out of ' num2str(length(hiter)) ' values'];
    if i==1; fprintf('%s',string{i}); else fprintf([repmat('\b',1,length(string{i-1})) '%s'],string{i}); end
    
    % Take values in edgemap less than current iteration
    maskmetrics = edgemetric<hiter(i);
    maskmetrics = maskmetrics & ~label>0 & ~watershedzone;
    
    maskpos = find(sum(maskmetrics,2)>0);
    %maskpos = maskpos(randperm(length(maskpos)));
    
    for m = 1:length(maskpos) %For all nodes at this threshold
        
        %get node neighbors and labels of those neighbors
        nodeneigh = neighbors(maskpos(m),:);
        nodeneigh(nodeneigh==0) = [];
        nodeneighlab = label(nodeneigh);
        
        %Find minimum value other than 0 among neighbors
        minfindnodeneighlab = nodeneighlab;
        minfindnodeneighlab(nodeneighlab==0) = 100000;
        minnodeneighlab = min(minfindnodeneighlab,[],1);
        
        %Find maximum value other than 0 among neighbors
        maxfindnodeneighlab = nodeneighlab;
        maxfindnodeneighlab(nodeneighlab==0) = -100000;
        maxnodeneighlab = max(maxfindnodeneighlab,[],1);
        
        %If min and max differ (i.e. two or more neighbor parcels), it's a
        %border
        maskinthismetric = maskmetrics(maskpos(m),:);
        watershed_nodes = (minnodeneighlab~=maxnodeneighlab) & (minnodeneighlab~=100000) & maskinthismetric;
        label(maskpos(m),watershed_nodes) = 0;
        watershedzone(maskpos(m),watershed_nodes) = 1;
        
        %If min and max are the same but different from 0, make the node
        %that value
        next_to_water = (minnodeneighlab==maxnodeneighlab) & (minnodeneighlab~=100000) & maskinthismetric;
        label(maskpos(m),next_to_water) = minnodeneighlab(next_to_water);
        
    end
end
clear i hiter


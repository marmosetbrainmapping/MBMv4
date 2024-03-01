function [vertexNbors_lh_value1_update] = spatial_smooth(labels_left_raw_average,vertexNbors_lh)


vertexNbors_lh_value = index_neighbor(labels_left_raw_average,vertexNbors_lh);
vertexNbors_lh_value = [labels_left_raw_average,vertexNbors_lh_value];
vertexNbors_lh_value_label_value = index(vertexNbors_lh_value);

% save(gifti(single(vertexNbors_lh_value_label_value(:,1))),'edge_map_no_Left_fisher_mooth03_mean.shape.gii','Base64Binary');
% save(gifti(single(vertexNbors_lh_value_label_value(:,2))),'edge_map_no_Left_fisher_mooth03_max.shape.gii','Base64Binary');
% save(gifti(single(vertexNbors_lh_value_label_value(:,3))),'edge_map_no_Left_fisher_mooth03_min.shape.gii','Base64Binary');

temp1 = vertexNbors_lh_value_label_value(:,1);
stop1 = 0;num = 0;
while 1 ~= 0
    
    if stop1 == 0
        num = num+1;
        vertexNbors_lh_value1 = index_neighbor(temp1(:,end),vertexNbors_lh);
        vertexNbors_lh_value1 = [temp1(:,end),vertexNbors_lh_value1];
        vertexNbors_lh_value_label_value1 = index(vertexNbors_lh_value1);vertexNbors_lh_value_label_value1 = vertexNbors_lh_value_label_value1(:,1);
        temp1 = [temp1,vertexNbors_lh_value_label_value1];
        clear vertexNbors_lh_value_label_value1 vertexNbors_lh_value1
    end
    
    if size(temp1,2)>1
        A = temp1(:,end)-temp1(:,end-1);
        if isempty(find(A~=0)) == 1 || num > 10
            stop1 = 1;
        end
        clear A
    end
    
    if  stop1 == 1
        break
    end

end
    
vertexNbors_lh_value1_update = temp1;%average

end


function vertexNbors_lh_value_label_value = index(vertexNbors_lh_value)

for i = 1:size(vertexNbors_lh_value,1)
    temp = vertexNbors_lh_value(i,:);temp(temp==0)=[];
    [~,I] = min(temp); 
    if I == 1
        vertexNbors_lh_value_label(i,:) = 1;
    else
        vertexNbors_lh_value_label(i,:) = 0;
    end
    clear I temp
end
clear i

for i = 1:size(vertexNbors_lh_value,1) 
    if vertexNbors_lh_value(i,1) ~= 0
        if vertexNbors_lh_value_label(i) == 1
            temp = vertexNbors_lh_value(i,:);temp(temp==0)=[];
            vertexNbors_lh_value_label_value(i,1) = mean(temp(2:end));
%             vertexNbors_lh_value_label_value(i,2) = max(temp(2:end));
%             vertexNbors_lh_value_label_value(i,3) = min(temp(2:end));
        else
            vertexNbors_lh_value_label_value(i,1) = vertexNbors_lh_value(i,1);
        end
        clear temp
    else
        vertexNbors_lh_value_label_value(i,1) = vertexNbors_lh_value(i,1);
    end
end
clear i
end

function vertexNbors_lh_value = index_neighbor(labels_left_raw_average,vertexNbors_lh)
vertexNbors_lh_value = double( vertexNbors_lh );
for i = 1:size(vertexNbors_lh_value,1)
    for j = 1:size(vertexNbors_lh_value,2)
        if vertexNbors_lh(i,j)~=0
            vertexNbors_lh_value(i,j) = labels_left_raw_average(vertexNbors_lh(i,j));
        end
    end
    clear j
end
clear i
end
function minimametric = metric_minima_all(metric,neighdist,neighbors)

         

minimametric = zeros(size(metric));
for i = 1:size(metric,1)
    
    nodeneigh = i;
    newneigh = nodeneigh;
    curneigh = newneigh;
    
    for n = 1:neighdist
        
        for t = 1:length(curneigh)
            try
%             newneigh = [newneigh neighbors(curneigh(t),2:end)];%for human
%             data
            newneigh = [newneigh neighbors(curneigh(t),1:end)];
            catch
                1;
            end
%             newneigh(isnan(newneigh)) = [];%for human data
            newneigh(newneigh==0) = [];
        end
        curneigh = setdiff(newneigh,nodeneigh);
        nodeneigh = union(nodeneigh,newneigh,'stable');
        
    end
    
    nodeneighval_all = metric(nodeneigh(2:end),:);
    origval_all = repmat(metric(nodeneigh(1),:),[length(nodeneigh)-1 1]);
    
    minval = sum(origval_all<nodeneighval_all,1);
    minimametric(i,:) = single(minval==(length(nodeneigh)-1));
    
end

    
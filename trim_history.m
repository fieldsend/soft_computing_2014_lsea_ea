function    [archive_modes] = trim_history(archive_modes,max_l)

% function  [archive_modes] = trim_history(archive_modes,max_l)
%
% for each mode in acrhive_modes, prune history to limit (max_l), 
% based on fitness

for i=1:length(archive_modes)
    
    %clean duplicate locations on merging
    l=size(archive_modes(i).local_region.history_locations,1);
    
    loc_dist=dist2(archive_modes(i).local_region.history_locations,...
        archive_modes(i).local_region.history_locations);
    loc_dist(1:l+1:l*l) = inf; % make distance to itself infinite
    [mn,ind]=min(min(loc_dist));
    
    while (mn==0)
        archive_modes(i).local_region.history_locations(ind,:)=[];
        archive_modes(i).local_region.history_values(ind,:)=[];
        l=l-1;
        loc_dist=dist2(archive_modes(i).local_region.history_locations,...
            archive_modes(i).local_region.history_locations);
        loc_dist(1:l+1:l*l) = inf; % make distance to itself infinite
        [mn,ind]=min(min(loc_dist));
    end
    
    if length(archive_modes(i).local_region.history_values)>max_l
        l=length(archive_modes(i).local_region.history_values);
        
        % remove least fit
        if (l>max_l)
            % order from best to worst
            [~,si] = sort(archive_modes(i).local_region.history_values,'descend');
            archive_modes(i).local_region.history_locations(si(max_l+1:end),:)=[];
            archive_modes(i).local_region.history_values(si(max_l+1:end))=[];
        end
    end
end
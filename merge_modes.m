function [active_modes,M_loc,V_loc,number_of_mid_evals,X,Y,index] = ...
    merge_modes(active_modes,active_modes_changed,problem_func,M_loc,V_loc,params,X,Y,index,tol_value,remaining_evals)

% function [active_modes,M_loc,V_loc,number_of_mid_evals,X,Y,index] = ...
%    merge_modes(active_modes,active_modes_changed,problem_func,M_loc,V_loc,params,X,Y,index,tol_value,remaining_evals)
%
% function merges the mode estimates if they are deemed to be
% optimising the same niche (via mid-point checking), or they are
% within some tiny tol_value. Only concerned with pairing of modes
% that are marked to have changed since the last merge check.

% only concern ourselves with modes that have actually shifted, or are new
% since the last generation, as no need to check others
I = find(active_modes_changed==1);

n = length(I);
to_compare = zeros(n,2);
to_compare(:,1) = I;
for i=1:n
    d = dist2(M_loc(I(i),:),M_loc);
    d(I(i)) = inf; % will be closest to itself, so need to get second closest
    [tmp, to_compare(i,2)] = min(d);
    active_modes(I(i)).local_region.dist = sqrt(tmp); % track Euc dist to nearest neighbour mode
end

% to_compare now contains the pairs of indices of closest modes, where at
% least one mode has shifted location/is new since last generation. However,
% there may be duplicated pairs (through reversals), so need to omit these.
to_compare = sort(to_compare,2);
% now sorted so that first column elements are always smaller than second
% column elemets on same row
[~,ind] = sort(to_compare(:,1));
to_compare = to_compare(ind,:); %now sorted so that first column is sorted smallest to highest

% remove_duplicates
for i=n:-1:2
    I = find(to_compare(:,1) == to_compare(i,1)); %get indices of all with first index element same
    if sum(sum(repmat(to_compare(i,:),length(I),1) == to_compare(I,:),2) == 2) > 1 % if more than one vector duplication
        to_compare(i,:) = [];
    end
end

% now check for merging
n = size(to_compare,1);
to_merge = [];
number_of_mid_evals=0;
for i = 1:n
    % merge if sufficiently close
    if (sqrt(dist2(active_modes(to_compare(i,1)).local_region.mode_location,active_modes(to_compare(i,2)).local_region.mode_location)) ...
            < tol_value)
        to_merge = vertcat(to_merge, i);
    else % otherwise merge if mid region is fitter
        % evaluate exact mid point between modes, and add to mode 2
        % history
        if (number_of_mid_evals < remaining_evals) % at very end of optimiser run, may exhaust evals before all pairs are compared
            mid_loc = 0.5*(active_modes(to_compare(i,1)).local_region.mode_location...
                -active_modes(to_compare(i,2)).local_region.mode_location)...
                +active_modes(to_compare(i,2)).local_region.mode_location;
            
            active_modes(to_compare(i,2)).local_region.new_location = mid_loc;
            [active_modes(to_compare(i,2)).local_region,mode_shift,y,X,Y,index] = ...
                evaluate(active_modes(to_compare(i,2)).local_region,problem_func,params,X,Y,index);
            if mode_shift==1 % better than mode 2 current mode, so merge
                M_loc(to_compare(i,2),:)= active_modes(to_compare(i,2)).local_region.mode_location;
                V_loc(to_compare(i,2))= active_modes(to_compare(i,2)).local_region.mode_value;
                to_merge = vertcat(to_merge, i);
            elseif (active_modes(to_compare(i,1)).local_region.mode_value < y) % better than mode 1 current mode, so merge
                to_merge = vertcat(to_merge, i);
            end
            number_of_mid_evals = number_of_mid_evals+1;
        end
    end
end
% merge those marked pairs, and flag the lower one for deletion
delete_index= zeros(size(to_merge));
for i=1:length(to_merge)
    
    kept_ind=1;
    % if peak of mode 1 is higher than mode 2, then replace
    if (active_modes(to_compare(to_merge(i),1)).local_region.mode_value > ...
            active_modes(to_compare(to_merge(i),2)).local_region.mode_value)
        active_modes(to_compare(to_merge(i),1)).local_region.history_locations = ...
            vertcat(active_modes(to_compare(to_merge(i),1)).local_region.history_locations,...
            active_modes(to_compare(to_merge(i),2)).local_region.history_locations);
        
        active_modes(to_compare(to_merge(i),1)).local_region.history_values = ...
            vertcat(active_modes(to_compare(to_merge(i),1)).local_region.history_values,...
            active_modes(to_compare(to_merge(i),2)).local_region.history_values);
        %'merge'
        delete_index(i) = to_compare(to_merge(i),2);
        
        if (active_modes(to_compare(to_merge(i),1)).local_region.width_multiplier > ...
                active_modes(to_compare(to_merge(i),2)).local_region.width_multiplier)
            active_modes(to_compare(to_merge(i),1)).local_region.width_multiplier = ...
                active_modes(to_compare(to_merge(i),2)).local_region.width_multiplier;
            active_modes(to_compare(to_merge(i),1)).local_region.less_fit_move = ...
                active_modes(to_compare(to_merge(i),2)).local_region.less_fit_move;  
        end
        
    else
        active_modes(to_compare(to_merge(i),2)).local_region.history_locations = ...
            vertcat(active_modes(to_compare(to_merge(i),2)).local_region.history_locations,...
            active_modes(to_compare(to_merge(i),1)).local_region.history_locations);
        
        active_modes(to_compare(to_merge(i),2)).local_region.history_values = ...
            vertcat(active_modes(to_compare(to_merge(i),2)).local_region.history_values,...
            active_modes(to_compare(to_merge(i),1)).local_region.history_values);
        
        %'merge'
        delete_index(i) = to_compare(to_merge(i),1);
        kept_ind=2;
        
        if (active_modes(to_compare(to_merge(i),2)).local_region.width_multiplier > ...
                active_modes(to_compare(to_merge(i),1)).local_region.width_multiplier)
            active_modes(to_compare(to_merge(i),2)).local_region.width_multiplier = ...
                active_modes(to_compare(to_merge(i),1)).local_region.width_multiplier;
            active_modes(to_compare(to_merge(i),2)).local_region.less_fit_move = ...
                active_modes(to_compare(to_merge(i),1)).local_region.less_fit_move;
        end
    end
    
    %clean duplicate locations on merging
    l=size(active_modes(to_compare(to_merge(i),kept_ind)).local_region.history_locations,1);
    
    loc_dist=dist2(active_modes(to_compare(to_merge(i),kept_ind)).local_region.history_locations,...
        active_modes(to_compare(to_merge(i),kept_ind)).local_region.history_locations);
    loc_dist(1:l+1:l*l) = inf; % make distance to itself infinite
    [mn,ind]=min(min(loc_dist));
    
    while (mn==0)
        active_modes(to_compare(to_merge(i),kept_ind)).local_region.history_locations(ind,:)=[];
        active_modes(to_compare(to_merge(i),kept_ind)).local_region.history_values(ind,:)=[];
        l=l-1;
        loc_dist=dist2(active_modes(to_compare(to_merge(i),kept_ind)).local_region.history_locations,...
            active_modes(to_compare(to_merge(i),kept_ind)).local_region.history_locations);
        loc_dist(1:l+1:l*l) = inf; % make distance to itself infinite
        [mn,ind]=min(min(loc_dist));
    end
    
end

% remove one of the merged pair, go for 1 as ordered lowest to highest on
% to_compare index vector
prev_merge=-1;
delete_index=sort(delete_index);
for i=length(delete_index):-1:1
    if (delete_index(i)~= prev_merge) % if not duplicated
        %delete_index(i)
        prev_merge=delete_index(i);
        %active_modes(delete_index(i)).local_region.mode_value
        active_modes(delete_index(i))=[];
        M_loc(delete_index(i),:)=[];
        V_loc(delete_index(i))=[];
    end
end
        
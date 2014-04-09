function [local_region] = select_new_location_close_to_current_mode(local_region,params,M_loc,tol_value)

% [local_region] =
% select_new_location_close_to_current_mode(local_region,params,M_loc,tol_value)
%
% function evolves the selected local EA, and generates a new
% location which is close to its current mode estimate -- stored as
% new_location in the local_region structure

new_location = params.minimum_values-1;

dt = dist2(local_region.mode_location,M_loc);
dt = sort(dt);
if length(dt)>1
	local_region.dist = sqrt(dt(2)); % will be second closest, as closest to itself
else
	% only one mode currently maintained -- so choose d as tenth of valid search space
	local_region.dist = sqrt(dist2(params.minimum_values,params.maximum_values));
end

d = max(local_region.dist, tol_value);
reject_count=0;
r=rand();
while (sum(new_location < params.minimum_values)> 0 || ...
        sum(new_location > params.maximum_values)>0 || ...
        sqrt(dist2(new_location,local_region.mode_location))>d/2)
    % reweight down by dimensionality
    if size(local_region.history_values,1)==1
        new_location = local_region.mode_location + randn(size(new_location))*d/10;
    else
        new_location = local_region.mode_location;
        [v,k]=max(std(local_region.history_locations));
        if r<0.5 && v>0
            new_location(k) = new_location(k) + randn()*v;
        else
            k = randperm(length(new_location));
            k=k(1);
            new_location(k) = local_region.mode_location(k) + ...
                randn()*max(sqrt(dist2(local_region.mode_location,local_region.history_locations)));
        end
    end
    reject_count =reject_count+1;
    if (reject_count>=10)
        new_location = local_region.mode_location + randn(size(new_location))*d/10;
    end
end
local_region.new_location = new_location;

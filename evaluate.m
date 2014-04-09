function [local_region,mode_shift,y,X,Y,index] = ...
    evaluate(local_region,problem_func,test_function_params,X,Y,index)

% evaluates the new_location member of the local_region argument,
% and updates the mode_estimate of the locla EA niche if it is
% better. Also updates histories, matrices of sampled locations and
% if the mode has moved, the mode_shift output is set to 1
% (otherwise will be 0)

y = feval(problem_func, local_region.new_location,test_function_params);
mode_shift=0;
if (y > local_region.mode_value)
    local_region.mode_location = local_region.new_location;
    local_region.mode_value = y;
    mode_shift = 1;
end

local_region.history_locations = [local_region.history_locations; local_region.new_location];
local_region.history_values = [local_region.history_values; y];

X(index,:) = local_region.new_location;
Y(index) = y;
index = index+1;
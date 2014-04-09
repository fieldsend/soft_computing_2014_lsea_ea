function [local_region,X,Y,index] = ...
    evaluate_first(local_region,problem_func,test_function_params,X,Y,index)

% function [local_region,X,Y,index] = ...
%    evaluate_first(local_region,problem_func,test_function_params,X,Y,index)
%
% new_location stored in the local region structure is the only
% solution thus far in mode region estimate, so by definition is
% also the mode estimate, and the only history thus far
% this function evalutes the new_location, and updates associated
% structures and matrices

y = feval(problem_func, local_region.new_location,test_function_params);
local_region.mode_location = local_region.new_location;
local_region.mode_value = y;

local_region.history_locations = local_region.new_location;
local_region.history_values = y;
local_region.less_fit_move=0; % initialise tracker of moves -- not
                              % used in current version
local_region.width_multiplier=1;
X(index,:) = local_region.new_location;
Y(index) = y;
index = index +1;
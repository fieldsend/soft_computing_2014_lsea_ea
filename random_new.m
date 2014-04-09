function [active_modes,active_modes_changed,number_rand_modes,M_loc,V_loc,XX,Y,index] = ...
    random_new(active_modes,active_modes_changed,problem_func,M_loc,V_loc,rng,params,test_function_params,XX,Y,index)

% function [active_modes,active_modes_changed,number_rand_modes,M_loc,V_loc,XX,Y,index] = ...
%    random_new(active_modes,active_modes_changed,problem_func,M_loc,V_loc,rng,params,test_function_params,XX,Y,index)
%
% function updates the active_modes structure with a single new
% mode, randomly generated in the legal search domain

number_rand_modes = 1;
X = feval(rng,number_rand_modes,params);
active_modes_changed = [active_modes_changed; ones(number_rand_modes,1)];
active_modes(end+1).local_region.new_location = X(1,:);

[active_modes(end).local_region,XX,Y,index] = evaluate_first(active_modes(end).local_region,problem_func,test_function_params,XX,Y,index);
M_loc = [M_loc; X];
V_loc = [V_loc; active_modes(end).local_region.mode_value];
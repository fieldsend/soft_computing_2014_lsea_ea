function [active_modes,active_modes_changed,number_of_new_modes,M_loc,V_loc,t,X,Y,index] = ...
    evolve(active_modes,active_modes_changed,problem_func,M_loc,V_loc,t,params,gpc,test_function_params,X,Y,index,evals_remaining)

% function [active_modes,active_modes_changed,number_of_new_modes,M_loc,V_loc,t,X,Y,index] = ...
%    evolve(active_modes,active_modes_changed,problem_func,M_loc,V_loc,t,params,gpc,test_function_params,X,Y,index,evals_remaining)
%
% function evolves modes via crossover (if the input t is equal to the input gpc)


min_vals = params.minimum_values;
max_vals = params.maximum_values;

n = length(active_modes);
number_of_new_modes =0;
if (t==gpc)
    if (n>1)
        II = randperm(n);
        if rem(n,2)==0
            R = zeros(n,length(active_modes(1).local_region.new_location));
            nn=n;
        else
            nn=n-1;
            R = zeros(nn,length(active_modes(1).local_region.new_location));
        end
        for i=1:2:nn;
            [R(i,:), R(i+1,:)] = SBX(active_modes(II(i)).local_region.mode_location,active_modes(II(i+1)).local_region.mode_location,20,min_vals,max_vals);
        end
        
        if (size(R,1)<=evals_remaining) % standard situation, evaluate all children
            M_loc = [M_loc; R];
            for i=1:size(R,1)
                local_region.new_location = R(i,:);
                [local_region,X,Y,index] = evaluate_first(local_region,problem_func,test_function_params,X,Y,index);
                V_loc = vertcat(V_loc, local_region.mode_value);
                active_modes(end+1).local_region = local_region;
            end
            % mark these as new
            active_modes_changed = [active_modes_changed; ones(size(R,1),1)];
            number_of_new_modes = size(R,1);
        else % special situation when in last iteration and do not have capacity to evaluate all children
            R = R(1:evals_remaining,:);
            M_loc = [M_loc; R];
            for i=1:size(R,1)
                local_region.new_location = R(i,:);
                [local_region,X,Y,index] = evaluate_first(local_region,problem_func,test_function_params,X,Y,index);
                V_loc = vertcat(V_loc, local_region.mode_value);
                active_modes(end+1).local_region = local_region;
            end
            % mark these as new
            active_modes_changed = [active_modes_changed; ones(size(R,1),1)];
            number_of_new_modes = size(R,1);
        end
    end
    t=0;
end
t =t+1;

function [active_modes,active_modes_changed,number_of_new_locations,M_loc,V_loc,X,Y,index] = ...
    evaluate_new_locations(active_modes,problem_func,M_loc,V_loc,test_function_params,X,Y,index,I)

% function [active_modes,active_modes_changed,number_of_new_locations,M_loc,V_loc,X,Y,index] = ...
%    evaluate_new_locations(active_modes,problem_func,M_loc,V_loc,test_function_params,X,Y,index,I)
%
% evalute all new locations generated in local mode regions, and
% track any changes to mode estimates

active_modes_changed = zeros(length(active_modes),1);
for i=1:length(I)
    [active_modes(I(i)).local_region,mode_shift,~,X,Y,index] = ...
        evaluate(active_modes(I(i)).local_region,problem_func,test_function_params,X,Y,index);
    if mode_shift==1
        active_modes_changed(I(i)) = 1;
        M_loc(I(i),:) = active_modes(I(i)).local_region.new_location;
        V_loc(I(i)) = active_modes(I(i)).local_region.mode_value;
        active_modes(I(i)).local_region.less_fit_move=0;
    else
        active_modes(I(i)).local_region.less_fit_move= ...
            active_modes(I(i)).local_region.less_fit_move+1; 
        % incremenet count on bad moves -- not currently expolited
        % in algorithm
    end
end
number_of_new_locations = length(I);
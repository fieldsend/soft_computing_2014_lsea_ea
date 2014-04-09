function [active_modes,active_modes_changed] = ...
    get_initial_locations(n,rsg,params)

% function [active_modes,active_modes_changed] = ...
%    get_initial_locations(n,rsg,params)
%
% function initialises population of niche estimates, with n
% niches, distributed in design space according to the 
% solution generator function 'rsg' and the legal bounds stored
% in the params structure

X = feval(rsg, n, params);
active_modes_changed = ones(n,1);
for i=1:n
    active_modes(i).local_region.new_location = X(i,:);
end
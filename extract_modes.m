function [RES,RES_Y] = extract_modes(active_modes)

% function [RES,RES_Y] = extract_modes(active_modes)
%
% function extracts the contents of the active_modes structure and
% places it into a matrix of locations (RES) and a matrix of
% responses (RES_Y)

RES = zeros(length(active_modes), length(active_modes(1).local_region.mode_location));
RES_Y = zeros(length(active_modes),1);

for i=1:length(active_modes)
    RES(i,:) = active_modes(i).local_region.mode_location;
    RES_Y(i) = active_modes(i).local_region.mode_value;
end
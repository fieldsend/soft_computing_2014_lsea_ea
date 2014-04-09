function [RES,RES_Y,active_modes,X,Y,state] = lsea_ea(problem_func,...
    problem_function_params,decision_variables,max_evaluations,...
    initial_solutions,gens_per_crossover,elite_type,elite_p_or_num,...
    rsg_params,max_hist,random_solution_generator)

% [RES,RES_Y,active_modes,X,Y,state] = lsea_ea(problem_func,
%         problem_function_params,decision_variables,max_evaluations,
%         initial_solutions,igens_per_crossover,elite_type,
%         elite_percent,random_solution_generator,rsg_params)
%
% NOTE: Algorithm assumes MAXIMISATION
%
% MANDATORY ARGUMENTS
%
% problem_func = string containing function to be optimised, expected to 
%     take two arguments - a vector containing a putative solution for 
%     evaluatione, and a second argument of problem specific meta-
%     parameters. This will be evaluated internally using feval()
% problem_function_params = scalar/vector/structure of any problem function meta-
%     parameters that are required
% decision_variables = number of decision variables
% max_evaluations = maximum number of evaluations to be taken through the
%     problem function
%
% OPTIONAL ARGUMENTS
%
% initial_soltions = number of initial random solutions, default 100
% gens_per_crossover = number of generations per peak population crossover, 
%    default 10
% elite_type = flag for elite process to use in alternate generations. 1 
%    means no elitism is used. 2 means the best `elite_p_or_num' % peaks by 
%    range are selected. 3 means the best `elite_p_or_num' % by rank, 
%    4 means the 'elite_p_or_num' value will be taken as max number to use
%    of ranked peaks. If not entered, or and invalid value entered, then 
%    elitism by ranking (number) will be selected (i.e. elite_type = 4).
% elite_p_or_num = percentage of `best' solutions to be selected 
%    every other generation for peak evolution (this percentage will be 
%    repeatedly selected until the number of function evaluations reached 
%    achieves that which would occur if *all* peaks had been evolved) if 
%    elite type is 1, 2 or 3. Otherwise it is the maximum number of best 
%    ranked individuals to repeatedly select. If a value is not entered 
%    this is set at 10 (subject to it effectively being fixed at 100 if 
%    elite_type==1) 
% rsg_params = meta-parameters required by the random solution generator. 
%    If not entered, then a structure with two components, maximum_values 
%    all ones and minimum_values all zeros (i.e. legal solutions defined as 
%    the unit hypercube) is used
% max_hist = the maximum population of a local EA. Default 20. 
% random_solution_generator = function to return legal random solution for
%    given problem -- must take two arguments, the number of decision 
%    variables and a structure holding any required meta parameters (e.g. 
%    maximum and minimum bounds for all decision variables). If this 
%    argument isn't provided the function a latin hypercube random solution
%    generator will be used  
%
% This algorithm forms part of a condititionally accepted work for 
% publication in Soft Computing (as part of the special issue on UKCI 2013)
%
%
%
% (c) Jonathan Fieldend 2013 & 2014, University of Exeter 

% set up any optional arguments not provided

if nargin()<4
    error('insufficient arguments provided to run optimiser');
end
if (exist('initial_solutions','var')==0)
    initial_solutions = 100;
    display('Number of initial solutions set at 100');
end
if (exist('gens_per_crossover','var')==0)
    gens_per_crossover = 10;
    display('Peak crossover set at every 10 generations');
end
if (exist('elite_type','var')==0)
    elite_type = 4;
    display('Elitism type not entered, ranking used, with max number');
end
if (elite_type~=1) && (elite_type~=2) && (elite_type~=3) && (elite_type~=4)
    elite_type = 4;
    display('Elitism type not valid value, ranking used, with max number');
end
if (exist('elite_p_or_num','var')==0)
    if elite_type~=1
        elite_p_or_num=10;
        display('Elistism percentage not entered, set as 10');
    end
end
if exist('max_hist','var')==0
    max_hist=20;
    display('No maximum local EA population size set, 20 used as default');
end
if exist('random_solution_generator','var')==0
    random_solution_generator = 'lhsg';
    display('No random solution generator specified, latin hypercube sampling used');
end
if (exist('rsg_params','var')==0)
    rsg_params.maximum_values = ones(decision_variables,1);
    rsg_params.minimum_values = zeros(decision_variables,1);
    display('Meta-parameters for random number generator not set. Bounds set as unit hypercube');
end
if (elite_p_or_num < 1)
    elite_p_or_num=1;
    display('Minimum value for elite percent or number term is 1, so has been set to 1');
elseif (elite_p_or_num > 100) && ((elite_type==2) || (elite_type==3))
    elite_p_or_num=100;
    display('Maximum value for elite percentage term is 100, so has been set to 100');
end
if max_hist<1
   max_hist=1;
   display('Cannot have the local EA population maximum less than 1, so set as 1');
end

tol_value = 10^-6; % small tolerance to merge close peaks together

elite_proportion = elite_p_or_num/100; % if percentage is being use, convert to proportion
elite_max_number = elite_p_or_num;% if absolute number is being used, store here

X = zeros(max_evaluations, decision_variables);
Y = zeros(max_evaluations,1);
index =1;
state=[];

% get initial locations
[active_modes,active_modes_changed] = ...
    get_initial_locations(initial_solutions,random_solution_generator,rsg_params);

% evaluation initial locations
for i=1:initial_solutions
    [active_modes(i).local_region,X,Y,index] = ...
        evaluate_first(active_modes(i).local_region,problem_func,problem_function_params,X,Y,index);
end
% track number of evaluations taken
evaluations = initial_solutions;

% keep modes in matrix for effciency on some computations
M_loc = zeros(initial_solutions, length(active_modes(1).local_region.new_location));
V_loc = zeros(initial_solutions,1);
for i=1:initial_solutions
    M_loc(i,:) =  active_modes(i).local_region.new_location;
    V_loc(i,:) = active_modes(i).local_region.mode_value;
end

t=0;
iteration = 1;
while (evaluations<max_evaluations)
    % first see if modes should be merged together
    [active_modes,M_loc,V_loc,number_of_mid_evals,X,Y,index] = ...
        merge_modes(active_modes,active_modes_changed,problem_func,M_loc,V_loc,problem_function_params,X,Y,index,tol_value,max_evaluations-evaluations);
    evaluations = evaluations+number_of_mid_evals;
    
    internal_count=0; % counter for how many times the elite subset are evolved this generation
    I = 1:length(active_modes);
    
    number_of_new_locations=0;
    if (rem(iteration,2) == 0) || (elite_type == 1)
        if length(I) > max_evaluations-evaluations
           % special case at end of optimiser run where there are fewer
           % function evaluations available than the total number of modes
           % As such, we preferentially evolve just the best remaining
           [~, fit_I] = sort(V_loc(I,:),'descend');
           I = fit_I(1:max_evaluations-evaluations);
        end
        % every even generation evolve *all* peaks
        for i=1:length(I) 
            [active_modes(I(i)).local_region] = ...
                select_new_location_close_to_current_mode(active_modes(I(i)).local_region,rsg_params,M_loc,tol_value);
        end
        [active_modes,active_modes_changed,number_of_new_locations,M_loc,V_loc,X,Y,index] = ...
            evaluate_new_locations(active_modes,problem_func,M_loc,V_loc,problem_function_params,X,Y,index,I);
    else % every odd generation bias peak selection to the elites
        LL = length(I);
        internal_its = 0;
        temp = zeros(length(active_modes),1);
        % first check on repititions amassed by elite function evals, the second 
        % condition covers the special case of the final generation
        while ((internal_its < LL) && (internal_its < max_evaluations-evaluations))
            if elite_type == 2 % top elite_percent by range
                range = max(V_loc)-min(V_loc);
                I2 = find(V_loc(I)>=min(V_loc)+range*(1.0-elite_proportion));
            elseif elite_type == 3 % top elite percent by rank
                [~, fit_I] = sort(V_loc(I,:),'descend');
                limit = ceil(length(I)*elite_proportion);
                I2 = I(fit_I(1:limit));
            else % top elite by number
                [~, fit_I] = sort(V_loc(I,:),'descend');
                limit = min(elite_max_number, length(I));
                I2 = I(fit_I(1:limit));
            end
            % special case at end of run
            if (length(I2) + number_of_new_locations + evaluations) > max_evaluations
                reduced_number = max_evaluations-(evaluations+number_of_new_locations);
                I2 = I2(1:reduced_number); 
            end
            
            for jj=1:length(I2)
                [active_modes(I2(jj)).local_region] = ...
                    select_new_location_close_to_current_mode(active_modes(I2(jj)).local_region,rsg_params,M_loc,tol_value);
            end
            [active_modes,active_modes_changed,temp_new_loc_num,M_loc,V_loc,X,Y,index] = ...
                evaluate_new_locations(active_modes,problem_func,M_loc,V_loc,problem_function_params,X,Y,index,I2);
            
            temp = temp+active_modes_changed;     
            I = 1:length(active_modes);
            internal_its = internal_its + length(I2);
            internal_count = internal_count+1;
            number_of_new_locations = number_of_new_locations+temp_new_loc_num;
        end
        active_modes_changed = temp;
    end
    evaluations = evaluations+number_of_new_locations;
    
    % evaluate the proposed locations
    % evolve new locations to look for modes
    number_of_new_modes=0;
    if (evaluations < max_evaluations)       
        [active_modes,active_modes_changed,number_of_new_modes,M_loc,V_loc,t,X,Y,index] = ...
            evolve(active_modes,active_modes_changed,problem_func,M_loc,V_loc,t,rsg_params,gens_per_crossover,problem_function_params,X,Y,index,max_evaluations-evaluations);
    end
    evaluations = evaluations+number_of_new_modes;
    % evolve random locations
    number_rand_modes=0;
    if (evaluations < max_evaluations)       
        [active_modes,active_modes_changed,number_rand_modes,M_loc,V_loc,X,Y,index] = ...
            random_new(active_modes,active_modes_changed,problem_func,M_loc,V_loc,random_solution_generator,rsg_params,problem_function_params,X,Y,index);
    end
    % keep track of how many calls have been made to the cost function
    evaluations = evaluations+number_rand_modes;
    %ss = sort(V_loc,'descend');
    %ss(1:min(4,length(ss)))'
    fprintf('pop size %d, iteration %d, evals %d, elite subset reps %d Best peak estimate %f\n',...
        length(active_modes),iteration,evaluations,internal_count, max(V_loc));
    [active_modes] = trim_history(active_modes,max_hist);
    %ob = max(V_loc);
    state = vertcat(state, length(active_modes));
    iteration = iteration+1;
end

[RES,RES_Y] = extract_modes(active_modes);




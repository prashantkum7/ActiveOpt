function [cost,time,status] = experiment_cost(plasmid,status)

t_plasmid = 4; % Time (in days) to build a plasmid, get primers and check/confirm the plasmid
t_expt = 4; % Time (in days) to insert a plasmid, perform the experiment and do the analysis
plasmid_cost = 33; % Swap RBS; primers ($25); to build a plasmid ($7) + $1 to confirmin $

% Status of the plasmid tells whether the plasmid has been already made or
% needs to be made

if(status == 0) % plasmid not constructed before
    time = t_plasmid + t_expt;
    cost = plasmid_cost;
    status = 1; % Updating the status because the plasmid will be constructed now
else % plasmid was already made before--> so only time is spent on performing experiments
    time = t_expt;
    cost = 0;
end

end
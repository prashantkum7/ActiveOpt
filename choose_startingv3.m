function [chosen_cases,chosen,cutoff,norm_best_yield,n_exp] = choose_startingv3(nsample,yield,strength,diff_yield,norm_best_yield,opt_sampling,n_exp,best_yield,run)

flag = 1;
chosen = 1; % Number of experiments required to get data from two different class labels


if(opt_sampling ~= 1)
chosen_cases(1:nsample/2,1) = randi([1 size(yield,1)]);
n_exp=n_exp+1; % 1 experiment chosen
norm_best_yield(n_exp,run) = max(yield(chosen_cases,1))*100/best_yield;
while(flag)
    a=zeros(size(yield,1),2);
    for i=1:size(yield,1)
        a(i,1) = i;
        for j = 1:chosen
            a(i,2) = a(i,2) + sqrt(sum((strength(chosen_cases(j,1),:)-strength(i,:)).^2));
        end
    end
    a = sortrows(a,-2);
    
    % Choosing the furthest point from the already chosen points
    chosen_cases=[chosen_cases;a(1,1)];
    n_exp=n_exp+1; % 1 experiment chosen
    norm_best_yield(n_exp,run) = max(yield(chosen_cases,1))*100/best_yield;
    chosen = size(chosen_cases,1);
    k1 = max(yield(chosen_cases(1:chosen),1));
    k2 = min(yield(chosen_cases(1:chosen),1));
    
    if(k1-k2>=diff_yield)
        flag = 0;
    end
    clearvars a;
end

else
    chosen_cases = randi(size(yield,1));
    n_exp=n_exp+1; % 1 experiment chosen
    norm_best_yield(n_exp,run) = max(yield(chosen_cases,1))*100/best_yield;
    chosen = size(chosen_cases,1);
    k1 = 0;
    k2 = 0;
end

cutoff = 0.5*(k1+k2);

% chosen_cases(nsample/2+1:nsample,1) = randi([high_cases+1 size(yield,1)]);


end
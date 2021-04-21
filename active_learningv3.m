% This code is same as v2 except that I have added one more feature for the
% truncated plasmids


function [nstart,cutoff_run,result_validation,yield,pseudo_sv,wt_features,max_fraction,starting_cases,starting_yield,strength,Ratio1,time_iter,itr_high_yield,norm_best_yield,plasmid_comb] = active_learningv3()

clear; clc; close all; stopflag = 1; run = 0;
opt_sampling = 2; % 1 - Random, 2 - furthest from the hyperplane, 3 - closest from the hyperplane
cutoff = 0.29; % Cut-off to label the data
diff_yield = 0.05; %  Minimum 5% difference in yield of the two classes
% Load Data
noise=0; % 1: Incorporating noise in the RBS strength value; 0: No noise
[plasmids,genes,RBS_strength,Valine_yield,std_rbs]=dataprocess(noise);
% Storing the valine yield and RBS strength
[yield,strength,plasmid_comb] = process_yield(Valine_yield,RBS_strength,plasmids,cutoff);
% yield = yield*13.32/11.4; % Conversion Factor

%% Incorporating the information for the truncated promoters
[prom_flag] = trunc_promoter();
runs=100;
norm_best_yield=zeros(93,runs);
%%
% Identify the top three high yield cases
for temp = 1:size(yield,1)
    a(temp,1)= temp;
    a(temp,2)=yield(temp,1);
end

a=sortrows(a,-2);
highest_yield=a(1:3,1); % Top three high yield cases
best_yield = a(1,2); % Best yield across all the experiments

for run = 1:runs
    % while(stopflag)
    run
    clearvars all_cases remaining_cases group training c_yd c_st validation c_prom c_feature;
    clearvars hyield lyield chosen_original chosen svmstruct;
    % Data===========================================================
    flag = 1;
    nsample = 2; % Number of experiments to start with
    accuracy = 0;
    iter = 0;
    n_exp=0;
    total_time = 0; total_cost = 0;
    nchoose = 1; % Number of experiments to choose in each iteration
    iter_highest=1000*ones(3,1);
    %================================================================
    RBS=RBS_strength; bias=0;
    yield_original=yield;
    
    high_cases=sum(yield(:,1)>cutoff);
    
    for temp = 1: size(yield,1)
        all_cases(temp,1) = temp;
    end
    
    chosen_cases=0; svmstruct=0;
    
    pseudo_sv=zeros(size(RBS,2)+size(prom_flag,2),3);
    k=[1,2,4,5,6,11]';
    for temp=1:size(RBS,2)+size(prom_flag,2)
        pseudo_sv(temp,1)=temp;
        %     pseudo_sv(temp,1)=k(temp,1);
    end
    
    while(flag)
        iter = iter + 1;
        if(iter==1)
            remaining_cases = all_cases;  % All cases are available
            
           [chosen_cases,chosen,cutoff,norm_best_yield,n_exp] = choose_startingv3(nsample,yield,strength,diff_yield,norm_best_yield,opt_sampling,n_exp,best_yield,run);
            nstart(run,1)=chosen;
            starting_cases(1:chosen,run) = chosen_cases(1:chosen,1);
            starting_yield(1:chosen,run) = yield(chosen_cases(1:chosen,1),1);
            
            % High and Low Yield groups
            cnt_h=0; cnt_l=0;
            for temp = 1:size(yield,1)
                if(yield(temp,1)>=cutoff)
                    cnt_h = cnt_h + 1; hyield(cnt_h,1) = temp;
                else
                    cnt_l = cnt_l + 1; lyield(cnt_l,1) = temp;
                end
            end
            
            % Time and cost to construct these plasmids
            for i = 1:size(nsample)
                plasmid1 = plasmid_comb(starting_cases(i,run),1); % Plasmid 1
                status = cell2mat(plasmids(strmatch(plasmid1,plasmids(:,1),'exact'),2));
                [cost1,time1,status] = experiment_cost(plasmid1,status);
                plasmids(strmatch(plasmid1,plasmids(:,1),'exact'),2) = num2cell(status); % Updating the plasmid status
                
                plasmid2 = plasmid_comb(starting_cases(i,run),2); % Plasmid 2
                status = cell2mat(plasmids(strmatch(plasmid2,plasmids(:,1),'exact'),2));
                [cost2,time2,status] = experiment_cost(plasmid2,status);
                plasmids(strmatch(plasmid2,plasmids(:,1),'exact'),2) = num2cell(status); % Updating the plasmid status
                
                total_cost = total_cost + cost1 + cost2; % Total cost to build and test the plasmid;do the experiments
                total_time = total_time + time1 + time2; % Total time to build and test the plasmid;do the experiments
            end
            chosen_original = chosen_cases;
            
            for temp = 1:size(chosen_cases,1)
                c_yd(temp,1)=yield(chosen_cases(temp,1),1);
                c_st(temp,:)=strength(chosen_cases(temp,1),:);
                c_prom(temp,1) = cell2mat(prom_flag(find(strcmp(plasmid_comb(chosen_cases(temp,1),1),prom_flag)==1),2));
                c_prom(temp,2) = cell2mat(prom_flag(find(strcmp(plasmid_comb(chosen_cases(temp,1),2),prom_flag)==1),2));
            end
            
            % Final Feature Matrix
            c_feature = zeros(size(chosen_cases,1),size(RBS,2)+size(prom_flag,2));
            
            c_feature = [c_st c_prom];
            
            %Classifying the valine yield as high (1)-low (0)
            group = zeros(size(chosen_cases,1),1);
            for temp=1:size(group,1)
                if(c_yd(temp,1)>=cutoff)
                    group(temp,1)=1;
                else
                    group(temp,1)=2;
                end
            end
            
            % Running SVM on the chosen cases
            if(opt_sampling ~= 1)
                svmstruct=svmtrain(c_feature,group,'kernel_function','linear','autoscale','false','ShowPlot',false);
                remaining_cases = setxor(all_cases,chosen_cases);
            else
                svmstruct = [];
                remaining_cases = setxor(all_cases,chosen_cases);
            end
            %            [chosen,total_high,flag] = choose_best(svmstruct,remaining_cases,strength,nchoose,yield,cutoff,flag); % Farthest from hyperplane
            %          [chosen,total_high,flag] = choose(svmstruct,remaining_cases,strength,nchoose,yield,cutoff,flag); % Closest from hyperplane
            
        end
        
        cutoff_run(iter,run)=cutoff;  % Updated Cutoff after each run
        
        clearvars c_st group c_yd;
        %Updating the chosen experiments, cost/time and the status of the plasmids
        [chosen_cases,remaining_cases,group,c_yd,c_st,c_prom,flag,total_cost,total_time,plasmids,chosen] = select_cases(all_cases,chosen_cases,remaining_cases,nchoose,yield,strength,iter,svmstruct,nsample,cutoff,flag,total_time,total_cost,plasmids,plasmid_comb,opt_sampling,prom_flag);
        n_exp=n_exp+1; % 1 experiment chosen
        norm_best_yield(n_exp,run) = max(yield(chosen_cases,1))*100/best_yield;
        % Final Feature Matrix
        
        c_feature = zeros(size(chosen_cases,1),size(RBS,2)+size(prom_flag,2));
        c_feature = [c_st c_prom];
        % Updating cutoff after each iteration
        k1 = max(yield(chosen_cases(1:size(chosen_cases,1),1),1));
        k2 = min(yield(chosen_cases(1:size(chosen_cases,1),1),1));
        cutoff = 0.5*(k1+k2);
        
        if(flag==1)
            % Running SVM on the chosen cases
            if(opt_sampling ~= 1)
                svmstruct=svmtrain(c_feature,group,'kernel_function','linear','autoscale','false','ShowPlot',false);
            else
                svmstruct = [];
            end
            
            for i = 1:nchoose
                if(yield(chosen(i,1))>=cutoff && iter > 1)
                    accuracy = accuracy + 1;
                end
            end
            
            validation(iter,1)=accuracy;
            
            if(opt_sampling ~= 1)
                pseudo_sv(:,2) = pseudo_sv(:,2) + (svmstruct.Alpha'*svmstruct.SupportVectors)';
                pseudo_sv(:,3) = (svmstruct.Alpha'*svmstruct.SupportVectors)';
                wt_features(:,run) = (svmstruct.Alpha'*svmstruct.SupportVectors)';
                bias = bias + svmstruct.Bias;
            else
                pseudo_sv(:,2) = 0;
                pseudo_sv(:,3) = 0;
                wt_features(:,run) = 0;
                bias = 0;
            end
            
        end
        %     classify_test=svmclassify(svmstruct,strength(test(1:ntest,1),:));
        
        if(or(iter==(high_cases-nsample/2),and(~flag,iter<(high_cases-nsample/2))))
            Ratio1(run,1) = trapz(validation)/((high_cases-nsample/2)*(high_cases-nsample/2)*0.5);
        end
        
        cost_iter(iter,run) = total_cost;
        time_iter(iter,run) = total_time;
        
        % Storing the normalized best yield after each iteration
        norm_best_yield(n_exp,run) = max(yield(chosen_cases,1))*100/best_yield;
        
    end
    
    %  ITER ENDS
    % Iterations required to reach the highest yield case(s)
    %         if(size(intersect(chosen_cases,highest_yield)))
    
    for k1 = 1:size(highest_yield,1)
        t = highest_yield(k1,1);
        row = find(chosen_cases==t);
        if(~isempty(row))
            iter_highest(k1,1)=row;
        end
        clearvars row;
    end
    %         end
    
    %     if(opt_sampling ~= 1)
    %         result = prediction(svmstruct,plasmids,RBS_strength,prom_flag,plasmid_comb);
    %         pseudo_sv(:,2)=pseudo_sv(:,2)/iter;
    %         bias=bias/iter;
    %     else
    %         result = 0;
    %         pseudo_sv(:,2)=0;
    %         bias=0;
    %     end
    
    for temp =1:iter-1
        %         result_validation(temp,run)=validation(temp,1)/total_high; %
        result_validation(temp,run)=validation(temp,1)/(high_cases-nsample/2); % nsample/2 high cases were chosen as the starting cases
    end
    
    %     max_fraction(run,1) = max(validation)/total_high;
    max_fraction(run,1) = max(validation)/(high_cases-nsample/2);
    max_fraction(run,2) = iter-1;
    itr_high_yield(:,run) = iter_highest;
    
end

% If the high yield strain is not found, then assign 0 iteration to it
itr_high_yield(itr_high_yield==1000)=0;
% figure;
% plot(result_validation,'*','MarkerSize',10);
% hold on;
% xlabel('Iteration','FontSize',12,'FontWeight','bold');
% ylabel('Fraction of Correct High Yield Identification','FontSize',12,'FontWeight','bold');
% figure;
% histogram(Ratio1)
% ylabel('Counts','FontSize',12,'FontWeight','bold');
% xlabel('Area Ratio','FontSize',12,'FontWeight','bold');
save('temp.mat');

end


function [left_plasmid,svmstruct,distance,train,yield,RBS_leftplasmid,Rank_pseudo_sv] = predict_instance()

clear;
clc;
close all;

% load('data_v2.mat');
noise=0; % 1: 10% noise incorporated in the RBS strength value; 0: No noise
[plasmids,genes,RBS_strength,Valine_yield,std_rbs]=dataprocess(noise);  % Standardizing the RBS strength

ntrain=80; %Number of data points in the training set
itermax=0; %Number of iterations
strength=0;
cutoff=0.29; %Cut-off to classify valine yield
option2=0;  % To correct for evaporation in the old samples
sampling='K-fold'; %K-fold : Choose the sampling type

RBS=RBS_strength;

tol=0.25;

%Generating a test case
%Parameters for post-processing
Rank_pseudo_sv=zeros(size(genes,1),2*itermax);

pseudo_sv=zeros(size(RBS,2),2);
k=[1,2,4,5,6,11]';
for temp=1:size(RBS,2)
    pseudo_sv(temp,1)=temp;
end

% Storing the valine yield and RBS strength
% [yield,strength]=process_yield(Valine_yield,RBS,yield_ev_marker);
[yield,strength,plasmid_comb] = process_yield(Valine_yield,RBS_strength,plasmids,cutoff);
% yield = yield*13.32/11.4;  % Conversion Factor
yield_original=yield;

% Simulating 1000 independent runs
for ct=1:1
    ct
    iter=0;
    final_accuracy=zeros(2,2);
    TP=0;
    TN=0;
    FP=0;
    FN=0;
    bias=0;
    
    %Classifying the valine yield as high (1)-low (0)
    for temp=1:size(yield,1)
        if(yield(temp,1) >= cutoff)
            group(temp,1) = 1;
        else
            group(temp,1) = 2;
        end
    end
    
    for temp=1:size(group,1)
        nsample(temp,1)=temp;
    end
    
    if(sampling == 'K-fold')
        [Indices] = crossvalind('Kfold', size(nsample,1),91);
    else
        [A,B]=crossvalind('HoldOut',size(nsample,1),0.8);
    end
    
    itermax = itermax + max(Indices);
    iter = 0;
    while(iter<1)
        iter=iter+1;
        count=0;
        clearvars train test svmstruct classify_test A B;
        
        % Randomly picking the yield for combinations not corrected for
        % evaporation
        if(option2 == 1)
            [yield,group]=evp_correction(yield_original,cutoff,yield);
        end
        
        %Dividing the data into training and testing set
        %Choosing a training and a testing set
        test = (Indices == iter); train =~ test;
        
        %         % Temporary adjustment
        %         tu = test; test = train; train = tu;
        
        A=find(train==1);
        B=find(test==1);
        
        ntrain=size(A,1);
        ntest=size(B,1);
        
        clearvars train test;
        for temp=1:ntrain
            train(temp,1)=nsample(A(temp,1),1);
        end
        
        for temp=1:ntest
            test(temp,1)=nsample(B(temp,1),1);
        end
        
        train = [train;test];
        ntrain = ntrain + 1;
        
        %Identify the experimental LOW and HIGH in the train set
        train_low=sum(yield(train(1:ntrain),1) < cutoff);
        train_high=sum(yield(train(1:ntrain),1) >= cutoff);
        
        %Classifier
        %         svmstruct=(svmtrain(strength(train(1:ntrain),:),group(train(1:ntrain),1),'Kernel_Function','rbf','rbf_sigma',1,'Autoscale','false'));
        svmstruct=svmtrain(strength(train(1:ntrain),:),group(train(1:ntrain),1),'kernel_function','linear','autoscale','false');
        bias=bias + svmstruct.Bias;
        
        % Creating a pseudo-support vector matrix
        % How much one RBS (gene) affect the classifier
        for j = 1:size(svmstruct.SupportVectors,2)
            pseudo_sv(j,2)= pseudo_sv(j,2) + svmstruct.Alpha'*svmstruct.SupportVectors(:,j); % Input is standardized
            % pseudo_sv(temp,2)= pseudo_sv(temp,2) + svmstruct.ScaleData.scaleFactor(1,temp)*(svmstruct.SupportVectors(temp1,temp)+svmstruct.ScaleData.shift(1,temp))*svmstruct.Alpha(temp1,1);
        end
    end
    
    avg_bias(ct,1)=bias/itermax;
    
    %Ranking the features
    pseudo_sv(:,2)=pseudo_sv(:,2)/itermax;
    Rank_pseudo_sv=sortrows(pseudo_sv,-2);
    
    count=0;
    for temp = 1:size(Rank_pseudo_sv,1)
        if(abs(Rank_pseudo_sv(temp,2))>0)
            count = count+1;
            Rank_pseudo_sv(count,1)=Rank_pseudo_sv(temp,1);
            Rank_pseudo_sv(count,2)=Rank_pseudo_sv(temp,2);
        end
    end
    
    Rank_pseudo_sv(count+1:size(Rank_pseudo_sv,1),1)=0;
    Rank_pseudo_sv(count+1:size(Rank_pseudo_sv,1),2)=0;
end

[distance,left_plasmid,RBS_leftplasmid]=ensemble_svm(svmstruct);

end



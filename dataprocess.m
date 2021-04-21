function [plasmids,genes,RBS_strength,yield,std_rbs]=dataprocess(noise)

% Reading the input file and storing plasmids | genes | RBS_strength
% filename = 'Matlab Input_v2.xlsx';
% Sheet = 'Input_v3';
% xlRange ='N3:Y46';

filename = 'Matlab Input_v3.xlsx';
Sheet = 'Input_v4';
xlRange ='O3:Z47';

[A,B,C]=xlsread(filename,Sheet,xlRange,'basic');

for temp = 1:size(C,1)-1
    plasmids(temp,1)= C(temp+1,1);
end
plasmids(:,2) = num2cell(0); % Status of the plasmid -> zero means the plasmid is yet to be made
% 1 means that the plasmid was already made before

for temp = 1:size(C,2)-1
    genes(temp,1)= C(1,temp+1);
end

for temp = 1:size(A,1)
    for temp1 = 1:size(A,2)
        if(~noise)
            RBS_strength(temp,temp1)=A(temp,temp1); % Value from the Salis Calculator
        else
            RBS_strength(temp,temp1)=0.9*A(temp,temp1)+rand(1,1)*(0.20*A(temp,temp1)); % Value randomly picked from RBS +- 10%
        end
    end
end

% Standardize the RBS Strength
for temp = 1:size(RBS_strength,2)
    count=0;
    for temp1 = 1:size(RBS_strength,1)
        if(isnan(RBS_strength(temp1,temp))==0)
            count=count+1;
            strength(count,1)=RBS_strength(temp1,temp);
        end
    end
    avg=mean(strength);
    sd=std(strength);
    std_rbs(temp,1)=sd;
    RBS_strength(:,temp) = (RBS_strength(:,temp) - avg)/sd;
    clearvars strength;
end

clearvars A B C;

% Reading the input file and storing Valine Yield
% xlRange ='B3:K37';
xlRange ='B3:K38';
[A,B,C]=xlsread(filename,Sheet,xlRange,'basic');

yield=zeros(size(plasmids,1),size(plasmids,1))*NaN;

for temp = 1:size(C,1)
    for temp1 = 1:size(C,2)
        [row1,col1]=find(strcmp(C(temp,size(C,2)),plasmids));
        [row2,col2]=find(strcmp(C(size(C,1),temp1),plasmids));
        if(isempty(row1)~=1 && isempty(row2)~=1)
            yield(row1,row2)=cell2mat(C(temp,temp1));
            clearvars row1 col1 row2 col2;
        end
    end
end
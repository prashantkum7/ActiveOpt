function [yield,strength,plasmid_comb] = process_yield(Valine_yield,RBS,plasmids,cutoff)

count=0;
for temp=1:size(Valine_yield,1)
    for temp1=1:size(Valine_yield,2)
        if(isnan(Valine_yield(temp,temp1))==0)
            count=count+1;
            %Storing the valine yield in the order of the plasmid
            %combinations
            t_yield(count,1)=Valine_yield(temp,temp1);
            t_plasmid_comb(count,1) = plasmids(temp,1);
            t_plasmid_comb(count,2) = plasmids(temp1,1);
            %RBS Strength
            for i = 1:size(RBS,2)
                a=0;
                if(isnan(RBS(temp,i))==0)
                    a=a+RBS(temp,i);
                end
                
                if(isnan(RBS(temp1,i))==0)
                    a=a+RBS(temp1,i);
                end
                %RBS strength per gene in the plasmid
                t_strength(count,i)=a;
            end
        end
    end
end

% Sorting and arranging the yields and plasmid_combs
c1 = 0;
c2 = 0;
for i=1:size(t_yield,1)
    if(t_yield(i,1)>=cutoff)
        c1=c1+1;
        yield1(c1,1)=t_yield(i,1);
        plasmid_comb1(c1,1) = t_plasmid_comb(i,1);
        plasmid_comb1(c1,2) = t_plasmid_comb(i,2);
        strength1(c1,:) = t_strength(i,:);
    else
        c2=c2+1;
        yield2(c2,1)=t_yield(i,1);
        plasmid_comb2(c2,1) = t_plasmid_comb(i,1);
        plasmid_comb2(c2,2) = t_plasmid_comb(i,2);
        strength2(c2,:) = t_strength(i,:);
    end
end
yield=[yield1;yield2];
plasmid_comb=[plasmid_comb1;plasmid_comb2];
strength = [strength1;strength2];

end
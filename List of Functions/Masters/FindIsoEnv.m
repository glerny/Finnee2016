function IE = FindIsoEnv(Clust, int)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


Clust = sortrows(Clust, 6);
IE = {};

while ~isempty(Clust)
    % Check 0.33 - 0.66
    
    Id = find(Clust(:,6)-Clust(1,6) >= int(1)/3 & Clust(:,6)-Clust(1,6) <= int(2)/3 | ...
        Clust(:,6)-Clust(1,6) >= 2*int(1)/3 & Clust(:,6)-Clust(1,6) <= 2*int(2)/3 );
    if any(Id)
        
        LId = 1;
        while any(Id)
            LId = [LId; Id];
            Id = find(Clust(:,6)-Clust(LId(end),6) >= int(1)/3 & Clust(:,6)-Clust(LId(end),6) <= int(2)/3 | ...
                Clust(:,6)-Clust(LId(end),6) >= 2*int(1)/3 & Clust(:,6)-Clust(LId(end),6) <= 2*int(2)/3 );
        end
        IE{end+1} = Clust(LId,:);
        Clust(LId,:) = [];
        
    else
        Id = find(Clust(:,6)-Clust(1,6) >= int(1)/2 & Clust(:,6)-Clust(1,6) <= int(2)/2);
        if any(Id)
            
            LId = 1;
            while any(Id)
                LId = [LId; Id];
                Id = find(Clust(:,6)-Clust(LId(end),6) >= int(1)/2 & Clust(:,6)-Clust(LId(end),6) <= int(2)/2);
            end
            
            IE{end+1} = Clust(LId,:);
            Clust(LId,:) = [];
            
        else
            Id = find(Clust(:,6)-Clust(1,6) >= int(1) & Clust(:,6)-Clust(1,6) <= int(2) | ...
                Clust(:,6)-Clust(1,6) >= 2*int(1) & Clust(:,6)-Clust(1,6) <= 2*int(2));
            if any(Id)
                
                LId = 1;
                while any(Id)
                    LId = [LId; Id];
                    Id = find(Clust(:,6)-Clust(LId(end),6) >= int(1) & Clust(:,6)-Clust(LId(end),6) <= int(2) | ...
                        Clust(:,6)-Clust(LId(end),6) >= 2*int(1) & Clust(:,6)-Clust(LId(end),6) <= 2*int(2) );
                end
                IE{end+1} = Clust(LId,:);
                Clust(LId,:) = [];
                
            else
                Clust(1,:) = [];
            end
        end
    end
end


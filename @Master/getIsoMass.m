function obj = getIsoMass(obj, polarity)

% 1- Initialisation and options
interval = [0.99 1.01];
maxCharge = 3;
Sum = obj.HACA.Summary;
IsoIons = {};
for ii = 1:max(Sum(:,2))
    cCluster = Sum(Sum(:,2) == ii, :);
    if size(cCluster, 1) > 1
        List = getDist(cCluster);
        
        for jj = [maxCharge:-1:1, 0.5]
            [l, ~] = find(List >= interval(1)/jj & ...
                List <= interval(2)/jj);
            
            while ~isempty(l)
                if isempty(cCluster), break; end
                cIons = [];
                ICharge = jj;
                [~, i2s] = max(cCluster(l, 8));
                cIons = cCluster(l(i2s), :);
                cCluster(l(i2s), :) = [];
                
                i2add = find(abs(cCluster(:, 6) - cIons(1,6))  >= interval(1)/jj & ...
                    abs(cCluster(:, 6) - cIons(1,6))  <= interval(2)/jj);
                while ~isempty(i2add)
                    cIons = [cIons; cCluster(i2add, :)];
                    cCluster(i2add, :) = [];
                    
                    i2add = [];
                    for kk = 1:size(cIons, 1)
                        i2add = [i2add; find(...
                            abs(cCluster(:, 6) - cIons(kk,6))  >= interval(1)/jj & ...
                            abs(cCluster(:, 6) - cIons(kk,6))  <= interval(2)/jj)];
                    end
                end
                
                IsoIons{end+1}.charge = ICharge;
                IsoIons{end}.IsotopicEnvelope = cIons;
                if isempty(cCluster), break; end
                List = getDist(cCluster);
                [l, ~] = find(List >= interval(1)/jj & ...
                    List <= interval(2)/jj);
            end
        end
    end
end


ListIons = [];
switch polarity
    case '+'
        MonoisotopicMasses.polarity = 'positive';
        for ii = 1:length(IsoIons)
            IE = IsoIons{ii}.IsotopicEnvelope;
            IE = sortrows(IE, -8);
            ListIons(ii,1) = ii;
            ListIons(ii,2) = IE(1,2);
            if IsoIons{ii}.charge == 0.5
                IsoIons{ii}.charge = 1;
            end
            ListIons(ii,3) = IsoIons{ii}.charge*(IE(1,6) - 1.007276);
            ListIons(ii,4) =  IsoIons{ii}.charge*IE(1,7);
            ListIons(ii,5) =  IE(1,3);
            ListIons(ii,6) =  IE(1,4);
            IE = [IE(:, 1:6), IsoIons{ii}.charge*(IE(:, 6)-1.007276), IE(:, 7:end)];
            IsoIons{ii}.IsotopicEnvelope = IE;
        end
    case '-'
        MonoisotopicMasses.polarity = 'negative';
        for ii = 1:length(IsoIons)
            IE = IsoIons{ii}.IsotopicEnvelope;
            IE = sortrows(IE, -8);
            ListIons(ii,1) = ii;
            ListIons(ii,2) = IE(1,2);
            if IsoIons{ii}.charge == 0.5
                IsoIons{ii}.charge = 1;
            end
            ListIons(ii,3) = IsoIons{ii}.charge*(IE(1,6) + 1.007276);
            ListIons(ii,4) =  IsoIons{ii}.charge*IE(1,7);
            ListIons(ii,5) =  IE(1,3);
            ListIons(ii,6) =  IE(1,4);
            IE = [IE(:, 1:6), IsoIons{ii}.charge*(IE(:, 6)+1.007276), IE(:, 7:end)];
            IsoIons{ii}.IsotopicEnvelope = IE;
        end
end
MonoisotopicMasses.ListIons = IsoIons;
MonoisotopicMasses.Summary = ListIons;
obj.HACA.MonoisotopicMasses = MonoisotopicMasses;
myMaster = obj;
try
    save(fullfile(obj.Path, obj.Name), 'myMaster')
catch
    uisave('myMaster')
end





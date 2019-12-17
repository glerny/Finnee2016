function List = getDiff(obj)

% 1- Initialisation and options
List = [];
Sum = obj.HACA.Summary;
for ii = 1:max(Sum(:,2))
    tgtCluster = Sum(Sum(:,2) == ii, :);
    if size(tgtCluster, 1) > 1
        l1 = getDist(tgtCluster);
        List = [List; l1(:)];
        List = unique(List);
    end
end


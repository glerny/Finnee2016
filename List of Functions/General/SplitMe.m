function [IdSplit, PearSplit] = SplitMe(myMat)

RMat = corrcoef(myMat);
PearSplit = min(RMat, [], 'all');
[l, m] = find(RMat == PearSplit);

[~, IdSplit] = min(pdist2(myMat(:, l)', myMat', 'correlation'), [], 1);
end


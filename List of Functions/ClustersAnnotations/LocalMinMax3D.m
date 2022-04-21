function [LocMax, LocMin]  = LocalMinMax3D(data, m, n)
thres = 0.05;

[Lin, Col] = size(data);
hyperData = nan(Lin+2*m, Col+2*n, m*n);

cnt = 1;
for ii = m:-1:-m
    for jj = n:-1:-n
        hyperData(m+1-ii:m+Lin-ii, n+1-jj:n+Col-jj, cnt) = data;
        cnt = cnt+1;
    end
end

hypMax = max(hyperData(m+1:m+Lin, n+1:n+Col, :), [], 3);
LocalMax = data == hypMax & data >= thres*(max(max(data)));
[lm, cm] = find(LocalMax);
LocMax = [lm, cm, data(LocalMax)];

hypMin = min(hyperData(m+1:m+Lin, n+1:n+Col, :), [], 3);
LocalMin = data == hypMin;
[lm, cm] = find(LocalMin);
LocMin = [lm, cm, data(LocalMin)];


end


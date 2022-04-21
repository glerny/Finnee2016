function DataOut = ConcatMe(DataIn, Bucket)
[m, n] = size(DataIn);

DataOut = [];
ip = 1;
for ii = Bucket:Bucket:n
    DataOut(:, end+1) = mean(DataIn(:, ip:ii-1), 2);
    ip = ii;
end

end


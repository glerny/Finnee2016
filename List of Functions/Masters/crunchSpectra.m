function XYout = crunchSpectra(XYin, int)

id = [1; find(diff(XYin(:,1)) > int)+1];

for ii = 1:length(id)-1
    XYout(ii,1) = sum(XYin(id(ii):id(ii+1)-1,1).*XYin(id(ii):id(ii+1)-1,2))...
        /sum((XYin(id(ii):id(ii+1)-1,2)));
    XYout(ii,2) = sum((XYin(id(ii):id(ii+1)-1,2)));
end


end


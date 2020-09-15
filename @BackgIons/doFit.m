function obj = doFit(obj)
thresPts = 5;

AxisX = obj.ROIs{1}.AxisTm.Data;

h = waitbar(0,'Calculating fitting coefficient, please wait');
for ii =1:length(AxisX)
    waitbar(ii/length(AxisX))
    Dt4c = zeros(length(obj.ROIs), 2);
    for jj = 1:length(obj.ROIs)
        XY = [obj.ROIs{jj}.AxisMZ.Data, obj.ROIs{jj}.StoredData(:, ii)];
        
        if nnz(XY(:,2)) > thresPts
            Dt4c(jj, 1) = obj.ROIs{jj}.TgtMz;
            Res = LocalMaxima(XY, 3, 0);
            if size(Res, 1) == 1
                Dt4c(jj, 2) = Res(1);
            end
        end
    end
    
    Dt4c(Dt4c(:,2) == 0, :) = [];
    Dt4c(:, 3) = (Dt4c(:,2)-Dt4c(:,1))./Dt4c(:,2);
    if size (Dt4c, 1) >= 2*obj.n+1
        [p, ~] = polyfitweighted(Dt4c(:,2),Dt4c(:,3), obj.n);
        P(ii,:) = p;
    else
        fprintf('\pscan %i not corrected', ii)
        P(ii,:) = [0 , 0];
    end
    
end
try close(h); catch, end %#ok<CTCH>

obj.P = P;
save(obj.path, 'obj');





function [p, S, mu, nIdO] = polyfit_NOTL(x, y, n, alfa, sim) 
% NOTE: note from polyval Standard error for prediction, returned as a
% vector of the same size as the query points x. Generally, an interval of
% y ± Δ corresponds to a roughly 68% prediction interval for future
% observations of large samples, and y ± 2Δ a roughly 95% prediction
% interval.

nIdO = isfinite(x) & isfinite(y);
[p,S,mu] = polyfit(x(nIdO),y(nIdO),n);
[y_, delta] = polyval(p, x, S, mu);

while 1
    
    switch sim
        case 'both'
            n2IdO =  nIdO & abs(y-y_) <= alfa*delta;
            
        case 'up'
            n2IdO =  nIdO & y-y_ <= alfa*delta;
            
        case 'down'
            n2IdO =  nIdO & y_-y <= alfa*delta;
            
    end
    
    if sum(n2IdO) == sum(nIdO)
        break
    else
        nIdO = n2IdO;
    end
    
    [p,S,mu] = polyfit(x(nIdO),y(nIdO),n);
    [y_, delta] = polyval(p, x, S, mu);
end

end


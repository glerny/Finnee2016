function profOut = extrapol2axis(profIn, axis)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

axis(:,3)   = 1;
axis(:,4)   = 1:length(axis(:,1));
profIn(:,3) = -1;
profIn(:,4) = 1:length(profIn(:,1));
Dt2wk = sortrows([axis; profIn], 1);
in2rem = Dt2wk(:,1) < axis(1,1) | Dt2wk(:,1) > axis(end,1);
Dt2wk(in2rem, :) = [];
Dt2wk = [zeros(5,4);Dt2wk; zeros(5, 4)];

ind2chk  = find(Dt2wk(2:end-1, 3) == 1 & (Dt2wk(1:end-2, 3) == -1 | Dt2wk(3:end, 3) == -1))+1;
spd      = Dt2wk(ind2chk, :);
spd      = [spd, Dt2wk(ind2chk-1, :)];
spd      = [spd, Dt2wk(ind2chk-2, :)];
spd      = [spd, Dt2wk(ind2chk-3, :)];
spd      = [spd, Dt2wk(ind2chk-4, :)];
spd      = [spd, Dt2wk(ind2chk-5, :)];
ind      = find(spd(:,7) == 1);
while ~isempty(ind)
    spd(ind, 5:end-4) = spd(ind, 9:end);
    spd(ind, end-3:end) = 0;
    ind      = find(spd(:,7) == 1);
end

ind      = find(spd(:,11) == 1);
while ~isempty(ind)
    spd(ind, 9:end-4) = spd(ind, 13:end);
    spd(ind, end-3:end) = 0;
    ind      = find(spd(:,11) == 1);
end
aMin = spd(:,[1:6, 9,10]);

spd      = Dt2wk(ind2chk, :);
spd      = [spd, Dt2wk(ind2chk+1, :)];
spd      = [spd, Dt2wk(ind2chk+2, :)];
spd      = [spd, Dt2wk(ind2chk+3, :)];
spd      = [spd, Dt2wk(ind2chk+4, :)];
spd      = [spd, Dt2wk(ind2chk+5, :)];

ind      = find(spd(:,7) == 1);
while ~isempty(ind)
    spd(ind, 5:end-4) = spd(ind, 9:end);
    spd(ind, end-3:end) = 0;
    ind      = find(spd(:,7) == 1);
end

ind      = find(spd(:,11) == 1);
while ~isempty(ind)
    spd(ind, 9:end-4) = spd(ind, 13:end);
    spd(ind, end-3:end) = 0;
    ind      = find(spd(:,11) == 1);
end
aPls = spd(:,[1:6, 9,10]);

ind2rem = aMin(:,5) == 0 | aPls(:,5) == 0;
aMin(ind2rem, :) = [];
aPls(ind2rem, :) = [];
TT = [aMin(:,4), aMin(:, 7:8), aMin(:, 5:6),aPls(:, [1,2, 5:8])];
a = (TT(:,3)-TT(:,5))./((TT(:,2)-TT(:,4)).*(TT(:,4)-TT(:,8)))...
    -(TT(:,3)-TT(:,9))./((TT(:,2)-TT(:,8)).*(TT(:,4)-TT(:,8)));
b = ((TT(:,3)-TT(:,5))-a.*(TT(:,2).^2- TT(:,4).^2))./(TT(:,2)- TT(:,4));
c = TT(:,3)-a.*TT(:,2).^2-b.*TT(:,2);
y_1 = a.*TT(:,6).^2+b.*TT(:,6)+c;
a = (TT(:,5)-TT(:,9))./((TT(:,4)-TT(:,8)).*(TT(:,8)-TT(:,10)))...
    -(TT(:,5)-TT(:,11))./((TT(:,4)-TT(:,10)).*(TT(:,8)-TT(:,10)));
b = ((TT(:,5)-TT(:,9))-a.*(TT(:,4).^2- TT(:,8).^2))./(TT(:,4)- TT(:,8));
c = TT(:,5)-a.*TT(:,4).^2-b.*TT(:,4);
y_2 = a.*TT(:,6).^2+b.*TT(:,6)+c;

profOut = axis(:,1);
profOut(TT(:,1), 2) = mean([y_1, y_2], 2);
end


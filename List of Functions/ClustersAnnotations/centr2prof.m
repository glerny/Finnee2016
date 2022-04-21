function XYout = centr2prof(XYin, Rp, axis)
XYout = zeros(size(axis));
for lp1 = 1: size(XYin, 1)
    XYout = XYout + XYin(lp1, 2)*...
        pdf('normal', axis, XYin(lp1, 1), (XYin(lp1, 1)/Rp)/2.354);
end
XYout = [axis; XYout];
end
    
function M = ChrMoment3D(X, Y, XY)


for ii = 1:length(X)
    M0(ii) = trapz(Y, XY(:, ii));
    M1(ii) = trapz(Y, Y.* XY(:, ii))/M0(ii);
    M2(ii) = trapz(Y, (Y - M1(ii)).^2.* XY(:, ii))/M0(ii);
end

M = ChrMoment([X, M0']);
M(5) = sum(M0.*M1, 'omitnan')/sum(M0, 'omitnan');
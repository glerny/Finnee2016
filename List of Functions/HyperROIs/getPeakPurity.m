function [Crit1, Crit2, DeconvSignal] = getPeakPurity(Model, Data)
Areas = Data'/Model';

for ii = 1:size(Data, 2)
    Filter = Model.*Areas(ii, :);
    Filter = Filter./sum(Filter, 2);
    Filter(isnan(Filter)) = 0;
    DeconvSignal(:, :, ii) = Data(:,ii).*Filter;
end
SDS = squeeze(mean(DeconvSignal, 3, 'omitnan'));

warning off
for ii = 1:size(DeconvSignal, 2)
    Crit1(ii) = pdist2(Model(:, ii)', SDS(:, ii)', 'correlation');
    Crit2(:, :, ii) = ...
        [sum(squeeze(DeconvSignal(:, ii, :)))' ...
         pdist2(SDS(:, ii)', squeeze(DeconvSignal(:, ii, :))', 'correlation')'];
end
warning on
end


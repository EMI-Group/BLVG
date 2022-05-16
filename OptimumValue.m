function OptVal = OptimumValue
global MPB  Benchmark
fit = NaN(1,Benchmark.MPBnumber);
for ii=1 : Benchmark.MPBnumber
    tmp = NaN(1,MPB{ii}.PeakNumber);
    for jj=1 : MPB{ii}.PeakNumber
        tmp(jj) = Basefunction(MPB{ii}.PeaksPosition(jj,:)',ii);
    end
    fit(ii) = max(tmp);
end
OptVal = fit;  
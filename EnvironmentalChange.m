function EnvironmentalChange
global TotalResults TotalSumResults SumResults EnvironmentNumber test MPB ChangeFlag ContexVector evals Results ChangeFrequency Benchmark CurrentEnvironmentNumber
rng(evals * CurrentEnvironmentNumber);
CurrentEnvironmentNumber = CurrentEnvironmentNumber + 1;
ContexVectorUpdater;
ContexVectorFitness = ComputeContextVectorFitness(ContexVector, Benchmark.MPBnumber);
optimumValue = OptimumValue;
for i = 1 : Benchmark.MPBnumber
    Results{i}(1,round(evals/ChangeFrequency)) = ContexVectorFitness(:, i);
    Results{i}(2,round(evals/ChangeFrequency)) = optimumValue(:, i);
end
for i = 1 : Benchmark.MPBnumber
    TotalResults{i}(1, EnvironmentNumber * (test-1) + round(evals/ChangeFrequency)) = ContexVectorFitness(:, i);
    TotalResults{i}(2, EnvironmentNumber * (test-1) + round(evals/ChangeFrequency)) = optimumValue(:, i);
end
SumResults(1,round(evals/ChangeFrequency)) = sum(ContexVectorFitness);
SumResults(2,round(evals/ChangeFrequency)) = sum(optimumValue);
TotalSumResults(1, EnvironmentNumber * (test-1) + round(evals/ChangeFrequency)) = sum(ContexVectorFitness);
TotalSumResults(2, EnvironmentNumber * (test-1) + round(evals/ChangeFrequency)) = sum(optimumValue);
for ii=1 : Benchmark.MPBnumber
    for jj = 1: MPB{ii}.PeakNumber
        R = rand(1,MPB{ii}.Dimension)-0.5;
        Shift = MPB{ii}.ShiftSeverity(jj)  * (R/pdist2(R,zeros(size(R))));
        for kk = 1 : MPB{ii}.Dimension
            if ((MPB{ii}.PeaksPosition(jj,kk) + Shift(kk)) < MPB{ii}.MinCoordinate)
                MPB{ii}.PeaksPosition(jj,kk) = (2*MPB{ii}.MinCoordinate) - MPB{ii}.PeaksPosition(jj,kk) - Shift(kk);
            elseif ((MPB{ii}.PeaksPosition(jj,kk) + Shift(kk)) > MPB{ii}.MaxCoordinate)
                MPB{ii}.PeaksPosition(jj,kk) = (2*MPB{ii}.MaxCoordinate) - MPB{ii}.PeaksPosition(jj,kk)	- Shift(kk);
            else
                MPB{ii}.PeaksPosition(jj,kk) = MPB{ii}.PeaksPosition(jj,kk) + Shift(kk);
            end
        end
        %/* change MPB{ii}.PeaksPosition width */
        offset = randn * MPB{ii}.WidthSeverity(jj);
        if ((MPB{ii}.PeaksWidth(jj) + offset) < MPB{ii}.MinWidth)
            MPB{ii}.PeaksWidth(jj) = 2.0 * MPB{ii}.MinWidth - MPB{ii}.PeaksWidth(jj) - offset;
        elseif ((MPB{ii}.PeaksWidth(jj) + offset) > MPB{ii}.MaxWidth)
            MPB{ii}.PeaksWidth(jj) = 2.0 * MPB{ii}.MaxWidth -MPB{ii}.PeaksWidth(jj) - offset;
        else
            MPB{ii}.PeaksWidth(jj) = MPB{ii}.PeaksWidth(jj) + offset;
        end
        %/* change MPB{ii}.PeaksPosition height */
        offset = MPB{ii}.HeightSeverity(jj) * randn;
        if ((MPB{ii}.PeaksHeight(jj) + offset) < MPB{ii}.MinHeight)
            MPB{ii}.PeaksHeight(jj) = 2.0 * MPB{ii}.MinHeight - MPB{ii}.PeaksHeight(jj) - offset;
        elseif ((MPB{ii}.PeaksHeight(jj) + offset) > MPB{ii}.MaxHeight)
            MPB{ii}.PeaksHeight(jj) = 2.0 * MPB{ii}.MaxHeight - MPB{ii}.PeaksHeight(jj) - offset;
        else
            MPB{ii}.PeaksHeight(jj) = MPB{ii}.PeaksHeight(jj) + offset;
        end
    end
end
ChangeFlag=1;
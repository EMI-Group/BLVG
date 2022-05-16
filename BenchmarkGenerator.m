global Benchmark MPB ScenarioNumber;
%% input benchmark setting
Benchmark = [];
if ScenarioNumber==1
    Benchmark.NonSeparableSubFunctionDimensions = [2,2,2,3,3,3];
    Benchmark.SeparableSubFunctionDimension = [10];
elseif ScenarioNumber==2
    Benchmark.NonSeparableSubFunctionDimensions = [];
    Benchmark.SeparableSubFunctionDimension = [25];
elseif ScenarioNumber==3
    Benchmark.NonSeparableSubFunctionDimensions = [25];
    Benchmark.SeparableSubFunctionDimension = [];
elseif ScenarioNumber==4
    Benchmark.NonSeparableSubFunctionDimensions = [2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4];
    Benchmark.SeparableSubFunctionDimension = [20];
elseif ScenarioNumber==5
    Benchmark.NonSeparableSubFunctionDimensions = [];
    Benchmark.SeparableSubFunctionDimension = [50];
elseif ScenarioNumber==6
    Benchmark.NonSeparableSubFunctionDimensions = [50];
    Benchmark.SeparableSubFunctionDimension = [];
elseif ScenarioNumber==7
    Benchmark.NonSeparableSubFunctionDimensions = [2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4];
    Benchmark.SeparableSubFunctionDimension = [40];
elseif ScenarioNumber==8
    Benchmark.NonSeparableSubFunctionDimensions = [2, 2, 2, 3, 3, 3, 5, 10, 15, 25, 30];
    Benchmark.SeparableSubFunctionDimension = [];
elseif ScenarioNumber==9
    Benchmark.NonSeparableSubFunctionDimensions = [2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 30, 40];
    Benchmark.SeparableSubFunctionDimension = [5];
elseif ScenarioNumber==10
    Benchmark.NonSeparableSubFunctionDimensions = [];
    Benchmark.SeparableSubFunctionDimension = [100];
elseif ScenarioNumber==11
    Benchmark.NonSeparableSubFunctionDimensions = [100];
    Benchmark.SeparableSubFunctionDimension = [];
elseif ScenarioNumber==12
    Benchmark.NonSeparableSubFunctionDimensions = [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4];
    Benchmark.SeparableSubFunctionDimension = [120];
elseif ScenarioNumber==13
    Benchmark.NonSeparableSubFunctionDimensions = [2,3,5,10,20,30];
    Benchmark.SeparableSubFunctionDimension = [130];
elseif ScenarioNumber==14
    Benchmark.NonSeparableSubFunctionDimensions = [2, 2, 2, 3, 5, 5, 5, 5, 5, 8, 8, 10, 10, 10, 20, 20, 30, 50];
    Benchmark.SeparableSubFunctionDimension = [];
elseif ScenarioNumber==15
    Benchmark.NonSeparableSubFunctionDimensions = [];
    Benchmark.SeparableSubFunctionDimension = [200];
elseif ScenarioNumber==16
    Benchmark.NonSeparableSubFunctionDimensions = [200];
    Benchmark.SeparableSubFunctionDimension = [];
elseif ScenarioNumber==17
    Benchmark.NonSeparableSubFunctionDimensions = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 8, 8, 8, 8, 8];
    Benchmark.SeparableSubFunctionDimension = [150];
elseif ScenarioNumber==18
    Benchmark.NonSeparableSubFunctionDimensions = [2, 2, 3, 5, 5, 6, 6, 8, 8, 10, 10, 15, 20, 20, 30, 30, 40, 60];
    Benchmark.SeparableSubFunctionDimension = [20];
elseif ScenarioNumber==19
    Benchmark.NonSeparableSubFunctionDimensions = [2, 2, 2, 2, 3, 3, 5, 5, 5, 5, 5, 5, 8, 8, 10, 10, 20, 20, 50, 60, 70];
    Benchmark.SeparableSubFunctionDimension = [];
elseif ScenarioNumber==20
    Benchmark.NonSeparableSubFunctionDimensions = [];
    Benchmark.SeparableSubFunctionDimension = [300];
elseif ScenarioNumber==21
    Benchmark.NonSeparableSubFunctionDimensions = [300];
    Benchmark.SeparableSubFunctionDimension = [];
end
Benchmark.SeparableBaseFunction = 1;   %%%%%% usage
Benchmark.Dimension = sum(Benchmark.NonSeparableSubFunctionDimensions) + sum(Benchmark.SeparableSubFunctionDimension);
Benchmark.TotalSubFunctionNumber = length(Benchmark.NonSeparableSubFunctionDimensions) + length(Benchmark.SeparableSubFunctionDimension);
Benchmark.PermutationMap = randperm(Benchmark.Dimension);
Benchmark.UB = 50 * ones(1,Benchmark.Dimension);
Benchmark.LB = -50 * ones(1,Benchmark.Dimension);
Benchmark.NonSeparableSubFunctionShiftSeverity = ones(1,length(Benchmark.NonSeparableSubFunctionDimensions));
Benchmark.MPBnumber = length(Benchmark.NonSeparableSubFunctionDimensions) + sum(Benchmark.SeparableSubFunctionDimension);
MPB = cell(1,Benchmark.MPBnumber);
AverageNonSeparableDimension = mean(Benchmark.NonSeparableSubFunctionDimensions);
SeparableWeight=1;   %%% SeparableWeight is a regulatory factor
if ~isnan(Benchmark.SeparableSubFunctionDimension)
    if ~isnan(AverageNonSeparableDimension)
        if Benchmark.SeparableSubFunctionDimension > AverageNonSeparableDimension
            SeparableWeight = 1/AverageNonSeparableDimension;
        else 
            SeparableWeight = 1/Benchmark.SeparableSubFunctionDimension;
        end
    end
end
counter = 0;
if ~isempty(Benchmark.NonSeparableSubFunctionDimensions)
    for ii=1 : length(Benchmark.NonSeparableSubFunctionDimensions)
        counter = counter + 1;
        InitializingPeaks(Benchmark.NonSeparableSubFunctionDimensions(ii),counter,SeparableWeight);
    end
end
if ~isempty(Benchmark.SeparableSubFunctionDimension)
    for ii=1 : Benchmark.SeparableSubFunctionDimension
        counter = counter + 1;
        InitializingPeaks(1,counter,SeparableWeight);
    end
end
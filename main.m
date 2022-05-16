clear all;close all;clc; 
global TotalResults TotalSumResults EnvironmentNumber test Benchmark evals ContexVector PSO ChangeFlag Results SumResults ChangeFrequency CurrentEnvironmentNumber ScenarioNumber MaxEpochSize UpdateEnvironmentStepSize; 
RunNumber = 20;
FinalResults = cell(1,40);
FinalSumResults = NaN(2,RunNumber);
TotalResults = cell(1,40);
TotalSumResults = NaN(2,100*RunNumber);
for test=1 : RunNumber
    rng(test);
    %% Initialzing problems
    ScenarioNumber = 8;
    CurrentEnvironmentNumber = 1;
    evals = 0;
    ChangeFlag = 0;
    BenchmarkGenerator;
    EnvironmentNumber = 100;
    ChangeFrequency = Benchmark.Dimension * 500;
    MaxEvals = ChangeFrequency * EnvironmentNumber;
    Results = cell(1,Benchmark.MPBnumber);
    SumResults = NaN(2,EnvironmentNumber);
     [delta, lambda, evaluations] = ism();  %%% DG2
     [nonseps, seps, theta, epsilon] = dsm(evaluations, lambda, Benchmark.Dimension);  %%% DG2
     disp("DG2");
    %% Genaral Parameters setting
    MaxDecompositionNumber = 25;  % Dimensions >= MaxDecompositionNumber are decomposition component.
    MaxCombinedNumber = 10;       % The Maximal size of a combination number = MaxCombinedNumber.
    MinCombinedNumber = 5;        % The size of the last combination component > MinCombinedNumber. 
    MinIndependentNumber = 10;    % A nonseparable component's dimensions >= MinIndependentNumber && < MaxDecompositionNumber is an independent component
    %% Parameters of decomposition parts 
    EpochSizeMaxPercent = 0.5;  % EpochSize means exploration, this is 80% of total fitness evaluations of decomposition part
    EpochSizeMinPercent = 0.5;
    MaxEpochSize = 50;     % get by experiments
    UpdateEnvironmentStepSize = 20;   % Update EpochSize every 10 environmental change
    DimensionK1 = 10;
    %% Preliminary work: problem division and combination
    NonsepsNumber = length(nonseps);
    SepNumber = length(seps);
    TotalProblemNum = NonsepsNumber + SepNumber;  % When DG2 is inaccurate, CombinedProblemNum is not equal to Benchmark.MPBnumber
    DecompositionParts = cell(1, NonsepsNumber);
    count1 = 1;
    IndependentParts = cell(1, NonsepsNumber);
    count2 = 1;
    CombinationParts = cell(1, TotalProblemNum);
    count3 = 1;
    for i = 1 : NonsepsNumber
        [~, dimension] = size(nonseps{:,i});
        if dimension >= MaxDecompositionNumber                                          % Decomposition parts
            DecompositionParts(1, count1) = nonseps(:,i);
            count1 = count1 + 1;
        elseif dimension < MaxDecompositionNumber && dimension >= MinIndependentNumber  % Independent parts
            IndependentParts(1, count2) = nonseps(:,i);
            count2 = count2 + 1;
        elseif dimension < MinIndependentNumber                                         % Combination parts
            CombinationParts(1, count3) = nonseps(:,i);
            count3 = count3 + 1;
        end
    end
    if count1 == 1
        DecompositionParts = [];
    else
        DecompositionParts(cellfun(@isempty,DecompositionParts)) = [];
    end
    if count2 == 1
        IndependentParts = [];
    else
        IndependentParts(cellfun(@isempty,IndependentParts)) = [];
    end
    [~, dimension] = size(seps);
    for i = 1 : dimension
        CombinationParts(1, count3) = num2cell(seps(:,i));
        count3 = count3 + 1;
    end
    if count3 == 1
        CombinationParts = [];
    else
        CombinationParts(cellfun(@isempty,CombinationParts)) = [];
    end
    CombinationParts1 = [];
    if ~isempty(CombinationParts)
        [~, TotalCombinationPartsNumber] = size(CombinationParts);
        % sort members in CombinationParts according to Benchmark
        % Permutation Map£º When DG2 is not accurate...
        tmp = cell(size(CombinationParts));
        count = 1;
        ldim = 1;
        for i = 1 : size(CombinationParts, 2)
            C(:, i).CombinationParts = CombinationParts{i};
            C(:, i).index = -1;
        end
        for i = 1 : Benchmark.MPBnumber 
            for j = 1 : TotalCombinationPartsNumber
%                 disp(C(:, j).CombinationParts);
                if isempty(setdiff(Benchmark.PermutationMap(ldim:ldim+MPB{i}.Dimension-1), C(:, j).CombinationParts))
                    if size(Benchmark.PermutationMap(ldim:ldim+MPB{i}.Dimension-1),2) == size(C(:, j).CombinationParts,2)
                        if(C(:, j).index == -1)
                            tmp{count} = Benchmark.PermutationMap(ldim:ldim+MPB{i}.Dimension-1);
                            count = count + 1;
                            break;
                        else
                            tmp{C(:, j).index} = [tmp{C(:, j).index}, Benchmark.PermutationMap(ldim:ldim+MPB{i}.Dimension-1)];
                            break;
                        end
                    else
                        if (C(:, j).index == -1)
                            tmp{count} = [tmp{count}, Benchmark.PermutationMap(ldim:ldim+MPB{i}.Dimension-1)];
                            C(:, j).index = count;
                            count = count + 1;
                        else
                            tmp{C(:, j).index} = [tmp{C(:, j).index}, Benchmark.PermutationMap(ldim:ldim+MPB{i}.Dimension-1)];
                        end
                        C(:, j).CombinationParts = setdiff(C(:, j).CombinationParts, Benchmark.PermutationMap(ldim:ldim+MPB{i}.Dimension-1));
                    end
                end
            end
            ldim = ldim + MPB{i}.Dimension;
        end
        CombinationParts = tmp;
        % Combine Combination parts
        CumulatedDimension = 0;
        CombinedProblemNum = 1;
        count = 1;
        CombinationParts1{count}.Variable = [];
        for i = 1 : TotalCombinationPartsNumber
            [~, dimension] = size(CombinationParts{:,i});
            if CumulatedDimension + dimension <= MaxCombinedNumber
                CombinationParts1{count}.Variable = [CombinationParts1{count}.Variable, CombinationParts(:,i)];
                CombinationParts1{count}.CombinedProblemNum = CombinedProblemNum;
                CumulatedDimension = CumulatedDimension + dimension;
                CombinedProblemNum = CombinedProblemNum + 1;
            else
                count = count + 1;
                CombinationParts1{count}.Variable = [];
                CumulatedDimension = 0;
                CombinedProblemNum = 1;
                CombinationParts1{count}.Variable = [CombinationParts1{count}.Variable, CombinationParts(:,i)];
                CombinationParts1{count}.CombinedProblemNum = CombinedProblemNum;
                CumulatedDimension = CumulatedDimension + dimension;
                CombinedProblemNum = CombinedProblemNum + 1;
            end
        end
        lastCombinationNumber = 0;
        for i = 1 : CombinationParts1{count}.CombinedProblemNum
            lastCombinationNumber = lastCombinationNumber + length(CombinationParts1{count}.Variable{i});
        end
        if lastCombinationNumber <= MinCombinedNumber
            CombinationParts1{count-1}.Variable = [CombinationParts1{count-1}.Variable, CombinationParts1{count}.Variable];
            CombinationParts1{count-1}.CombinedProblemNum = CombinationParts1{count-1}.CombinedProblemNum + CombinationParts1{count}.CombinedProblemNum;
            CombinationParts1{count} = [];
        end 
        CombinationParts1(cellfun(@isempty,CombinationParts1)) = [];
    end
    %% Initialing PSOs 
    ContexVector = zeros(1,Benchmark.Dimension);
    [~, DecompositionPartsNumber] = size(DecompositionParts);
    [~, IndependentPartsNumber] = size(IndependentParts);
    [~, CombinationPartsNumber] = size(CombinationParts1);
    TotalSwarmNumber = IndependentPartsNumber + DecompositionPartsNumber + CombinationPartsNumber;
    PSO = cell(1, TotalSwarmNumber);
    for i = 1 : IndependentPartsNumber
        PSO{i}.ProblemType = 'independent parts';
        PSO{i}.Variable = IndependentParts{i};
    end
    count = 1;
    for i = IndependentPartsNumber+1 : IndependentPartsNumber + DecompositionPartsNumber
        PSO{i}.ProblemType = 'decomposition parts';
        PSO{i}.Variable = DecompositionParts{count};
        count = count + 1;
    end
    count = 1;
    for i = IndependentPartsNumber+DecompositionPartsNumber+1 : TotalSwarmNumber
        PSO{i}.ProblemType = 'combination parts';
        PSO{i}.Variable = CombinationParts1{count}.Variable;
        PSO{i}.CombinedProblemNum = CombinationParts1{count}.CombinedProblemNum;
        count = count + 1;
    end
    for i = 1 : TotalSwarmNumber
        if strcmp(PSO{i}.ProblemType, 'independent parts')
            PSO{i}.PopulationSize = 10;
            PSO{i}.Dimension = length(PSO{i}.Variable);
            PSO{i}.MaxFitnessEvaluations = PSO{i}.Dimension / Benchmark.Dimension * ChangeFrequency;
            PSO{i}.FitnessEvaluations = 0;
            PSO{i}.MaxCoordinate = Benchmark.UB(1);
            PSO{i}.MinCoordinate = Benchmark.LB(1);
            PSO{i}.SwarmNumber = 1;
            PSO{i}.FreeSwarmID = 1;
            PSO{i}.MinExclusionLimit = 0.5 * ((PSO{i}.MaxCoordinate-PSO{i}.MinCoordinate) / ((10) ^ (1 / PSO{i}.Dimension)));
            PSO{i}.ExclusionLimit = PSO{i}.MinExclusionLimit;
            PSO{i}.ConvergenceLimit = PSO{i}.ExclusionLimit;
            PSO{i}.ShiftSeverity = 1;
            PSO{i}.x = 0.729843788;
            PSO{i}.c1 = 2.05;
            PSO{i}.c2 = 2.05;
            PSO{i}.DeactivationLimit = 0.1;
            [PSO{i}.particle(PSO{i}.SwarmNumber), PSO{i}.FitnessEvaluations] = InitializingPSO(PSO{i});
        elseif strcmp(PSO{i}.ProblemType, 'decomposition parts')
            PSO{i}.Dimension = length(PSO{i}.Variable);
            PSO{i}.MaxFitnessEvaluations = PSO{i}.Dimension / Benchmark.Dimension * ChangeFrequency;
            PSO{i}.FitnessEvaluations = 0;
            PSO{i}.MaxCoordinate = Benchmark.UB(1);
            PSO{i}.MinCoordinate = Benchmark.LB(1);
            PSO{i}.DimensionK1 = DimensionK1;
            PSO{i}.K1 = floor(PSO{i}.Dimension / PSO{i}.DimensionK1);
            if rem(PSO{i}.Dimension, PSO{i}.DimensionK1) ~= 0
                PSO{i}.K2 = 1;
            else
                PSO{i}.K2 = 0;
            end
            PSO{i}.K = PSO{i}.K1 + PSO{i}.K2;
            PSO{i}.DimensionK2 = PSO{i}.Dimension - PSO{i}.DimensionK1 * PSO{i}.K1;
            PSO{i}.DecompContext.DecompContextMapIndex = 1 : PSO{i}.Dimension;
            PSO{i}.DecompContext.DecompContextVector = PSO{i}.MinCoordinate + (PSO{i}.MaxCoordinate - PSO{i}.MinCoordinate).*rand(1, PSO{i}.Dimension);
            DecompContextVector_ = NaN(1, size(ContexVector, 2));
            DecompContextVector_(:, PSO{i}.Variable) =  PSO{i}.DecompContext.DecompContextVector;
            [PSO{i}.DecompContext.DecompContextVectorFitness, PSO{i}.FitnessEvaluations] = benchmark_func_decomp(DecompContextVector_, PSO{i}.FitnessEvaluations);  
            PSO{i}.VariableMapIndex = randperm(PSO{i}.Dimension);
            PSO{i}.SearchSwarms = cell(1, PSO{i}.K);
            PSO{i}.TrackerSwarms.PopulationSize = 10;
            PSO{i}.TrackerSwarms.Dimension = PSO{i}.Dimension;
            PSO{i}.TrackerSwarms.x = 0.729843788;
            PSO{i}.TrackerSwarms.c1 = 2.05;
            PSO{i}.TrackerSwarms.c2 = 2.05; 
            PSO{i}.TrackerSwarms.MaxCoordinate = PSO{i}.MaxCoordinate;
            PSO{i}.TrackerSwarms.MinCoordinate = PSO{i}.MinCoordinate;
            PSO{i}.TrackerSwarms.ShiftSeverity = 1;
            PSO{i}.TrackerSwarms.MinExclusionLimit = 0.5 * ((PSO{i}.TrackerSwarms.MaxCoordinate-PSO{i}.TrackerSwarms.MinCoordinate) / ((10) ^ (1 / PSO{i}.TrackerSwarms.Dimension)));
            PSO{i}.TrackerSwarms.ExclusionLimit = PSO{i}.TrackerSwarms.MinExclusionLimit;
            PSO{i}.TrackerSwarms.SwarmNumber = 0;
            % Initializing K component SearchParticles
            beginIndex = 1;
            endIndex = beginIndex;
            for j = 1 : PSO{i}.K
                if j <= PSO{i}.K1
                    PSO{i}.SearchSwarms{j}.Dimension = PSO{i}.DimensionK1;
                else
                    PSO{i}.SearchSwarms{j}.Dimension = PSO{i}.DimensionK2;
                end
                PSO{i}.SearchSwarms{j}.PopulationSize = 10;
                PSO{i}.SearchSwarms{j}.MaxCoordinate = PSO{i}.MaxCoordinate;
                PSO{i}.SearchSwarms{j}.MinCoordinate = PSO{i}.MinCoordinate;
                PSO{i}.SearchSwarms{j}.x = 0.729843788;
                PSO{i}.SearchSwarms{j}.c1 = 2.05;
                PSO{i}.SearchSwarms{j}.c2 = 2.05;
                endIndex = beginIndex + PSO{i}.SearchSwarms{j}.Dimension - 1;   % Begin and end index of each component in complete decision.
                PSO{i}.SearchSwarms{j}.BeginIndex = beginIndex;
                PSO{i}.SearchSwarms{j}.EndIndex = endIndex;
                beginIndex = endIndex + 1;
                PSO{i}.SearchSwarms{j}.MapIndex = PSO{i}.VariableMapIndex(:, PSO{i}.SearchSwarms{j}.BeginIndex: PSO{i}.SearchSwarms{j}.EndIndex);
                [PSO{i}.SearchSwarms{j}.particle, PSO{i}.FitnessEvaluations] = InitializingPSO(PSO{i}, j);
                PSO{i}.EpochSizePercent = EpochSizeMaxPercent;
                PSO{i}.EpochSize = (PSO{i}.MaxFitnessEvaluations * PSO{i}.EpochSizePercent) / (PSO{i}.SearchSwarms{1}.PopulationSize * PSO{i}.K);
                if PSO{i}.EpochSize > MaxEpochSize
                    PSO{i}.EpochSize = MaxEpochSize;  % EpochSize always >= 50, because the EpochSizeMinPercent = 0.2
                end
                PSO{i}.EpochSizePercentStep = (EpochSizeMaxPercent - EpochSizeMinPercent) / floor(EnvironmentNumber / UpdateEnvironmentStepSize);
                PSO{i}.EpochSize2Percent = 1 - PSO{i}.EpochSizePercent;
                PSO{i}.EpochSize2 = (PSO{i}.MaxFitnessEvaluations * PSO{i}.EpochSize2Percent) / (PSO{i}.SearchSwarms{1}.PopulationSize * PSO{i}.TrackerSwarms.SwarmNumber);
            end
        elseif strcmp(PSO{i}.ProblemType, 'combination parts')
            PSO{i}.PopulationSize = 10;
            PSO{i}.Dimension = 0;
            for j = 1 : PSO{i}.CombinedProblemNum
                PSO{i}.Dimension = PSO{i}.Dimension + length(PSO{i}.Variable{j});
            end
            PSO{i}.MaxFitnessEvaluations = PSO{i}.Dimension / Benchmark.Dimension * ChangeFrequency;
            PSO{i}.FitnessEvaluations = 0;
            PSO{i}.MaxCoordinate = Benchmark.UB(1);
            PSO{i}.MinCoordinate = Benchmark.LB(1);
            for j = 1 : PSO{i}.CombinedProblemNum
                PSO{i}.MinExclusionLimit(j) = 0.5 * ((PSO{i}.MaxCoordinate-PSO{i}.MinCoordinate) / ((10) ^ (1 / length(PSO{i}.Variable{j}))));
                PSO{i}.ExclusionLimit(j) = PSO{i}.MinExclusionLimit(j);
                PSO{i}.ConvergenceLimit(j) = PSO{i}.ExclusionLimit(j);
                PSO{i}.ShiftSeverity(j) = 1;  
            end
            PSO{i}.x = 0.729843788;
            PSO{i}.c1 = 2.05;
            PSO{i}.c2 = 2.05;
            PSO{i}.SwarmNumber = 1;
            PSO{i}.FreeSwarmID = 1;
            PSO{i}.DeactivationLimit = 0.1;
            [PSO{i}.particle(PSO{i}.SwarmNumber), PSO{i}.FitnessEvaluations] = InitializingPSO(PSO{i});
        end
    end
    %% Main loop
    while 1
        mmPSO;
        if ChangeFlag == 1
            ChangeFlag = 0;  
            Reaction; 
            disp([test,CurrentEnvironmentNumber]);
        end
        if  evals > MaxEvals
            break;
        end
    end
    for i = 1 : TotalProblemNum
        FinalResults{i}(1,test) = mean(Results{i}(1,:));
        FinalResults{i}(2,test) = mean(Results{i}(2,:));
    end
    FinalSumResults(1, test) = mean(SumResults(1,:));
    FinalSumResults(2, test) = mean(SumResults(2,:));
end
OutputSA = FinalSumResults(2,:)-FinalSumResults(1,:);
x="a2_f20";
tmp1 = fopen(strcat(x,'_SA.txt'),'wt');
for ii = 1:RunNumber
    fprintf(tmp1,'%g',OutputSA(ii));
    fprintf(tmp1,'\n');
end

fclose(tmp1);
Output(1) = median(OutputSA);
Output(2) = mean(OutputSA);
Output(3) = std(OutputSA)/sqrt(RunNumber);
tmp2 = fopen(strcat(x,'.txt'),'wt');
for ii = 1:3
    fprintf(tmp2,'%g',Output(ii));
    fprintf(tmp2,'\n');
end
fclose(tmp2);
function PSO = Reaction_decomposition(PSO)
global ContexVector MaxEpochSize ChangeFlag evals CurrentEnvironmentNumber UpdateEnvironmentStepSize
tmpEvals = evals;
PSO.FitnessEvaluations = 0;
%% Context vector become a new gbest && Exclusion 
PSO.TrackerSwarms.SwarmNumber = PSO.TrackerSwarms.SwarmNumber + 1;
PSO.TrackerSwarms.particle(PSO.TrackerSwarms.SwarmNumber).GbestPosition = PSO.DecompContext.DecompContextVector;
PSO.TrackerSwarms.particle(PSO.TrackerSwarms.SwarmNumber).GbestValue = PSO.DecompContext.DecompContextVectorFitness;
PSO.TrackerSwarms.particle(PSO.TrackerSwarms.SwarmNumber).Gbest_past_environment = PSO.TrackerSwarms.particle(PSO.TrackerSwarms.SwarmNumber).GbestPosition;
PSO.TrackerSwarms.particle(PSO.TrackerSwarms.SwarmNumber).Shifts = [];
PSO.TrackerSwarms.particle(PSO.TrackerSwarms.SwarmNumber).Velocity = zeros(PSO.TrackerSwarms.PopulationSize,PSO.TrackerSwarms.Dimension);
ExclusionFlag = 0;
for ii = 1 : PSO.TrackerSwarms.SwarmNumber
    for jj = 1 : PSO.TrackerSwarms.SwarmNumber
        if ii ~= jj && ExclusionFlag==0 && pdist2(PSO.TrackerSwarms.particle(ii).GbestPosition,PSO.TrackerSwarms.particle(jj).GbestPosition) < PSO.TrackerSwarms.ExclusionLimit
            if PSO.TrackerSwarms.particle(ii).GbestValue < PSO.TrackerSwarms.particle(jj).GbestValue
                PSO.TrackerSwarms.particle(ii) = [];
                PSO.TrackerSwarms.SwarmNumber = PSO.TrackerSwarms.SwarmNumber-1;
            else
                PSO.TrackerSwarms.particle(jj) = [];
                PSO.TrackerSwarms.SwarmNumber = PSO.TrackerSwarms.SwarmNumber-1;
            end
            ExclusionFlag = 1;
        end
    end
end
%% Introduce diversity around GbestPosition
for i = 1 : PSO.TrackerSwarms.SwarmNumber
    PSO.TrackerSwarms.particle(i).X = repmat(PSO.TrackerSwarms.particle(i).GbestPosition,PSO.TrackerSwarms.PopulationSize,1)+ (rands(PSO.TrackerSwarms.PopulationSize,PSO.TrackerSwarms.Dimension)*PSO.TrackerSwarms.ShiftSeverity);
    PSO.TrackerSwarms.particle(i).X(1,:) = PSO.TrackerSwarms.particle(i).GbestPosition; 
end
%% Update memory for PSO.TrackerSwarms
for i=1 : PSO.TrackerSwarms.SwarmNumber
    tmp = NaN(size(ContexVector));
    Y = repmat(tmp,PSO.TrackerSwarms.PopulationSize,1);
    Y(:, PSO.Variable) = PSO.TrackerSwarms.particle(i).X;
    [PSO.TrackerSwarms.particle(i).FitnessValue, PSO.FitnessEvaluations] = benchmark_func_decomp(Y, PSO.FitnessEvaluations); 
    PSO.TrackerSwarms.particle(i).PbestValue = PSO.TrackerSwarms.particle(i).FitnessValue;
    PSO.TrackerSwarms.particle(i).PbestPosition = PSO.TrackerSwarms.particle(i).X;
    PSO.TrackerSwarms.particle(i).Gbest_past_environment = PSO.TrackerSwarms.particle(i).GbestPosition;
    [PSO.TrackerSwarms.particle(i).GbestValue,BestPbestID] = max(PSO.TrackerSwarms.particle(i).PbestValue);
    PSO.TrackerSwarms.particle(i).GbestPosition = PSO.TrackerSwarms.particle(i).PbestPosition(BestPbestID,:);
end
%% Update memory for PSO.SearchSwarm
for i = 1 : PSO.K
    [PSO.SearchSwarms{i}.particle.FitnessValue, PSO.FitnessEvaluations] = Fitness(PSO.SearchSwarms{i}.particle.X, PSO.Variable, PSO.ProblemType, PSO.SearchSwarms{i}.MapIndex, PSO.DecompContext, PSO.FitnessEvaluations);
    PSO.SearchSwarms{i}.particle.PbestValue = PSO.SearchSwarms{i}.particle.FitnessValue;
    PSO.SearchSwarms{i}.particle.PbestPosition = PSO.SearchSwarms{i}.particle.X;
    [PSO.SearchSwarms{i}.particle.GbestValue, BestIndex] = max(PSO.SearchSwarms{i}.particle.PbestValue);
    PSO.SearchSwarms{i}.particle.GbestPosition = PSO.SearchSwarms{i}.particle.X(BestIndex,:);
end
%% Acclerate velocity for particle stagnation
for i = 1 : PSO.K
    [PopulationSize, Dimension]  = size(PSO.SearchSwarms{i}.particle.X);
    for ii=1 : PopulationSize
        for jj=1 : Dimension
            if PSO.SearchSwarms{i}.particle.Velocity(ii,jj) < 1 && PSO.SearchSwarms{i}.particle.Velocity(ii,jj) > -1
                PSO.SearchSwarms{i}.particle.Velocity(ii,jj) = randn;
            end
        end
    end
end
%% Reinitialize Context Vector
PSO.DecompContext.DecompContextVector = PSO.MinCoordinate + (PSO.MaxCoordinate - PSO.MinCoordinate).*rand(1, PSO.Dimension);
DecompContextVector_ = NaN(1, size(ContexVector, 2));
DecompContextVector_(:, PSO.Variable) =  PSO.DecompContext.DecompContextVector;
[PSO.DecompContext.DecompContextVectorFitness, PSO.FitnessEvaluations] = benchmark_func_decomp(DecompContextVector_, PSO.FitnessEvaluations);  
%% Update EpochSize & EpochSize2
if rem(CurrentEnvironmentNumber , UpdateEnvironmentStepSize) == 0
    PSO.EpochSizePercent = PSO.EpochSizePercent - PSO.EpochSizePercentStep;
    PSO.EpochSize = (PSO.MaxFitnessEvaluations * PSO.EpochSizePercent) / (PSO.SearchSwarms{1}.PopulationSize * PSO.K);
    if PSO.EpochSize > MaxEpochSize
        PSO.EpochSize = MaxEpochSize;  
    end
    PSO.EpochSize2Percent = 1 - PSO.EpochSizePercent;
    PSO.EpochSize2 = (PSO.MaxFitnessEvaluations * PSO.EpochSize2Percent) / (PSO.SearchSwarms{1}.PopulationSize * PSO.TrackerSwarms.SwarmNumber);
else
    PSO.EpochSize2 = (PSO.MaxFitnessEvaluations * PSO.EpochSize2Percent) / (PSO.SearchSwarms{1}.PopulationSize * PSO.TrackerSwarms.SwarmNumber);
end
%% Move mPSO for EpochSize2 times
if PSO.TrackerSwarms.SwarmNumber ~= 0
    for i = 1 : PSO.EpochSize2
        for j = 1 : PSO.TrackerSwarms.SwarmNumber
            PSO.TrackerSwarms.particle(j).Velocity = PSO.TrackerSwarms.x * (PSO.TrackerSwarms.particle(j).Velocity ...
                + (PSO.TrackerSwarms.c1 * rand(PSO.TrackerSwarms.PopulationSize, PSO.TrackerSwarms.Dimension).*(PSO.TrackerSwarms.particle(j).PbestPosition - PSO.TrackerSwarms.particle(j).X))...
                + (PSO.TrackerSwarms.c2 * rand(PSO.TrackerSwarms.PopulationSize, PSO.TrackerSwarms.Dimension).*(repmat(PSO.TrackerSwarms.particle(j).GbestPosition,PSO.TrackerSwarms.PopulationSize,1) - PSO.TrackerSwarms.particle(j).X)));
            PSO.TrackerSwarms.particle(j).X = PSO.TrackerSwarms.particle(j).X + PSO.TrackerSwarms.particle(j).Velocity;
            for jj=1 : PSO.TrackerSwarms.PopulationSize
                for kk=1 : PSO.TrackerSwarms.Dimension
                    if PSO.TrackerSwarms.particle(j).X(jj,kk) > PSO.TrackerSwarms.MaxCoordinate
                        PSO.TrackerSwarms.particle(j).X(jj,kk) = PSO.TrackerSwarms.MaxCoordinate;
                        PSO.TrackerSwarms.particle(j).Velocity(jj,kk) = 0;
                    elseif PSO.TrackerSwarms.particle(j).X(jj,kk) < PSO.TrackerSwarms.MinCoordinate
                        PSO.TrackerSwarms.particle(j).X(jj,kk) = PSO.TrackerSwarms.MinCoordinate;
                        PSO.TrackerSwarms.particle(j).Velocity(jj,kk) = 0;
                    end
                end
            end
            tmp = NaN(size(ContexVector));
            Y = repmat(tmp,PSO.TrackerSwarms.PopulationSize,1);
            Y(:, PSO.Variable) = PSO.TrackerSwarms.particle(j).X;
            [PSO.TrackerSwarms.particle(j).FitnessValue, PSO.FitnessEvaluations] = benchmark_func_decomp(Y, PSO.FitnessEvaluations); 
            if ChangeFlag ~= 0
                return;
            end
            tmp = PSO.TrackerSwarms.particle(j).FitnessValue > PSO.TrackerSwarms.particle(j).PbestValue ;
            PSO.TrackerSwarms.particle(j).PbestValue(tmp)  = PSO.TrackerSwarms.particle(j).FitnessValue(tmp);
            PSO.TrackerSwarms.particle(j).PbestPosition(tmp, :) = PSO.TrackerSwarms.particle(j).X(tmp, :);
            [PSO.TrackerSwarms.particle(j).GbestValue, BestIndex] = max(PSO.TrackerSwarms.particle(j).PbestValue);
            PSO.TrackerSwarms.particle(j).GbestPosition = PSO.TrackerSwarms.particle(j).PbestPosition(BestIndex,:);
        end
    end
end
%% Disrupt the order of subproblem
PSO.VariableMapIndex = randperm(PSO.Dimension);
for i = 1 : PSO.K
    if ChangeFlag == 0
        PSO.SearchSwarms{i}.MapIndex = PSO.VariableMapIndex(:, PSO.SearchSwarms{i}.BeginIndex: PSO.SearchSwarms{i}.EndIndex);
        [PSO.SearchSwarms{i}.particle, PSO.FitnessEvaluations] = InitializingPSO(PSO, i);
    end
    if ChangeFlag ~= 0
        return;
    end
end
% disp("evals: ");
% disp(evals - tmpEvals);
end
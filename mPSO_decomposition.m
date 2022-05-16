function Pso = mPSO_decomposition(Pso)
global ChangeFlag ContexVector;
%% Each Pso component run EpochSize times
for i = 1 : Pso.K  
    for j = 1 : Pso.EpochSize
        if ChangeFlag == 0
            if Pso.FitnessEvaluations >= Pso.MaxFitnessEvaluations
                return;
            end
            Pso.SearchSwarms{i}.particle.Velocity =  Pso.SearchSwarms{i}.x * (Pso.SearchSwarms{i}.particle.Velocity + (Pso.SearchSwarms{i}.c1 * rand(Pso.SearchSwarms{i}.PopulationSize,Pso.SearchSwarms{i}.Dimension).*(Pso.SearchSwarms{i}.particle.PbestPosition - Pso.SearchSwarms{i}.particle.X)) + (Pso.SearchSwarms{i}.c2 * rand(Pso.SearchSwarms{i}.PopulationSize,Pso.SearchSwarms{i}.Dimension).*(repmat(Pso.SearchSwarms{i}.particle.GbestPosition,Pso.SearchSwarms{i}.PopulationSize,1)-Pso.SearchSwarms{i}.particle.X)));
            [Pso.SearchSwarms{i}.particle.X, Pso.SearchSwarms{i}.particle.Velocity] = reflectMethodForParticles(repmat(Pso.SearchSwarms{i}.MinCoordinate, Pso.SearchSwarms{i}.PopulationSize, Pso.SearchSwarms{i}.Dimension), repmat(Pso.SearchSwarms{i}.MaxCoordinate, Pso.SearchSwarms{i}.PopulationSize, Pso.SearchSwarms{i}.Dimension), Pso.SearchSwarms{i}.particle.X, Pso.SearchSwarms{i}.particle.Velocity); 
            [Pso.SearchSwarms{i}.particle.FitnessValue,Pso.FitnessEvaluations] = Fitness(Pso.SearchSwarms{i}.particle.X, Pso.Variable, Pso.ProblemType, Pso.SearchSwarms{i}.MapIndex, Pso.DecompContext, Pso.FitnessEvaluations);
            if ChangeFlag ~= 0
                return;
            end
            tmp = Pso.SearchSwarms{i}.particle.FitnessValue > Pso.SearchSwarms{i}.particle.PbestValue;
            Pso.SearchSwarms{i}.particle.PbestValue(tmp) = Pso.SearchSwarms{i}.particle.FitnessValue(tmp);
            Pso.SearchSwarms{i}.particle.PbestPosition(tmp, :) = Pso.SearchSwarms{i}.particle.X(tmp, :);
            [Pso.SearchSwarms{i}.particle.GbestValue, GbestPositionIndex] = max(Pso.SearchSwarms{i}.particle.PbestValue);
            Pso.SearchSwarms{i}.particle.GbestPosition = Pso.SearchSwarms{i}.particle.PbestPosition(GbestPositionIndex,:); 
        end
    end
    if Pso.SearchSwarms{i}.particle.GbestValue > Pso.DecompContext.DecompContextVectorFitness
        Pso.DecompContext.DecompContextVector(:, Pso.SearchSwarms{i}.MapIndex) = Pso.SearchSwarms{i}.particle.GbestPosition;
        DecompContextVector_ = NaN(1, size(ContexVector, 2));
        DecompContextVector_(:, Pso.Variable) =  Pso.DecompContext.DecompContextVector;
        [Pso.DecompContext.DecompContextVectorFitness, Pso.FitnessEvaluations] = benchmark_func_decomp(DecompContextVector_, Pso.FitnessEvaluations); 
        if ChangeFlag ~= 0
            return;
        end
    end
end
%% Disrupt the order of subproblem
Pso.VariableMapIndex = randperm(Pso.Dimension);
for i = 1 : Pso.K
    if ChangeFlag == 0
        Pso.SearchSwarms{i}.MapIndex = Pso.VariableMapIndex(:, Pso.SearchSwarms{i}.BeginIndex: Pso.SearchSwarms{i}.EndIndex);
        if Pso.FitnessEvaluations >= Pso.MaxFitnessEvaluations
            return;
        end
        [Pso.SearchSwarms{i}.particle, Pso.FitnessEvaluations] = InitializingPSO(Pso, i);
    end
    if ChangeFlag ~= 0
        return;
    end
end
end

%% Adjust particles overstep the boundary 
function [reflectResult, step]  = reflectMethodForParticles(min, max, oldValue, step)
    tmp = oldValue + step;
    reflectResult = tmp;
    stepTmp = zeros(size(reflectResult));
    minJugde = tmp < min;
    maxJudge = tmp > max;
    tmp1 = min;
    tmp2 = max;
    reflectResult(minJugde) = tmp1(minJugde);
    step(minJugde) = stepTmp(minJugde);
    reflectResult(maxJudge) = tmp2(maxJudge);
    step(maxJudge) = stepTmp(maxJudge);
end

function [particle, FitnessEvaluations] = InitializingPSO(varargin)
global Benchmark MPB;
narginchk(1,2);
if nargin == 1
    PSO = varargin{1};
elseif nargin == 2
    PSO = varargin{1};
    SearchSwarmsIndex = varargin{2};
end
if strcmp(PSO.ProblemType, 'decomposition parts')
    particle.X = PSO.SearchSwarms{SearchSwarmsIndex}.MinCoordinate + ((PSO.SearchSwarms{SearchSwarmsIndex}.MaxCoordinate - PSO.SearchSwarms{SearchSwarmsIndex}.MinCoordinate) .* rand(PSO.SearchSwarms{SearchSwarmsIndex}.PopulationSize, PSO.SearchSwarms{SearchSwarmsIndex}.Dimension));
    [particle.FitnessValue, FitnessEvaluations] = Fitness(particle.X, PSO.Variable, PSO.ProblemType, PSO.SearchSwarms{SearchSwarmsIndex}.MapIndex, PSO.DecompContext, PSO.FitnessEvaluations);
    particle.Velocity = zeros(PSO.SearchSwarms{SearchSwarmsIndex}.PopulationSize,PSO.SearchSwarms{SearchSwarmsIndex}.Dimension);
else
    particle.X = PSO.MinCoordinate + ((PSO.MaxCoordinate - PSO.MinCoordinate) .* rand(PSO.PopulationSize, PSO.Dimension));
    [particle.FitnessValue, FitnessEvaluations] = Fitness(particle.X, PSO.Variable, PSO.ProblemType, PSO.FitnessEvaluations);
    particle.Velocity = zeros(PSO.PopulationSize,PSO.Dimension);
end
particle.Shifts = [];
particle.ShiftSeverity = 1;
particle.Gbest_past_environment = NaN(1, PSO.Dimension);
particle.DeactivationFlag = 0;
particle.PbestPosition = particle.X;
particle.PbestValue = particle.FitnessValue;
particle.PastPbestValue = particle.PbestValue;
[particle.GbestValue, GbestID] = max(particle.PbestValue);
if strcmp(PSO.ProblemType, 'combination parts')
    ldim = 1;
    for i = 1 : PSO.CombinedProblemNum
        particle.GbestPosition(:, ldim:ldim + length(PSO.Variable{i}) - 1) = particle.PbestPosition(GbestID(:,i),ldim:ldim + length(PSO.Variable{i}) - 1);
        particle.ShiftSeverity(i) = 1;
        ldim = ldim + length(PSO.Variable{i});
    end
else
    particle.GbestPosition = particle.PbestPosition(GbestID, :);
end


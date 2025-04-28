classdef GSEA < ALGORITHM
% <multi> <real> <multimodal>
% A Grid Self-adaptive Exploration-based Algorithm for Multimodal Multiobjective Optimization

    methods
        function main(Algorithm,Problem)
            T_stage = Algorithm.ParameterSet(0.2);
            
            %% Generate random population
            Population = Problem.Initialization();
            
            %% Initiate the grid and table
            V = ceil((Problem.maxFE/Problem.D)^(1/Problem.D)).*ones(1,Problem.D);   % the number of segments along each dimension in decision space
            w = (Problem.upper - Problem.lower)./V;                 % the width of each segment
            State = zeros(V,'int8');
            SearchCount = zeros(V,'int8');
            
            %% update the table
            Archive = [];
            [Archive,State,SearchCount] = updateTable(State,SearchCount,Archive,Population,w,Problem.lower);
            [~,FrontNo,CrowdDisDec] = EnvironmentalSelection5(Population,Problem.N,V,w,Problem.lower,1);
      
            %% Optimization
            while Algorithm.NotTerminated(Population)
				MatingPool = TournamentSelection(2,Problem.N,CrowdDisDec,FrontNo);
                parents = Population(MatingPool);
                OffDec  = OperatorGA(Problem,parents.decs);
                
                OffLocation = getLocation(OffDec,V,w,Problem.lower);
				ColIndex = getIndex(OffLocation,V(1));
                Status = State(ColIndex);
                
                % State==-1, possible to reexplore
                if ~isempty(find(Status==-1,1))
                    temp = find(Status==-1);
                    Gen = Problem.FE/Problem.N;
                    Possibility = min(exp((double(SearchCount(ColIndex(temp)))-min(Problem.D,2))/Gen),1);
                    temp = temp(rand(length(temp),1)<=Possibility);
                    if ~isempty(temp)
                        % DE operator generating a new individual, x_i and
                        % x_r2 are solutions with good convergence quality
                        % in population, x_r1 is randomly selected from
                        % archive
                        FirstParent = Population(MatingPool(temp));
                        MatingPool = TournamentSelection(2,2*length(temp),FrontNo,CrowdDisDec);
						SecondParent = Population(MatingPool(1:length(temp)));
                        ThirdParent = Population(MatingPool(length(temp)+1:end));
                        OffDecDE = OperatorDE(Problem,FirstParent.decs,SecondParent.decs,ThirdParent.decs);
                        OffDec(temp,:) = OffDecDE;
                    end
                end
                
                Offspring = Problem.Evaluation(OffDec);	% evaluate offspring
                [Archive,State,SearchCount] = updateTable(State,SearchCount,Archive,Offspring,w,Problem.lower);
                if Problem.FE<=Problem.maxFE*T_stage
                    Mod = 1;	% stage 1 for environmental selection
                else
                    Mod = 2;	% stage 2
                end
                [Population,FrontNo,CrowdDisDec] = EnvironmentalSelection5([Population,Offspring],Problem.N,V,w,Problem.lower,Mod);
            end
        end
    end
end
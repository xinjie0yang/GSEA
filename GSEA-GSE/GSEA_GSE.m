classdef GSEA_GSE < ALGORITHM
% <multi> <real> <multimodal>
% GSEA variation without GSE method

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
                
                % 没有GSE
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
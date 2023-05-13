classdef MMEA_SaGE < ALGORITHM
% <multi> <real> <multimodal>
% Grid Search Based Multimodal Multiobjective Evolutionary Algorithm

    methods
        function main(Algorithm,Problem)
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
                    Possibility = min(exp((double(SearchCount(ColIndex(temp)))-Problem.N*0.1)/Gen),1);
                    temp = temp(rand(length(temp),1)<=Possibility);
                    if ~isempty(temp)
                        % DE operator generating a new individual, x_r1 is one of the parents, x_r2 and x_r3 are randomly selected
                        FirstP = MatingPool(temp);
                        FirstParent = Population(FirstP);
                        SecondParent = Population(randi(Problem.N,1,length(temp)));
                        ThirdParent = Population(randi(Problem.N,1,length(temp)));
                        OffDecDE = OperatorDE(Problem,FirstParent.decs,SecondParent.decs,ThirdParent.decs);
                        OffDec(temp,:) = OffDecDE;
                    end
                end
                
                Offspring = Problem.Evaluation(OffDec);	% evaluate offspring
                [Archive,State,SearchCount] = updateTable(State,SearchCount,Archive,Offspring,w,Problem.lower);
                if Problem.FE<=Problem.maxFE*0.2
                    Mod = 1;	% stage 1 for environmental selection
                else
                    Mod = 2;	% stage 2
                end
                [Population,FrontNo,CrowdDisDec] = EnvironmentalSelection5([Population,Offspring],Problem.N,V,w,Problem.lower,Mod);
            end
        end
    end
end
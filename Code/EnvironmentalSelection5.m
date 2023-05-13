function [Population,FrontNo,CrowdDisDec] = EnvironmentalSelection5(Population,N,V,w,lower,Mod)
% Sort the population according to non-dominated relationship and special crowding
    %% Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,inf);
    [n,D] = size(Population.decs);
    
    if n==N
        CrowdDisDec = CrowdingGrid(Population.decs,Population.decs,V,w,lower);
        return;
    end
    
    if Mod==1           % stage 1：diversity-first selection
        %% 每一层保持一个小格中只有一个个体
        
        % 每一个小格中保留非支配排序值最小的一个个体
        Select = false(1,n);
        if MaxFNo==1
            Select = true(1,n);
        else
            [FrontNo,srt] = sort(FrontNo);
            Population = Population(srt);
            Location = getLocation(Population.decs,V,w,lower);
            [~,ia,~] = unique(Location,'rows','stable');
            numf = 1;
            while length(ia)<N
                ib = find(FrontNo==numf);
                ia = union(ia,ib);
                numf = numf+1;
            end
            
            % convergence-first selection
            Fia = setdiff((1:n)',ia);
            FrontNo(Fia) = inf;
            Pop = Population(ia);
            [FrontNo(ia),MaxFNo] = NDSort(Pop.objs,N);
            Select(FrontNo<=MaxFNo) = true;
        end
        front = FrontNo==MaxFNo;
        front = find(front);
    else                % stage 2: convergence-first selection
        Select = false(1,n);
        k = 1;
        while 1
            front = FrontNo==k;
            front = find(front);
            Select(front) = true;
            if sum(Select)>=N
                break;
            end
            k = k+1;
        end
        MaxFNo = k;         % the critical front
    end
    
    %% truncation on critical front
    if sum(Select)>N
        NumSel = sum(Select);
        TempPop = Population(front);
        if MaxFNo>1          % critical front is not R1, find two nearest neighbors in R1
            FirstPop = Population(FrontNo==1&Select);
            DecNom1 = TempPop.decs;
            DecNom2 = FirstPop.decs;
            MaxDec = max([DecNom1;DecNom2],[],1);
            MinDec = min([DecNom1;DecNom2],[],1);
            DecNom1 = (DecNom1-repmat(MinDec,size(DecNom1,1),1))./repmat(MaxDec-MinDec,size(DecNom1,1),1);
            DecNom2 = (DecNom2-repmat(MinDec,size(DecNom2,1),1))./repmat(MaxDec-MinDec,size(DecNom2,1),1);
            dist = pdist2(DecNom1,DecNom2,'euclidean');
            [~,rank] = sort(dist,2);
            nn = min(2,size(dist,2));
            rank = rank(:,1:nn);
            GCD = CrowdingGrid(FirstPop.decs,Population(FrontNo<MaxFNo&Select).decs,V,w,lower);
            GCDSum = sum(GCD(rank),2);	% the bigger GCDSum value, the more crowded neighbors are
            [~,index] = sort(GCDSum,'descend');
            Select(front(index(1:NumSel-N))) = false;
        else            % critical front is R1, use f_DKN
			Choose = true(1,NumSel);
            Decs = TempPop.decs;
            Objs = TempPop.objs;
            while sum(Choose)>N
                CD = Crowding(Decs,Objs);
                x = find(CD==max(CD));
                index = randperm(length(x));
                y = x(index(1));
                Del = find(Choose);
                Choose(Del(y)) = false;
                Decs(y,:) = [];
                Objs(y,:) = [];
            end
            Select(front) = Choose;
        end
    end
    Population = Population(Select);
    FrontNo = FrontNo(Select);
    CrowdDisDec = CrowdingGrid(Population.decs,Population.decs,V,w,lower);
end

function GCD = CrowdingGrid(Dec1,Dec2,V,w,lower)
% grid-based crowding distance in decision space
	[n,D] = size(Dec1);
	Location1 = getLocation(Dec1,V,w,lower);
	Location2 = getLocation(Dec2,V,w,lower);
	DistCity = pdist2(Location1,Location2,'cityblock');
	GCD = sum(max(D-DistCity,0),2);
	GCD = GCD-D;
end

function [CrowdDis]=Crowding(PopDec,PopObj)
% double K-nearest neighbor method referring to CPDEA
	K = 3;
    [N, ~] = size(PopDec);
    Z = min(PopDec,[],1);
    Zmax = max(PopDec,[],1);
    popdec = (PopDec-repmat(Z,N,1))./repmat(Zmax-Z,N,1);
    distance = pdist2(popdec,popdec);
    distance(logical(eye(N))) = inf;
    distance = sort(distance,2);
    CrowdDis = sum(distance(:,1:K),2)./(mean(sum(distance(:,1:K),2))*K);
    
    Z = min(PopObj,[],1);
    Zmax = max(PopObj,[],1);
    popobj = (PopObj-repmat(Z,N,1))./repmat(Zmax-Z,N,1);
    distance = pdist2(popobj,popobj);
    distance(logical(eye(N))) = inf;
    distance = sort(distance,2);
    CrowdDis = CrowdDis+sum(distance(:,1:K),2)./(mean(sum(distance(:,1:K),2))*K);
    
    CrowdDis = 1./(1+CrowdDis);
end
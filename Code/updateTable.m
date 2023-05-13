function [Arc,S,SCount] = updateTable(State,SearchCount,Archive,Population,w,lower)
% update Table according to Archive and present population, the status of
% one subgrid is -1,0,or 1
    Decs = Population.decs;
    R = [Archive Population];
    [~,ia,~] = unique(R.decs,'rows','stable');
    R = R(ia);
    [FrontNo,~] = NDSort(R.objs,1);
    Archive = R(FrontNo==1);
    V = size(State);
    Location1 = getLocation(Archive.decs,V,w,lower);
    [Location1,ia,~] = unique(Location1,'rows','stable');
    Arc = Archive(ia);
    
    State(State==1) = -1;       % identify subgrids that are not unexplored areas
    ColIndex1 = getIndex(Location1,V(1));
    State(ColIndex1) = 1;       % identify nondominance areas
    Location = getLocation(Decs,V,w,lower);
    ColIndex = getIndex(Location,V(1));
    [m n] = hist(ColIndex,unique(ColIndex));
    SearchCount(n) = SearchCount(n)+int8(m');
    Status = State(n);
    tmp = Status==0;
    State(n(tmp)) = -1;         % identify dominance areas
    
    S = State;
    SCount = SearchCount;
end
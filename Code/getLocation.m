function Location = getLocation(Dec,V,w,lower)
% find the corresponding subgrid from Grid, according to decs of solution
    [N,~] = size(Dec);
    Loc = floor((Dec-repmat(lower,N,1))./repmat(w,N,1))+1;
    Location = max(1,min(Loc,repmat(V,N,1)));
end
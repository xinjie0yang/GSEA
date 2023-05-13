function Index = getIndex(Location,V)
% find the index of Location in table State and SearchCount
    [N,D] = size(Location);
    Location = fliplr(Location);
    bArr = repmat([1,cumprod((V-1)*ones(1,D-1),2)],N,1);
    Index = sum(Location.*bArr,2);
end
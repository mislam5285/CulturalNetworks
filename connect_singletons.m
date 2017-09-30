function adjmat=connect_singletons(adjmat,numagents)
%if the agent has no ties, create a single tie to a random agent:
for i=1:numagents;
    if isempty(find(adjmat(i,:),1))
       agidx=1:1:numagents;
       agidx(agidx==i)=[];
       k=randsample(agidx,1);
       adjmat(i,k)=1;
       adjmat(k,i)=1;
    end

end
function adjmat=rewire_no_bias(i,agents,adjmat)

if agents.last_outcome(i,1)==0 && ~isempty(find(adjmat(i,:),1))
% if the last interaction was a failure, cut ties with the previous
% interaction partner
   adjmat(i,agents.last_outcome(i,2))=nan;
   adjmat(agents.last_outcome(i,2),i)=nan;
% if it was a success then randomly pick someone that the agent is not connected to
% yet, and creat a link to them
elseif (agents.last_outcome(i)==1 || isempty(find(adjmat(i,:),1)))
    idx=find(adjmat(i,:)==0);
    new_neighbor=idx(randi(length(idx),1));
   while ~isempty(find(adjmat(i,:)==new_neighbor,1)) || new_neighbor==i
      new_neighbor=randi(length(adjmat),1); 
   end
     adjmat(i,new_neighbor)=1;
     adjmat(new_neighbor,i)=1; 
end

end
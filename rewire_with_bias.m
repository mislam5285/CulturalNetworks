function [adjmat,agents]=rewire_with_bias(i,agents,adjmat,k,pref_locus)
%Function for adjusting the social ties (neighborhoods) of agents based on
%their preferences for tags.

%Calculate the agent's base-level activations for each tag from the
%positive experiences:
pos_bla=log(agents.pref_tlast(i,:,pref_locus).^(-0.5)+2*(agents.pref_n(i,:,pref_locus)-1)./(sqrt(agents.pref_tlast(i,:,pref_locus))+sqrt(agents.pref_tfirst(i,:,pref_locus))));
%neg_bla=log(agents.bad_tlast(i,:,pref_locus).^(-0.5)+2*(agents.bad_n(i,:,pref_locus)-1)./(sqrt(agents.bad_tlast(i,:,pref_locus))+sqrt(agents.bad_tfirst(i,:,pref_locus))));


bad_inds=find(pos_bla<agents.negt(i)); %mark tags belog the negative threshold as bad
%cut bad ties:
nbhd=find(adjmat(i,:));
nbhd_inds=agents.indicators(nbhd);
for j=1:length(nbhd)
   if ~isempty(find(bad_inds==nbhd_inds(j),1))
     adjmat(i,nbhd(j))=0;
     adjmat(nbhd(j),i)=0;
   end
end


good=find(pos_bla>agents.post(i));%mark tags above the positive threshold as good
%find potential new neighbors with those tags:
potential_neighbors=[];
for j=1:length(good)
   potential_neighbors=[potential_neighbors; find(agents.indicators==good(j))];
end

potential_neighbors(potential_neighbors==i)=[];

%pick some to connect to with own dangling ties
nbhd_size=length(find(adjmat(i,:)));
if length(potential_neighbors)>k-nbhd_size && ~isempty(potential_neighbors) 
    potential_neighbors=randsample(potential_neighbors,k-nbhd_size); %select a random subset
end

%each chosen candidate has to decide whether it will accept and connect:
for j=1:length(potential_neighbors)
    %calculate the candidate's BLA for the ego's tag:
        neighbor_pref_bla=log(agents.pref_tlast(potential_neighbors(j),:,pref_locus).^(-0.5)+2*(agents.pref_n(potential_neighbors(j),:,pref_locus)-1)./(sqrt(agents.pref_tlast(potential_neighbors(j),:,pref_locus))+sqrt(agents.pref_tfirst(potential_neighbors(j),:,pref_locus))));
        %neighbor_neg_bla=log(agents.bad_tlast(potential_neighbors(j),:,pref_locus).^(-0.5)+2*(agents.bad_n(potential_neighbors(j),:,pref_locus)-1)/(sqrt(agents.bad_tlast(potential_neighbors(j),:,pref_locus))+sqrt(agents.bad_tfirst(potential_neighbors(j),:,pref_locus))));

         %if it accepts, creat tie: 
        if neighbor_pref_bla(agents.indicators(i))>agents.post(potential_neighbors(j))
            adjmat(i,potential_neighbors(j))=1;
            adjmat(potential_neighbors(j),i)=1;             
        end
end




end
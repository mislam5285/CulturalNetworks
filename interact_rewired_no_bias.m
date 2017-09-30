function [agents,pref_n,pref_t]=interact_rewired_no_bias(i,agents,partner,trait_locus,pref_n,pref_t)
%Function for interaction between two agents in unbiased version of the cultural evolution
%model

if agents.traits(partner,trait_locus)==agents.traits(i,trait_locus) %if they match boost both agents' fitness
     agents.fitness(partner)=agents.fitness(partner)+1;
     agents.fitness(i)=agents.fitness(i)+1;
     %mark the outcome as positive:
     agents.last_outcome(i,1)=1;
     agents.last_outcome(partner,1)=1;
          agents.last_outcome(i,2)=partner;
     agents.last_outcome(partner,2)=i;

else %otherwise penalize both agents
     agents.fitness(partner)=agents.fitness(partner)-1;
      agents.fitness(i)=agents.fitness(i)-1;
      %mark the outcome as negative:
     agents.last_outcome(i,1)=0;
     agents.last_outcome(partner,1)=0;
     agents.last_outcome(i,2)=partner;
     agents.last_outcome(partner,2)=i;
end

end


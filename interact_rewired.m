function agents = interact_rewired(i,agents,partner,ind_locus,partner_ind_locus,trait_locus,pref_locus,preferred_ind_variant,config)
%Function for interaction between two agents in the cultural evolution
%model
%   This function controls the interaction between two agents and its
%   results based on the active tags, traits and preferences of the two
%   actors.
agent_ind_variant=agents.indicators(i);

if agents.traits(partner,trait_locus)==agents.traits(i,trait_locus) %if they match boost both agents' fitness
      agents.fitness(partner)=agents.fitness(partner)+1;
      agents.fitness(i)=agents.fitness(i)+1;
      %updatees the preferences in the Lamarckian ACT-R scenario:
      if strcmp(config(2),'1') && strcmp(config(9),'0') 
         agents.pref_n(i,preferred_ind_variant,pref_locus)=agents.pref_n(i,preferred_ind_variant,pref_locus)+1;
         agents.pref_n(partner,agent_ind_variant,pref_locus)=agents.pref_n(partner,agent_ind_variant,pref_locus)+1;
         agents.pref_tlast(i,preferred_ind_variant,pref_locus)=0;
         agents.pref_tlast(partner,agent_ind_variant,pref_locus)=0; %reset the clock on the the last activation of this preferences
        
      elseif strcmp(config(2),'0') && strcmp(config(9),'0')
          %updates memory in the Lamarckian naive scenario:
          free_spot=find(isnan(agents.prefs(i,:,pref_locus)),1);%find empty slots in the memory vectors
          if ~isempty(free_spot)
             agents.prefs(i,free_spot,pref_locus)=agents.indicators(partner,ind_locus);%put the partner's tag in the free slot
             agents.preflens(i,free_spot,pref_locus)=0;%set its age to zero
          else
             %if there are no empty slots just replace the oldest memory
              %with the current experience
             oldest=find(agents.prefs(i,:,pref_locus)==max(agents.prefs(i,:,pref_locus)),1);
             agents.prefs(i,oldest,pref_locus)=agents.indicators(partner,ind_locus);
             agents.preflens(i,oldest,pref_locus)=0;
          end
          
          %same procedure for the partner:
          free_spot=find(isnan(agents.prefs(partner,:,pref_locus)),1);
          if ~isempty(free_spot)
             agents.prefs(partner,free_spot,pref_locus)=agents.indicators(i,partner_ind_locus);
             agents.preflens(partner,free_spot,pref_locus)=0;
          else
             oldest=find(agents.prefs(partner,:,pref_locus)==max(agents.prefs(partner,:,pref_locus)),1);
             agents.prefs(partner,oldest,pref_locus)=agents.indicators(i,partner_ind_locus);
             agents.preflens(partner,oldest,pref_locus)=0;
          end
          
          
          agents.prefs(i,:,pref_locus)=[agents.prefs(i,2:end,pref_locus) agents.indicators(partner,ind_locus)]; %FLAG
          agents.prefs(partner,:,pref_locus)=[agents.prefs(partner,2:end,pref_locus) agents.indicators(i,partner_ind_locus)];
      end
else %otherwise penalize both agents
      agents.fitness(partner)=agents.fitness(partner)-1;
      agents.fitness(i)=agents.fitness(i)-1;
      %and update memories of negative experiences in the Lamarckian ACT-R
      %scenario:
      if strcmp(config(2),'1') && strcmp(config(9),'0')
          agents.bad_n(i,preferred_ind_variant,pref_locus)=agents.bad_n(i,preferred_ind_variant,pref_locus)+1;
          agents.bad_n(partner,agent_ind_variant,pref_locus)=agents.bad_n(partner,agent_ind_variant,pref_locus)+1;
          agents.bad_tlast(i,preferred_ind_variant,pref_locus)=0;
          agents.bad_tlast(partner,agent_ind_variant,pref_locus)=0;
      end
end

end


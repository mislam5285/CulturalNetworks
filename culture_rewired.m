function [ind_shares,trait_shares,ind_skew,trait_skew,fitness_skew, ...
avg_ind_entropy,avg_trait_entropy,trait_modularity,...
ind_modularity,degree_distr,cc,avg_path_length] ...
 =culture_rewired(config,range,numagents,numgenerations,extra,r,c,plots,adjmat)

%CONFIG is a number from 0 to 512 that translates to binary where each position
%represents presence or absence of  certain features
%1st bit: 0=no 'learning'/1='learning'
%2nd: Naive/ACT memory        
%3rd: no space/space
%4th: singe/multiple indicators
%5th: single/multiple prefs
%6th: single/multiple traits
%7th: b
%8th: no network/ network
%9th: lamarckian  preferences/ genetic preferences
% Most of the bits were fixed throghout experimentation for Revay & Cioffi
% (2018). The following configurations were tested: 010000110, 010000010,
% 010000011.
% Furthermore, many of the bits are deprecated, i.e. they do not affect the
%anymore in any way, and are simply remnants from older versions of the
%model
%OTHER PARAMETERS:
%RANGE: number of possible traits and tags 
%NUMAGENTS: number of agents in the simulation
%NUMGENERATIONS: number of generations in the simulation
%EXTRA: in the ACT-R case this refers to the BLA threshold. In the naive
%case this refers to the memory 'length'.
%R: ratio of interaction/evolution (horizontal to vertical transmission)
%PLOTS: if 'plots' then charts are plotted at the end of the simulations
%ADJMAT: adjacency matrix pf the network on which the population lives
%DIFFS: a matrix derived from the ADJMAT which gives the difference between
%the expected and actual link counts. Used for modularity calculations.
%%% The model outputs a large number of variables, however, there are
%%% other variables on top of thes that the model keeps track of, that the 
%%%user can output if they wish.


%initialize these as the sizes for matrices holding info on the
%indicators/traits/prefernces
if strcmp(config(4),'1')
 ni=2;
else
 ni=1;  
end

if strcmp(config(6),'1')
 nt=2;
else
 nt=1;
end

if strcmp(config(5),'1')
 np=2;
else
 np=1;
end

%these are all book-keeping vector for in-run and output stats:
cc=nan(numgenerations,1);
degree_distr=nan(numgenerations,1);
trait_modularity=nan(numgenerations,1);
ind_modularity=nan(numgenerations,1);
two_m=length(find(adjmat));
avg_degree=round(mean(sum(adjmat)));
avg_path_length=nan(numgenerations,1);
avg_post=zeros(numgenerations,1);
avg_negt=zeros(numgenerations,1);
avg_pos_bla=zeros(numgenerations,range);
avg_neg_bla=zeros(numgenerations,range);
ind_changes=[];
trait_changes=[];
neighborhood=[];
indicator_entropy=nan(range,numgenerations,ni,nt);
trait_entropy=nan(range,numgenerations,nt,ni);
local_entropy=zeros(numgenerations,numagents,ni,nt);
indicators_interacted_avg=zeros(numgenerations,ni);
good_interactions_avg=zeros(numgenerations,1);
interactions_avg=zeros(numgenerations,1);
unique_indicators=zeros(numgenerations,ni);
unique_traits=zeros(numgenerations,nt);
fitness_skew=zeros(numgenerations,1);
shares_per_inds=zeros(numgenerations,range,ni,nt);
trait_counts_per_inds=cell(range,ni);
ind_counts=zeros(numgenerations,range,ni);
trait_counts=zeros(numgenerations,range,nt);
pref_counts=zeros(numgenerations,range,np);
ind_shares=zeros(ni,numgenerations);
trait_shares=zeros(nt,numgenerations);
ind_skew=zeros(ni,numgenerations);
trait_skew=zeros(nt,numgenerations);
pref_shares=zeros(np,numgenerations);
traits_converged=0;
traits_time_converged=2001;
inds_converged=0;
inds_time_converged=2001;
two_ind_shares=zeros(1,numgenerations);
two_inds_time_converged=2001;
two_inds_converged=0;
for i=1:range
    for k=1:ni
        trait_counts_per_inds{i,k}=zeros(numgenerations,range,nt);
    end
end

%%initialize agent variables:
indicators=randi(range,numagents,ni);
traits=randi(range,numagents,nt);
fitness=zeros(numagents,1);
good_interactions=zeros(numagents,1);
interactions=zeros(numagents,1);
indicators_interacted=zeros(numagents,range,ni);

%initialize the agent population as a structure from the previously
%initiliazed variables:
%minimal agent genotype:
agents=struct('fitness',fitness,'interactions',interactions,...
    'good_interactions',good_interactions,'indicators',indicators,...
    'traits',traits,'indicators_interacted',indicators_interacted);

%%% initialize the old 'learning' variables for th agents
if strcmp(config(1),'1') 
    learnfreq_avg=zeros(numgenerations,1);
    learnfreq=rand(numagents,1);
    agents.learnfreq=learnfreq;
    if strcmp(config(2),'0')
        learned_traits=randi(range,numagents,nt,extra);
        agents.learned_traits=learned_traits;  
    else
        learned_n=zeros(numagents,range,nt);
        learned_t=zeros(numagents,range,nt);
        agents.learned_n=learned_n;
        agents.learned_t=learned_t;
    end
end

%initialize the variables needed for ACT-R memory:
if strcmp(config(2),'1') && strcmp(config(7),'0')
    pref_n=zeros(numagents,range,np);
    pref_t=inf(numagents,range,np);
    bad_n=zeros(numagents,range,np);
    bad_t=inf(numagents,range,np);
    pref_bla_avgs=zeros(100,range,np);
     agents.pref_n=pref_n;
     agents.pref_tlast=pref_t;
     agents.pref_tfirst=pref_t;
     agents.bad_n=bad_n;
     agents.bad_tlast=bad_t;
     agents.bad_tfirst=bad_t;
%unbiased version:
elseif strcmp(config(2),'1') && strcmp(config(7),'1')
    agents.last_outcome=nan(numagents,2); % just tracks the outcome of the last interaciton
%initialize variables for naive memory:    
else
    prefs=randi(range,numagents,extra,np);
    preflens=randi(floor(r*numagents),numagents,extra,np);
    agents.prefs=prefs;
    agents.preflens=preflens;
end

%initiliaze neighoborrhoods and book-keeping variables for the moore
%neighborhood version of the model:
if strcmp(config(3),'1')
   neighborhood=make_moore_nbhd(numagents);
    num_ind_clusters=zeros(numgenerations/10,ni);
    num_trait_clusters=zeros(numgenerations/10,nt);
    ind_cluster_sizes=cell(numgenerations/10,ni);
    trait_cluster_sizes=cell(numgenerations/10,nt);
    trait_changes=zeros(numagents,nt);
    ind_changes=zeros(numagents,ni);
    agents.nbhd=neighborhood;
  
end

% extra book-keeping variables for the multiple-tag scenario
if strcmp(config(4),'1')
    inds_time_converged=[2001 2001];
    if strcmp(config(5),'1')
       p2i_map=randi(2,numagents,2); 
       agents.p2i=p2i_map;
    end
    
end
%extra book-keeping variables for the multiple preferences scenario 
if strcmp(config(5),'1')
   if strcmp(config(2),'1')
     pref_n=zeros(numagents,range,2);
     pref_t=inf(numagents,range,2);
     bad_n=zeros(numagents,range,2);
     bad_t=inf(numagents,range,2);
     pref_bla_avgs=zeros(numgenerations,range,2);
     agents.pref_n=pref_n;
     agents.pref_t=pref_t;
     agents.bad_n=bad_n;
     agents.bad_t=bad_t;
   else
     prefs=randi(range,numagents,extra,2);
     agents.prefs=prefs;
   end
end

%visualization details in case the user wants to plot results:
if strcmp(plots,'plots')
     % jsut some neat colormaps for the visualizatoins:
                trait_map=[ 0    0.0794    0.9603;
         0    0.1746    0.9127;
         0    0.2698    0.8651;
         0    0.3651    0.8175;
         0    0.4603    0.7698;
         0    0.5556    0.7222;
         0    0.6508    0.6746;
         0    0.7460    0.6270;
         0    0.8413    0.5794;
         0    0.9365    0.5317];
     
    ind_map=[1.0000    0.0794    0.9206;
    1.0000    0.1746    0.8254;
    1.0000    0.2698    0.7302;
    1.0000    0.3651    0.6349;
    1.0000    0.4603    0.5397;
    1.0000    0.5556    0.4444;
    1.0000    0.6508    0.3492;
    1.0000    0.7460    0.2540;
    1.0000    0.8413    0.1587;
    1.0000    0.9365    0.0635];

end

%initialize the BLA thresholds reasonably (based on the initial
%distribution of BLAs) in the biased scenario:
if strcmp (config(7),'0')
    %initialze the threshold to be within 1 sd of the mean BLA for each agent
    avg=mean(log(agents.pref_tlast(:,:,1).^(-0.5)+2*(agents.pref_n(:,:,1)-1)./(sqrt(agents.pref_tlast(:,:,1))+sqrt(agents.pref_tfirst(:,:,1)))),2);
    sd=std(log(agents.pref_tlast(:,:,1).^(-0.5)+2*(agents.pref_n(:,:,1)-1)./(sqrt(agents.pref_tlast(:,:,1))+sqrt(agents.pref_tfirst(:,:,1)))),0,2);
    agents.post=rand(numagents,1).*sd*2+avg-sd; 
    avg=mean(log(agents.bad_tlast(:,:,1).^(-0.5)+2*(agents.bad_n(:,:,1)-1)./(sqrt(agents.bad_tlast(:,:,1))+sqrt(agents.bad_tfirst(:,:,1)))),2);
    sd=std(log(agents.bad_tlast(:,:,1).^(-0.5)+2*(agents.bad_n(:,:,1)-1)./(sqrt(agents.bad_tlast(:,:,1))+sqrt(agents.bad_tfirst(:,:,1)))),0,2);
    agents.negt=rand(numagents,1).*sd*2+avg-sd;
end

%bits 8, and 9 do not require any special preallocation
%times at which the network data is written out:
%gens=[5 10 15 20 30 50 75 100];
%begin generation:
for n=1:numgenerations
    
    tic;
   sprintf('generation: %d', n)  
   
% write out network data at the chosen times: 
%       if ~isempty(find(gens==n,1))
%        fid = fopen(sprintf('rewiring-gen%d-nodes.csv',n),'w');
%        fprintf(fid,'Id,Trait,Indicator\n');
%        fclose(fid);
%        ids=1:1:numagents;
%        nodes=ids';
%        nodes=[nodes agents.traits agents.indicators];
%        dlmwrite(sprintf('rewiring-gen%d-nodes.csv',n),nodes,'delimiter',',','-append');
%        list=edgelist(adjmat,numagents,'unweighted','zeros');
%        csvwrite(sprintf('rewiring-edgelist-gen%d',n),list);
%       end
    
   for k=1:r
       
       %REWIRING
       order=randperm(numagents); %randomly permute the order in which the agents will be activated
       for j=1:numagents 

           i=order(j); %grab the current agents

           if strcmp(config(5),'1') && strcmp(config(6),'1')
              pref_locus=trait_locus;   %if two sets of preferences and traits, ...
          %get the one associated with the selected trait locus          
           elseif strcmp(config(5),'1') && strcmp(config(6),'0')
              pref_locus=randi(2,1);  %if two sets of preferneces but only one ...
          %trait, choose the preference locus randomly
           else 
              pref_locus=1;
           end  

           %have the agent select potential partners based on their
           %preferences:
           if strcmp(config(7),'0') %biased case
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                adjmat=rewire_with_bias(i,agents,adjmat,avg_degree,pref_locus);
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%
           else %unbiased case
                adjmat=rewire_no_bias(i,agents,adjmat,avg_degree);
           end   
       end
      
        %connect singleton nodes (to prevent the network from collapsing to a totally unconnected one)
       if strcmp(config(7),'0')
        adjmat=connect_singletons(adjmat,numagents);
       else
        adjmat=connect_singletons_no_bias(adjmat,numagents);   
       end
       
       %adjust everyone's fitness by the link maintenance costs:
       for i=1:numagents
           nbhd=find(adjmat(i,:));
           agents.fitness(i)=agents.fitness(i)-c*length(nbhd);      
       end
       
       %INTERACTION
       for j=1:numagents

           i=randi(numagents); %pick a random agent

           if strcmp(config(6),'1')
              trait_locus=randi(2,1); %randomly select cultural trait
           else
              trait_locus=1;
           end

           if strcmp(config(5),'1') && strcmp(config(6),'1')
              pref_locus=trait_locus;  %if two sets of preferences and traits, ...
          %get the one associated with the selected trait locus        
           elseif strcmp(config(5),'1') && strcmp(config(6),'0')
              pref_locus=randi(2,1); %if two sets of preferneces but only one ...
          %trait, choose the preference locus randomly
           else 
              pref_locus=1;
           end    


           if strcmp(config(5),'1') && strcmp(config(4),'1')
              ind_locus=agents.p2i(i,pref_locus); %if two sets of tags and ...
          %preferences assign the tag locus to the preference locuse based...
          %on the agent's mapping
           else    
              ind_locus=1;
           end



           nbhd=find(adjmat(i,:)); %get the agent's nbhd
           if ~isempty(nbhd)
               partner=nbhd(randsample(length(nbhd),1)); %randomly select connection
               partner_ind_variant=agents.indicators(partner); %get its tag
               %interact:
               if strcmp(config(7),'0') %baised case
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                agents=interact_rewired(i,agents,partner,ind_locus,1,trait_locus,pref_locus,partner_ind_variant,config);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
               else %unbiased case
                agents=interact_rewired_no_bias(i,agents,partner,trait_locus);
               end
           end

           %ACT-R memory update:
           if strcmp(config(2),'1') && strcmp(config(7),'0')
               agents.pref_tlast=agents.pref_tlast+1; %tick the time since last activation
               agents.bad_tlast=agents.bad_tlast+1;
               agents.pref_tfirst=agents.pref_tfirst+1; %tick the time since last activation
               agents.bad_tfirst=agents.bad_tfirst+1;
%            elseif strcmp(config(2),'1') && strcmp(config(7),'1')
%                [row,col]=find(pref_t);
%                preft_vals=nonzeros(pref_t);
%                preft_vals(isnan(preft_vals))=0;
%                preft_vals=preft_vals+1;
%                pref_t=sparse(row,col,preft_vals);
           end
           
           %record interaction stats
           agents.interactions(i)=agents.interactions(i)+1;
           agents.interactions(partner)=agents.interactions(partner)+1;
           agents.indicators_interacted(i,agents.indicators(partner),ind_locus)=1;
           agents.indicators_interacted(partner,agents.indicators(i),1)=1;
       end
   end
   
    %once all the rounds are done, create new generation of agents:
    [offspring,trait_changes,ind_changes]=create_offspring_rewired(range,numagents,ni,nt,np,extra,r,agents,ind_changes,trait_changes,neighborhood,config,n,adjmat);
    agents=offspring;
    agents.last_outcome=nan(numagents,2); %reset the last outcomes 
    
    %mark generation level statistics:
    good_interactions_avg(n)=mean(agents.good_interactions);
    interactions_avg(n)=mean(agents.interactions);
    this_indicators_interacted=zeros(numagents,ni);
    for j=1:ni
        for i=1:numagents
           this_indicators_interacted(i,j)=length(find(agents.indicators_interacted(i,:,j)));
        end
        indicators_interacted_avg(n,j)=mean(this_indicators_interacted(:,j));
        unique_indicators(n,j)=length(unique(agents.indicators(:,j)));
    end
    for j=1:nt
        unique_traits(n,j)=length(unique(agents.traits(:,j)));
    end
    
    for j=1:ni
        for i=1:range
            indices=find(agents.indicators(:,j)==i);
            for k=1:nt
                trait_counts_per_inds{i,j}(n,:,k)=histc(agents.traits(indices,k),1:1:range);
                shares_per_inds(n,i,j,k)=max(trait_counts_per_inds{i,j}(n,:,k))/sum(trait_counts_per_inds{i,j}(n,:,k));
            end
        end   
    end
    
    fitness_skew(n)=skewness(agents.fitness);
      
    avg_post(n)=mean(agents.post);
    avg_negt(n)=mean(agents.negt);
    avg_pos_bla(n,:)=mean(log(agents.pref_tlast(:,:,1).^(-0.5)+2*(agents.pref_n(:,:,1)-1)./(sqrt(agents.pref_tlast(:,:,1))+sqrt(agents.pref_tfirst(:,:,1)))));
    avg_neg_bla(n,:)=mean(log(agents.bad_tlast(:,:,1).^(-0.5)+2*(agents.bad_n(:,:,1)-1)./(sqrt(agents.bad_tlast(:,:,1))+sqrt(agents.bad_tfirst(:,:,1)))));


  
        %keep track of the distributions of tags, traits, and preferences and
    %various derived statistics:
    for j=1:ni
        ind_counts(n,:,j)=histc(agents.indicators(:,j),1:1:range);
        ind_skew(j,n)=skewness(ind_counts(n,:,j));
        ind_counts_now=ind_counts(n,:,j);
        max_ind=max(ind_counts_now);
        ind_shares(j,n)=max_ind/sum(ind_counts(n,:,j));

    end
    for j=1:nt
        trait_counts(n,:,j)=histc(agents.traits(:,j),1:1:range);
        trait_shares(j,n)=max(trait_counts(n,:,j))/sum(trait_counts(n,:,j));
        if (n==10)
           sprintf('yo') 
        end
        trait_skew(j,n)=skewness(trait_counts(n,:,j));
    end
  
    %DEPRECATED
%     for j=1:np
%         if strcmp(config(2),'1')
%             pref_bla_avgs(n,:,j)=mean(log(2*agents.pref_n(:,j).*agents.pref_t(:,j).^(-0.5)));
%             if min(pref_bla_avgs(n,:,j))<0
%                 pref_bla_avgs_plus=pref_bla_avgs(n,:,j)-min(pref_bla_avgs(n,:,j));
%             else 
%                 pref_bla_avgs_plus=pref_bla_avgs(n,:,j);
% 
%             end
%             pref_shares(j,n)=max(pref_bla_avgs_plus)/sum(pref_bla_avgs_plus);
%         else
%             prefs=nan(numagents*extra,1);
%             for k=1:extra
%                prefs((k-1)*numagents+1:k*numagents)=agents.prefs(:,k,j);
%             end
% 
% 
%             pref_counts(n,:,j)=histc(prefs,1:1:range);
%             pref_shares(j,n)=max(pref_counts(n,:,j))/sum(pref_counts(n,:,j));            
%         end
%     end
 
   %calculate avg shortest paths 
if mod(n,10)==0 %this is expensive, so we only do it every 10 generations
    d = zeros(numagents);
    for i=1:numagents 
        [d(:, i),~,~] =graphshortestpath(sparse(adjmat), i);
    end
    %path_lengths = d;
    d = d(:);
    d = d(d ~= Inf & d ~= 0); % not unreachable and not self
    avg_path_length(n) = nanmean(d); % sum and average across everything but the diagonal
end

%%%% MODULARITY CALCULATION %%%%
two_m=length(find(adjmat));
degrees=sum(adjmat);
len=length(degrees);
exp_edges=nan(numagents,numagents);
degrees_rotated=nan(numagents,numagents);
for i=1:numagents
    exp_edges(i,:)=(degrees.*(degrees(i)*ones(1,len)))/two_m;
end
diffs=adjmat-exp_edges;

   coeff=clustering2(adjmat,trait_counts(n,:,1),agents,'graph','traits',diffs);
   trait_modularity(n)=coeff/two_m;
   coeff=clustering2(adjmat,ind_counts(n,:,1),agents,'graph','inds',diffs);
   ind_modularity(n)=coeff/two_m;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate entropy values
    for i=1:range
        for j=1:ni
            for k=1:nt
                indicator_entropy(i,n,j,k)=shannon_entropy(agents.traits(:,k),agents.indicators(:,j),i,range,10);
                trait_entropy(i,n,k,j)=shannon_entropy(agents.indicators(:,j),agents.traits(:,k),i,range,10);
            end
        end    
    end
    
    inds=permute(ind_counts,[2 1 3]);
    traits=permute(trait_counts,[2 1 3]);
    avg_ind_entropy=zeros(numgenerations,ni,nt);
    avg_trait_entropy=zeros(numgenerations,nt,ni);
    for i=1:numgenerations
        for j=1:ni
            for k=1:nt
                avg_ind_entropy(i,j,k)=sum(inds(~isnan(indicator_entropy(:,i,j,k)),i).*indicator_entropy(~isnan(indicator_entropy(:,i,j,k)),i,j))/numagents;
                avg_trait_entropy(i,k,j)=sum(traits(~isnan(trait_entropy(:,i,k,j)),i).*trait_entropy(~isnan(trait_entropy(:,i,k,j)),i,k))/numagents; 
            end
        end
    end
 
 %calculate degree and clustering coefficients
    degree_distr(n)=skewness(sum(adjmat));
    cc(n)=clustering_coeff(adjmat,'undirected','unweighted');


%     if strcmp(plots,'plots')
%         for i=1:numagents
%             if strcmp(config(2),'1')
%                 [~,best_pref]=max(log(2*agents.pref_n(i,:).*agents.pref_t(i,:).^(-0.5)));
%                 if best_pref<=extra
%                     agents.prefmapping(i)=3*(best_pref-1)+1;
%                 elseif best_pref<=extra+1
%                     agents.prefmapping(i)=3*(best_pref-1)+2;
%                 else
%                     agents.prefmapping(i)=3*(best_pref-1)+3;
%                 end 
%             else
%                 mode_pref=mode(agents.prefs(i,:));
%                 mode_pref_share=length(find(agents.prefs(i,:)==mode_pref))/extra;
%                 if mode_pref_share<=0.33
%                     agents.prefmapping(i)=3*(mode_pref-1)+1;
%                 elseif mode_pref_share<=0.67
%                     agents.prefmapping(i)=3*(mode_pref-1)+2;
%                 else
%                     agents.prefmapping(i)=3*(mode_pref-1)+3;
%                 end
%             end
%         end
%     end
end


%dlmwrite(sprintf('rewiring-n%d-r%d-nodes.csv',numagents,r),nodes,'delimiter',',','-append');


%plot stuff at the end:
     if strcmp(plots,'plots')


        figure
        for k=1:ni
            subplot(1,ni,k)
            area(ind_counts(:,:,k));
            colormap(gca,ind_map)
            title('indicators')
            xlabel('time')
        end
        
        figure
        for k=1:nt
            subplot(1,nt,k)
            area(trait_counts(:,:,k));
            colormap(gca,trait_map)
            title('traits')
            xlabel('time')
        end

        figure
        for k=1:ni
            subplot(1,ni,k)
            plot(ind_shares(k,:));
            colormap(gca,ind_map)
            title('indicators')
            xlabel('time')
        end
        ylim([0 1])
        
        figure
        for k=1:nt
            subplot(1,nt,k)
            plot(trait_shares(k,:));
            colormap(gca,trait_map)
            title('traits')
            xlabel('time')
        end
        ylim([0 1])
        
        for j=1:ni
            for k=1:nt
                figure
                for i=1:range
                   subplot(range,1,i)
                   area(trait_counts_per_inds{i,j}(:,:,k))
                   colormap(trait_map)
                end
                title(sprintf('Traits %d in inds %d',k,j))
            end
        end
        
        figure
        plot(degree_distr)
        xlabel('Degree distributio skewness')
        
        figure
        plot(cc)
        xlabel('clustering coefficient')
        
        figure
        plot(ind_modularity)
        xlabel('indicator modularity')
        
        figure
        plot(trait_modularity)
        xlabel('trait modularity')
        
        figure
        plot(squeeze(avg_ind_entropy))
        xlabel('indicator entropy')
        
        figure
        plot(squeeze(avg_trait_entropy))
        xlabel('trait entropy')
        
        figure
        plot(avg_path_length)
        xlabel('avg path length')
        
        figure
        plot(avg_post)
        xlabel('avg post')
        
        figure
        plot(avg_negt)
        xlabel('avg negt')
        
        figure
        plot(avg_pos_bla)
        xlabel('avg pos bla')
        
        figure
        plot(avg_neg_bla)
        xlabel('avg neg bla')

            figure
            plot(avg_deg)
            xlabel('avg degree')
            
             figure
            plot(good_ints)
            xlabel('good ints')
            
                       figure
            plot(bad_ints)
            xlabel('bad ints')
            

    end
    
    toc;
end

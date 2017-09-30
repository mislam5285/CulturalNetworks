
numagents=1024; %number of agents
for i=1:30 %make 100 of each
   for k=[50 100 150 200] %average degrees
       for beta=[0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1] %rewiring coefficients
           adjmat=makesw(numagents,k,beta); %make the small-world matrix
           %list=edgelist(adjmat,1024);
           %cd('../networks')
           csvwrite(sprintf('small-world-adjmat-n1024-k%d-beta%d-%d.csv',k,floor(beta*100),i),adjmat);
           %csvwrite(sprintf('small-world-edgelist-n1024-k%d-beta%d-%d.csv',k,floor(beta*100),i),list);
           
           %make the difference matrices:
           two_m=length(find(adjmat));
           exp_edges=nan(numagents,numagents);

           for m=1:numagents;
               for j=1:numagents;
                    exp_edges(m,j)=(sum(adjmat(m,:))*sum(adjmat(j,:)))/two_m;
               end
           end

           diffs=adjmat-exp_edges;

           csvwrite(sprintf('small-world-diffs-n1024-k%d-beta%d-%d.csv',k,floor(beta*100),i),diffs);
           
       end
   end
end
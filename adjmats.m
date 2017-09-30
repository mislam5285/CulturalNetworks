c_int=[70]; %the different avg degrees
numagents=1000; %number of agents in each network

for m=1:length(c_int)
    for k=1:100 %make 100 of each
       rng(k)
       k
       tic;
       adjmat=zeros(numagents,numagents);
       for i=1:numagents
           for j=1:i
               if i~=j && rand<(c_int(m)/numagents)
                  adjmat(i,j)=1; 
                  adjmat(j,i)=1;
               end
           end
       end

       %here we produce the difference matrices:
       two_m=length(find(adjmat));
       exp_edges=nan(numagents,numagents);

       for i=1:numagents;
           for j=1:numagents;
                exp_edges(i,j)=(sum(adjmat(i,:))*sum(adjmat(j,:)))/two_m;
           end
       end

       diffs=adjmat-exp_edges;

       %write both matrices into a csv file:
       csvwrite(sprintf('diffs-c%d-%d.csv',c_int(m),k),diffs);
       csvwrite(sprintf('adjmat-n%d-c%d-%d.csv',numagents,c_int(m),k),adjmat);

       toc

    end
end
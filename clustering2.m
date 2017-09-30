function c=clustering2(mat,counts,agents,mode,label,diffs)
%returns the sum of the diffs matrix at coordinates where pairs of connected agents in the 'mat' matrix has the same label
%used for the network modularity calculation

%grid layout version
if strcmp(mode,'grid')

      pctexp=counts./length(mat)^2;

      coeffs=nan(length(mat)^2,1);
      for i=1:length(mat)^2
          nbhd=get_moore(i,length(mat)^2);
          nbhd_size=length(nbhd);
          num_same=length(find(mat(nbhd)==mat(i)));
          pct=num_same/nbhd_size;
          coeffs(i)=pct/pctexp(mat(i));
      end

      c=mean(coeffs);

      %graph layout version:
elseif strcmp(mode,'graph')

           %network modularity
    %define what is the label (traits or tags)
    if strcmp(label,'traits')
        labels=agents.traits;
    elseif strcmp(label,'inds')
        labels=agents.indicators;
    end
    %tic;
       idx=false(length(mat),length(mat));
       for i=1:length(mat)
           for j=1:length(mat)
               if labels(i)==labels(j)
                idx(i,j)=true;
               end
           end
       end

       c=sum(diffs(idx));
    %toc
end


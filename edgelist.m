function list=edgelist(m,n,mode,missvals)
%Function for converting and adjacency matrix 'm' into an edgelist
% 'n' - numebr of nodes
% 'mode' - 'weighted' or 'unweighted'
% 'missvals' - 'zeros' if 0 stands for missing value, otherwise 'nana'
if strcmp(missvals,'zeros')
    num_edges=length(find(m))/2;
elseif strcmp(missvals,'nan')
   num_edges=length(find(~isnan(m)))/2; 
end
if strcmp(mode,'weighted')
   list=nan(floor(num_edges),3); 
else
    list=nan(floor(num_edges),2);
end
    idx=1;
for i=1:n
    for j=1:i
      if (strcmp(missvals,'zeros') && m(i,j)~=0) || strcmp(missvals,'nan') && ~isnan(m(i,j));
         list(idx,1)=i;
         list(idx,2)=j;
         if strcmp(mode,'weighted')
             list(idx,3)=m(i,j);
         end
         idx=idx+1;
      end
    end
end

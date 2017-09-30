function entropy=shannon_entropy(source,stream,value,range,threshold)
%calculate the shannon entropy of a encoded stream of values given a soouce
%alphabet
       p=histc(source(stream==value),1:range)/length(source(stream==value)); 
       p=p(p~=0);
       summands=p(~isnan(p)).*log(p(~isnan(p)))/log(2);
       summands(isnan(summands))=0;
       if length(source(stream==value))>threshold
        entropy=-sum(summands)/length(source(stream==value));
       else
        entropy=nan;
       end       
       

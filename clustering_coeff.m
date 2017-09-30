function avg_cc=clustering_coeff(adjmat,dir,w)

%num_connected=0;
%num_closed=0;
cc=nan(length(adjmat),1);

if strcmp(w,'weighted')
  if strcmp(dir,'directed')
     for i=1:length(adjmat)
       out_nbhd=find(adjmat(i,:));
       %in_nbhd=find(adjmat(:,i));
       %nbhd=union(out_nbhd,in_nbhd);
       %if length(nbhd)<2
       if length(out_nbhd)<2
           cc(i)=0;
       else    
           num_connected=(sum(adjmat(i,out_nbhd)))*(length(out_nbhd)-1);
           %num_connected=(sum(adjmat(i,out_nbhd))+sum(adjmat(in_nbhd,i)))*(length(nbhd)-1);
           closed_weights=0;
           %for j=1:length(nbhd)
           for j=1:length(out_nbhd)
               %others=nbhd;
               others=out_nbhd;
               others(j)=[];
               for k=1:length(others)
                   %if adjmat(nbhd(j),others(k))>0 | adjmat(others(k),nbhd(j))>0 ;
                   if adjmat(out_nbhd(j),others(k))>0 | adjmat(others(k),out_nbhd(j))>0 ;
                      closed_weights=closed_weights+(adjmat(i,out_nbhd(j))+adjmat(i,others(k)))/2;
                      %closed_weights=closed_weights+(adjmat(i,nbhd(j))+adjmat(nbhd(j),i)+adjmat(others(k),i)+adjmat(i,others(k)))/4; 
                   end
               end       
           end

           cc(i)=closed_weights/num_connected;
       end
     end
  else
     for i=1:length(adjmat)
       nbhd=find(adjmat(i,:));
       if length(nbhd)<2
           cc(i)=0;
       else    
           num_connected=(sum(adjmat(i,:)))*(length(nbhd)-1);
           closed_weights=0;
           for j=1:length(nbhd)
               others=nbhd;
               others(j)=[];
               for k=1:length(others)
                   if adjmat(nbhd(j),others(k))>0
                      closed_weights=closed_weights+(adjmat(i,nbhd(j))+adjmat(i,others(k)))/2; 
                   end
               end       
           end

           cc(i)=2*closed_weights/num_connected;
       end
     end
  end
else
    if strcmp(dir,'directed')
        for i=1:length(adjmat)
           out_nbhd=find(adjmat(i,:));
           %in_nbhd=find(adjmat(:,i));
           %nbhd=union(out_nbhd,in_nbhd);
           %if length(nbhd)<2
           if length(out_nbhd)<2
               cc(i)=0;
           else    
               %num_connected=length(nbhd)*(length(nbhd)-1);
               num_connected=length(out_nbhd)*(length(out_nbhd)-1);
               num_closed=0;
               %for j=1:length(nbhd)
                for j=1:length(out_nbhd)    
                   %others=nbhd;
                   others=out_nbhd;
                   %others(j)=[];
                   for k=1:length(others)
%                        if adjmat(nbhd(j),others(k))>0
%                           num_closed=num_closed+1; 
%                        end
%                        if adjmat(others(k),nbhd(j))>0
%                          num_closed=num_closed+1; 
%                        end  
                       if adjmat(out_nbhd(j),others(k))>0 | adjmat(others(k),out_nbhd(j))>0
                          num_closed=num_closed+1; 
                       end
                   end
               end

               
           end
           
           cc(i)=num_closed/num_connected;
           
        end
    else
         for i=1:length(adjmat)
           nbhd=find(adjmat(i,:));
           if length(nbhd)<2
               cc(i)=0;
           else    
               num_connected=length(nbhd)*(length(nbhd)-1);
               num_closed=0;
               for j=1:length(nbhd)
                   others=nbhd;
                   others(j)=[];
                   for k=1:length(others)
                       if adjmat(nbhd(j),others(k))>0;
                          num_closed=num_closed+1; 
                       end
                   end       
               end

               cc(i)=2*num_closed/num_connected;
           end
         end 
        
    end
end



avg_cc=mean(cc);
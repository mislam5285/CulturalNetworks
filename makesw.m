function adjmat=makesw(n,k,beta)
%make a small-world network with N nodes, avg degree K and rewiring
%probability BETA

k=k/2;
lattice=zeros(n);
%create regular lattice
for i=1:n
    for j=i+1:i+k
        if j<=n;
            lattice(i,j)=1;
        end
    end
end

%rewire with probability beta:
[row,col]=find(lattice);
adjmat=zeros(n);
for i=1:length(row)
    if rand<beta
      old_col=col(i);
      new_col=randi(n,1);
      while new_col==i || new_col==old_col;
          new_col=randi(n,1);
      end
      adjmat(row(i),new_col)=1;
      adjmat(new_col,row(i))=1;

    else
      adjmat(row(i),col(i))=1;
      adjmat(col(i),row(i))=1;
    end
    
end
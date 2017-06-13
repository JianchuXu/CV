function [cost,Inline]=calcost(error)
cost=zeros(size(error,1),1);
Inline=zeros(size(error));
alpha=0.95;
m=1;
for i=1:size(error,1)
    tolerance=chi2inv(alpha,m);
    for j=1:size(error,2)
        if error(i,j)<=tolerance
            cost(i)=cost(i)+error(i,j);
            Inline(i,j)=1;
        else
            cost(i)=cost(i)+tolerance;
        end
    end
end

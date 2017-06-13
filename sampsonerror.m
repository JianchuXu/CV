function error=sampsonerror(data1,data2,F_in)
%Use all given 2D inhomogeneous data to calculate error
for m=1:size(F_in,3)
    F=F_in(:,:,m);
    for i=1:size(data1,2)
        da1=homo(data1(:,i));
        da2=homo(data2(:,i));
        error(m,i)=(da2'*F*da1)^2/...
            ((da2'*F(:,1))^2+(da2'*F(:,2))^2+(F(1,:)*da1)^2+(F(2,:)*da1)^2);
    end
end

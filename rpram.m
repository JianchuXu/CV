function omega=rpram(R)
I=eye(size(R));
[~,~,temp]=svd(R-I);
V=temp(:,end);
Vhat=[R(3,2)-R(2,3);
    R(1,3)-R(3,1);
    R(2,1)-R(1,2)];
sin=V'*Vhat/2;
cos=(trace(R)-1)/2;
theta=atan2(sin,cos);
omega=theta*V/norm(V);
if theta<0
    omega=omega*(1-2*pi/theta*ceil((theta-pi)/(2*pi)));
end

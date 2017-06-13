function [P1, P2]=projectF(F)
P1=[eye(3), zeros(3,1)];
[U,D,V]=svd(F);
if (det(U)<0)
    U=-U;
end
if (det(V)<0)
    V=-V;
end
W=[0, 1, 0;
    -1, 0, 0;
    0, 0, 0];
Z=[0, -1, 0;
    1, 0, 0;
    0, 0, 1];
D(3,3)=(D(1,1)+D(2,2))/2;
S=U*W*U';
M=U*Z*D*V';
P2=[M, deskew(S)];
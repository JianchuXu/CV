function out=paramderiv(in)
omega=rpram(in);
I=eye(3);
theta=norm(omega);
dtheta_domega=1/theta*omega';
dvec_domega=[0, 0, 0;
            0, 0, -1;
            0, 1, 0;
            0, 0, 1;
            0, 0, 0;
            -1, 0, 0;
            0, -1, 0;
            1, 0, 0;
            0, 0, 0];
s=(1-cos(theta))/(theta^2);
M=omega*omega';
m=reshape(M',9,1);
dm_domega=kron(omega,I)+kron(I,omega);
ds_domega=(theta*sin(theta)-2*(1-cos(theta)))/(theta^3)*dtheta_domega;
out=-reshape(I,9,1)*sin(theta)*dtheta_domega+...
    csinc(theta)*dvec_domega+...
    reshape(skew(omega)',9,1)*dcsinc(theta)*dtheta_domega+...
    s*dm_domega+...
    m*ds_domega;
if (abs(theta)<1e-5)
    out=dvec_domega;
end
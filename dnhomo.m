function [trans, result] = dnhomo(data)
%Data normalization for 2D & 3D homogeneous data
nd=size(data,1);
if nd==3
    mux=mean(data(1,:));
    muy=mean(data(2,:));
    vari=var(data(1,:))+var(data(2,:));
    s=sqrt(2/vari);
    trans=[s  0  -mux*s;
           0  s  -muy*s;
           0  0       1];
    result=trans*data;
end
if nd==4
    mux=mean(data(1,:));
    muy=mean(data(2,:));
    muz=mean(data(3,:));
    vari=var(data(1,:))+var(data(2,:))+var(data(3,:));
    s=sqrt(3/vari);
    trans=[s  0  0 -mux*s;
           0  s  0 -muy*s;
           0  0  s -muz*s;
           0  0  0      1];
    result=trans*data;
end
function cost = getcost(data3d, data2d, P)
data2d=dehomo(data2d);
data2d_est=dehomo(P*data3d);
cost=0;
for i=1:size(data3d,2)
    dif=data2d(:,i)-data2d_est(:,i);
    cost=cost+dif'*dif;
end
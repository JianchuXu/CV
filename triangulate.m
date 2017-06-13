function X=triangulate(data1, data2, F)
n=size(data1,2);
for i=1:n
    da1=homo(data1(:,i));
    da2=homo(data2(:,i));
    T1=[da1(3), 0, -da1(1);
        0, da1(3), -da1(2);
        0, 0, da1(3)];
    T2=[da2(3), 0, -da2(1);
        0, da2(3), -da2(2);
        0, 0, da2(3)];
    Fs=inv(T2)'*F*inv(T1);
    %Calculate epipoles of Fs
    [V2,~,V1]=svd(Fs);
    e1=V1(:,end);
    e2=V2(:,end);
    e1=sqrt(1/(e1(1)^2+e1(2)^2)).*e1;
    e2=sqrt(1/(e2(1)^2+e2(2)^2)).*e2;
    %Form rotation matrix
    R1=[e1(1), e1(2), 0;
        -e1(2), e1(1), 0;
        0, 0, 1];
    R2=[e2(1), e2(2), 0;
        -e2(2), e2(1), 0;
        0, 0, 1];
    %Calculate t
    Fs=R2*Fs*R1';
    f1=e1(3);f2=e2(3);
    a=Fs(2,2);  b=Fs(2,3);  c=Fs(3,2);  d=Fs(3,3);
    syms t;
    gt=t*((a*t+b)^2+f2^2*(c*t+d)^2)^2-...
        (a*d-b*c)*(1+f1^2*t^2)^2*(a*t+b)*(c*t+d);
    t_res=real(double(solve(gt==0)));
    %Take real part of t give least cost
    for j=1:size(t_res,1)
        %cost function
        s(j)=t_res(j)^2/(1+f1^2*t_res(j)^2)+...
            (c*t_res(j)+d)^2/((a*t_res(j)+b)^2+f2^2*(c*t_res(j)+d)^2);
    end
    [~,index]=min(s);
    t=t_res(index);
    %deter points as cloest points on l1 and l2
    l1=[t*f1, 1, -t]';
    l2=[-f2*(c*t+d), a*t+b, c*t+d]';
    x1l=[-l1(1)*l1(3);
        -l1(2)*l1(3);
        l1(1)^2+l1(2)^2];
    x2l=[-l2(1)*l2(3);
        -l2(2)*l2(3);
        l2(1)^2+l2(2)^2];
    %correct points back to original coordinates
    x1_cor(:,i)=homo(dehomo(inv(T1)*R1'*x1l));
    x2_cor(:,i)=homo(dehomo(inv(T2)*R2'*x2l));
    %map x to line under F
    line=F*x1_cor(:,i);
    line_or=[-line(2)*x2_cor(3,i), line(1)*x2_cor(3,i), line(2)*x2_cor(1,i)-line(1)*x2_cor(2,i)]';
    [P1, P2]=projectF(F);
    pi_plane=P2'*line_or;
    n_plane=pi_plane(1:3);
    d_plane=pi_plane(4);
    X(:,i)=[d_plane*x1_cor(:,i);
        -n_plane'*x1_cor(:,i)];
end






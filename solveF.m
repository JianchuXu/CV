function F=solveF(data1_norm, data2_norm, T1, T2)
n=size(data1_norm,2);
A=zeros(n,9);
for i=1:n
    A(i,:)=[data2_norm(1,i).*data1_norm(:,i)'...
        data2_norm(2,i).*data1_norm(:,i)'...
        data2_norm(3,i).*data1_norm(:,i)'];
end
[~,~,V]=svd(A);
a=V(:,end); %a corresponds to F1, last column of V
b=V(:,end-1); %b corresponds to F2, second from last column of V
F1=reshape(a,3,3)';
F2=reshape(b,3,3)';
syms alpha_a
detF=det(alpha_a*F1+F2);
alpha=single(solve(detF==0));
% a1=a(1);a2=a(2);a3=a(3);a4=a(4);a5=a(5);a6=a(6);a7=a(7);a8=a(8);a9=a(9);
% b1=b(1);b2=b(2);b3=b(3);b4=b(4);b5=b(5);b6=b(6);b7=b(7);b8=b(8);b9=b(9);
% coeff=[a1*a5*a9 - a1*a6*a8 - a2*a4*a9 + a2*a6*a7 + a3*a4*a8 - a3*a5*a7 ...
%     a1*a5*b9 - a1*a6*b8 - a1*a8*b6 + a1*a9*b5 - a2*a4*b9 + a2*a6*b7 + a2*a7*b6 - a2*a9*b4 + a3*a4*b8 - a3*a5*b7 - a3*a7*b5 + a3*a8*b4 + a4*a8*b3 - a4*a9*b2 - a5*a7*b3 + a5*a9*b1 + a6*a7*b2 - a6*a8*b1...
%     a1*b5*b9 - a1*b6*b8 - a2*b4*b9 + a2*b6*b7 + a3*b4*b8 - a3*b5*b7 - a4*b2*b9 + a4*b3*b8 + a5*b1*b9 - a5*b3*b7 - a6*b1*b8 + a6*b2*b7 + a7*b2*b6 - a7*b3*b5 - a8*b1*b6 + a8*b3*b4 + a9*b1*b5 - a9*b2*b4...
%     b1*b5*b9 - b1*b6*b8 - b2*b4*b9 + b2*b6*b7 + b3*b4*b8 - b3*b5*b7];
% alpha=roots(coeff);
count=1;
for i=1:size(alpha,1)
    if (isreal(alpha(i)))
        F_temp=alpha(i).*F1+F2;
        F(:,:,count)=T2'*F_temp*T1;
        count=count+1;
    end
end





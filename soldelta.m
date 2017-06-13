function [newX,  newP2, newomegau, newomegav, news]=soldelta(data1,data2,X,omegau,omegav,s,lambda)
%data1 and data2 are 2D inhomogeneous points
%X is 3D homogeneous points
n=size(data1,2);
%F parameterize then deparamterize
U_F=deRpram(omegau);
V_F=deRpram(omegav);
sigma=deparam(s);
Z=[0, -1, 0;
    1, 0, 0;
    0, 0, 1];
m=U_F*Z*diag([sigma(1),sigma(2),(sigma(1)+sigma(2))/2])*V_F';
P1=[eye(3),zeros(3,1)];
P2=[m, -U_F(:,3)];
v1=V_F(:,1); v2=V_F(:,2); v3=V_F(:,3);
sig1=sigma(1); sig2=sigma(2);
A12=paramderiv(U_F);
A22=paramderiv(V_F);
A32=derivaparam(sigma);
data1_est=dehomo(P1*X);
dif1=data1-data1_est;
data2_est=dehomo(P2*X);
dif2=data2-data2_est;
%% Jacobian and V,W
for i=1:n
    data1now=data1(:,i);
    data2now=data2(:,i);
    X_now=X(:,i);
    w1=P1(3,:)*X_now;
    w2=P2(3,:)*X_now;
    %A
    A_left=(1/w2).*[X_now', 0, 0, 0, 0, -data2now(1)*X_now';
                    0, 0, 0, 0, X_now', -data2now(2)*X_now'];
    temp=[-sig2*v2, sig1*v1, (sig1+sig2)/2*v3; 0, 0, -1];
    A11=kron(eye(3),temp);
    A1=A11*A12;
    A21=[kron(eye(3), [sig1*U_F(1,2), -sig2*U_F(1,1), (sig1+sig2)/2*U_F(1,3)]);
        zeros(1, 9);
        kron(eye(3), [sig1*U_F(2,2), -sig2*U_F(2,1), (sig1+sig2)/2*U_F(2,3)]);
        zeros(1, 9);
        kron(eye(3), [sig1*U_F(3,2), -sig2*U_F(3,1), (sig1+sig2)/2*U_F(3,3)]);
        zeros(1, 9)];
    A2=A21*A22;
    A31=[U_F(1,2)*v1+U_F(1,3)/2*v3, U_F(1,3)/2*v3-U_F(1,1)*v2;
        0, 0;
        U_F(2,2)*v1+U_F(2,3)/2*v3, U_F(2,3)/2*v3-U_F(2,1)*v2;
        0, 0;
        U_F(3,2)*v1+U_F(3,3)/2*v3, U_F(3,3)/2*v3-U_F(3,1)*v2;
        0, 0];
    A3=A31*A32;
    A_right=[A1,A2,A3];
    A(:,:,i)=A_left*A_right;
    %B'
    Bcommon=derivaparam(X_now);
    B11=(1/w1).*[P1(1,:)-data1now(1)*P1(3,:);
             P1(2,:)-data1now(2)*P1(3,:)];
    B1(:,:,i)=B11*Bcommon;   
    %B''
    B21=(1/w2).*[P2(1,:)-data2now(1)*P2(3,:);
             P2(2,:)-data2now(2)*P2(3,:)];
    B2(:,:,i)=B21*Bcommon; 
    V(:,:,i)=B1(:,:,i)'*inv(eye(2))*B1(:,:,i)...
        +B2(:,:,i)'*inv(eye(2))*B2(:,:,i);
    W(:,:,i)=A(:,:,i)'*inv(eye(2))*B2(:,:,i);
end
%% e
ea=zeros(7,1);
for i=1:n
    ea=ea+A(:,:,i)'*inv(eye(2))*dif2(:,i);
    eb(:,i)=B1(:,:,i)'*inv(eye(2))*dif1(:,i)+...
        B2(:,:,i)'*inv(eye(2))*dif2(:,i);
end

sume_temp=zeros(7,1);
for i=1:n
    sume_temp=sume_temp+W(:,:,i)*inv(V(:,:,i)+lambda*eye(3))*eb(:,i);
end
e=ea-sume_temp;
%% U S
U=zeros(7);
for i=1:n
    U=U+A(:,:,i)'*inv(eye(2))*A(:,:,i);
end

tempS=zeros(7);
for i=1:n
    tempS=tempS+W(:,:,i)*inv(V(:,:,i)+lambda*eye(3))*W(:,:,i)';
end
S=U+lambda*eye(7)-tempS;
deltaa=S\e;
for i=1:n
    deltab=inv(V(:,:,i)+lambda*eye(3))*(eb(:,i)...
        -W(:,:,i)'*deltaa);
    [Xhat,~,~]=paramn(X(:,i));
    newXhat=Xhat+deltab;
    newX(:,i)=deparam(newXhat);
end
newomegau=omegau+deltaa(1:3);
newomegav=omegav+deltaa(4:6);
news=s+deltaa(7);

newU=deRpram(newomegau);
newV=deRpram(newomegav);
newS=deparam(news);
newF=newS(1)*newU(:,1)*newV(:,1)'+newS(2)*newU(:,2)*newV(:,2)';
Z=[0, -1, 0; 1, 0, 0; 0, 0, 1];
newm=newU*Z*diag([newS(1), newS(2), (newS(1)+newS(2))/2])*newV';
newP2=[newm, -newU(:,3)];

%% feature mapping
% partab;
%% outlier rejection
clear;
% partc;
%% DLT estimate
clear;
load('data_inliers.mat');
data1=data1_next;
data2=data2_next;
[T1,data1_norm]=dnhomo(homo(data1));
[T2,data2_norm]=dnhomo(homo(data2));
n=size(data1_norm,2);
for i=1:n
    A_DLT(i,:)=kron(data2_norm(:,i)', data1_norm(:,i)');
end
[~,~,V_DLT]=svd(A_DLT);
F_temp=reshape(V_DLT(:,end),3,3)';
[U,D,Vf]=svd(F_temp);
D(3,3)=0;
F_temp=U*D*Vf';
F=T2'*F_temp*T1;
F=F./norm(F,'fro');
format longg
display(F);

%% LM
%%%%Initial estimate
lambda=0.001;
terminate=0.001;
count=1;
flag1=true;
flag2=true;
X_orig=triangulate(data1, data2, F);
[P1, P2]=projectF(F);
[U_F,sigma_temp,V_F]=svd(F);
if (det(U_F)<0)
    U_F=-U_F;
end
if (det(V_F)<0)
    V_F=-V_F;
end
omegau=rpram(U_F);
omegav=rpram(V_F);
s=paramn([sigma_temp(1,1), sigma_temp(2,2)]);
for i=1:n
    X(:,i)=deparam(paramn(X_orig(:,i)));
end
%% 
cost(count)=getcost(X,homo(data1), P1)+getcost(X,homo(data2),P2);
while(flag1)

    while(flag2)
        [newX,  newP2, newomegau, newomegav, news]=soldelta(data1,data2,X,omegau,omegav,s,lambda);
        cost0(count)=getcost(newX,homo(data1),P1)+getcost(newX,homo(data2),newP2);
        if cost(count)-cost0(count)>0
            flag2=false;
        else
            lambda=lambda*10;
        end
    end
    flag2=true;
    if cost(count)-cost0(count)<terminate
        cost(count+1)=cost0(count);
        flag1=false;
        modelu=newomegau;
        modelv=newomegav;
        models=news;
    else
        cost(count+1)=cost0(count);
        count=count+1;
        X=newX;
        omegau=newomegau;
        omegav=newomegav;
        s=news;
        lambda=0.1*lambda;
    end
end

U=deRpram(omegau);
V=deRpram(omegav);
S=deparam(s);
model=S(1)*U(:,1)*V(:,1)'+S(2)*U(:,2)*V(:,2)';

modeldisplay=model./norm(model,'fro');
disp('Eventual F');
format longg
disp(modeldisplay);
disp('cost are:')
disp(cost);
disp(cost0(end));

%%
clear;
close all;
load('Flm.mat');
F=modeldisplay;
load('data_inliers.mat');
rng(7);
num=randperm(size(data1_next,2),3);
data1=homo(data1_next(:,num));
data2=homo(data2_next(:,num));
figure;
imshow('IMG_5030.jpg');
hold on
for i=1:3
    corner=[data1(1,i)-6, data1(2,i)-6];
    rectangle('Position',[corner, 12, 12],...
        'Curvature',[1 1],'FaceColor','r','EdgeColor','r');
end
hold off;
figure;
imshow('IMG_5031.jpg');
hold on
for i=1:3
    corner=[data2(1,i)-6, data2(2,i)-6];
    rectangle('Position',[corner, 12, 12],...
        'Curvature',[1 1],'FaceColor','b','EdgeColor','b');
    linecor=F*data1(:,i);
    linecor=linecor./(-linecor(2));
    h=refline(linecor(1),linecor(3));
    h.Color='r';
    h.LineWidth=2;
end
hold off;

load('data_outliers.mat');
rng(7);
num=randperm(size(data1_out,2),3);
data1=homo(data1_out(:,num));
data2=homo(data2_out(:,num));
figure;
imshow('IMG_5030.jpg');
hold on
for i=1:3
    corner=[data1(1,i)-6, data1(2,i)-6];
    rectangle('Position',[corner, 12, 12],...
        'Curvature',[1 1],'FaceColor','r','EdgeColor','r');
end
hold off;
figure;
imshow('IMG_5031.jpg');
hold on
for i=1:3
    corner=[data2(1,i)-6, data2(2,i)-6];
    rectangle('Position',[corner, 12, 12],...
        'Curvature',[1 1],'FaceColor','b','EdgeColor','b');
    linecor=F*data1(:,i);
    linecor=linecor./(-linecor(2));
    h=refline(linecor(1),linecor(3));
    h.Color='r';
    h.LineWidth=2;
end
hold off;

clear;
load('dataconstrain.mat');
data1=flipud(data1constrain);
data2=flipud(data2constrain);

rng(7);

min_cost=Inf;
trial=1;
p=0.99;
while (1)
    num=randperm(size(data1,2),7);
    [T1,data1_norm]=dnhomo(homo(data1(:,num)));
    [T2,data2_norm]=dnhomo(homo(data2(:,num)));
    F=solveF(data1_norm, data2_norm,T1,T2);
    error=sampsonerror(data1,data2,F);
    for m=1:size(error,1)
        [cost_temp(m,:),aaa(m,:)]=calcost(error(m,:));
    end
    [cost,index]=min(cost_temp);
    inlier_temp=aaa(index,:);
    if cost<min_cost
        min_cost=cost;
        Inlier=inlier_temp;
    end
    w=size(find(Inlier==1),2)/size(data1,2);
    max_trial=log(1-p)/log(1-w^7);
    if trial>=max_trial
        break;
    end
    trial=trial+1;
end

data1_next=data1(:,Inlier==1);
data2_next=data2(:,Inlier==1);
data1_out=data1(:,Inlier~=1);
data2_out=data2(:,Inlier~=1);

figure;
imshow('IMG_5030.JPG');
hold on
for i=1:size(data1_next,2)
    plot([data1_next(1,i) data2_next(1,i)],...
        [data1_next(2,i) data2_next(2,i)],'r');
    rectangle('Position',[data1_next(1,i)-4, data1_next(2,i)-4, 9, 9],...
        'EdgeColor','r');
end
hold off
title('Feature mapping from Image1 to Image2');
set(gca,'FontSize',20);
figure;
imshow('IMG_5031.JPG');
hold on
for i=1:size(data1_next,2)
    plot([data2_next(1,i) data1_next(1,i)],...
        [data2_next(2,i) data1_next(2,i)],'r');
    rectangle('Position',[data2_next(1,i)-4, data2_next(2,i)-4, 9, 9],...
        'EdgeColor','r');
end
hold off
title('Feature mapping from Image2 to Image1');
set(gca,'FontSize',20);
disp('Number of trials:');
disp(trial-1);
disp('Number of Inliners');
disp(size(find(Inlier==1),2));
save('data_inliers.mat','data1_next', 'data2_next');
save('data_outliers.mat','data1_out', 'data2_out');


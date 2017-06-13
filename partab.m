clear ;close all;clc;
im1.data = double(rgb2gray(imread('IMG_5030.JPG')));
im2.data = double(rgb2gray(imread('IMG_5031.JPG')));
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 1/12.*[-1 8 0 -8 1];% conv kernel
% gradient images
im1.Xfiltered = imfilter(im1.data,h);im1.Yfiltered = imfilter(im1.data,h');
im2.Xfiltered = imfilter(im2.data,h);im2.Yfiltered = imfilter(im2.data,h');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thres.minor=4.5;    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
im1.minor = minor(im1.Xfiltered,im1.Yfiltered,9,thres.minor);
im2.minor = minor(im2.Xfiltered,im2.Yfiltered,9,thres.minor);
j1 = nonmax(im1.minor,9);j2 = nonmax(im2.minor,9);
[im1.Xorig_temp, im1.Yorig_temp] = find(j1~=0);
[im2.Xorig_temp, im2.Yorig_temp] = find(j2~=0);
% exclude features at corner
count = 1;
for i = 1:1:size(im1.Xorig_temp,1)
    if ( 20 < im1.Xorig_temp(i) && im1.Xorig_temp(i)< 740 && 20 < im1.Yorig_temp(i) && im1.Yorig_temp(i) < 1000)
        im1.Xorig(count,1) = im1.Xorig_temp(i);
        im1.Yorig(count,1) = im1.Yorig_temp(i);
        count = count +1;
    end
end

count = 1;
for i = 1:1:size(im2.Xorig_temp,1)
    if ( 20 < im2.Xorig_temp(i) && im2.Xorig_temp(i)< 740 && 20 < im2.Yorig_temp(i) && im2.Yorig_temp(i) < 1000)
        im2.Xorig(count,1) = im2.Xorig_temp(i);
        im2.Yorig(count,1) = im2.Yorig_temp(i);
        count = count +1;
    end
end

[im1.Xsub, im1.Ysub] = sub(im1.Xfiltered,im1.Yfiltered,...
    im1.Xorig,im1.Yorig,9);
[im2.Xsub, im2.Ysub] = sub(im2.Xfiltered,im2.Yfiltered,...
    im2.Xorig,im2.Yorig,9);

%% draw feature images
figure;imshow('IMG_5030.jpg');hold on
for i = 1:size(im1.Xsub,1)
    rectangle('Position',[im1.Ysub(i), im1.Xsub(i), 9, 9], 'EdgeColor','r');
end
hold off;title('Detected Features for Image1');
set(gca,'FontSize',15);

figure;imshow('IMG_5031.jpg');hold on
for i = 1:size(im2.Xsub,1)
    rectangle('Position',[im2.Ysub(i), im2.Xsub(i), 9, 9], 'EdgeColor','r');
end
hold off;title('Detected Features for Image2');
set(gca,'FontSize',15);
disp(size(im1.Xsub,1));
disp(size(im2.Xsub,1));
%%
% calculate the correlation coefficient with interpolation
corr = zeros(size(im1.Xorig,1),size(im2.Xorig,1));
for i = 1:size(im1.Xorig,1)
    [temp1x,temp1y]=meshgrid(im1.Ysub(i)-4:im1.Ysub(i)+4,...
        im1.Xsub(i)-4:im1.Xsub(i)+4);
    im1interp(:,:,i)=interp2(im1.data, temp1x, temp1y);
end

for j = 1:size(im2.Xorig,1)
    [temp2x,temp2y]=meshgrid(im2.Ysub(j)-4:im2.Ysub(j)+4,...
        im2.Xsub(j)-4:im2.Xsub(j)+4);
    im2interp(:,:,j)=interp2(im2.data, temp2x, temp2y);
end

for i = 1:size(corr,1)
    for j = 1:size(corr,2)
        corr(i,j) = corr2(im1interp(:,:,i),im2interp(:,:,j));
    end
end
% one to one feature matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thres.cc = 0.81;    %%%%%%%%
thres.d = 18;       %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kkk = 1;

for i = 1:size(corr,1)
    for j = 1:size(corr,2)
        if corr(i,j) > thres.cc
            curt = corr(i,j);
            corr(i,j) = -1;
            next = max(max(corr(i,:)),max(corr(:,j)));
            if (1-curt) < (1-next)*thres.d
                feature1(kkk) = i;
                feature2(kkk) = j;
                kkk = kkk +1;
            end
            corr(i,:) = -1;
            corr(:,j) = -1;
        end
    end
end
%%
disp('Features befor constrain:');
disp(size(feature1,2));
figure;
imshow('IMG_5030.jpg');
hold on
im1.Xdraw = im1.Xorig(feature1)-4;
im1.Ydraw = im1.Yorig(feature1)-4;
count=1;
for i = 1:size(feature1,2)
    data1full(:,i)=[im1.Xorig(feature1(i));im1.Yorig(feature1(i))];
    if sqrt(abs(im1.Yorig(feature1(i)) - im2.Yorig(feature2(i))).^2+abs(...
            abs(im1.Xorig(feature1(i)) - im2.Xorig(feature2(i)))).^2) < 150
        data1constrain(:,count)=[im1.Xorig(feature1(i));im1.Yorig(feature1(i))];
        rectangle('Position',[im1.Ydraw(i), im1.Xdraw(i), 9, 9],...
            'EdgeColor','r');
        plot([im1.Yorig(feature1(i)) im2.Yorig(feature2(i))],...
            [im1.Xorig(feature1(i)) im2.Xorig(feature2(i))],'Color','r');
        count=count+1;
    end
end
hold off;
title('Feature mapping from Image1 to Image2');
set(gca,'FontSize',20);

figure;
imshow('IMG_5031.jpg');
hold on
im2.Xdraw = im2.Xorig(feature2)-4;
im2.Ydraw = im2.Yorig(feature2)-4;
count=1;
for i = 1:size(feature2,2)
    data2full(:,i)=[im2.Xorig(feature2(i));im2.Yorig(feature2(i))];
    if sqrt(abs(im2.Yorig(feature2(i)) - im1.Yorig(feature1(i))).^2+abs(...
            abs(im2.Xorig(feature2(i)) - im1.Xorig(feature1(i)))).^2) < 150
        data2constrain(:,count)=[im2.Xorig(feature2(i));im2.Yorig(feature2(i))];
        rectangle('Position',[im2.Ydraw(i), im2.Xdraw(i), 9, 9],...
            'EdgeColor','r');
        plot([im2.Yorig(feature2(i)) im1.Yorig(feature1(i))],...
            [im2.Xorig(feature2(i)) im1.Xorig(feature1(i))],'Color','r');
        count=count+1;
    end
end
hold off;
title('Feature mapping from Image2 to Image1');
set(gca,'FontSize',20);
disp('Part b features after constrain of \sqrt{(x_1^2-x_2^2)+(y_1^2-y_2^2))} \leq 125: ')
disp(count-1);

save('data.mat','data1full','data2full');
save('dataconstrain.mat','data1constrain','data2constrain');


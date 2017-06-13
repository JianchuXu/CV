clear all;
close all;
clc;

load('TrainingSamplesDCT_8_new.mat');

[BG_row, BG_col] = size(TrainsampleDCT_BG);
[FG_row, FG_col] = size(TrainsampleDCT_FG);

u_BG = zeros(64,1);
u_FG = zeros(64,1);

% find the mean of the 64D gaussian, which is a vector of 64 dimentions
for i = 1:1:64
   u_BG(i,:) = mean(TrainsampleDCT_BG(:,i));
   u_FG(i,:) = mean(TrainsampleDCT_FG(:,i));
end

U_BG = zeros(64,BG_row);
for i = 1:1:BG_row
   U_BG(:,i) = u_BG;  
end

U_FG = zeros(64,FG_row);
for i = 1:1:FG_row
   U_FG(:,i) = u_FG; 
end

% find the covariance matrix of the 64D gaussian
cov_BG = (TrainsampleDCT_BG' - U_BG)*((TrainsampleDCT_BG' - U_BG)')./BG_row;
cov_FG = (TrainsampleDCT_FG' - U_FG)*((TrainsampleDCT_FG' - U_FG)')./FG_row;

% find the variance of each dimension
var_BG = zeros(1,64);
var_FG = zeros(1,64);

for i = 1:1:64
    var_BG(i) = cov_BG(i,i);
    var_FG(i) = cov_FG(i,i);   
end

% find the CCD of each feature(64 in total) of each class(2 in total)

Train_BG_sort = sort(TrainsampleDCT_BG,1);
Train_BG_sort = Train_BG_sort';
BG_difference = Train_BG_sort - U_BG;

Train_FG_sort = sort(TrainsampleDCT_FG,1);
Train_FG_sort = Train_FG_sort';
FG_difference = Train_FG_sort - U_FG;

for i =1:1:64
    p_BG = exp((BG_difference(i,:)./var_BG(i)) .* (BG_difference(i,:)).*(-0.5))./(sqrt(2*pi)*var_BG(i));
    p_FG = exp((FG_difference(i,:)./var_BG(i)) .* (FG_difference(i,:)).*(-0.5))./(sqrt(2*pi)*var_FG(i));
    
    figure,
    plot(Train_BG_sort(i,:),p_BG,Train_FG_sort(i,:),p_FG,':');
end

Zig_Zag = textread('Zig-Zag Pattern.txt');
Zig_Zag = Zig_Zag +1;
img = double(imread('cheetah.bmp'))/255;
[row, col] = size(img);
img2 = zeros(row-7,col-7);
jj = 1;
ii = 1;
feature = zeros(64,1);
%%  
p_grass = 0.8081;
p_cheet = 0.1919;


for i = 4:1:(col-4)
    for j = 4:1:(row-4)
        A = img(j-3:j+4,i-3:i+4);
        B = dct2(A);
        
        for k = 1:1:64
          [row2, col2] = find(Zig_Zag == k);
          feature(k) = B(row2,col2);
        end
        
       P_BG = exp((-0.5)*(feature - u_BG)'*(inv(cov_BG))*(feature - u_BG)) / ((2*pi)^16 * sqrt(det(cov_BG)));
       P_FG = exp((-0.5)*(feature - u_FG)'*(inv(cov_FG))*(feature - u_FG)) / ((2*pi)^16 * sqrt(det(cov_FG)));
       
       if P_BG * p_grass < P_FG * p_cheet
          img2(jj,ii) = 1;
       else
          img2(jj,ii) = 0;
       end
        jj = jj+1;
    end
    jj = 1;
    ii = ii+1;
end

figure,
imagesc(img2);
colormap(gray(255));
%% 

img3 = double(imread('cheetah_mask.bmp'))/255;
[row,col] = size(img3);
img4 = img3(4:row-4,4:col-4);
[row,col] = size(img4);
img5 = ones(row,col);
img5 = img5-img4;

P1 = (sum(sum(img4)) - sum(sum(img4.*img2)))/sum(sum(img4));
P3 = P1*p_cheet;

P2 = sum(sum(img5.*img2))/sum(sum(img5));
P4 = P2*p_grass;

error = P3+P4;

%%

% 8-D gaussian
clear all;
close all;
clc;

load('TrainingSamplesDCT_8_new.mat');

[BG_row, BG_col] = size(TrainsampleDCT_BG);
[FG_row, FG_col] = size(TrainsampleDCT_FG);

u_BG = zeros(8,1);
u_FG = zeros(8,1);

TrainsampleDCT_BG = TrainsampleDCT_BG(:,1:8);
TrainsampleDCT_FG = TrainsampleDCT_FG(:,1:8);


% find the mean of the 64D gaussian, which is a vector of 64 dimentions
for i = 1:1:8
   u_BG(i,:) = mean(TrainsampleDCT_BG(:,i));
   u_FG(i,:) = mean(TrainsampleDCT_FG(:,i));
end

U_BG = zeros(8,BG_row);
for i = 1:1:BG_row
   U_BG(:,i) = u_BG;  
end

U_FG = zeros(8,FG_row);
for i = 1:1:FG_row
   U_FG(:,i) = u_FG; 
end

% find the covariance matrix of the 8D gaussian
cov_BG = (TrainsampleDCT_BG' - U_BG)*((TrainsampleDCT_BG' - U_BG)')./BG_row;
cov_FG = (TrainsampleDCT_FG' - U_FG)*((TrainsampleDCT_FG' - U_FG)')./FG_row;

% find the variance of each dimension
var_BG = zeros(1,8);
var_FG = zeros(1,8);

for i = 1:1:8
    var_BG(i) = cov_BG(i,i);
    var_FG(i) = cov_FG(i,i);   
end

% find the CCD of each feature(8 in total) of each class(2 in total)

Train_BG_sort = sort(TrainsampleDCT_BG,1);
Train_BG_sort = Train_BG_sort';
BG_difference = Train_BG_sort - U_BG;

Train_FG_sort = sort(TrainsampleDCT_FG,1);
Train_FG_sort = Train_FG_sort';
FG_difference = Train_FG_sort - U_FG;

for i =1:1:8
    p_BG = exp((BG_difference(i,:)./var_BG(i)) .* (BG_difference(i,:)).*(-0.5))./(sqrt(2*pi)*var_BG(i));
    p_FG = exp((FG_difference(i,:)./var_BG(i)) .* (FG_difference(i,:)).*(-0.5))./(sqrt(2*pi)*var_FG(i));
    
    figure,
    plot(Train_BG_sort(i,:),p_BG,Train_FG_sort(i,:),p_FG,':');
end

Zig_Zag = textread('Zig-Zag Pattern.txt');
Zig_Zag = Zig_Zag +1;
img = double(imread('cheetah.bmp'))/255;
[row, col] = size(img);
img2 = zeros(row-7,col-7);
jj = 1;
ii = 1;
feature = zeros(8,1);

p_grass = 0.8081;
p_cheet = 0.1919;
img2=zeros(size(img));


for i = 4:1:(col-4)
    for j = 4:1:(row-4)
        A = img(j-3:j+4,i-3:i+4);
        B = dct2(A);
        
        for k = 1:1:8
          [row2, col2] = find(Zig_Zag == k);
          feature(k) = B(row2,col2);
        end
        
       P_BG = exp((-0.5)*(feature - u_BG)'*(inv(cov_BG))*(feature - u_BG)) / ((2*pi)^4 * sqrt(det(cov_BG)));
       P_FG = exp((-0.5)*(feature - u_FG)'*(inv(cov_FG))*(feature - u_FG)) / ((2*pi)^4 * sqrt(det(cov_FG)));
       
       if P_BG * p_grass < P_FG * p_cheet
          img2(i,j) = 1;
       else
          img2(i,j) = 0;
       end
    end
end

figure,
imagesc(img2);
colormap(gray(255));


img3 = double(imread('cheetah_mask.bmp'))/255;
[row,col] = size(img3);
img4 = img3(4:row-4,4:col-4);
[row,col] = size(img4);
img5 = ones(row,col);
img5 = img5-img4;

P1 = (sum(sum(img4)) - sum(sum(img4.*img2)))/sum(sum(img4));
P3 = P1*p_cheet;

P2 = sum(sum(img5.*img2))/sum(sum(img5));
P4 = P2*p_grass;

error = P3+P4;





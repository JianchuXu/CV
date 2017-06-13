clearvars -except error;
close all;
clc;

load('TrainingSamplesDCT_subsets_8.mat');

[BG_row, BG_col] = size(D3_BG);
[FG_row, FG_col] = size(D3_FG);

u_BG = zeros(64,1);
u_FG = zeros(64,1);

% find the mean of the 64D gaussian, which is a vector of 64 dimentions
for i = 1:1:64
    u_BG(i,:) = mean(D3_BG(:,i));
    u_FG(i,:) = mean(D3_FG(:,i));
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
cov_BG = (D3_BG' - U_BG)*((D3_BG' - U_BG)')./BG_row;
cov_FG = (D3_FG' - U_FG)*((D3_FG' - U_FG)')./FG_row;

% find the variance of each dimension
var_BG = zeros(1,64);
var_FG = zeros(1,64);

for i = 1:1:64
    var_BG(i) = cov_BG(i,i);
    var_FG(i) = cov_FG(i,i);
end

% find the CCD of each feature(64 in total) of each class(2 in total)

Train_BG_sort = sort(D3_BG,1);
Train_BG_sort = Train_BG_sort';
BG_difference = Train_BG_sort - U_BG;

Train_FG_sort = sort(D3_FG,1);
Train_FG_sort = Train_FG_sort';
FG_difference = Train_FG_sort - U_FG;

Zig_Zag = textread('Zig-Zag Pattern.txt');
Zig_Zag = Zig_Zag +1;
img = double(imread('cheetah.bmp'))/255;
[row, col] = size(img);
img2 = zeros(row,col);
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
            img2(j,i) = 1;
        else
            img2(j,i) = 0;
        end
    end
end
%%
pixelnum=255*270;
imagemask=im2double(imread('cheetah_mask.bmp'));
errorindex=find(imagemask~=img2);
errorch=0;
errorgr=0;
for i=1:size(errorindex)
    if imagemask(i)==1
        errorch=errorch+1;
    else
        errorgr=errorgr+1;
    end
end
errorch=errorch/size(find(imagemask==1),1);
errorgr=errorgr/(pixelnum-size(find(imagemask==1),1));
ppch=size(find(imagemask==1),1)/(pixelnum);
ppgr=1-ppch;
error.ML.D3(1)=errorch*ppch+errorgr*ppgr;
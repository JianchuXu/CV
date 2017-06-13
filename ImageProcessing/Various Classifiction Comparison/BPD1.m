for kkkkk=1:9
    close all
    clearvars -except error kkkkk
    clc;
    load ('TrainingSamplesDCT_subsets_8.mat');
    load ('Prior_1.mat');
    load ('Alpha.mat');
    p_grass = 0.8081;
    p_cheet = 0.1919;
    Zig_Zag =load('Zig-Zag Pattern.txt');
    Zig_Zag = Zig_Zag +1;
    img = double(imread('cheetah.bmp'))/255;
    [row, col] = size(img);
    img3 = double(imread('cheetah_mask.bmp'))/255;
    mu0_FG=mu0_FG';
    mu0_BG=mu0_BG';
    img2=zeros(size(img));
    
    %% Predictive distribution
    %%%%%%%%%%%%%%%%%% For BG
    feature = zeros(64,1);
    mu_BG=mean(D1_BG);
    mu_BG=mu_BG';
    
    covar_BG=cov(D1_BG,1);
    covar0_BG=diag(W0*alpha(kkkkk));
    n=size(D1_BG,1);
    mu1_BG=covar0_BG*inv(covar0_BG+1/n*covar_BG)*mu_BG+1/n*covar_BG*inv(covar0_BG+1/n*covar_BG)*mu0_BG;
    %Calculate U1
    covar1_BG=covar0_BG*inv(covar0_BG+1/n*covar_BG)*1/n*covar_BG;
    %Calculate E1
    %IN predictive distribution
    mun_BG=mu1_BG;
    covarn_BG=covar_BG+covar1_BG;
    
    %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% For FG
    mu_FG=mean(D1_FG);
    mu_FG=mu_FG';
    
    covar_FG=cov(D1_FG,1);
    covar0_FG=diag(W0*alpha(kkkkk));
    n=size(D1_FG,1);
    mu1_FG=covar0_FG*inv(covar0_FG+1/n*covar_FG)*mu_FG+1/n*covar_FG*inv(covar0_FG+1/n*covar_FG)*mu0_FG;
    %Calculate U1
    covar1_FG=covar0_FG*inv(covar0_FG+1/n*covar_FG)*1/n*covar_FG;
    %Calculate E1
    %IN predictive distribution
    mun_FG=mu1_FG;
    covarn_FG=covar_FG+covar1_FG;
    
    %% ML classification
    for i = 4:1:(col-4)
        for j = 4:1:(row-4)
            A = img(j-3:j+4,i-3:i+4);
            B = dct2(A);
            
            for k = 1:64
                [row2, col2] = find(Zig_Zag == k);
                feature(k) = B(row2,col2);
            end
            
            P_BG = exp((-0.5)*(feature - mun_BG)'*(inv(covarn_BG))*(feature - mun_BG)) / ((2*pi)^16 * sqrt(det(covarn_BG)));
            P_FG = exp((-0.5)*(feature - mun_FG)'*(inv(covarn_FG))*(feature - mun_FG)) / ((2*pi)^16 * sqrt(det(covarn_FG)));
            
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
    error.BP.D1(kkkkk)=errorch*ppch+errorgr*ppgr;
    
end



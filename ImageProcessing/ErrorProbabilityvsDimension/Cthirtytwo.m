
range=[1 2 4 8 16 24 32 40 48 56 64];
for dimension = range;
    clearvars -except dimension range error errorb;
    load ('TrainingSamplesDCT_8_new.mat');
    p_BG = 0.8081;
    p_FG = 0.1919;
    count=1;
    C=32;
    
    updateMean_BG=zeros(dimension,C)';
    updateMean_FG=zeros(dimension,C)';
    er=0.1*dimension;
    %%%%%%%%%%%%%%%%%%%%%%    Initilization      %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Mean_FG=rand(dimension,C)';
    Mean_BG=rand(dimension,C)';                 %Mean Initilization
    
    for i=1:C
        Variance_BG(i,:,:)=diag(rand(1,dimension));
    end
    for i=1:C
        Variance_FG(i,:,:)=diag(rand(1,dimension));
    end                                         %Variance Initilization
    
    Weight_BG(1:C)=1/C;
    Weight_FG(1:C)=1/C;                         %Weighting Initilization
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    Background;
    Frontground;
    
    %% BDR
    Zig_Zag = load('Zig-Zag Pattern.txt');
    Zig_Zag = Zig_Zag +1;
    img = double(imread('cheetah.bmp'))/255;
    [row, col] = size(img);
    img2 = zeros(row,col);
    feature = zeros(1,64);
    for i = 4:1:(col-4)
        for j = 4:1:(row-4)
            A = img(j-3:j+4,i-3:i+4);
            B = dct2(A);
            
            for k = 1:1:64
                [row2, col2] = find(Zig_Zag == k);
                feature(k) = B(row2,col2);
            end
            %         feature=feature';
            for m=1:C
                if (dimension ~= 1)
                    finalBDR_BG_temp(m)=Weight_BG(m)*exp((-0.5)*(feature(1:dimension) - Mean_BG(m,:))...
                        *(inv(squeeze(Variance_BG(m,:,:))))*(feature(1:dimension) - Mean_BG(m,:))')...
                        / ( sqrt((2*pi)^dimension*det(squeeze(Variance_BG(m,:,:)))));
                    finalBDR_FG_temp(m)=Weight_FG(m)*exp((-0.5)*(feature(1:dimension) - Mean_FG(m,:))...
                        *(inv(squeeze(Variance_FG(m,:,:))))*(feature(1:dimension) - Mean_FG(m,:))') ...
                        / ( sqrt((2*pi)^dimension*det(squeeze(Variance_FG(m,:,:)))));
                else
                    finalBDR_BG_temp(m)=exp((-0.5)*((feature(1:dimension) - Mean_BG(m))...
                        /Variance_BG(m))^2) / (Variance_BG(m)*sqrt(2*pi));
                    finalBDR_FG_temp(m)=exp((-0.5)*((feature(1:dimension) - Mean_FG(m))...
                        /Variance_FG(m))^2) / (Variance_FG(m)*sqrt(2*pi));
                end
            end
            finalBDR_BG=sum(finalBDR_BG_temp)*p_BG;
            finalBDR_FG=sum(finalBDR_FG_temp)*p_FG;
            if finalBDR_BG < finalBDR_FG
                img2(j,i) = 1;
            else
                img2(j,i) = 0;
            end
        end
    end
    
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
    errorb.thirtytwo{count}=errorch*ppch+errorgr*ppgr;
    count=count+1;
end



%% BackGround

while (abs(sum((updateMean_BG-Mean_BG)))>0.001)
    clearvars Class Index
    updateMean_BG=Mean_BG;
    for i=1:size(TrainsampleDCT_BG,1)
        for j=1:C
            if (dimension ~= 1)
                BDR_BG(j)=Weight_BG(j)*exp((-0.5)*(TrainsampleDCT_BG(i,1:dimension) - Mean_BG(j,:))...
                    *(inv(squeeze(Variance_BG(j,:,:))))*(TrainsampleDCT_BG(i,1:dimension) - Mean_BG(j,:))')...
                    / ( sqrt((2*pi)^dimension*det(squeeze(Variance_BG(j,:,:)))));
            else
                BDR_BG(j)=exp((-0.5)*((TrainsampleDCT_BG(i,dimension) - Mean_BG(j))/Variance_BG(j))^2)...
                    / (Variance_BG(j)*sqrt(2*pi));
            end
        end
        Class(i)=find(BDR_BG==max(BDR_BG));
    end
    
    %%%%%%%%%%%%%%%%%%%%%    M-step        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:C
        Index=find(Class==j);
        if ((~isempty(Index)) && (size(Index,2)>1))
            Mean_BG(j,:)=mean(TrainsampleDCT_BG(Index,1:dimension))';
            Variance_BG_temp=cov(TrainsampleDCT_BG(Index,1:dimension),1);
            Variance_BG(j,:,:)=diag(max(diag(Variance_BG_temp),er));
            Weight_BG(j)=size(Index,2)/size(TrainsampleDCT_BG,1);
        end
    end
    Weight_BG=Weight_BG/sum(Weight_BG);
%     iteration=iteration+1;
end
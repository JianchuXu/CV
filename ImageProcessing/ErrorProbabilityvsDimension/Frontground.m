%% FrondGround
clearvars Class Index
while (abs(sum((updateMean_FG-Mean_FG)))>0.001)
    clearvars Class Index
    updateMean_FG=Mean_FG;
    for i=1:size(TrainsampleDCT_FG,1)
        for j=1:C
            if (dimension ~= 1)
                BDR_FG(j)=Weight_FG(j)*exp((-0.5)*(TrainsampleDCT_FG(i,1:dimension) - Mean_FG(j,:))...
                    *(inv(squeeze(Variance_FG(j,:,:))))*(TrainsampleDCT_FG(i,1:dimension) - Mean_FG(j,:))')...
                    / ( sqrt((2*pi)^dimension*det(squeeze(Variance_FG(j,:,:)))));
            else
                BDR_FG(j)=exp((-0.5)*((TrainsampleDCT_FG(i,dimension) - Mean_FG(j))/Variance_FG(j))^2)...
                    / (Variance_FG(j)*sqrt(2*pi));
            end
        end
        Class(i)=find(BDR_FG==max(BDR_FG));
    end
    
    %%%%%%%%%%%%%%%%%%%%%    M-step        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=1:C
        Index=find(Class==j);
        if ((~isempty(Index)) && (size(Index,2)>1))
            Mean_FG(j,:)=mean(TrainsampleDCT_FG(Index,1:dimension))';
            Variance_FG_temp=cov(TrainsampleDCT_FG(Index,1:dimension),1);
            Variance_FG(j,:,:)=diag(max(diag(Variance_FG_temp),er));
            Weight_FG(j)=size(Index,2)/size(TrainsampleDCT_FG,1);
        end
    end
    Weight_FG=Weight_FG/sum(Weight_FG);
%     iteration=iteration+1;
end

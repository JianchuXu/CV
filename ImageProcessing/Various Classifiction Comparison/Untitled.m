% error.ML.D1(1:9)=error.ML.D1(1);
% error.ML.D2(1:9)=error.ML.D2(1);
% error.ML.D3(1:9)=error.ML.D3(1);
% error.ML.D4(1:9)=error.ML.D4(1);
plot(alpha,error.BP.D1,'LineWidth',4)
hold on
plot(alpha,error.BP.D2,'LineWidth',4)
plot(alpha,error.BP.D3,'LineWidth',4)
plot(alpha,error.BP.D4,'LineWidth',4)
set(gca,'XScale','log')
title('Baysian Parameters Stratgy 1')
legend('Data1','Data2','Data3','Data4')
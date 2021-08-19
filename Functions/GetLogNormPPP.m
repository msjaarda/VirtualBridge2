function [ECDF,ECDFRank,PPx,PPy,Fity,FitR2] = GetLogNormPPP(Data,Plot)
%GETLogNormPPP Fits, and optionally plots, BlockMaximumData
%   Data    - simply the block maximum data (max moment, shear, etc. during period)
%   Plot    - boolean

% Get ECDF and ECDFRank
[ECDF, ECDFRank] = ecdf(Data);
ECDFRank = ECDFRank'; ECDF(1) = []; ECDFRank(1) = [];
mdl = fitlm(norminv((1:length(ECDFRank))/(length(ECDFRank) + 1)),log(ECDFRank),'linear');
PPx = norminv((1:length(ECDFRank))/(length(ECDFRank) + 1));
PPy = log(ECDFRank);
Fity = mdl.Fitted;
FitR2 = mdl.Rsquared.Ordinary;
    
if Plot
    
    figure
   
    
    scatter(PPx,PPy,7,[0.6 0.6 0.6],'filled');
    hold on
    plot(PPx,mdl.Fitted,'r--','LineWidth',1);
    
    box on
    y1 = ylim;
    x1 = xlim;
    text(x1(1)+(x1(2)-x1(1))*0.25,y1(1)+(y1(2)-y1(1))*0.6,sprintf('R^2:  %.1f%%',FitR2*100),"Color",'k')
    
    ylabel('log ( Bridge Action Effect )')
    xlabel('Standard Normal Percentile')
    
end


end


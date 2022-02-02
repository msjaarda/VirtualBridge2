function [pd,x_values,y_PDF_Fit,y_CDF_Fit] = GetBlockMaxFit_60t(Data,Dist,Plot)
%GETBLOCKMAXFIT Fits, and optionally plots, BlockMaximumData
%   Data    - simply the block maximum data (max moment, shear, etc. during period)
%   Dist    - string, 'Nomral, 'Lognormal'
%   Plot    - boolean

pd = fitdist(Data,Dist);
Top = ceil(max(Data)/10)*10;
Bot = floor(min(Data)/10)*10;
TBDiff = Top-Bot;

x_values = linspace(max(0,Bot-TBDiff*.1),Top+TBDiff*.1);
y_PDF_Fit = pdf(pd,x_values);
y_CDF_Fit = cdf(pd,x_values);

if Plot
    % X is for the plot
    X = x_values;
    % x is for the bar
    x = X(1:end-1) + diff(X);
    
    y = histcounts(Data,'BinEdges',X,'normalization','pdf');
    
    figure
    bar(x,y,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5)
    hold on
    plot(X,y_PDF_Fit,'r--','LineWidth',1)
    
    % We can also put the fit type and R^2 value on the plot
    
    [MaxECDF, MaxECDFRank] = ecdf(Data); MaxECDFRank = MaxECDFRank'; MaxECDF(1) = []; MaxECDFRank(1) = [];
    if strcmp(Dist,'Normal')
        mdl = fitlm(norminv((1:length(MaxECDFRank))/(length(MaxECDFRank) + 1)),MaxECDFRank,'linear');
    elseif strcmp(Dist,'Lognormal')
        mdl = fitlm(norminv((1:length(MaxECDFRank))/(length(MaxECDFRank) + 1)),log(MaxECDFRank),'linear');
    end

    y1 = ylim;
    text(x_values(70),y1(1)+(y1(2)-y1(1))*.75,sprintf('Dist:  %s',Dist),"Color",'k')
    if strcmp(Dist,'Normal') || strcmp(Dist,'Lognormal')
        text(x_values(70),y1(1)+(y1(2)-y1(1))*.70,sprintf('R^2:    %.1f%%',mdl.Rsquared.Ordinary*100),"Color",'k')
%         y_CDF_Fit =  mdl.Rsquared.Ordinary*100; % temp!
    end
    
    set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
    ylabel('Normalized Histograms (NTS)')
    xlabel('Bridge Action Effect')

end

% [MaxECDF, MaxECDFRank] = ecdf(Data); MaxECDFRank = MaxECDFRank'; MaxECDF(1) = []; MaxECDFRank(1) = [];
%     if strcmp(Dist,'Normal')
%         mdl = fitlm(norminv((1:length(MaxECDFRank))/(length(MaxECDFRank) + 1)),MaxECDFRank,'linear');
%     elseif strcmp(Dist,'Lognormal')
%         mdl = fitlm(norminv((1:length(MaxECDFRank))/(length(MaxECDFRank) + 1)),log(MaxECDFRank),'linear');
%     end
%     
%         if strcmp(Dist,'Normal') || strcmp(Dist,'Lognormal')
%         %text(x_values(70),y1(1)+(y1(2)-y1(1))*.70,sprintf('R^2:    %.1f%%',mdl.Rsquared.Ordinary*100),"Color",'k')
%         y_CDF_Fit =  mdl.Rsquared.Ordinary*100; % temp!
%     end
    
    
end
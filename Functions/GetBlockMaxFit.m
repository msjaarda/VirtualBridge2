function [pd,x_values,y_PDF_Fit,y_CDF_Fit,pdm,x_valuesm,y_PDF_Fitm,y_CDF_Fitm] = GetBlockMaxFit(Data,Dist,Plot,PropTruck,x_values)
%GETBLOCKMAXFIT Fits, and optionally plots, BlockMaximumData
%   Data      - simply the block maximum data (max moment, shear, etc. during period)
%   Dist      - string, 'Normal, 'Lognormal'
%   Plot      - boolean
%   PropTruck - put 0 if you don't know!
%   x_values  - put 0 if you don't know!

% Data = Max.Max for Matt [old code] Data = Max for Lucas

% For Lucas Code adaption
try
    Types = Data.m;
    Data = Data.Max;
    m1 = Types == 1;
    m2 = Types == 2;
    m3 = Types == 3;
    oldcode = false;
catch
    oldcode = true;
end

%{%
if contains(Dist,'Zero') && PropTruck ~= 0 % Consider the zero if needed
    Data(end+1:round(length(Data)/PropTruck)) = 0;
    Dist = erase(Dist,'Zero');
end
%}

if oldcode
    
    % New fit
    Prop = 0.95;
    if strcmp(Dist,'Normal')
        %mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),sort(Data),'linear','Weights',[0.1*ones(round(length(Data)*Prop),1);1*ones(length(Data)-round(length(Data)*(Prop)),1)]);
        mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),sort(Data),'linear');
        pd = makedist('normal',mdlx.Coefficients.Estimate(1),mdlx.Coefficients.Estimate(2));
    elseif strcmp(Dist,'Lognormal')
        %mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear','Weights',[0.1*ones(round(length(Data)*Prop),1);1*ones(length(Data)-round(length(Data)*(Prop)),1)]);
        mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear');
        pd = makedist('lognormal',mdlx.Coefficients.Estimate(1),mdlx.Coefficients.Estimate(2));
    end
    
    % Old fit
    %pd = fitdist(Data,Dist);
    
    Top = ceil(max(Data)/10)*10;
    Bot = floor(min(Data)/10)*10;
    TBDiff = Top-Bot;
    
    % Must add try catch to be compatible with versions before specifying
    % x_values was offered
    try
        if x_values == 0
            x_values = linspace(max(0,Bot-TBDiff*.1),Top+TBDiff*.1);
        end
    catch
        x_values = linspace(max(0,Bot-TBDiff*.1),Top+TBDiff*.1);
    end
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
        
        %figure
        %bar(x,ym1*sum(m1)/height(m1)+ym2*sum(m2)/height(m2)+ym3*sum(m3)/height(m3),'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5)
        
        % We can also put the fit type and R^2 value on the plot
        
        %     [MaxECDF, MaxECDFRank] = ecdf(Data); MaxECDFRank = MaxECDFRank'; MaxECDF(1) = []; MaxECDFRank(1) = [];
        %     if strcmp(Dist,'Normal')
        %         %mdl = fitlm(norminv((1:length(MaxECDFRank))/(length(MaxECDFRank) + 1)),MaxECDFRank,'linear');
        %         mdl = fitlm(norminv((1:length(Data))/(length(Data) + 1)),sort(Data),'linear');
        %     elseif strcmp(Dist,'Lognormal')
        %         %mdl = fitlm(norminv((1:length(MaxECDFRank))/(length(MaxECDFRank) + 1)),log(MaxECDFRank),'linear');
        %         mdl = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear');
        %     end
        
        y1 = ylim;
        text(x_values(70),y1(1)+(y1(2)-y1(1))*.75,sprintf('Dist:  %s',Dist),"Color",'k')
        if strcmp(Dist,'Normal') || strcmp(Dist,'Lognormal')
            text(x_values(70),y1(1)+(y1(2)-y1(1))*.70,sprintf('R^2:    %.1f%%',mdlx.Rsquared.Ordinary*100),"Color",'k')
            %         y_CDF_Fit =  mdl.Rsquared.Ordinary*100; % temp!
        end
        
        set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
        ylabel('Normalized Histograms (NTS)')
        xlabel('Bridge Action Effect')
        
    end
else
    
    % New fit
    Prop = 0.95;
    if strcmp(Dist,'Normal')
        mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),sort(Data),'linear','Weights',[0.1*ones(round(length(Data)*Prop),1);1*ones(length(Data)-round(length(Data)*(Prop)),1)]);
        pd = makedist('normal',mdlx.Coefficients.Estimate(1),mdlx.Coefficients.Estimate(2));
    elseif strcmp(Dist,'Lognormal')
        mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear','Weights',[0.1*ones(round(length(Data)*Prop),1);1*ones(length(Data)-round(length(Data)*(Prop)),1)]);
        pd = makedist('lognormal',mdlx.Coefficients.Estimate(1),mdlx.Coefficients.Estimate(2));
    end
    
    for i=1:3
        DataTemp = eval(append('Data(m',int2str(i),')'));
        if height(DataTemp)>=35
            if strcmp(Dist,'Normal')
                eval(append('mdlxm',int2str(i),' = fitlm(norminv((1:length(DataTemp))/(length(DataTemp) + 1)),sort(DataTemp),''linear'',''Weights'',[0.1*ones(round(length(DataTemp)*Prop),1);1*ones(length(DataTemp)-round(length(DataTemp)*(Prop)),1)]);'));
                eval(append('pdm',int2str(i),' = makedist(''normal'',mdlxm',int2str(i),'.Coefficients.Estimate(1),mdlxm',int2str(i),'.Coefficients.Estimate(2));'));
            elseif strcmp(Dist,'Lognormal')
                eval(append('mdlxm',int2str(i),' = fitlm(norminv((1:length(DataTemp))/(length(DataTemp) + 1)),log(sort(DataTemp)),''linear'',''Weights'',[0.1*ones(round(length(DataTemp)*Prop),1);1*ones(length(DataTemp)-round(length(DataTemp)*(Prop)),1)]);'));
                eval(append('pdm',int2str(i),' = makedist(''lognormal'',mdlxm',int2str(i),'.Coefficients.Estimate(1),mdlxm',int2str(i),'.Coefficients.Estimate(2));'));
            end
        else
            eval(append('mdlxm',int2str(i),'.Rsquared.Ordinary = 0;'));
            eval(append('pdm',int2str(i),' = 0;'));
        end
        eval(append('Topm',int2str(i),' = ceil(max(DataTemp)/10)*10;'));
        eval(append('Botm',int2str(i),' = floor(min(DataTemp)/10)*10;'));
        eval(append('TBDiffm',int2str(i),' = Topm',int2str(i),'-Botm',int2str(i),';'));
    end
    
    % Old fit
    %pd = fitdist(Data,Dist);
    
    Top = ceil(max(Data)/10)*10;
    Bot = floor(min(Data)/10)*10;
    TBDiff = Top-Bot;
    
    x_values = linspace(max(0,Bot-TBDiff*.1),Top+TBDiff*.1);
    y_PDF_Fit = pdf(pd,x_values);
    y_CDF_Fit = cdf(pd,x_values);
    
    try
        x_valuesm1 = linspace(max(0,Botm1-TBDiffm1*.1),Topm1+TBDiffm1*.1);
        y_PDF_Fitm1 = pdf(pdm1,x_valuesm1);
        y_CDF_Fitm1 = cdf(pdm1,x_valuesm1);
    catch
        x_valuesm1 = 0; y_PDF_Fitm1 = 0; y_CDF_Fitm1 = 0;
    end
    
    try
        x_valuesm2 = linspace(max(0,Botm2-TBDiffm2*.1),Topm2+TBDiffm2*.1);
        y_PDF_Fitm2 = pdf(pdm2,x_valuesm2);
        y_CDF_Fitm2 = cdf(pdm2,x_valuesm2);
    catch
        x_valuesm2 = 0; y_PDF_Fitm2 = 0; y_CDF_Fitm2 = 0;
    end
    
    try
        x_valuesm3 = linspace(max(0,Botm3-TBDiffm3*.1),Topm3+TBDiffm3*.1);
        y_PDF_Fitm3 = pdf(pdm3,x_valuesm3);
        y_CDF_Fitm3 = cdf(pdm3,x_valuesm3);
    catch
        x_valuesm3 = 0; y_PDF_Fitm3 = 0; y_CDF_Fitm3 = 0;
    end
    
    if Plot
        % X is for the plot
        X = x_values;
        % x is for the bar
        x = X(1:end-1) + diff(X);
        
        y = histcounts(Data,'BinEdges',X,'normalization','pdf');
        ym1 = histcounts(Data(m1),'BinEdges',X,'normalization','pdf'); ym1(isnan(ym1)) = 0;
        ym2 = histcounts(Data(m2),'BinEdges',X,'normalization','pdf'); ym2(isnan(ym2)) = 0;
        ym3 = histcounts(Data(m3),'BinEdges',X,'normalization','pdf'); ym3(isnan(ym3)) = 0;
        figure
        bar(x,ym1*sum(m1)/height(m1),'EdgeColor','none','FaceColor',[0 .4470 .7410],'FaceAlpha',0.5)
        hold on
        bar(x,ym2*sum(m2)/height(m2),'EdgeColor','none','FaceColor',[.85 .3250 .098],'FaceAlpha',0.5)
        bar(x,ym3*sum(m3)/height(m3),'EdgeColor','none','FaceColor',[.929 .694 .125],'FaceAlpha',0.5)
        plot(x_valuesm2,y_PDF_Fitm2*sum(m2)/height(m2),'Color',[.85 .3250 .098],'LineStyle','--','LineWidth',1)
        plot(x_valuesm3,y_PDF_Fitm3*sum(m3)/height(m3),'Color',[.929 .694 .125],'LineStyle','--','LineWidth',1)
        plot(x_valuesm1,y_PDF_Fitm1*sum(m1)/height(m1),'Color',[0 .4470 .7410],'LineStyle','--','LineWidth',1)
        
        y1 = ylim;
        text(x_values(70),y1(1)+(y1(2)-y1(1))*.75,sprintf('Dist:  %s',Dist),"Color",'k')
        if strcmp(Dist,'Normal') || strcmp(Dist,'Lognormal')
            
            text(x_values(70),y1(1)+(y1(2)-y1(1))*.70,sprintf('R^2:    %.1f%%',mdlxm1.Rsquared.Ordinary*100),"Color",[0 .4470 .7410])
            text(x_values(70),y1(1)+(y1(2)-y1(1))*.65,sprintf('R^2:    %.1f%%',mdlxm2.Rsquared.Ordinary*100),"Color",[.85 .3250 .098])
            text(x_values(70),y1(1)+(y1(2)-y1(1))*.60,sprintf('R^2:    %.1f%%',mdlxm3.Rsquared.Ordinary*100),"Color",[.929 .694 .125])
            %         y_CDF_Fit =  mdl.Rsquared.Ordinary*100; % temp!
        end
        
        set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
        ylabel('Normalized Histograms (NTS)')
        xlabel('Bridge Action Effect')
        
        figure
        bar(x,y,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5)
        hold on
        plot(X,y_PDF_Fit,'r--','LineWidth',1)
        
        %figure
        %bar(x,ym1*sum(m1)/height(m1)+ym2*sum(m2)/height(m2)+ym3*sum(m3)/height(m3),'EdgeColor','none','FaceColor',[.929 .694 .125],'FaceAlpha',0.5)
        
        y1 = ylim;
        text(x_values(70),y1(1)+(y1(2)-y1(1))*.75,sprintf('Dist:  %s',Dist),"Color",'k')
        if strcmp(Dist,'Normal') || strcmp(Dist,'Lognormal')
            text(x_values(70),y1(1)+(y1(2)-y1(1))*.70,sprintf('R^2:    %.1f%%',mdlx.Rsquared.Ordinary*100),"Color",'k')
            %         y_CDF_Fit =  mdl.Rsquared.Ordinary*100; % temp!
        end
        
        set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
        ylabel('Normalized Histograms (NTS)')
        xlabel('Bridge Action Effect')
        
        pdm.pdm1 = pdm1; pdm.pdm2 = pdm2; pdm.pdm3 = pdm3;
        x_valuesm.x_valuesm1 = x_valuesm1; x_valuesm.x_valuesm2 = x_valuesm2; x_valuesm.x_valuesm3 = x_valuesm3;
        y_PDF_Fitm.y_PDF_Fitm1 = y_PDF_Fitm1; y_PDF_Fitm.y_PDF_Fitm2 = y_PDF_Fitm2; y_PDF_Fitm.y_PDF_Fitm3 = y_PDF_Fitm3;
        y_CDF_Fitm.y_CDF_Fitm1 = y_CDF_Fitm1; y_CDF_Fitm.y_CDF_Fitm2 = y_CDF_Fitm2; y_CDF_Fitm.y_CDF_Fitm3 = y_CDF_Fitm3;
        
    end
end
end
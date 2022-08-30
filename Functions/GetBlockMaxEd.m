function [Ed, AQ, Aq] = GetBlockMaxEd(Data,BlockM,Dist,ESIAT,ESIAEQ,ESIAEq,AQ1,AQ2,PropTruck,FitType)
%GETBLOCKMAXEd Fits, and optionally plots, BlockMaximumData
%   Data    - simply the block maximum data (max moment, shear, etc. during period)
%   BlockM  - string, 'Daily', 'Weekly', 'Monthly', 'Yearly', or 'Lifetime'
%   Dist    - string, 'Nomral, 'Lognormal'
%   PropTruck - double, proportion of special transport yearly maxima with
%   accompaniment, = 1 if not needed
%   FitType - double, 1 : original fit method, 2 : tail fitting method
%   If 'Zero' added then the zero values will be added


% Data = Max.Max for Matt [old code] or Data = Max for Lucas

% --- Calculate Real Design Value Ed ---

if strcmp(BlockM,'Yearly')
    n = 1;
elseif strcmp(BlockM,'Weekly')
    n = 1/50;
elseif strcmp(BlockM,'Daily')
    n = 1/(5*50);
elseif strcmp(BlockM,'Monthly')
    n = 1/12;
elseif strcmp(BlockM,'Lifetime')
    n = 50;
else
    n = BlockM;
end
    
Beta = norminv(1-n*0.0000013/PropTruck);
Alpha = 0.7;

if contains(Dist,'Zero') && PropTruck ~= 0 % Consider the zero if needed
    Data(end+1:round(length(Data)/PropTruck)) = 0;
    Beta = norminv(1-n*0.0000013);
    Dist = erase(Dist,'Zero');
elseif contains(Dist,'Zero')
    Dist = erase(Dist,'Zero');
end

if contains(Dist,'Plot')
    Plotter = 1;
    Dist = erase(Dist,'Plot');
else
    Plotter = 0;
end

Prop = 0.95;
if FitType == 2
    if strcmp(Dist,'Normal')
        %mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear');
        mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),sort(Data),'linear','Weights',[0.1*ones(round(length(Data)*Prop),1);1*ones(length(Data)-round(length(Data)*(Prop)),1)]);
        pd = makedist('normal',mdlx.Coefficients.Estimate(1),mdlx.Coefficients.Estimate(2));
        Em = mean(Data);
        Stdev = std(Data);
        COV = Stdev/Em;
        Delta2 = log(COV^2+1);
    elseif strcmp(Dist,'Lognormal')
        mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear','Weights',[0.1*ones(round(length(Data)*Prop),1);1*ones(length(Data)-round(length(Data)*(Prop)),1)]);
        %mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear');
        muu = mdlx.Coefficients.Estimate(1);
        sig = mdlx.Coefficients.Estimate(2);
        pd = makedist('lognormal',mdlx.Coefficients.Estimate(1),mdlx.Coefficients.Estimate(2));
        Em = exp(muu+sig^2/2);
        %Stdev = sqrt(exp(2*muu+sig^2)*(exp(sig^2)-1));
        %COV = Stdev/Em;
        Delta2 = sig^2;
    end
elseif FitType == 1
    Em = mean(Data);
    Stdev = std(Data);
    COV = Stdev/Em;
    Delta2 = log(COV^2+1);
end

if strcmp(Dist,'Normal')
    Ed = Em*(1+Alpha*Beta*COV);
    % FYI ESIAT has the 1.5 in it already...
    AQ = Ed/(ESIAT);
    Aq = ((Ed/1.5)-AQ1*ESIAEQ(1)-AQ2*ESIAEQ(2))/(sum(ESIAEq));
    if FitType == 1
        pd = makedist('normal',Em,Stdev);
        mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),sort(Data),'linear');
    end
elseif strcmp(Dist,'Lognormal')
    Ed = Em*exp(Alpha*Beta*sqrt(Delta2)-0.5*Delta2);
    AQ = Ed/(ESIAT);
    Aq = ((Ed/1.5)-AQ1*ESIAEQ(1)-AQ2*ESIAEQ(2))/(sum(ESIAEq));
    if FitType == 1
        pd = makedist('lognormal',mean(log(Data)),std(log(Data)));
        mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear');
    end
elseif strcmp(Dist,'Extreme Value')
    Ed = Em*(1 + COV*(0.45 + 0.78*log(-log(normpdf(Alpha*Beta)))));    %            exp(Alpha*Beta*sqrt(Delta2)-0.5*Delta2);
    AQ = Ed/(ESIAT);
    Aq = 1;
end

%Plot
if Plotter
    Top = ceil(max(Data)/10)*10;
    Bot = floor(min(Data)/10)*10;
    TBDiff = Top-Bot;
    x_values = linspace(max(0,Bot-TBDiff*.1),Top+TBDiff*.1);
    y_PDF_Fit = pdf(pd,x_values);
    % X is for the plot
    X = x_values;
    % x is for the bar
    x = X(1:end-1) + diff(X);
    y = histcounts(Data,'BinEdges',X,'normalization','pdf');
    
    figure
    bar(x,y,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5)
    hold on
    plot(X,y_PDF_Fit,'r--','LineWidth',1)
    
    y1 = ylim;
    plot([Ed Ed],[0 y1(1)+(y1(2)-y1(1))*.20],'k--','LineWidth',1)
    text(X(70),y1(1)+(y1(2)-y1(1))*.80,sprintf('Block:       %s',BlockM),"Color",'k')
    text(x_values(70),y1(1)+(y1(2)-y1(1))*.75,sprintf('Dist:         %s',Dist),"Color",'k')
    Rsquared = sum((y_PDF_Fit(2:end)-mean(y)).^2)/sum((y-mean(y)).^2);
    if strcmp(Dist,'Normal') || strcmp(Dist,'Lognormal')
        text(x_values(70),y1(1)+(y1(2)-y1(1))*.70,sprintf('R^2:           %.1f%%',mdlx.Rsquared.Ordinary*100),"Color",'k')
        %         y_CDF_Fit =  mdl.Rsquared.Ordinary*100; % temp!
    end
    text(X(70),y1(1)+(y1(2)-y1(1))*.65,sprintf('TailFit:      %s',mat2str(FitType == 2)),"Color",'k')
    text(X(70),y1(1)+(y1(2)-y1(1))*.60,sprintf('NoZeros:  %s',mat2str(contains(Dist,'Zero') && PropTruck ~= 0)),"Color",'k')
    text(Ed,y1(1)+(y1(2)-y1(1))*.25,sprintf('Ed = %.1f',Ed),"Color",'k','HorizontalAlignment','center')
    text(X(5),y1(1)+(y1(2)-y1(1))*.85,sprintf('%.1f%% Accompaniment Rate',100*PropTruck),"Color",'k')
    text(X(5),y1(1)+(y1(2)-y1(1))*.80,sprintf('%i Total Events',length(Data)),"Color",'k')
    
    set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
    ylabel('Normalized Histograms (NTS)')
    xlabel('Bridge Action Effect')
    
end

end
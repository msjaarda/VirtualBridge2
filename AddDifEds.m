clear, clc%, close all

% 1. Add ability to modify parent file (OutInfo) by adding 'Ed'
% Rows are fit details (type, coefficients, tail status, zeros status)
% 2. Add convergence score

% Load File
load('Output\WIM9\Mar10-22 170621.mat'); OI = OutInfo;

% In the end we can put this stiff in a loop (including file name above)
BlockM{1} = 'Yearly';  
Class{1} = 'ClassOW';
AE = [9];
Dist = 'Normal';
TailFit = false;
RemZeros = false;

for r = 1:length(AE)
    for k = 1:length(Class)
        for j = 1:length(BlockM)

% Get a subset of the maximum results for the three selections made
Data = OI.Max(AE(r)).(Class{k}).(BlockM{j}).Max;

% Set influence line you want - take a look at OutInfo.ILData to know which
ILName = OI.ILData(AE(r)).Name(7:end); ILName = strrep(ILName,'.',' ');

if strcmp(BlockM{j},'Yearly')
    n = 1;
elseif strcmp(BlockM{j},'Weekly')
    n = 1/50;
elseif strcmp(BlockM{j},'Monthly')
    n = 1/12;
end

if RemZeros
    Beta = norminv(1-n*0.0000013/OI.PropTrucks.(Class{k}).(BlockM{j})(AE(r)));
    NumData = round(length(Data)/OI.PropTrucks.(Class{k}).(BlockM{j})(AE(r)));
else
    % Ask Luca why this don't match
    Data = [Data; zeros(round(length(Data)/OI.PropTrucks.(Class{k}).(BlockM{j})(AE(r)))-length(Data),1)];
    NumData = length(Data);
    Beta = norminv(1-n*0.0000013);
end
Alpha = 0.7;

if TailFit
    Prop = 0.95;
    Weights = [0.1*ones(round(length(Data)*Prop),1);1*ones(length(Data)-round(length(Data)*(Prop)),1)];
else
    Weights = ones(length(Data),1);
end

if strcmp(Dist,'Normal')
    mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),sort(Data),'linear','Weights',Weights);
    if TailFit
        Em = mdlx.Coefficients.Estimate(1);
        Stdev = mdlx.Coefficients.Estimate(2);
    else
        Em = mean(Data);
        Stdev = std(Data);
    end
elseif strcmp(Dist,'Lognormal')
    mdlx = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear','Weights',Weights);
    muu = mdlx.Coefficients.Estimate(1);
    sig = mdlx.Coefficients.Estimate(2);
    Em = exp(muu+sig^2/2);
    Stdev = sqrt(exp(2*muu+sig^2)*(exp(sig^2)-1));
end

pd = makedist(Dist,mdlx.Coefficients.Estimate(1),mdlx.Coefficients.Estimate(2));  

COV = Stdev/Em;
Delta2 = log(COV^2+1);

if strcmp(Dist,'Normal')
    Ed = Em*(1+Alpha*Beta*COV);
elseif strcmp(Dist,'Lognormal')
    Ed = Em*exp(Alpha*Beta*sqrt(Delta2)-0.5*Delta2);
end

% Check if matches (it don't right now...)
Edx = OI.EdLN.(Class{k}).(BlockM{j})(AE(r));

% Plot
Top = ceil(max(Data)/10)*10;
Bot = floor(min(Data)/10)*10;
TBDiff = Top-Bot;

%X = linspace(max(0,Bot-TBDiff*.1),Top+TBDiff*.1);   % X is for the plot     
%X = linspace(Bot-TBDiff*.1,Top+TBDiff*.1);   % X is for the plot   
X = linspace(Bot-TBDiff*.1,max(Ed*1.05,Top+TBDiff*.1));   % X is for the plot
XB = linspace(Bot-TBDiff*.1,max(Ed*1.05,Top+TBDiff*.1),25);   % X is for the plot
x = XB(1:end-1) + diff(XB);                           % x is for the bar
y = histcounts(Data,'BinEdges',XB,'normalization','pdf');

figure
bar(x,y,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5)
hold on
plot(X,pdf(pd,X),'r--','LineWidth',1)
y1 = ylim; ylim([0 y1(2)*1.2]); y1 = ylim;
plot([Ed Ed],[0 y1(1)+(y1(2)-y1(1))*.20],'k--','LineWidth',1)

% Add Text
title(ILName)
text(X(70),y1(1)+(y1(2)-y1(1))*.80,sprintf('Block:       %s',BlockM{j}),"Color",'k')
text(X(70),y1(1)+(y1(2)-y1(1))*.75,sprintf('Dist:         %s',Dist),"Color",'k')
text(X(70),y1(1)+(y1(2)-y1(1))*.70,sprintf('R^2:           %.1f%%',mdlx.Rsquared.Ordinary*100),"Color",'k')
text(X(70),y1(1)+(y1(2)-y1(1))*.65,sprintf('TailFit:      %s',mat2str(TailFit)),"Color",'k')
text(X(70),y1(1)+(y1(2)-y1(1))*.60,sprintf('NoZeros:  %s',mat2str(RemZeros)),"Color",'k')
text(Ed,y1(1)+(y1(2)-y1(1))*.25,sprintf('Ed = %.1f',Ed),"Color",'k','HorizontalAlignment','center')
text(X(5),y1(1)+(y1(2)-y1(1))*.85,sprintf('%.1f%% Accompaniment Rate',100*OI.PropTrucks.(Class{k}).(BlockM{j})(AE(r))),"Color",'k')
text(X(5),y1(1)+(y1(2)-y1(1))*.80,sprintf('%i Total Events',NumData),"Color",'k')

set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
ylabel('Normalized Histograms (NTS)')
xlabel('Bridge Action Effect')

z = 1;
% Convergence
for p = 5:5:100
    r = randperm(height(Data),round(height(Data)*p/100));
    Data2 = Data(r);
    
    if strcmp(Dist,'Normal')
        mdlx = fitlm(norminv((1:length(Data2))/(length(Data2) + 1)),sort(Data2));
        if TailFit
            Em = mdlx.Coefficients.Estimate(1);
            Stdev = mdlx.Coefficients.Estimate(2);
        else
            Em = mean(Data2);
            Stdev = std(Data2);
        end
    elseif strcmp(Dist,'Lognormal')
        mdlx = fitlm(norminv((1:length(Data2))/(length(Data2) + 1)),log(sort(Data2)));
        muu = mdlx.Coefficients.Estimate(1);
        sig = mdlx.Coefficients.Estimate(2);
        Em = exp(muu+sig^2/2);
        Stdev = sqrt(exp(2*muu+sig^2)*(exp(sig^2)-1));
    end
    
    pd = makedist(Dist,mdlx.Coefficients.Estimate(1),mdlx.Coefficients.Estimate(2));
    
    COV = Stdev/Em;
    Delta2 = log(COV^2+1);
    
    if strcmp(Dist,'Normal')
        Ed2(z) = Em*(1+Alpha*Beta*COV);
    elseif strcmp(Dist,'Lognormal')
        Ed2(z) = Em*exp(Alpha*Beta*sqrt(Delta2)-0.5*Delta2);
    end
    z = z+1;

end

figure
plot(5:5:100,Ed2)
ylim([0 mean(Ed2)*1.3])
ylabel('E_{d}')
xlabel('% of all Data')
title('Sample Size Verification for Extreme Value')
yline(Ed2(end))

        end
    end
end

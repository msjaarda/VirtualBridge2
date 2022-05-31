function pd = GetFit(Data,BlockM,DistTypes,Plot,IncZ)
% This function will give back pd and gof (goodness of fit)
% pd will be a structure pd.(Dist)
% pd.(Dist).pd
% pd.(Dist).IncZ
% pd.(Dist).PropZ
% pd.(Dist).gof
% pd.(Dist).Ed
% pd.(Dist).R2 (if linear model)
% pd.edcf, pd.ecdfx
% pd.Best = name of Dist with lowest gof
% can use pd.(pd.Best) to access best fit

% Turn off fit warning
warning('off','stats:gevfit:IterLimit')

% Consider adding convergence score

% Gives Ed as well, although a more complete GetEd includes support for
% different alphas?

% If you give Dist == 'All', all fits will be performed and given

% Dist can be:
% 'Normal' fitdist
    % 'NormalLM' linearmodel (has r2)
% 'Lognormal' fitdist
    % 'LognormalLM' linearmodel (has r2)
% 'Lognormal' linearmodel w/ tailfit
% 'gev' Generalized Extreme Value fitdist
% 'gev Gumbel' Generalized Extreme Value k == 0 Type 1 Max equation
% 'All'


% In future turn off: msg =

%     'Maximum likelihood estimation did not converge.  Iteration limit exceeded.'
% 
% 
% warnID =
% 
%     'stats:gevfit:IterLimit'
    
    % but save events in pd

if strcmp(DistTypes,'All')
    if length(Data) < 30
        DistTypes = {'Normal', 'Lognormal', 'LognormalTF'};
    else
        DistTypes = {'Normal', 'Lognormal', 'LognormalTF', 'gev', 'gevGumbel'};
    end
elseif ~iscell(DistTypes) % Turn Dist into cell if it is individual
    DistTypes = cellstr(DistTypes);
else % if manually selected, remove gev for less than 30 data
    if length(Data) < 30 || length(Data(Data>0)) < 10
        DistTypes(contains(DistTypes,'gev')) = [];
    end
end

n = GetnBlockM(BlockM);

if ~IncZ
    % Delete zeros from Data, but save proportion
    PropZ = sum(Data == 0)/length(Data);
    Data(Data == 0) = [];
    Beta = norminv(1-n*0.0000013/PropZ);
    
    if isempty(Data) % do a normal fit if Data is empty
    Data = 0;
    Dist = "Normal";    
    warning('off','stats:LinearModel:RankDefDesignMat')
    mdl = fitlm(norminv((1:length(Data))/(length(Data) + 1)),sort(Data),'linear');
    pd.(Dist).pd = makedist('normal',mdl.Coefficients.Estimate(1),mdl.Coefficients.Estimate(2));
    pd.(Dist).R2 = mdl.Rsquared.Ordinary*100;
    pd.Best = Dist;
    pd.(Dist).Ed = Data;
    return
    end
else
    Beta = norminv(1-n*0.0000013);
end
Alpha = 0.7;

% Gather a few stats from Data
Em = mean(Data);
Stdev = std(Data);
COV = Stdev/Em;
Delta2 = log(COV^2+1);

% Get EdSIA values... only works for Normal!
EdSIA(1) = Em*(1+Alpha*Beta*COV);
EdSIA(2) = Em*exp(Alpha*Beta*sqrt(Delta2)-0.5*Delta2);
EdSIA(3) = Em*(1 + COV*(0.45 + 0.78*log(-log(normpdf(Alpha*Beta)))));
% See Table C3 Eurocode 1990
EdEC(1) = Em*(1+Alpha*Beta*COV);
EdEC(2) = Em*exp(Alpha*Beta*COV);
EdEC(3) = (Em+0.577/(pi/(Stdev*sqrt(6))))+(1/(pi/(Stdev*sqrt(6))))*log(-log(normpdf(Alpha*Beta)));

for k = 1:length(DistTypes)
    Dist = DistTypes{k};
    if strcmp(Dist,'gevGumbel')
        % So we can use GEV afterall with fixed k = 0? Method of moments
        gamma = 0.5772;
        sigmaHat = sqrt(6)*std(Data)/pi;
        muHat = mean(Data) - gamma*sigmaHat;
        pd.(Dist).pd = makedist('GeneralizedExtremeValue','k',0,'sigma',sigmaHat,'mu',muHat);
    elseif strcmp(Dist,'LognormalTF')
        Prop = 0.95;
        Weight = [0.1*ones(round(length(Data)*Prop),1);1*ones(length(Data)-round(length(Data)*(Prop)),1)];
        mdl = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear','Weights',Weight);
        pd.(Dist).pd = makedist('lognormal',mdl.Coefficients.Estimate(1),mdl.Coefficients.Estimate(2));
        pd.(Dist).R2 = mdl.Rsquared.Ordinary*100;
    elseif  strcmp(Dist,'NormalLM')
        mdl = fitlm(norminv((1:length(Data))/(length(Data) + 1)),sort(Data),'linear');
        pd.(Dist).pd = makedist('normal',mdl.Coefficients.Estimate(1),mdl.Coefficients.Estimate(2));
        pd.(Dist).R2 = mdl.Rsquared.Ordinary*100;
    elseif  strcmp(Dist,'LognormalLM') % Min Sum ^2 Error
        mdl = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear');
        pd.(Dist).pd = makedist('normal',mdl.Coefficients.Estimate(1),mdl.Coefficients.Estimate(2));
        pd.(Dist).R2 = mdl.Rsquared.Ordinary*100;
    else
        pd.(Dist).pd = fitdist(Data,Dist); % Maximum Likelihood
    end
end
[pd.ecdf, pd.ecdfx] = ecdf(Data);

% CheckSIA(1,:) = [-norminv(1-cdf(pd.Normal.pd,EdSIA(1))), Beta*Alpha];
% CheckSIA(2,:) = [-norminv(1-cdf(pd.Lognormal.pd,EdSIA(2))), Beta*Alpha];
% CheckSIA(3,:) = [-norminv(1-cdf(pd.gevGumbel.pd,EdSIA(3))), Beta*Alpha];

if Plot
    % Method #1 fitdist
    Top = ceil(max(Data)/10)*10;
    Bot = floor(min(Data)/10)*10;
    TBDiff = Top-Bot;
    x_values = linspace(max(0,Bot-TBDiff*.1),Top+TBDiff*.1);
    
    figure('Name','Dists','NumberTitle','off'), hold on
    X = x_values;               % X is for the plot
    x = X(1:end-1) + diff(X);   % x is for the bar
    y = histcounts(Data,'BinEdges',X,'normalization','pdf');
    bar(x,y,'EdgeColor','none','FaceColor',[.6 .6 .6],'FaceAlpha',0.5,'DisplayName','Data')
    C = linspecer(length(DistTypes));
    
    for k = 1:length(DistTypes)
        Dist = DistTypes{k};
        plot(X,pdf(pd.(Dist).pd,x_values),'-','Color',C(k,:),'LineWidth',1,'DisplayName',Dist)
    end
    
    set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
    ylabel('Normalized Histogram and Fit')
    xlabel('Bridge Action Effect')
    title('Fits'); legend('location','best'); box on

end

for k = 1:length(DistTypes)
    Dist = DistTypes{k};
    pd.(Dist).Ed = icdf(pd.(Dist).pd,1-normcdf(-Beta*Alpha));
end
pd.ecdfEd = interp1(pd.ecdf,pd.ecdfx,1-normcdf(-Beta*Alpha),'linear','extrap');

if Plot
    % Probability Paper
    fGumbelPP = figure('Name','Probability Paper','NumberTitle','off'); hold on
    scatter(sort(Data),-log(-log((1:length(Data))/(length(Data) + 1))),7,'k','filled','DisplayName','Daily Max Data');
    C = linspecer(length(DistTypes));
    
    for k = 1:length(DistTypes)
        Dist = DistTypes{k};
        plot(sort(Data),-log(-log(cdf(pd.(Dist).pd,sort(Data)))),'--','Color',C(k,:),'LineWidth',1,'DisplayName',Dist)
    end
    
    xlabel('X - Load Effect'); ylabel('-log(-log(Probability of non-exceedance))')
    title('Gumbel Probability Paper'); legend('location','best'); box on
end

% Goodness of Fit -Tail
Vals = [0.9:0.01:0.99 0.999 0.9999 0.99993];
EDataTail = interp1(pd.ecdf,pd.ecdfx,Vals,'linear','extrap');
% Initialize pd.Best
pd.Best = DistTypes{1};
for k = 1:length(DistTypes)
    Dist = DistTypes{k};
    pd.(Dist).gof = sqrt(sum((EDataTail-icdf(pd.(Dist).pd,Vals)).^2));
    if pd.(Dist).gof < pd.(pd.Best).gof
        pd.Best = Dist;
    end
end

end

% Demo

% clear, clc, close all
% 
% load('C:\Users\mjsja\Desktop\SwissTraffic2\Misc\BMAxles.mat')
% Data = BMAxles.CH.MaxAx.Tandem.Class.Weekly.Max; clear BMAxles;
% Dist = 'All';
% Plot = 1;
% IncZ = 1;
% BlockM = 'Weekly';
% 
% pd = GetFit(Data,BlockM,Dist,Plot,IncZ);
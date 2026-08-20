function pd = GetFit_alphaConsistent(Data,BlockM,DistTypes,Plot,IncZ)
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

% 2024.04.15 Qianhui Yu
% The treatment for alphaE is changed for VIM simulation. AlphaE is
% considered to be corresponding to yearly target relibility index. (in older versions, AlphaE is multiplied to yearly target index for JAM simulation and it is multiplied to weekly target for VIM. This is not consistent)
% ! Note: not Beta in this function refers to the target reliability for
% design action effect for a given reference period already including the
% influence of alphaE. No additional alphaE is needed to be multiplied to
% Beta. This is a siginificant big change comparing with prvious versions.

%2026.07.23 Lucas time of the event 1 --> 100%, for exemple jammed can take
%place only 2% of the time so EventDuration = 0.02; The value below is a
%default value, meaning that if it has been well defined in the Input File
%it should do it automatically
EventDuration = 1;

if length(IncZ) > 1
    BETATarget = IncZ(2);
    IncZ(2) = [];
    if length(IncZ) > 1
        EventDuration = IncZ(2);
        IncZ(2) = [];
    end
else
    BETATarget = 4.2;
end

% 2024.04.15 Qianhui Yu alpha is applied to the yearly target reliability
Alpha = 0.7;

%BETATarget = 4.2; %4.7 CHANGE BACK!!
PFTarget = 1-normcdf(BETATarget*Alpha);

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
        %DistTypes = {'Normal', 'Lognormal', 'LognormalTF', 'gev', 'gevGumbel','GeneralizedPareto'};
    end
elseif ~iscell(DistTypes) % Turn Dist into cell if it is individual
    DistTypes = cellstr(DistTypes);
else % if manually selected, remove gev for less than 30 data
    if length(Data) < 30 || length(Data(Data>0)) < 10
        DistTypes(contains(DistTypes,'gev')) = [];
    end
end

n = GetnBlockM(BlockM);
n = n/EventDuration;

if ~IncZ
    if sum(Data == 0) == 0
        Beta = norminv(1-n*PFTarget);
    else
    % Delete zeros from Data, but save proportion
    PropZ = sum(Data == 0)/length(Data);
    Data(Data == 0) = [];
    Beta = norminv(1-n*PFTarget/PropZ);
    
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
    end
else
    Beta = norminv(1-n*PFTarget);
    
end



% Gather a few stats from Data
Em = mean(Data);
Stdev = std(Data);
COV = Stdev/Em;
Delta2 = log(COV^2+1);

% Get EdSIA values... only works for Normal!
EdSIA(1) = Em*(1+1*Beta*COV);
EdSIA(2) = Em*exp(1*Beta*sqrt(Delta2)-0.5*Delta2);
EdSIA(3) = Em*(1 + COV*(0.45 + 0.78*log(-log(normpdf(1*Beta)))));
% See Table C3 Eurocode 1990
EdEC(1) = Em*(1+1*Beta*COV);
EdEC(2) = Em*exp(1*Beta*COV);
EdEC(3) = (Em+0.577/(pi/(Stdev*sqrt(6))))+(1/(pi/(Stdev*sqrt(6))))*log(-log(normpdf(1*Beta)));

[pd.ecdf, pd.ecdfx] = ecdf(Data);
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
        pd.(Dist).pd = makedist('lognormal',mdl.Coefficients.Estimate(1),mdl.Coefficients.Estimate(2));
        pd.(Dist).R2 = mdl.Rsquared.Ordinary*100;
    elseif  strcmp(Dist,'GeneralizedPareto')
        % --- 1. Threshold ---
        Prop = 0.95;
        u = quantile(Data, Prop);
        % --- 2. Excesses ---
        Excess = Data(Data > u) - u;
        Nu = length(Excess);
        nTot = length(Data);
        if Nu < 10
            warning('Few exceedances for GPD fit')
        end
        % --- 3. Fit GPD ---
        pd.(Dist).pd = fitdist(Excess,'GeneralizedPareto');
        % --- 4. Store POT info ---
        pd.(Dist).u = u;
        pd.(Dist).p_u = Nu / nTot;
        % --- 5. Custom quantile function ---
        xi = pd.(Dist).pd.k;
        betaGP = pd.(Dist).pd.sigma;
        p_u = pd.(Dist).p_u;
        if abs(xi) < 1e-6
        pd.(Dist).icdf = @(p) u - betaGP * log((1-p)./p_u);
        else
        pd.(Dist).icdf = @(p) u + (betaGP/xi)*((( (1-p)./p_u ).^(-xi)) - 1);
        end
        pd.(Dist).pdf = @(x) (x > u) .* (p_u .* pdf(pd.(Dist).pd, x - u));
        [Xuniq, ia] = unique(pd.ecdfx);
        ECDFuniq = pd.ecdf(ia);
        pd.(Dist).cdf = @(x) (x <= u) .* interp1(Xuniq, ECDFuniq, x, 'linear', 'extrap') + (x > u)  .* (1 - p_u .* (1 - cdf(pd.(Dist).pd, x - u)));
    else
        pd.(Dist).pd = fitdist(Data,Dist); % Maximum Likelihood
    end
end

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
        if strcmp(Dist,'GeneralizedPareto')
            plot(X,pd.(Dist).pdf(x_values),'-','Color',C(k,:),'LineWidth',1,'DisplayName',Dist)
            %xline(pd.(Dist).u,'k--','Threshold u')
        else
            plot(X,pdf(pd.(Dist).pd,x_values),'-','Color',C(k,:),'LineWidth',1,'DisplayName',Dist)
        end
    end
    
    set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
    ax = gca;
    text(ax.XLim(1)+0.62*(ax.XLim(2)-ax.XLim(1)),ax.YLim(1)+0.30*(ax.YLim(2)-ax.YLim(1)),append(int2str(size(Data,1)),' ',BlockM,' Maxima Values'));
    ylabel('Normalized Histogram and Fit')
    xlabel('Bridge Action Effect')
    title('Fits'); legend('location','best'); box on

end

for k = 1:length(DistTypes)
    Dist = DistTypes{k};
    if strcmp(Dist,'GeneralizedPareto')
    pd.(Dist).Ed = pd.(Dist).icdf(1-normcdf(-Beta));
    else
    pd.(Dist).Ed = icdf(pd.(Dist).pd,1-normcdf(-Beta*1));
    end
    
    % Probability of failure
    PF = n/1000;
    if strcmp(Dist,'GeneralizedPareto')
    pd.(Dist).E1000 = pd.(Dist).icdf(1-PF);
    else
    pd.(Dist).E1000 = icdf(pd.(Dist).pd,1-PF);
    end
    
end
% The reason we can't take ecdfEd is because it can never be larger than
% your largest data value... the last data value is defined with the cdf
% value of 1! Therefore... we have to approximate...
% Note that the top value won't show on the scatter - it is infinity!
% Neither will zero... but this is added by matlab into the ecdfx for the
% purpose of plotting stairs.
pd.ecdfEd = interp1(pd.ecdf,pd.ecdfx,1-normcdf(-Beta*1),'linear','extrap');
% % XXOld = sort(Data);
% % YYOld = ((1:length(Data))/(length(Data) + 1))';
% % [XXXOld, ia, ~] = unique(XXOld);
% % YYYOld = YYOld(ia);
% XX = [min(Data); sort(Data)];
% YY = ((0:length(Data))/length(Data))';
% [XXX, ia, ~] = unique(XX);
% YYY = YY(ia);
% pd.ecdfEd = interp1(YYYOld,XXXOld,1-normcdf(-Beta*Alpha),'linear','extrap');

if Plot
    % Probability Paper
    fGumbelPP = figure('Name','Probability Paper','NumberTitle','off'); hold on
    %scatter(sort(Data),-log(-log((1:length(Data))/(length(Data) + 1))),7,'k','filled','DisplayName','Max Data');
    %scatter([min(Data); sort(Data)],-log(-log((0:length(Data))/(length(Data)+0.00001))),7,'k','filled','DisplayName','Max Data');
    scatter([min(Data); sort(Data)],-log(-log((0:length(Data))/(length(Data)))),7,'k','filled','DisplayName','Max Data');
    C = linspecer(length(DistTypes));
    
    %TempX = 220:480;
    for k = 1:length(DistTypes)
        Dist = DistTypes{k};
        if strcmp(Dist,'GeneralizedPareto')
            Xs = sort(Data);
            FX = pd.(Dist).cdf(Xs);
            plot(Xs,-log(-log(FX)),'--','Color',C(k,:),'LineWidth',1,'DisplayName',Dist)
        else
        plot(sort(Data),-log(-log(cdf(pd.(Dist).pd,sort(Data)))),'--','Color',C(k,:),'LineWidth',1,'DisplayName',Dist)
        end
        %plot(TempX,-log(-log(cdf(pd.(Dist).pd,TempX))),'--','Color',C(k,:),'LineWidth',1,'DisplayName',Dist)
    end
    
    xlabel('X - Load Effect'); ylabel('-log(-log(Probability of non-exceedance))')
    title('Gumbel Probability Paper'); legend('location','northwest'); box on
end

% Goodness of Fit -Tail - SPECIAL ATTENTION PAID! THIS IS VERY TRICKY...
% ONGOING MONITOR.. BASICALLY LINSPACE COULD BE GOOD BECAUSE WHEN THERE IS
% LITTLE DATA, THE .99, .999, .9999. .99993 ALL COMPARE TO A THE MAX
% EMPIRICAL DIST

% New algorithm attempt... say we want to pay attention to sample size...


Vals = [0.9:0.01:0.99 0.999 0.9999 0.99993];
% Try this next
%Vals = linspace(0.9,(1-1/length(Data)),12);

%Vals = [0.9:0.01:(1-1/length(Data))];
EDataTail = interp1(pd.ecdf,pd.ecdfx,Vals,'linear','extrap');
% Initialize pd.Best
pd.Best = DistTypes{1};
% We need a smarter selection method - must take into account how much data
% we have...
for k = 1:length(DistTypes)
    Dist = DistTypes{k};
        if strcmp(Dist,'GeneralizedPareto')
        pd.(Dist).gof = sqrt(sum((EDataTail-pd.(Dist).icdf(Vals)).^2));
        else
        pd.(Dist).gof = sqrt(sum((EDataTail-icdf(pd.(Dist).pd,Vals)).^2));
        end
    if pd.(Dist).gof < pd.(pd.Best).gof
        pd.Best = Dist;
    end
end

if Plot
    %correction Lucas pour text out of the fit 30.07.2026
    ax = gca;
    text(ax.XLim(1)+0.65*(ax.XLim(2)-ax.XLim(1)),ax.YLim(1)+0.25*(ax.YLim(2)-ax.YLim(1)),append('Best Fit: ',pd.Best));
    text(ax.XLim(1)+0.65*(ax.XLim(2)-ax.XLim(1)),ax.YLim(1)+0.20*(ax.YLim(2)-ax.YLim(1)),append('Ed: ',num2str(pd.(pd.Best).Ed,5)));
    text(ax.XLim(1)+0.65*(ax.XLim(2)-ax.XLim(1)),ax.YLim(1)+0.15*(ax.YLim(2)-ax.YLim(1)),append('Gof: ',num2str(pd.(pd.Best).gof,5)));
    text(ax.XLim(1)+0.65*(ax.XLim(2)-ax.XLim(1)),ax.YLim(1)+0.10*(ax.YLim(2)-ax.YLim(1)),append('\beta_{corrected}=',num2str(Beta,3)));
    %xplot= [min(Data); sort(Data)]; xplot = xplot(end-1);
    %yplot= -log(-log((0:length(Data))/(length(Data)))); yplot = yplot(end-1);
    %scatter(xplot,yplot,10,'r','filled','DisplayName','Data of interest');
    %text(xplot*0.8,yplot*1.05,append('AE: ',num2str(xplot,5)),'Color','red');
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
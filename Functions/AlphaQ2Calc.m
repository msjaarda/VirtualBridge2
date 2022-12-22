function [AQ,AQ2,Summary] = AlphaQ2Calc(Country,FitPlots,AdvPlots,DispTab,BETATarget)
% ALPHAQ2CALC
% Calculates AlphaQ using:
% Country (or BMAxles Group Name) such as "CH", "DE", "US" or "AT"
% FitPlots/AdvPlots/DispTab Booleans 
% Returns AQ (value for Weekly BlockMax, with 1) All 2) Class+ 3) Class
 
% VBAxleStats    >> BMAxles       >> AlphaQ1Calc
% Q1Q2Stats      >> BMQ1Q2        >> AlphaQ2Calc
% VBWIMq         >> qMaxEvents    >> qInvestigation

load('BMQ1Q2.mat')

BM = {'Daily', 'Weekly', 'Yearly'};             % j
ClassType = {'All', 'ClassOW', 'Class'};        % i
ClassT = {'All', 'Classified+', 'Classified'};
DistTypes = {'All'};                      % k

Max = BMQ1Q2.(Country).Max;

% Initialize Fig Positions
Ri = 0; Up = 0;

Dist = DistTypes{1};

% FIGURE OF ALL DISTS
if FitPlots
    
figure('Position',[Ri Up 1200 400],'Name',[Country ' Q1Q2' ' Fits'],'NumberTitle','off');  Ri = Ri+25; Up = Ri+25;

% Set colours
C = linspecer(9);
% ScaleDown
ScaleDown = [1 2.5 5];
% X Stuff
Step = 5;
LimitL = 100;
LimitR = 700;
X = LimitL:Step:LimitR;
x = X(1:end-1) + diff(X);
    
for i = 1:length(ClassType)
    Class = ClassType{i};
    
    subplot(1,3,i)
    hold on
    
    for j = 1:length(BM)
        BlockM = BM{j};
        
        
        y = histcounts(Max.(Class).(BlockM).Max,'BinEdges',X,'normalization','pdf');
        bar(x,y/ScaleDown(j),1,'EdgeColor',C(j,:),'FaceColor',[.8 .8 .8],'FaceAlpha',0.5,'DisplayName',BlockM)
        
        Dist = DistTypes{1};
        %pd = GetFit(Max.(Class).(BlockM).Max,BlockM,[Dist 'TF'],1,[1 BETATarget]);
        pd = GetFit(Max.(Class).(BlockM).Max,BlockM,[Dist],1,[1 BETATarget]);
        plot(x,pdf(pd.(pd.Best).pd,x)/ScaleDown(j),'k-','DisplayName',[Dist 'Fit'])
        
    end
    
    % Set Plot Details
    set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
    box on
    xlim([LimitL+25 LimitR-25])
    if i == 1
        ylabel('Normalized Histograms (NTS)')
    end
    xlabel('Total Strip Load (kN)')
    title([ClassT{i}])
    legend('location','northeast')
    
end

% FIGURE OF ALL DISTS SEPARATED
figure('Position',[50 50 1200 400],'Name',[Country ' Q1Q2' ' Fits Separated'],'NumberTitle','off');  Ri = Ri+25; Up = Ri+25;

% X Stuff
Step = 2.5;
LimitL = 100; LimitR = 700;
Xp = LimitL-Step/2:Step:LimitR-Step/2;

for j = 1:length(BM)
    BlockM = BM{j};
    
    subplot(1,3,j)
    hold on
    k = 1;
    
    Class = 'Class';
    y = histcounts(Max.(Class).(BlockM).Max,'BinEdges',X,'normalization','pdf');
    bar(x,y/ScaleDown(j),1,'EdgeColor','k','FaceColor',[.8 .8 .8],'FaceAlpha',0.8,'DisplayName',Class)
    
    Dist = DistTypes{k};
    pd = GetFit(Max.(Class).(BlockM).Max,BlockM,[Dist],0,[1 BETATarget]);
    plot(x,pdf(pd.(pd.Best).pd,x)/ScaleDown(j),'k-','DisplayName',[Dist 'Fit'])    
    
    Class = 'All';
    y = histcounts(Max.(Class).(BlockM).Max,'BinEdges',X,'normalization','pdf');
    bar(x,y/ScaleDown(j),1,'EdgeColor',C(j,:),'FaceColor',[.8 .8 .8],'FaceAlpha',0.2,'DisplayName',Class)
    
    Dist = DistTypes{k};
    pd = GetFit(Max.(Class).(BlockM).Max,BlockM,[Dist],0,[1 BETATarget]);
    plot(x,pdf(pd.(pd.Best).pd,x)/ScaleDown(j),'k-','DisplayName',[Dist 'Fit'])
    
    % Set Plot Details
    a = ylim;
    ylim([a(1) a(2)*j])
    box on
    set(gca,'ytick',[],'yticklabel',[])
    if j == 1
        ylabel('Normalized Histograms (NTS)')
    end
    xlabel('Total Strip Load (kN)')
    title(BlockM)
    legend('location','best')
    xlim([LimitL+25 LimitR-25])
    
end

end

% PROBABILITY PAPER PLOT

if FitPlots
    figure('Position',[Ri Up 1200 400],'Name',[Country ' Q1Q2' ' Probability Paper'],'NumberTitle','off'); Ri = Ri+25; Up = Ri+25;
end
   
for i = 1:length(ClassType)
    Class = ClassType{i};
    
    for j = 1:length(BM)
        BlockM = BM{j};
        
        Data = Max.(Class).(BlockM).Max;
        Prop = 0.95;
        mdl = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear','Weights',[0.1*ones(round(length(Data)*Prop),1);1*ones(length(Data)-round(length(Data)*(Prop)),1)]);
        %mdl = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear');
        Fit(i,j) = mdl.Rsquared.Ordinary*100;
        
        if FitPlots
            subplot(1,3,i)
            hold on
            scatter(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),7,C(j,:),'filled','DisplayName','Max Data');
            plot(norminv((1:length(Data))/(length(Data) + 1)),mdl.Fitted,'-','Color',C(j,:),'DisplayName',['Fitted ' num2str(mdl.Rsquared.Ordinary,3)]);
        
            box on
            xlim([-4 4]); ylim([3.5 7]);  y1 = ylim;
            text(1,y1(1)+(y1(2)-y1(1))*(j/6),sprintf('R^2:  %.1f%%',mdl.Rsquared.Ordinary*100),"Color",C(j,:))
            if i == 1
                ylabel('log[ Total Strip Load (kN) ]')
            end
            xlabel('Standard Normal Percentile')
            title([ClassT{i}])
        end
    end
end

% TABLE GENERATION

Summary = array2table(zeros(4,length(BM)*length(ClassType)));
Summary.Properties.VariableNames = {[BM{3} ' ' ClassType{1}], [BM{3} ' ' ClassType{2}], [BM{3} ' ' ClassType{3}], [BM{2} ' ' ClassType{1}],...
    [BM{2} ' ' ClassType{2}], [BM{2} ' ' ClassType{3}], [BM{1} ' ' ClassType{1}], [BM{1} ' ' ClassType{2}], [BM{1} ' ' ClassType{3}]};
Summary.Properties.RowNames = {'Ed', 'LNFit (%)', 'AlphaQLN', 'AlphaQ2LN'};
AQ1 = 0.55;

for i = 1:length(ClassType)
    Class = ClassType{i};

    for j = 1:length(BM)
        BlockM = BM{j};
        
        if j == 2
            pd = GetFit(Max.(Class).(BlockM).Max,BlockM,'All',1,[1 BETATarget]);
            set(gcf,'Name',[get(gcf,'name') ' ' sprintf('%s %s %s \n',Class,BlockM,pd.Best)],'NumberTitle','off')
        end
        Data = Max.(Class).(BlockM).Max;
        pd = GetFit(Max.(Class).(BlockM).Max,BlockM,[Dist],0,[1 BETATarget]);
        %EdLN = pd.([Dist]).Ed;
        EdLN = pd.(pd.Best).Ed;
        AQLN = EdLN/(1000*1.5);
        AQ2LN = (EdLN/(1.5)-AQ1*600)/400;
        
        Summary.([BlockM ' ' Class]) = [EdLN; Fit(i,j); AQLN; AQ2LN];
        
        if strcmp(BlockM,'Weekly')
            AQ(i) = AQLN; Ed(i) = EdLN; AQ2(i) = AQ2LN; 
        end
    end
end

if DispTab
     format bank
%     fprintf('\n')
%     disp(Summary)

    figure('Position',[Ri Up 1300 125],'Name',[Country ' Q1Q2' ' Summary'],'NumberTitle','off'); Ri = Ri+25; Up = Ri+25;
    
    % Get the table in string form.
    TString = evalc('disp(Summary)');
    % Use TeX Markup for bold formatting and underscores.
    TString = strrep(TString,'<strong>','\bf');
    TString = strrep(TString,'</strong>','\rm');
    TString = strrep(TString,'_','\_');
    % Get a fixed-width font.
    FixedWidth = get(0,'FixedWidthFontName');
    % Output the table using the annotation command.
    annotation(gcf,'Textbox','String',TString,'Interpreter','Tex',...
        'FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);

end

if AdvPlots
    
[TrTyps,TN] = VBTypes2Names;

% % What vehicles were involved?
% TrTyps = [11; 12; 22; 23; 111; 11117; 1127; 12117; 122; 11127; 1128; 1138; 1238; 41; 42; 43; 44; 45; 46; 47; 48; 49; 0; 99];
% TN = ["11" "12" "22" "23" "111" "1111r" "112r" "1211r" "122" "11112r" "112a" "113a" "123a"...
%     "60t Crane" "6ax 60t" "7ax 72t" "8ax 84t" "9ax 96t" "96t Crane" "Libherr 132" "Libherr 15" "84t AT7" "No Class" "Empty"]';

% For each ClassType (or in this case just All and Class)
for i = 1:length(ClassType)
    Class = ClassType{i};
    
    % Select Daily, Weekly, or Yearly
    for j = 1:length(BM)
        BlockM = BM{j};
        
        figure('Position',[0 0 1500 400]);
        hold on
        
        % Get Primary Lane and Accompanying Lane Vehicles
        Primary = zeros(height(Max.(Class).(BlockM)),1); Accomp = Primary;
        [~, Ind] = max([Max.(Class).(BlockM).L1Load'; Max.(Class).(BlockM).L2Load']); Ind = Ind';
        Primary(Ind == 1) = Max.(Class).(BlockM).L1Veh(Ind == 1); Primary(Ind == 2) = Max.(Class).(BlockM).L2Veh(Ind == 2);
        Accomp(Ind == 1) = Max.(Class).(BlockM).L2Veh(Ind == 1); Accomp(Ind == 2) = Max.(Class).(BlockM).L1Veh(Ind == 2);
        
        PrimaryC = categorical(Primary); AccompC = categorical(Accomp);
        [N, ~] = histcounts(PrimaryC,categorical(TrTyps)); [Nx] = histcounts(AccompC,categorical(TrTyps));
        
        Label = cell(length(TN),1); Labelx = cell(length(TN),1);
        Perc = 100*N/sum(N); Percx = 100*Nx/sum(Nx);
        PercS = string(round(Perc,1))'; PercSx = string(round(Percx,1))';
        
        for p = 1:length(TN)
            if Perc(p) > 4
                Label{p} = [TN{p} ' - ' PercS{p} '%'];  else;  Label{p} = '';
            end
            if Percx(p) > 4
                Labelx{p} = [TN{p} ' - ' PercSx{p} '%'];  else;  Labelx{p} = '';
            end
        end
        
        % We need colors to be the same for each TrTyp
        % Labels should be strings that include the percentages
        
        subplot(1,2,1)
        sgtitle([Class ' ' BlockM],'fontweight','bold','fontsize',12);
        h = pie(N,Label);
        newColors = linspecer(numel(N));
        % Isolate the patch handles
        patchHand = findobj(h, 'Type', 'Patch');
        % Set the color of all patches using the nx3 newColors matrix
        set(patchHand, {'FaceColor'}, mat2cell(newColors, ones(size(newColors,1),1), 3))
        title('Primary Lane')
        subplot(1,2,2)
        h = pie(Nx,Labelx);
        % Isolate the patch handles
        patchHand = findobj(h, 'Type', 'Patch');
        set(patchHand, {'FaceColor'}, mat2cell(newColors, ones(size(newColors,1),1), 3))
        title('Accompanying Lane')
        
        if strcmp(Class,'ClassOW') && strcmp(BlockM,'Daily')
            NAD = N + Nx; NADnorm = NAD/sum(NAD);
        elseif  strcmp(Class,'ClassOW') && strcmp(BlockM,'Yearly')
            NAY = N + Nx; NAYnorm = NAY/sum(NAY);
        elseif strcmp(Class,'ClassOW') && strcmp(BlockM,'Weekly')
            NAW = N + Nx; NAWnorm = NAW/sum(NAW);
        end
    end
end

% OPTIONAL (MUST BE FIXED)

% % Lets do the comparison between AllDaily and AllYearly
% figure('Position',[0 0 1200 400]);
% b = bar(categorical(TN),NAW);
% b.FaceColor = 'flat';
% b.CData(:,:) = newColors;
% ylabel('Number #')
% xlabel('Truck Type')
% title('Total Involved in Weekly ClassOW')
% 
% figure('Position',[0 0 1200 400]);
% b = bar(categorical(TN),100*NAYnorm./NADnorm);
% b.FaceColor = 'flat';
% b.CData(:,:) = newColors;
% ylabel('Prevalance Increase (%)')
% xlabel('Truck Type')
% title('Prevalance Increase (%) From Daily to Yearly')

% I should add SPEED
figure('Position',[Ri Up 1300 400],'Name',[Country ' Q1Q2' ' Speeds'],'NumberTitle','off'); Ri = Ri+25; Up = Ri+25;

% X Stuff
Step = 1;
LimitL = 0;
LimitR = 150;
X = LimitL:Step:LimitR;
%Xp = LimitL-Step/2:Step:LimitR-Step/2;
x = X(1:end-1) + diff(X);

for i = 1:length(ClassType)
    Class = ClassType{i};
    
    subplot(1,3,i)
    hold on
    
    for j = 1:length(BM)
        BlockM = BM{j};
        
        y = histcounts(Max.(Class).(BlockM).L1Sp(Max.(Class).(BlockM).L1Sp>0),'BinEdges',X,'normalization','pdf');
        bar(x,y/ScaleDown(j),1,'EdgeColor',C(j,:),'FaceColor',[.8 .8 .8],'FaceAlpha',0.5,'DisplayName',BlockM)
        
    end
    
    % Set Plot Details
    set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
    box on
    xlim([LimitL+.25 LimitR-.25])
    if i == 1
        ylabel('Normalized Histograms (NTS)')
    end
    xlabel('Speed (kph)')
    title([ClassT{i}])   
    legend('location','northwest')
    
end

end


end
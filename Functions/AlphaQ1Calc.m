function [AQ1,Summary] = AlphaQ1Calc(Country,Type,FitPlots,AdvPlots,DispTab)
% ALPHAQCALC
% Calculates AlphaQ1 using:
% Country (or BMAxles Group Name) such as "SWISS", "GER", or "AUST"
% FitPlots/AdvPlots/DispTab Booleans 
% Returns AQ1 (value for Weekly BlockMax, with 1) All 2) Class+ 3) Class

% Make work for Single, Tridem...

% Load
load('Misc/BMAxles.mat')
Max = BMAxles.(Country).MaxAx.(Type);
% Set Divisor
if strcmp(Type,"Single")
    Div = 1;
elseif strcmp(Type,"Tridem")
    Div = 3;
else
    Div = 2;
end

% CURVE FITTING

% Block Maxima (always use j)
BM = {'Daily', 'Weekly', 'Yearly'};

% Select which Classification you want to analyze (always use i)
ClassType = {'All', 'ClassOW', 'Class'};
% Must be in the order 'All' 'ClassOW' 'Class' due to deletions
ClassT = {'All', 'Classified+', 'Classified'};

% Set Distribution Types
DistTypes = {'Normal', 'Lognormal'};

% Set x values to be use globally
x_values = 0:1:800;

% Fit Block Maxima to Normal Curve
for i = 1:length(ClassType)
    Class = ClassType{i};
    for j = 1:length(BM)
        BlockM = BM{j};
        for k = 1:length(DistTypes)
            Dist = DistTypes{k};
            Data = Max.(Class).(BlockM).Max/Div;
            [pd.(Class).(BlockM).(Dist),~,y_values.(Class).(BlockM).(Dist).PDF_Fit,y_values.(Class).(BlockM).(Dist).CDF_Fit] = GetBlockMaxFit(Data,Dist,0,0,x_values);
        end
    end
end


% Set colours
C = linspecer(9);
% ScaleDown
ScaleDown = [1 2.5 5];
% X Stuff
Step = 2.5;
LimitL = 10; LimitR = 300;
X = LimitL:Step:LimitR;
x = X(1:end-1) + diff(X);

% FIGURE OF ALL DISTS
if FitPlots
    
figure('Position',[0 0 1500 400]);
    
for i = 1:length(ClassType)
    Class = ClassType{i};
    
    subplot(1,3,i)
    hold on
    
    for j = 1:length(BM)
        BlockM = BM{j};
        
        y = histcounts(Max.(Class).(BlockM).Max/Div,'BinEdges',X,'normalization','pdf');
        bar(x,y/ScaleDown(j),1,'EdgeColor',C(j,:),'FaceColor',[.8 .8 .8],'FaceAlpha',0.5,'DisplayName',BlockM)
        
        for k = 2
            Dist = DistTypes{k};
            plot(x_values,y_values.(Class).(BlockM).(Dist).PDF_Fit/ScaleDown(j),'k-','DisplayName',[Dist 'Fit'])
        end
    end
    
    % Set Plot Details
    set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
    box on
    xlim([LimitL+25 LimitR-25])
    if i == 1
        ylabel('Normalized Histograms (NTS)')
    end
    xlabel(['Total ' Type ' Load (kN) /' num2str(Div)])
    title([ClassT{i}])
    legend('location','northeast')
end

% FIGURE OF ALL DISTS SEPARATED
figure('Position',[0 0 1500 400]);

% X Stuff
Step = 2.5;
LimitL = 10; LimitR = 300;

for j = 1:length(BM)
    BlockM = BM{j};
    
    subplot(1,3,j)
    hold on
    k = 2;
    
    Class = 'Class';
    y = histcounts(Max.(Class).(BlockM).Max/Div,'BinEdges',X,'normalization','pdf');
    bar(x,y/ScaleDown(j),1,'EdgeColor','k','FaceColor',[.8 .8 .8],'FaceAlpha',0.8,'DisplayName',Class)
    
    Dist = DistTypes{k};
    plot(x_values,y_values.(Class).(BlockM).(Dist).PDF_Fit/ScaleDown(j),'k-','DisplayName',[Dist 'Fit'])    
    
    Class = 'All';
    y = histcounts(Max.(Class).(BlockM).Max/Div,'BinEdges',X,'normalization','pdf');
    bar(x,y/ScaleDown(j),1,'EdgeColor',C(j,:),'FaceColor',[.8 .8 .8],'FaceAlpha',0.2,'DisplayName',Class)
    
    Dist = DistTypes{k};
    plot(x_values,y_values.(Class).(BlockM).(Dist).PDF_Fit/ScaleDown(j),'k-','DisplayName',[Dist 'Fit'])

    % Set Plot Details
    a = ylim;
    ylim([a(1) a(2)*j])
    box on
    set(gca,'ytick',[],'yticklabel',[])
    if j == 1
        ylabel('Normalized Histograms (NTS)')
    end
    xlabel(['Total ' Type ' Load (kN) /' num2str(Div)])
    title(['Maximum ' Type ' Axles ' BlockM])
    legend('location','northeast')
    xlim([LimitL+25 LimitR-25])
    
end

end

% PROBABILITY PAPER PLOT

if FitPlots
    figure('Position',[0 0 1500 400]);
end

for i = 1:length(ClassType)
    Class = ClassType{i};
    
    for j = 1:length(BM)
        BlockM = BM{j};
        
        
        Data = Max.(Class).(BlockM).Max/Div;
        mdl = fitlm(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),'linear');
        Fit(i,j) = mdl.Rsquared.Ordinary*100;
        
        if FitPlots
            subplot(1,3,i)
            hold on
            scatter(norminv((1:length(Data))/(length(Data) + 1)),log(sort(Data)),7,C(j,:),'filled','DisplayName','Max Data');
            plot(norminv((1:length(Data))/(length(Data) + 1)),mdl.Fitted,'-','Color',C(j,:),'DisplayName',['Fitted ' num2str(mdl.Rsquared.Ordinary,3)]);
            box on
            xlim([-4 4]); ylim([1.5 6])
            y1 = ylim;
            text(1,y1(1)+(y1(2)-y1(1))*(j/6),sprintf('R^2:  %.1f%%',mdl.Rsquared.Ordinary*100),"Color",C(j,:))
            
            if i == 1
                ylabel(['log [ Total ' Type ' Load (kN) /2 ]'])
            end
            xlabel('Standard Normal Percentile')
            title([ClassT{i}])
        end
    end
end

% Give one more figure which shoes Weekly Class with detailed annotations
if FitPlots


end

% TABLE GENERATION

Summary = array2table(zeros(4,length(BM)*length(ClassType)));
Summary.Properties.VariableNames = {[BM{3} ' ' ClassType{1}], [BM{3} ' ' ClassType{2}], [BM{3} ' ' ClassType{3}], [BM{2} ' ' ClassType{1}],...
    [BM{2} ' ' ClassType{2}], [BM{2} ' ' ClassType{3}], [BM{1} ' ' ClassType{1}], [BM{1} ' ' ClassType{2}], [BM{1} ' ' ClassType{3}]};
Summary.Properties.RowNames = {'Ed', 'LNFit (%)', 'AlphaQ1N', 'AlphaQ1LN'};

for i = 1:length(ClassType)
    Class = ClassType{i};

    for j = 1:length(BM)
        BlockM = BM{j};

        Data = Max.(Class).(BlockM).Max/Div;
        [~, AQN, ~] = GetBlockMaxEd(Data,BlockM,'Normal',300*1.5,[1 1],1,1,1,1,1);
        [EdLN, AQLN, ~] = GetBlockMaxEd(Data,BlockM,'Lognormal',300*1.5,[1 1],1,1,1,1,1);
        
        Summary.([BlockM ' ' Class]) = [EdLN; Fit(i,j); AQN; AQLN];
        
        if strcmp(BlockM,'Weekly')
            AQ1(i) = AQLN;
        end
        
    end
end

if DispTab
    format bank
    fprintf('\n')
    disp(Summary)
end

if AdvPlots
    
    % What vehicles were involved?
    [TrTyps,TN] = VBTypes2Names;
    
    % For each ClassType (or in this case just All and Class)
    for i = 1:length(ClassType)
        Class = ClassType{i};
        
        figure('Position',[0 0 1500 400]);
        hold on
        
        % Select Daily, Weekly, or Yearly
        for j = 1:length(BM)
            BlockM = BM{j};
            
            [N, Cats] = histcounts(categorical(Max.(Class).(BlockM).CLASS),categorical(TrTyps));
            
            Label = cell(length(TN),1);
            Perc = 100*N/sum(N);
            PercS = string(round(Perc,1))';
            
            for p = 1:length(TN)
                if Perc(p) > 2
                    Label{p} = [TN{p} ' - ' PercS{p} '%'];  else;  Label{p} = '';
                end
            end
            
            % We need colors to be the same for each TrTyp
            % Labels should be strings that include the percentages
            subplot(1,3,j)
            title([ClassT{i} ' ' BlockM],'fontweight','bold','fontsize',10);
            h = pie(N,Label);
            newColors = linspecer(numel(N));
            % Isolate the patch handles
            patchHand = findobj(h, 'Type', 'Patch');
            % Set the color of all patches using the nx3 newColors matrix
            set(patchHand, {'FaceColor'}, mat2cell(newColors, ones(size(newColors,1),1), 3))
            title([ClassT{i} ' ' BlockM])
            
        end
    end

    % What were the characteristics of those vehicles?
    % Axle Spacing
    
    if strcmp(Type,"Tandem")
    figure('Position',[0 0 1500 400]);
    
    % X Stuff
    Step = 0.025;
    LimitL = Step/2; LimitR = 3.15;
    X = LimitL:Step:LimitR;
    Xp = LimitL-Step/2:Step:LimitR-Step/2;
    x = X(1:end-1) + diff(X);
    
    for i = 1:length(ClassType)
        Class = ClassType{i};
        
        subplot(1,3,i)
        hold on
        
        for j = 1:length(BM)
            BlockM = BM{j};
            
            y = histcounts(Max.(Class).(BlockM).W1_2M,'BinEdges',X,'normalization','pdf');
%             if i == 2 && strcmp(BlockM,"Weekly")
%                 ystar = Max.(Class).(BlockM).W1_2M;
%             end
            bar(x,y/ScaleDown(j),1,'EdgeColor',C(j,:),'FaceColor',[.8 .8 .8],'FaceAlpha',0.5,'DisplayName',BlockM)
            
        end
        
        % Set Plot Details
        set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
        box on
        xlim([LimitL+.25 LimitR-.25])
        if i == 1
            ylabel('Normalized Histograms (NTS)')
        end
        xlabel('m')
        title([ClassT{i}])
        legend('location','best')
        %xticks(0.4:0.2:2.6)
    end
    end
    
    
end

end
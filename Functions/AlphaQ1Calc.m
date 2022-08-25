function [AQ1,Summary] = AlphaQ1Calc(Country,Type,FitPlots,AdvPlots,DispTab)
% ALPHAQCALC
% Calculates AlphaQ1 using:
% Country (or BMAxles Group Name) such as "CH", "DE", or "AT"
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

% Initialize Fig Positions
Ri = 0; Up = 0;

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
    
figure('Position',[Ri Up 1200 350],'Name',[Country ' ' Type ' Fits'],'NumberTitle','off'); Ri = Ri+25; Up = Ri+25;
    
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
            pd = GetFit(Max.(Class).(BlockM).Max/Div,BlockM,Dist,0,1);
            plot(x_values,pdf(pd.(Dist).pd,x_values)/ScaleDown(j),'k-','DisplayName',[Dist 'Fit'])
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
figure('Position',[Ri Up 1200 350],'Name',[Country ' ' Type ' Fits Separated'],'NumberTitle','off'); Ri = Ri+25; Up = Ri+25;

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
    pd = GetFit(Max.(Class).(BlockM).Max/Div,BlockM,Dist,0,1);
    plot(x_values,pdf(pd.(Dist).pd,x_values)/ScaleDown(j),'k-','DisplayName',[Dist 'Fit'])    
    
    Class = 'All';
    y = histcounts(Max.(Class).(BlockM).Max/Div,'BinEdges',X,'normalization','pdf');
    bar(x,y/ScaleDown(j),1,'EdgeColor',C(j,:),'FaceColor',[.8 .8 .8],'FaceAlpha',0.2,'DisplayName',Class)
    
    Dist = DistTypes{k};
    pd = GetFit(Max.(Class).(BlockM).Max/Div,BlockM,Dist,0,1);
    plot(x_values,pdf(pd.(Dist).pd,x_values)/ScaleDown(j),'k-','DisplayName',[Dist 'Fit'])

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
    figure('Position',[Ri Up 1200 350],'Name',[Country ' ' Type ' Probability Paper'],'NumberTitle','off'); Ri = Ri+25; Up = Ri+25;
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



% TABLE GENERATION

Summary = array2table(zeros(4,length(BM)*length(ClassType)));
Summary.Properties.VariableNames = {[BM{3} ' ' ClassType{1}], [BM{3} ' ' ClassType{2}], [BM{3} ' ' ClassType{3}], [BM{2} ' ' ClassType{1}],...
    [BM{2} ' ' ClassType{2}], [BM{2} ' ' ClassType{3}], [BM{1} ' ' ClassType{1}], [BM{1} ' ' ClassType{2}], [BM{1} ' ' ClassType{3}]};
Summary.Properties.RowNames = {'Ed', 'LNFit (%)', 'AlphaQ1Best', 'AlphaQ1LN'};

for i = 1:length(ClassType)
    Class = ClassType{i};

    for j = 1:length(BM)
        BlockM = BM{j};
        
        if j == 2
            pd = GetFit(Max.(Class).(BlockM).Max/Div,BlockM,'All',1,1);
            set(gcf,'Name',[get(gcf,'name') ' ' sprintf('%s %s %s \n',Class,BlockM,pd.Best)],'NumberTitle','off')
        end
        pd = GetFit(Max.(Class).(BlockM).Max/Div,BlockM,'All',0,1);
        EdBest = pd.(pd.Best).Ed;
        %fprintf('%s %s %s \n',Class,BlockM,pd.Best)
        EdLN = pd.Lognormal.Ed;
        AQBest = EdBest/(300*1.5);
        AQLN = EdLN/(300*1.5);
        
        Summary.([BlockM ' ' Class]) = [EdLN; Fit(i,j); AQBest; AQLN];
        
        if strcmp(BlockM,'Weekly')
            AQ1(i) = AQLN; Ed(i) = EdLN;
        end
    end
end

% Give one more figure which shoes Weekly Class with detailed annotations

if FitPlots
    Ri = 300; Up = 200;
    figure('Position',[Ri Up 600 450],'Name',[Country ' ' Type ' Extreme Value'],'NumberTitle','off'); Ri = Ri+50; Up = Ri+50; 
    j = 2; i = 3; k = 2;
    BlockM = BM{2};
    Class = ClassType{1};
    Dist = DistTypes{k}; Data = Max.(Class).(BlockM).Max/Div;
    pd = GetFit(Max.(Class).(BlockM).Max/Div,BlockM,Dist,0,1);
    
    y = histcounts(Max.(Class).(BlockM).Max/Div,'BinEdges',X,'normalization','pdf');
    bar(x,y/ScaleDown(j),1,'EdgeColor','k','FaceColor',[.8 .8 .8],'FaceAlpha',0.8,'DisplayName',Class)
    hold on
    plot(x_values,pdf(pd.(Dist).pd,x_values)/ScaleDown(j),'k-','DisplayName',[Dist 'Fit'])
    

    % Set Plot Details
    a = ylim;
    ylim([a(1) a(2)*1.5]); a = ylim;
    box on
    set(gca,'ytick',[],'yticklabel',[])
    ylabel('Normalized Histogram (NTS)')
    xlabel(['Total ' Type ' Load (kN) /' num2str(Div)])
    title(['Maximum ' Type ' Axles ' BlockM])
    legend('location','northeast')
    xlim([LimitL+75 LimitR-70])
    
    % Make two lines, one at Em and one at Ed > show that with 0.7BCOV
    line([mean(Data),mean(Data)],[0,a(2)*0.7],'Color','b','LineStyle','-','HandleVisibility','off')
    text(mean(Data),a(2)*0.83,'Mean Value','HorizontalAlignment','center')
    text(mean(Data),a(2)*0.77,'E_m','HorizontalAlignment','center')
    line([Ed(i),Ed(i)],[0,a(2)*0.35],'Color','b','LineStyle','-','HandleVisibility','off')
    text(Ed(i),a(2)*0.48,'Design Value','HorizontalAlignment','center')
    text(Ed(i),a(2)*0.42,'E_d','HorizontalAlignment','center')
    line([mean(Data),Ed(i)],[a(2)*0.30,a(2)*0.30],'Color','r','LineStyle','--','HandleVisibility','off')
    line([Ed(i)-3,Ed(i)],[a(2)*0.31,a(2)*0.30],'Color','r','LineStyle','--','HandleVisibility','off')
    line([Ed(i)-3,Ed(i)],[a(2)*0.29,a(2)*0.30],'Color','r','LineStyle','--','HandleVisibility','off')
    %annotation('textarrow',[(mean(Data)-LimitL)/(LimitR-70-LimitL-75),(Ed(i)-LimitL-150)/(LimitR-70-LimitL-75)],[0.30, 0.30],'String',' Increase ','FontSize',13,'Linewidth',2)
    text(mean([mean(Data),Ed(i)]),a(2)*0.34,'{\it increase}','HorizontalAlignment','center')
    text(mean([mean(Data),Ed(i)]),a(2)*0.26,'\alpha\beta\nu','HorizontalAlignment','center','Interpreter','tex')
    
end


if DispTab
     format bank
%     fprintf('\n')
%     disp(Summary)
    
    figure('Position',[Ri Up 1300 125],'Name',[Country ' ' Type ' Summary'],'NumberTitle','off'); Ri = Ri+25; Up = Ri+25;
    
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
    
    % What vehicles were involved?
    %[TrTyps,TN] = VBTypes2Names;
    VType = VBVTypes;
    
    % For each ClassType (or in this case just All and Class)
    for i = 1:length(ClassType)
        Class = ClassType{i};
        
        figure('Position',[Ri Up 1200 400],'Name',[Country ' ' Type ' ' Class ' Responsible'],'NumberTitle','off'); Ri = Ri+50; Up = Ri+50;
        hold on
        
        % Select Daily, Weekly, or Yearly
        for j = 1:length(BM)
            BlockM = BM{j};
            
            N = histcounts(categorical(Max.(Class).(BlockM).CLASS),categorical(VType.CLASS));
            
            Label = cell(height(VType),1);
            Perc = 100*N/sum(N);
            PercS = string(round(Perc,1))';
            
            for p = 1:length(VType.Name)
                if Perc(p) > 2
                    Label{p} = strcat(VType.Name(p), ' - ', PercS{p}, '%');  else;  Label{p} = '';
                end
            end
            
            % We need colors to be the same for each TrTyp
            % Labels should be strings that include the percentages
            subplot(1,3,j)
            title([ClassT{i} ' ' BlockM],'fontweight','bold','fontsize',10);
            h = pie(N,Label);
            newColors = linspecer(numel(N));
            newColors(1,:) = [1 1 1];
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
    figure('Position',[Ri Up 1200 350],'Name',[Country ' ' Type ' Axle Spacing'],'NumberTitle','off'); Ri = Ri+50; Up = Ri+50;
    
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
            
            y = histcounts(Max.(Class).(BlockM).W1_2,'BinEdges',X,'normalization','pdf');
%             if i == 2 && strcmp(BlockM,"Weekly")
%                 ystar = Max.(Class).(BlockM).W1_2;
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
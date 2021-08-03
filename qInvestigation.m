% qInvestigation
% Find alpha small q, given alpha Q1 Q2. From Q1 and Q1/Q2 Investigation we get aQ1 = 0.60 and aQ2 = 0.50.
% Produced using vars in qInvestigation folder (Misc) which comes from VBWIMq.
% This is part of an initial effort to estimate alpha q using OBSERVED WIM.

% Initial commands
clear, clc, close all

% We've got problems... this needs customization for the below options...

% Analysis Type Option
%AType = 'BoxTwin';    % - Box and Standard Twin
%AType = 'SlabShort';
AType = 'SlabLong';

% I want to process all at the same time... how can i do this?

% Which IL to Highlight?
HL = 5;

% Possible to load variables and skip first section
Initial = false;

% Perform initial analysis if necessary (necessary vars loaded either way)

% Inputs
BM = {'Daily', 'Weekly', 'Yearly'};             % j
ClassType = {'All', 'ClassOW', 'Class'};        % i
ClassT = {'All', 'Classified+', 'Classified'};
DistTypes = {'Normal', 'Lognormal'};%, 'GeneralizedExtremeValue'};
    
if Initial
    % --- Step 1: Input desired filters/parameters ---
    % Select which year to analyze %Years = 2011:2019;
    Years = 0;
    % Select stations (locations) to analyze %Stations = [408 409 405 406 402 415 416];
    Stations = 0;
    
    % Get Initial Vars from Steps 2 to 4
    [Max,pd,x_values,y_values,ILData,ESIA,BaseData] = qInvestInitial(AType,Years,Stations,BM,ClassType,ClassT,DistTypes);
    
else

    load(strcat(AType,'MaxEventsqMax'))
    load(strcat(AType,'ESIAq'))
    load(strcat(AType,'ILDataq'))
    load(strcat(AType,'BaseDataq'))
    load(strcat(AType,'MaxEventsqpd.mat'))
    load(strcat(AType,'MaxEventsqy_values.mat'))
    load(strcat(AType,'MaxEventsqx_values.mat'))

end

HLName = ILData(HL).Name;

% --- Step 5: Plot Block Maxima ---

% Set colours
C = linspecer(9);
% ScaleDown
ScaleDown = [1 2.5 5];

for r = [HL];
     
    figure('Position',[300 300 1200 350]);
    xLims = [0 max(x_values(r).All.Yearly)+100*max(x_values(r).All.Yearly)/1500];
    
    for i = 1:length(ClassType)
        Class = ClassType{i};

        subplot(1,3,i)
        hold on

        for j = 1:length(BM)
            BlockM = BM{j};
            
            % X is for the plot
            X = x_values(r).(Class).(BlockM);
            % x is for the bar
            x = X(1:end-1) + diff(X);

            y = histcounts(Max(r).(Class).(BlockM).Max,'BinEdges',X,'normalization','pdf');
            bar(x,y/ScaleDown(j),1,'EdgeColor',C(j,:),'FaceColor',[.8 .8 .8],'FaceAlpha',0.5,'DisplayName',BlockM)

            for k = 2
                Dist = DistTypes{k};
                plot(X,y_values(r).(Class).(BlockM).(Dist).PDF_Fit/ScaleDown(j),'k-','DisplayName',[Dist 'Fit'])

                if strcmp(BlockM,'Yearly') & (~strcmp(Dist,'ExtremeValue') | ~strcmp(Dist,'GeneralizedExtremeValue'))
                    plot(X,pdf(pd(r).(Class).D2Y.(Dist),X)/ScaleDown(j),'k--','DisplayName',[Dist 'D2Y'])
                    plot(X,pdf(pd(r).(Class).W2Y.(Dist),X)/ScaleDown(j),'k:','DisplayName',[Dist 'W2Y'])
                end
            end             
        end

        % Set Plot Details
        set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
        box on
        xlim(xLims)
        if i == 1
            ylabel('Normalized Histograms (NTS)')
        end
        xlabel('Bridge Action Effect (kNm)')
        title([ClassT{i}])
        legend('location','northeast')

    end
end


% --- Step 6: Probability Paper Plotting & Tail Fitting Methods ---
% When we apply the above equations, we are assuming the maximums are normally or log-normally distributed. 
% We should evaluate the fit with a probability paper plot.
PlotB = zeros(1,length(ILData)); PlotB(HL) = 1;

% For each Influence Case
for r = 1:length(ILData);
    if PlotB(r) == 1
        figure('Position',[300 300 1200 350]);
    end
    
    for i = 1:length(ClassType)
        Class = ClassType{i};
        
        for j = 1:length(BM)
            BlockM = BM{j};

            [MaxECDF, MaxECDFRank] = ecdf(Max(r).(Class).(BlockM).Max); MaxECDFRank = MaxECDFRank'; MaxECDF(1) = []; MaxECDFRank(1) = [];
            mdl = fitlm(norminv((1:length(MaxECDFRank))/(length(MaxECDFRank) + 1)),log(MaxECDFRank),'linear');
            
            x2 = norminv((1:length(MaxECDFRank))/(length(MaxECDFRank) + 1));
            y2 = log(MaxECDFRank);
            
            md2 = fitlm(x2(round(length(x2)*3/4):length(x2)),y2(round(length(x2)*3/4):length(x2)),'linear');
            
            if PlotB(r) == 1
                subplot(1,3,i)
                hold on
                scatter(norminv((1:length(MaxECDFRank))/(length(MaxECDFRank) + 1)),log(MaxECDFRank),7,C(j,:),'filled','DisplayName','Max Data');
                plot(norminv((1:length(MaxECDFRank))/(length(MaxECDFRank) + 1)),mdl.Fitted,'-','Color',C(j,:),'DisplayName',['Fitted ' num2str(mdl.Rsquared.Ordinary,3)]);
                box on
                y1 = ylim;
                text(1,y1(1)+(y1(2)-y1(1))*(j/6),sprintf('R^2:  %.1f%%',mdl.Rsquared.Ordinary*100),"Color",C(j,:))
            end

            Fit(i,j,r) = mdl.Rsquared.Ordinary*100;
            
        end
        
        if PlotB(r) == 1
            if i == 1
                ylabel('log[ Total Load Effect (kNm) ]')
            end
            xlabel('Standard Normal Percentile')
            %legend('Location', 'southeast')
            title([ClassT{i}])
        end
    end
end


% --- Step 7 : Summarize  Estimates ---
% The alpha we are finding and putting in the table right now is some kind of composite alpha...
% We need to add AQ1 and AQ2 to get the real Aq...
% According to Annex C of SIA 269
Beta.Yearly = 4.7; % Here we use the Beta annual - so we should use annual max effects
Beta.Weekly = 5.44422; % See Tail Fitting > Beta Conversion
Beta.Daily = 5.72397;
Alpha = 0.7;
AQ1 = 0.6;
AQ2 = 0.5;

% For Plot -- Modify!
ClassP = 'Class';
BlockMP = 'Weekly';

% For each Influence Case
for r = 1:length(ILData)
    
    % Initialize SumTable
    ST(r).SumTableLogNorm = array2table(zeros(5,length(BM)*length(ClassType)));
    ST(r).SumTableLogNorm.Properties.VariableNames = {[BM{3} ' ' ClassType{1}], [BM{3} ' ' ClassType{2}], [BM{3} ' ' ClassType{3}], [BM{2} ' ' ClassType{1}],...
        [BM{2} ' ' ClassType{2}], [BM{2} ' ' ClassType{3}], [BM{1} ' ' ClassType{1}], [BM{1} ' ' ClassType{2}], [BM{1} ' ' ClassType{3}]};
    ST(r).SumTableLogNorm.Properties.RowNames = {'Em', 'Ed', 'LNFit (%)', 'AlphaQLN', 'AlphaqLN'};

    ST(r).SumTableLogNormComp = array2table(zeros(2,(length(BM)-1)*length(ClassType)));
    ST(r).SumTableLogNormComp.Properties.RowNames = {'AlphaQLN', 'AlphaqLN'};
    ST(r).SumTableLogNormComp.Properties.VariableNames = {['W2Y ' ClassType{1}],['W2Y ' ClassType{2}],['W2Y ' ClassType{3}],['D2Y ' ClassType{1}],['D2Y ' ClassType{2}],['D2Y ' ClassType{3}]};
    
    for i = 1:length(ClassType)
        Class = ClassType{i};
        
        for j = 1:length(BM)
            BlockM = BM{j};
            
            Em = mean(Max(r).(Class).(BlockM).Max);
            Stdev = std(Max(r).(Class).(BlockM).Max);
            COV = Stdev/Em;
            Delta2 = log(COV^2+1);
            N5 = prctile(Max(r).(Class).(BlockM).Max,95);
            Dist95N = norminv(0.95,pd(r).(Class).(BlockM).Normal.mu,pd(r).(Class).(BlockM).Normal.sigma);
            Dist95LN = logninv(0.95,pd(r).(Class).(BlockM).Lognormal.mu,pd(r).(Class).(BlockM).Lognormal.sigma);
            
            EdN = Em*(1+Alpha*Beta.(BlockM)*COV);
            AQN = EdN/(ESIA.Total(r));
            AqN = ((EdN/1.5)-AQ1*ESIA.EQ(1,r)-AQ2*ESIA.EQ(2,r))/(ESIA.Eq(r));
            
            EdLN = Em*exp(Alpha*Beta.(BlockM)*sqrt(Delta2)-0.5*Delta2);
            AQLN = EdLN/(ESIA.Total(r));
            AqLN = ((EdLN/1.5)-AQ1*ESIA.EQ(1,r)-AQ2*ESIA.EQ(2,r))/(sum(ESIA.Eq(:,r)));
            
            % Find %ile corresponding to Ed
            N5 = prctile(Max(r).(Class).(BlockM).Max,95);
            I95 = invprctile(Max(r).(Class).(BlockM).Max,EdLN);
            
            % ST stands for Summary Table
            ST(r).SumTableLogNorm.([BlockM ' ' Class]) = [Em; EdLN; Fit(i,j,r); AQLN; AqLN];
            
            if strcmp(BlockM,'Weekly')
                pdx = pd(r).(Class).W2Y.Lognormal;
                Delta2 = pdx.sigma^2;
                Em = exp(pdx.mu+0.5*Delta2);
                ST(r).SumTableLogNormComp.(['W2Y ' Class]) = [Em*exp(Alpha*Beta.Yearly*sqrt(Delta2)-0.5*Delta2)/(ESIA.Total(r)); (((Em*exp(Alpha*Beta.Yearly*sqrt(Delta2)-0.5*Delta2))/1.5)-AQ1*ESIA.EQ(1,r)-AQ2*ESIA.EQ(2,r))/(ESIA.Eq(r))];
            elseif strcmp(BlockM,'Daily')
                pdx = pd(r).(Class).D2Y.Lognormal;
                Delta2 = pdx.sigma^2;
                Em = exp(pdx.mu+0.5*Delta2);
                ST(r).SumTableLogNormComp.(['D2Y ' Class]) = [Em*exp(Alpha*Beta.Yearly*sqrt(Delta2)-0.5*Delta2)/(ESIA.Total(r)); (((Em*exp(Alpha*Beta.Yearly*sqrt(Delta2)-0.5*Delta2))/1.5)-AQ1*ESIA.EQ(1,r)-AQ2*ESIA.EQ(2,r))/(ESIA.Eq(r))];
            end
            
            PlotAq.(Class).(BlockM)(r) = ST(r).SumTableLogNorm{'AlphaqLN',[BlockM ' ' Class]};
            PlotAQ.(Class).(BlockM)(r) = ST(r).SumTableLogNorm{'AlphaQLN',[BlockM ' ' Class]};
        end
    end
end

format bank
fprintf('\n')
disp(ST(HL).SumTableLogNorm)



% Prepare plot data for VBTriPlot
% Get xdata
xdata = 2:2:10;
xdata = 10:5:30;
% Titles
%Title{1} = 'M+ p1'; Title{2} = 'M+ p2'; Title{3} = 'M+ p3'; FigNum = 0;
Title{1} = 'M- p1'; Title{2} = 'M- p2'; Title{3} = 'M- p3'; FigNum = 0;
% Set starting points for ydata
%SP = [1 11 21];
SP = [6 16 26];
% Weekly data only
j = 2; BlockM = BM{j};

% Prep 3 subplots
for m = 1:3
    % All ClassTypes
    for i = 1:length(ClassType)
        Class = ClassType{i};
        % Get ydata
        ydata{m}(i,:) = PlotAQ.(Class).(BlockM)(SP(m):SP(m)+4);
        % DisplayName
        PDets.DN{m,i} = strcat(Class,BlockM);
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDets.MEC{m}(i,:) = [0 0 0]; PDets.MFC{m}(i,:) = C(i,:); 
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum);




% % Collect all Alphaq Estimates for each IL and display
% figure('Position',[200 200 730 400]);
% Title{1} = 'V'; Title{2} = 'M+'; Title{3} = 'M-';
% for m = 1:3
%     subplot(1,3,m)
%     hold on
%     
%     for i = 1:length(ClassType)
%         Class = ClassType{i};
%         
%         for j = 2%1:length(BM)
%             BlockM = BM{j};
%             
%             if m == 1
%                 plot(10:10:80,PlotAQ.(Class).(BlockM)(25:32),'-s','Color','k','MarkerEdgeColor','k','MarkerFaceColor',C(j,:),'MarkerSize',5,'DisplayName',strcat(Class,BlockM))
%             elseif m == 2
%                 plot(10:10:80,PlotAQ.(Class).(BlockM)(33:40),'-s','Color','k','MarkerEdgeColor','k','MarkerFaceColor',C(j,:),'MarkerSize',5)
%             else
%                 plot(10:10:80,PlotAQ.(Class).(BlockM)(41:48),'-s','Color','k','MarkerEdgeColor','k','MarkerFaceColor',C(j,:),'MarkerSize',5)
%                 %plot([2 4 6 8 10 12 15 20 25 30],PlotAQ(7:end),'-s','Color','k','MarkerEdgeColor','k','MarkerFaceColor','none','MarkerSize',7,'DisplayName','AlphaQ')
%             end
%         end
%     end
%     %plot(5:5:30,PlotAq(((i-1)*6+1):((i)*6)),'-s','Color','k','MarkerEdgeColor','none','MarkerFaceColor','k','MarkerSize',7,'DisplayName','Alphaq')
%     %plot(5:5:30,PlotAQ(((i-1)*6+1):((i)*6)),'-s','Color','k','MarkerEdgeColor','k','MarkerFaceColor','none','MarkerSize',7,'DisplayName','AlphaQ')
%     ytickformat('%.2f'); yticks(0:0.1:1); xticks(10:10:80); set(gca,'TickDir','out'); set(gca,'YGrid','on')
%     xlabel('Span (m)')
%     a = gca; set(a,'box','off','color','none');
%     % create new, empty axes with box but without ticks
%     b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
%     % set original axes as active, and link axes in case of zooming
%     axes(a); linkaxes([a b]);
% 
%     % Set axis limits
%     ylim([0 1]); xlim([10 80])
%     % Y-axis only for the leftmost (first) plot
%     if m == 1
%         ylabel('E_{WIM}/E_{SIA}')
%         yh = get(gca,'ylabel'); % handle to the label object
%         p = get(yh,'position'); % get the current position property
%         p(1) = 0.9*p(1);          % double the distance,
%         % negative values put the label below the axis
%         set(yh,'position',p)   % set the new position
%         legend('location','northeast')
%     end
%     title([Title{m}])
% end

% --- Step 8: Populate WIM, a comparison with AGB and MAT Results Variables ---
% First choose which result we will use (try Lognormal Weekly and Yearly)

Section = 'SlabShort';

% Distinguish between p1, p2, and p3
% Borrow from Output2Struct
% Don't go past 6 entries
% Load WIM first

%UniqNames = unique(Infl.Names,'stable');
ColumnNames = {'EQ1','EQ2','Eq','YClass','YClassOW','YAll','WClass','WClassOW','WAll','E'};

% Initialize the final matrix/table
X = NaN(8,length(ColumnNames));
XT = array2table(X,'VariableNames',ColumnNames);
SpanDiv = 2; % Set to 5 for Slab


% Must add p1 p2 p3 for SlabShort WIM
% for r = 1:30
%     % Get InfName (call it Temp)
%     %Temp = UniqNames{r}; %Fixed 26.5.20
%     Temp = ILData(r).Name;
%     Temp2 = strfind(Temp,'S');
%     Temp3 = strfind(Temp,'.');
%     % Get Span from InfName
%     %Span = str2num(Temp(end-1:end));
%     Span = str2num(Temp(Temp2(end)+1:end));
%     % Get Action Effect from InfName
%     %AE = Temp(1:end-2);
%     AE = Temp(Temp3(end-1)+1:Temp3(end)-1);
%     
%     try
%         if istable(WIM.(Section).(AE))
%             WIM.(Section).(AE).YClass(Span/SpanDiv) = ST(r).SumTableLogNorm.('Yearly Class')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).YClassOW(Span/SpanDiv) = ST(r).SumTableLogNorm.('Yearly ClassOW')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).YAll(Span/SpanDiv) = ST(r).SumTableLogNorm.('Yearly All')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).WClass(Span/SpanDiv) = ST(r).SumTableLogNorm.('Weekly Class')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).WClassOW(Span/SpanDiv) = ST(r).SumTableLogNorm.('Weekly ClassOW')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).WAll(Span/SpanDiv) = ST(r).SumTableLogNorm.('Weekly All')(end-1)*(ESIA.Total(r));
%             
%             WIM.(Section).(AE).E(Span/SpanDiv) = ESIA.Total(r);
%             WIM.(Section).(AE).Eq(Span/SpanDiv) = ESIA.Eq(r);
%             WIM.(Section).(AE).EQ1(Span/SpanDiv) = ESIA.EQ(1,r);
%             WIM.(Section).(AE).EQ2(Span/SpanDiv) = ESIA.EQ(2,r);
%             
%         end
%     catch
%         WIM.(Section).(AE) = XT;
%         
%         WIM.(Section).(AE).YClass(Span/SpanDiv) = ST(r).SumTableLogNorm.('Yearly Class')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).YClassOW(Span/SpanDiv) = ST(r).SumTableLogNorm.('Yearly ClassOW')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).YAll(Span/SpanDiv) = ST(r).SumTableLogNorm.('Yearly All')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).WClass(Span/SpanDiv) = ST(r).SumTableLogNorm.('Weekly Class')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).WClassOW(Span/SpanDiv) = ST(r).SumTableLogNorm.('Weekly ClassOW')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).WAll(Span/SpanDiv) = ST(r).SumTableLogNorm.('Weekly All')(end-1)*(ESIA.Total(r));
%         
%         WIM.(Section).(AE).E(Span/SpanDiv) = ESIA.Total(r);
%         WIM.(Section).(AE).Eq(Span/SpanDiv) = ESIA.Eq(r);
%         WIM.(Section).(AE).EQ1(Span/SpanDiv) = ESIA.EQ(1,r);
%         WIM.(Section).(AE).EQ2(Span/SpanDiv) = ESIA.EQ(2,r);
%         
%     end
% end

% Section = 'Twin';
% 
% for r = 25:length(ILData)
%     % Get InfName (call it Temp)
%     %Temp = UniqNames{r}; %Fixed 26.5.20
%     Temp = ILData(r).Name;
%     Temp2 = strfind(Temp,'S');
%     Temp3 = strfind(Temp,'.');
%     % Get Span from InfName
%     %Span = str2num(Temp(end-1:end));
%     Span = str2num(Temp(Temp2(end)+1:end));
%     % Get Action Effect from InfName
%     %AE = Temp(1:end-2);
%     AE = Temp(Temp3(end-1)+1:Temp3(end)-1);
%     
%     try
%         if istable(WIM.(Section).(AE))
%             WIM.(Section).(AE).YClass(Span/SpanDiv) = ST(r).SumTableLogNorm.('Yearly Class')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).YClassOW(Span/SpanDiv) = ST(r).SumTableLogNorm.('Yearly ClassOW')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).YAll(Span/SpanDiv) = ST(r).SumTableLogNorm.('Yearly All')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).WClass(Span/SpanDiv) = ST(r).SumTableLogNorm.('Weekly Class')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).WClassOW(Span/SpanDiv) = ST(r).SumTableLogNorm.('Weekly ClassOW')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).WAll(Span/SpanDiv) = ST(r).SumTableLogNorm.('Weekly All')(end-1)*(ESIA.Total(r));
%             
%             WIM.(Section).(AE).E(Span/SpanDiv) = ESIA.Total(r);
%             WIM.(Section).(AE).Eq(Span/SpanDiv) = ESIA.Eq(r);
%             WIM.(Section).(AE).EQ1(Span/SpanDiv) = ESIA.EQ(1,r);
%             WIM.(Section).(AE).EQ2(Span/SpanDiv) = ESIA.EQ(2,r);
%             
%         end
%     catch
%         WIM.(Section).(AE) = XT;
%         
%         WIM.(Section).(AE).YClass(Span/SpanDiv) = ST(r).SumTableLogNorm.('Yearly Class')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).YClassOW(Span/SpanDiv) = ST(r).SumTableLogNorm.('Yearly ClassOW')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).YAll(Span/SpanDiv) = ST(r).SumTableLogNorm.('Yearly All')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).WClass(Span/SpanDiv) = ST(r).SumTableLogNorm.('Weekly Class')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).WClassOW(Span/SpanDiv) = ST(r).SumTableLogNorm.('Weekly ClassOW')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).WAll(Span/SpanDiv) = ST(r).SumTableLogNorm.('Weekly All')(end-1)*(ESIA.Total(r));
%         
%         WIM.(Section).(AE).E(Span/SpanDiv) = ESIA.Total(r);
%         WIM.(Section).(AE).Eq(Span/SpanDiv) = ESIA.Eq(r);
%         WIM.(Section).(AE).EQ1(Span/SpanDiv) = ESIA.EQ(1,r);
%         WIM.(Section).(AE).EQ2(Span/SpanDiv) = ESIA.EQ(2,r);
%         
%     end
% end
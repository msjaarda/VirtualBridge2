% --------- qInvestigation ---------

% Find alpha small q, given alpha Q1 Q2. From Q1 and Q1/Q2 Investigation we get aQ1 = 0.60 and aQ2 = 0.50.
% Produced using vars in qInvestigation folder (Misc) which comes from VBWIMq.
% This is part of an initial effort to estimate alpha q using OBSERVED WIM.

% COULD IMPROVE COLORS... 

% Initial commands
clear, clc, close all, FigNum = 0;

% ---- INPUT ----
% Analysis Type Option
AType = 'BoxTwin';  % BoxTwin, SlabShort, or SlabLong
% Which IL to Highlight?
HL = [];
% -----  END  ----

% Auto-Inputs
BM = {'Daily', 'Weekly', 'Yearly'};             % j
ClassType = {'All', 'ClassOW', 'Class'};        % i
ClassT = {'All', 'Classified+', 'Classified'};
DistTypes = {'Lognormal'};%, 'GeneralizedExtremeValue'};
C = linspecer(9); D = linspecer(18);      % Set colours
ScaleDown = [1 2.5 5]; % ScaleDown

% Get ANames, xdata, and Title corresponding to AType input
if strcmp(AType,'BoxTwin')
    AName = 'Box Girder Bridges'; AName2 = 'Twin Girder Bridges';
    xdata = 10:10:80; Title{1} = 'V'; Title{2} = 'M+'; Title{3} = 'M-'; SP = [1, 9, 17]; SP2 = [25, 33, 41]; Title2 = Title;
else
    Title{1} = 'M+ p1'; Title{2} = 'M+ p2'; Title{3} = 'M+ p3'; SP = [1 11 21];
    Title2{1} = 'M- p1'; Title2{2} = 'M- p2'; Title2{3} = 'M- p3'; SP2 = [6 16 26];
    if strcmp(AType,'SlabShort')
        xdata = 2:2:10; AName = 'Short Slab Bridges'; AName2 = AName;
    elseif strcmp(AType,'SlabLong')
        xdata = 10:5:30; AName = 'Long Slab Bridges'; AName2 = AName;
    end
end

% Possible to load variables and skip first section
Initial = false;
    
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

% --- Step 5: Plot Block Maxima ---
if ~isempty(HL)
    HLName = ILData(HL).Name; HLName = HLName(7:end);
end

for r = HL;
     
    figure('Position',[300 300 1200 350]);
    xLims = [0 max(x_values(r).All.Yearly)+100*max(x_values(r).All.Yearly)/1500];
    
    for i = 1:length(ClassType)
        Class = ClassType{i};

        subplot(1,3,i)
        hold on

        for j = 1:length(BM)
            BlockM = BM{j};
            
          % X is for the plot                   x is for the bar
            X = x_values(r).(Class).(BlockM);   x = X(1:end-1) + diff(X);

            y = histcounts(Max(r).(Class).(BlockM).Max,'BinEdges',X,'normalization','pdf');
            bar(x,y/ScaleDown(j),1,'EdgeColor',C(j,:),'FaceColor',[.8 .8 .8],'FaceAlpha',0.5,'DisplayName',BlockM)

            for k = 1
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
        xlabel(HLName)
        title([ClassT{i}])
        legend('location','northeast')

    end
end


% --- Step 6: Probability Paper Plotting & Tail Fitting Methods ---
% We are assuming the maximums are normally or log-normally distributed. 
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
            
            [~,~,PPx,PPy,Fity,Fit(i,j,r)] = GetLogNormPPP(Max(r).(Class).(BlockM).Max,false);
            
            if PlotB(r) == 1
                subplot(1,3,i);  hold on
                scatter(PPx,PPy,7,C(j,:),'filled','DisplayName','Max Data');
                plot(PPx,Fity,'-','Color',C(j,:),'DisplayName',['Fitted ' num2str(Fit(i,j,r),3)]);
                box on
                y1 = ylim;
                text(1,y1(1)+(y1(2)-y1(1))*(j/6),sprintf('R^2:  %.1f%%',Fit(i,j,r)*100),"Color",C(j,:))
            end
        end
        
        if PlotB(r) == 1
            if i == 1
                ylabel(strcat('log ( ', HLName, ')'))
            end
            xlabel('Standard Normal Percentile')
            %legend('Location', 'southeast')
            title([ClassT{i}])
        end
    end
end


% --- Step 7 : Summarize  Estimates ---
% The alpha we are finding and putting in the table right is composite alpha...
% We can add AQ1 and AQ2 to get the real Aq
Beta.Yearly = 4.700;
Alpha = 0.7;
AQ1 = 0.6; AQ2 = 0.5;

% For each Influence Case
for r = 1:length(ILData)
    
    % Initialize SumTable, ST
    ST(r).SumTable = array2table(zeros(4,length(BM)*length(ClassType)));
    ST(r).SumTable.Properties.VariableNames = {[BM{3} ' ' ClassType{1}], [BM{3} ' ' ClassType{2}], [BM{3} ' ' ClassType{3}], [BM{2} ' ' ClassType{1}],...
        [BM{2} ' ' ClassType{2}], [BM{2} ' ' ClassType{3}], [BM{1} ' ' ClassType{1}], [BM{1} ' ' ClassType{2}], [BM{1} ' ' ClassType{3}]};
    ST(r).SumTable.Properties.RowNames = {'Ed', 'LNFit (%)', 'AlphaQLN', 'AlphaqLN'};

    ST(r).SumTableComp = array2table(zeros(2,(length(BM)-1)*length(ClassType)));
    ST(r).SumTableComp.Properties.RowNames = {'AlphaQLN', 'AlphaqLN'};
    ST(r).SumTableComp.Properties.VariableNames = {['W2Y ' ClassType{1}],['W2Y ' ClassType{2}],['W2Y ' ClassType{3}],['D2Y ' ClassType{1}],['D2Y ' ClassType{2}],['D2Y ' ClassType{3}]};
    
    for i = 1:length(ClassType)
        Class = ClassType{i};
        
        for j = 1:length(BM)
            BlockM = BM{j};
            
            [EdLN, AQLN, AqLN] = GetBlockMaxEd(Max(r).(Class).(BlockM).Max,BlockM,'Lognormal',...
                ESIA.Total(r),ESIA.EQ(:,r),ESIA.Eq(:,r),AQ1,AQ2);
            
            % ST stands for Summary Table
            ST(r).SumTable.([BlockM ' ' Class]) = [EdLN; Fit(i,j,r); AQLN; AqLN];
            
            if strcmp(BlockM,'Weekly')
                pdx = pd(r).(Class).W2Y.Lognormal;
                Delta2 = pdx.sigma^2;
                Em = exp(pdx.mu+0.5*Delta2);
                ST(r).SumTableComp.(['W2Y ' Class]) = [Em*exp(Alpha*Beta.Yearly*sqrt(Delta2)-0.5*Delta2)/(ESIA.Total(r)); (((Em*exp(Alpha*Beta.Yearly*sqrt(Delta2)-0.5*Delta2))/1.5)-AQ1*ESIA.EQ(1,r)-AQ2*ESIA.EQ(2,r))/(ESIA.Eq(r))];
            elseif strcmp(BlockM,'Daily')
                pdx = pd(r).(Class).D2Y.Lognormal;
                Delta2 = pdx.sigma^2;
                Em = exp(pdx.mu+0.5*Delta2);
                ST(r).SumTableComp.(['D2Y ' Class]) = [Em*exp(Alpha*Beta.Yearly*sqrt(Delta2)-0.5*Delta2)/(ESIA.Total(r)); (((Em*exp(Alpha*Beta.Yearly*sqrt(Delta2)-0.5*Delta2))/1.5)-AQ1*ESIA.EQ(1,r)-AQ2*ESIA.EQ(2,r))/(ESIA.Eq(r))];
            end
            
            PlotAq.(Class).(BlockM)(r) = ST(r).SumTable{'AlphaqLN',[BlockM ' ' Class]};
            PlotAQ.(Class).(BlockM)(r) = ST(r).SumTable{'AlphaQLN',[BlockM ' ' Class]};
        end
    end
end

format bank, fprintf('\n')
if ~isempty(HL), disp(ST(HL).SumTable), end

% Weekly data only
j = 2; BlockM = BM{j};

% Prepare plot data for VBTriPlot
% Titles & Set starting points for ydata are at TOP
FigTitle = ['LEs from WIM on ' AName];
% Prep 3 subplots
for m = 1:3
    % All ClassTypes
    for i = 1:length(ClassType)
        Class = ClassType{i};
        % Get ydata
        ydata{m}(i,:) = PlotAQ.(Class).(BlockM)(SP(m):SP(m)+length(xdata)-1);
        % DisplayName
        PDets.DN{m,i} = strcat(Class,BlockM);
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDets.MEC{m}(i,:) = [0 0 0]; PDets.MFC{m}(i,:) = D(i*6-4,:); 
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum,FigTitle);

% Plot #2
FigTitle = ['LEs from WIM on ' AName2];
% Prep 3 subplots
for m = 1:3
    % All ClassTypes
    for i = 1:length(ClassType)
        Class = ClassType{i};
        % Get ydata
        ydata{m}(i,:) = PlotAQ.(Class).(BlockM)(SP2(m):SP2(m)+length(xdata)-1);
        % DisplayName
        PDets.DN{m,i} = strcat(Class,BlockM);
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5; PDets.MEC{m}(i,:) = [0 0 0]; PDets.MFC{m}(i,:) = D(i*6-4,:); 
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title2,'WIM',FigNum,FigTitle);

% --- Step 8: Populate WIM, a comparison with AGB and MAT Results Variables ---
% First choose which result we will use (try Lognormal Weekly and Yearly)

% ------ NOT FINISHED YET... WAIT FOR NEW ILLIB FINALIZATION ------

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
%             WIM.(Section).(AE).YClass(Span/SpanDiv) = ST(r).SumTable.('Yearly Class')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).YClassOW(Span/SpanDiv) = ST(r).SumTable.('Yearly ClassOW')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).YAll(Span/SpanDiv) = ST(r).SumTable.('Yearly All')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).WClass(Span/SpanDiv) = ST(r).SumTable.('Weekly Class')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).WClassOW(Span/SpanDiv) = ST(r).SumTable.('Weekly ClassOW')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).WAll(Span/SpanDiv) = ST(r).SumTable.('Weekly All')(end-1)*(ESIA.Total(r));
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
%         WIM.(Section).(AE).YClass(Span/SpanDiv) = ST(r).SumTable.('Yearly Class')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).YClassOW(Span/SpanDiv) = ST(r).SumTable.('Yearly ClassOW')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).YAll(Span/SpanDiv) = ST(r).SumTable.('Yearly All')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).WClass(Span/SpanDiv) = ST(r).SumTable.('Weekly Class')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).WClassOW(Span/SpanDiv) = ST(r).SumTable.('Weekly ClassOW')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).WAll(Span/SpanDiv) = ST(r).SumTable.('Weekly All')(end-1)*(ESIA.Total(r));
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
%             WIM.(Section).(AE).YClass(Span/SpanDiv) = ST(r).SumTable.('Yearly Class')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).YClassOW(Span/SpanDiv) = ST(r).SumTable.('Yearly ClassOW')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).YAll(Span/SpanDiv) = ST(r).SumTable.('Yearly All')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).WClass(Span/SpanDiv) = ST(r).SumTable.('Weekly Class')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).WClassOW(Span/SpanDiv) = ST(r).SumTable.('Weekly ClassOW')(end-1)*(ESIA.Total(r));
%             WIM.(Section).(AE).WAll(Span/SpanDiv) = ST(r).SumTable.('Weekly All')(end-1)*(ESIA.Total(r));
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
%         WIM.(Section).(AE).YClass(Span/SpanDiv) = ST(r).SumTable.('Yearly Class')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).YClassOW(Span/SpanDiv) = ST(r).SumTable.('Yearly ClassOW')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).YAll(Span/SpanDiv) = ST(r).SumTable.('Yearly All')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).WClass(Span/SpanDiv) = ST(r).SumTable.('Weekly Class')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).WClassOW(Span/SpanDiv) = ST(r).SumTable.('Weekly ClassOW')(end-1)*(ESIA.Total(r));
%         WIM.(Section).(AE).WAll(Span/SpanDiv) = ST(r).SumTable.('Weekly All')(end-1)*(ESIA.Total(r));
%         
%         WIM.(Section).(AE).E(Span/SpanDiv) = ESIA.Total(r);
%         WIM.(Section).(AE).Eq(Span/SpanDiv) = ESIA.Eq(r);
%         WIM.(Section).(AE).EQ1(Span/SpanDiv) = ESIA.EQ(1,r);
%         WIM.(Section).(AE).EQ2(Span/SpanDiv) = ESIA.EQ(2,r);
%         
%     end
% end
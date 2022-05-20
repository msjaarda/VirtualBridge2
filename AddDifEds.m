clear, clc, close all

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
IncZ = false;

for r = 1:length(AE)
    for k = 1:length(Class)
        for j = 1:length(BlockM)
            
            % Get a subset of the maximum results for the three selections made
            Data = OI.Max(AE(r)).(Class{k}).(BlockM{j}).Max;
            Data = [Data; zeros(round(length(Data)/OI.PropTrucks.(Class{k}).(BlockM{j})(AE(r)))-length(Data),1)];
            pdNEW = GetFit(Data,BlockM{j},[Dist 'LM'],1,IncZ);
            
            % Set influence line you want - take a look at OutInfo.ILData to know which
            ILName = OI.ILData(AE(r)).Name(7:end); ILName = strrep(ILName,'.',' ');
            
            pd = pdNEW.(pdNEW.Best).pd;
            Ed = pdNEW.(pdNEW.Best).Ed;
            
            % Plot
            Top = ceil(max(Data)/10)*10;
            Bot = floor(min(Data)/10)*10;
            TBDiff = Top-Bot;
            
            X = linspace(Bot-TBDiff*.1,max(Ed*1.05,Top+TBDiff*.1));   % X is for the plot
            XB = linspace(Bot-TBDiff*.1,max(Ed*1.05,Top+TBDiff*.1),25);   % X is for the plot
            x = XB(1:end-1) + diff(XB);                           % x is for the bar
            if IncZ
                y = histcounts(Data,'BinEdges',XB,'normalization','pdf');
            else
                y = histcounts(OI.Max(AE(r)).(Class{k}).(BlockM{j}).Max,'BinEdges',XB,'normalization','pdf');
            end
            
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
            text(X(70),y1(1)+(y1(2)-y1(1))*.70,sprintf('R^2:           %.1f%%',pdNEW.(pdNEW.Best).R2),"Color",'k')
            text(X(70),y1(1)+(y1(2)-y1(1))*.65,sprintf('TailFit:      %s',mat2str(TailFit)),"Color",'k')
            text(X(70),y1(1)+(y1(2)-y1(1))*.60,sprintf('IncZeros:  %s',mat2str(IncZ)),"Color",'k')
            text(Ed,y1(1)+(y1(2)-y1(1))*.25,sprintf('Ed = %.1f',Ed),"Color",'k','HorizontalAlignment','center')
            text(X(5),y1(1)+(y1(2)-y1(1))*.85,sprintf('%.1f%% Accompaniment Rate',100*OI.PropTrucks.(Class{k}).(BlockM{j})(AE(r))),"Color",'k')
            text(X(5),y1(1)+(y1(2)-y1(1))*.80,sprintf('%i Total Events',length(Data)),"Color",'k')
            
            set(gca,'ytick',[],'yticklabel',[],'ycolor','k')
            ylabel('Normalized Histograms (NTS)')
            xlabel('Bridge Action Effect')
            
            z = 1;
            % Convergence
            for p = 5:5:100
                r = randperm(height(Data),round(height(Data)*p/100));
                Data2 = Data(r);
                pdNEW = GetFit(Data2,BlockM{j},Dist,0,IncZ);
                pd = pdNEW.(pdNEW.Best).pd;
                Ed2(z) = pdNEW.(pdNEW.Best).Ed;
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


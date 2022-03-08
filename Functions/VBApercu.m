function T = VBApercu(PDC,Title,ILData,BrStInd,TrLineUp,PEsia,DLF,Lane,ILRes)
% Plot a Series of WIM or VWIM Vehicles on a Bridge

% Take only the influence lines that apply to the current InfCase
Infv = ILData.v;
Infx = 0:ILRes:(length(ILData.v)-1)*ILRes;
Span = max(Infx);

% if isWIM we will have Lane.Sites and Lane.Station
if any(ismember(PDC.Properties.VariableNames,'DTS'))
    isWIM = true;
else
    isWIM = false;
end

% Initialize Figure
figure

% Just for displaying class, see Classify...
[TrTyps,TrNames] = VBTypes2Names;

% Get number and name of lanes
LaneswTr = unique(PDC.LANE); 
% Changed 26/08/21 to fix a bug when PDC has only lane 2 for ex.
LaneswTr = [1:LaneswTr(end)]';
NumLaneswTr = length(LaneswTr);
% Get TotalLanes and get Aperçu lanes (ALane)

TotalLanes = Lane.Sites.NumLanes;
%ALane = Lane.Details.ALANE;
% We need to sort Lane.Details by LANE first... then do this.
Lane.Details = sortrows(Lane.Details,3);
ALane = Lane.Details.ALANE;

% Get total number of subplots
NumSubplots = TotalLanes + 2;

% Define Plot Colors RBGY
Col{1} = [.94 .28 .18]; Col{2} = [0 .447 .74]; Col{3} = [.184 .8 .086];       
Col{4} = [.99 .67 0]; Col{5} = Col{2}; Col{6} = Col{3};

% Plot Influence Line - Open up subplot and choose the last subplot
sp(NumSubplots) = subplot(NumSubplots,1,NumSubplots);
% Note that trucks go the other way than is plotted... must flip ILs
if size(Infv,2) == 1 | all(all(Infv == Infv(:,1),2))
    plot(Infx,flip(-Infv),'k','LineWidth',1.5)
else
    for i = 1:size(Infv,2)
        hold on
        plot(Infx,flip(-Infv(:,i)),'Color',Col{i},'LineWidth',1.5)
    end
end
xlabel('Distance Along Bridge (m)'); ylabel('Inf Line'); xlim([-0.5 Span+0.5]); box on
text(1,-max(max(Infv))+max(max(Infv))/5,sprintf('%.1f%% of E_{SIA}',PEsia*100),'FontSize',11,'FontWeight','bold','Color','k')
PerI = find(ILData.Name == '.');
text(Span-1,-max(max(Infv))+max(max(Infv))/5,sprintf('%s',strrep(ILData.Name(PerI(1)+1:PerI(end)-1),'.',' ')),'FontSize',11,'FontWeight','bold','Color','k','HorizontalAlignment','right')
set(gca,'ytick',[]); set(gca,'yticklabel',[])

% Define overall plot title
sgtitle(Title);

% Q / q    [     1            2         3         4          5     ]
%           AllTrAxIndex  AxleValue   Truck#    LaneID   Station(m)

% Define Q, an excerpt of TrLineUp with just those vehicles on the bridge
Q = TrLineUp(TrLineUp(:,1) >= BrStInd & TrLineUp(:,1) <= BrStInd+length(Infx)-1,:); % added equals...
% Define T, an excerpt of WIM/VWIM PDC with just vehicles on the bridge
T = PDC(unique(Q(:,3)),:);

if height(T) == 0
posi = find(TrLineUp(:,1) <= BrStInd);
posi = posi(end);
T = PDC(TrLineUp(posi,3)+1,:);
end
% si T c'est zero alors utiliser BrStInd pour retrouver le dernier TrLineUp
% et donc le dernier PDC puis faire PDC+1

% Gather variables for the plots of each lane (vehicle corners and axle vals)
for i = 1:NumLaneswTr
  
    % If there is no traffic on a given lane
    if i > NumLaneswTr
        % Set vehicle corners and accummulated axles to zero
        ac{i} = 0; vc{i} = 0;
        NoVeh(i) = 1;
        continue
    end
    % q is a subset of Q, t is a subset of T
    q{i} = Q(Q(:,4) == LaneswTr(i),:); t{i} = T(T.LANE == LaneswTr(i),:);
    % normalize q values for start of the bridge at zero
    q{i}(:,1) = round((q{i}(:,1) - BrStInd)); q{i}(:,5) = q{i}(:,5) - BrStInd*ILRes;
    [a, b] = unique(q{i}(:,3));
    % vc stands for vehicle corners, ac for accummulated
    if ~isempty(b)
        for j = 1:length(b)
            locations = q{i}(q{i}(:,3) == a(j),1);
            axlevalues = q{i}(q{i}(:,3) == a(j),2);

            % Assign initial if necessary
            if sum(locations == 0) > 0
                % Added sum 1/6/21
                initial(i,j) = sum(axlevalues(locations == 0));
            else
                initial(i,j) = 0;
            end
            % Assign ac if necessary
            if sum(locations) > 0
                ac{i}(:,j) = accumarray(locations(locations ~= 0),axlevalues(locations ~= 0),[Infx(end)/ILRes 1]);
            else
                ac{i}(:,j) = zeros([Infx(end)/ILRes 1]);
            end
            
            Temp = TrLineUp(TrLineUp(:,3) == a(j),5);
            %Temp = sort(Temp,'descend');
            if Lane.Dir(i) == 1
                Temp = sort(Temp,'ascend');
                vc{i}(j,1) = Temp(1)-ILRes*BrStInd-1;
                vc{i}(j,2) = Temp(end)-ILRes*BrStInd+1;
            else
                 Temp = sort(Temp,'descend');
                 vc{i}(j,2) = Temp(1)-ILRes*BrStInd+1;
                 vc{i}(j,1) = Temp(end)-ILRes*BrStInd-1;
            end
        end
        barp(:,i) = [sum(initial(i,:),2); sum(ac{i},2)];
        NoVeh(i) = 0;
    else
        barp(1:Infx(end)/ILRes+1,i) = 0;
        ac{i} = 0; vc{i} = 0;
        NoVeh(i) = 1;
    end
end

% Need at least a value for all lanes
if size(barp,2) < TotalLanes
    barp(1:size(barp,1), size(barp,2)+1:TotalLanes) = 0.001;
end

% Plot axle loads - Open up subplot and choose the second last subplot
sp(TotalLanes+1) = subplot(NumSubplots,1,TotalLanes+1);
% Use handle, h, to assign colors later
h = bar(0:ILRes:Infx(end),barp/9.81,1.2/ILRes,'grouped','EdgeColor','k');
xlim([-0.5 Span+0.5]);
ylim([0 max(ceil(max(max(barp/9.81))/5)*5,1)])
ylabel('Axle Loads (t)')
set(gca,'xticklabel',[])

% Show text of DLA and Total Weight
text(1,ceil(max(max(barp/9.81))/5)*5-3,sprintf('Total: %.0f (DLF = %.2f)',sum(sum(barp)),DLF),'FontSize',11,'FontWeight','bold','Color','k')

if ILRes ~= 1
    xtemp = 0:ILRes:Infx(end);
    Infvtemp = interp1(Infx,Infv,xtemp);
    if size(Infvtemp,1) == 1
        Infvtemp = Infvtemp';
    end
else
    Infvtemp = Infv;
end
if size(Infvtemp,2) < TotalLanes
    Infvtemp(1:size(Infvtemp,1), size(Infvtemp,2)+1:TotalLanes) = 0;
end

% Show text of Load Effect and % of ESIA
text(Span-1,ceil(max(max(barp/9.81))/5)*5-3,sprintf('Load Effect: %.0f',DLF*sum(sum(barp(1:end,:).*flip(Infvtemp(:,1:end))))),'FontSize',11,'FontWeight','bold','Color','k','HorizontalAlignment','right')

% Set Line colors
for i = 1:NumLaneswTr
    if NoVeh(i) ~= i
        h(i).FaceColor = Col{i};
    end
end


% j counts up in order (1, 2, 3)
% ALanes is Aperçu lanes
% sp(1:TotalLanes) are in the order top to bottom

% We go in order j 1:TotalLanes, but actual we use ALane(j), so we don't
% plot in order downwards always ALane(j) is the plot we are on

% Code order, Actual Names, ALANE [[1:length(Lanes)]' Lane.Details.ALANE]
for j = 1:TotalLanes
    
    % Use ALANES for selecting subplot number
    sp(ALane(j)) = subplot(NumSubplots,1,ALane(j));
    
    % Write extra data on first subplot if WIM
    if ALane(j) == 1 && isWIM
        text(0,9,datestr(T.DTS(1),'mmm-dd-yyyy'),'FontSize',9,'FontWeight','bold','HorizontalAlignment','left','Color','w')
        text(0,1.25,datestr(T.DTS(1),'HH:MM:SS'),'FontSize',9,'FontWeight','bold','HorizontalAlignment','left','Color','w')
        % North Arrow and Canton/Hwy
        if isfield(Lane,'Details')
            if Lane.Details.NSEW == 4 | Lane.Details.NSEW == 3
                text(Span,9,['NORTH ' char(8593)],'FontSize',9,'FontWeight','bold','HorizontalAlignment','right','Color','w')
            else
                text(Span,9,['NORTH ' char(8594)],'FontSize',9,'FontWeight','bold','HorizontalAlignment','right','Color','w')
            end
            text(Span,1.25,Lane.Sites.CANTON + " (" + Lane.Sites.HWY + ")",'FontSize',9,'FontWeight','bold','HorizontalAlignment','right','Color','w')
        end
    end
    % Write extra data on second subplot if we have Lane Details
    %if ALane(j) == 2 && isfield(Lane,'Details')
    if ALane(j) == 2 && isWIM
        if Lane.Details.NSEW(j) == 2 | Lane.Details.NSEW(j) == 4
            text(0,1.25,char(8592) + " " + Lane.Details.DIR(j),'FontSize',9,'FontWeight','bold','HorizontalAlignment','left','Color','w')
            text(Span,1.25,Lane.Details.FROM(j) + " " + char(8594),'FontSize',9,'FontWeight','bold','HorizontalAlignment','right','Color','w')
        else
            text(0,1.25,char(8592) + " " + Lane.Details.FROM(j),'FontSize',9,'FontWeight','bold','HorizontalAlignment','left','Color','w')
            text(Span,1.25,Lane.Details.DIR(j) + " " + char(8594),'FontSize',9,'FontWeight','bold','HorizontalAlignment','right','Color','w')
        end
    end
    
    % If we are on a plot that has Tr, draw the trucks
    if j <= NumLaneswTr
        % Assign y scale for trucks
        Lo = 2; Hi = 8;
        % Adjust scale for middle lanes
        if ALane(j) > 1 & TotalLanes > 2 & ALane(j) < TotalLanes
            Lo = Lo + 0.20; Hi = Hi - 0.20;
        end
        % Adjust scale for Inf Length
        if Span > 50;
            Lo = min(Lo + (Span-50)/35,3.5); Hi = max(Hi - (Span-50)/35,6.5);
        end
        % Difference
        DiF = Hi - Lo; % Normally 6
        
        % For each truck... draw it
        for i = 1:numel((vc{j}))/2
            hold on
            EdgCol = 'k';
            % Truck Outline
            fill([vc{j}(i,1) vc{j}(i,1) vc{j}(i,2) vc{j}(i,2)],[Lo Hi Hi Lo],Col{j},'EdgeColor',EdgCol,'LineWidth',1.5);
            if Lane.Dir(j) == 1
                % Back Bumper
                fill([vc{j}(i,2)-0.1 vc{j}(i,2)-0.1 vc{j}(i,2)+0.1 vc{j}(i,2)+0.1],[Lo+0.5 Hi-0.5 Hi-0.5 Lo+0.5],'w','EdgeColor',EdgCol,'LineWidth',1.5);
                % Front of Truck
                fill([vc{j}(i,1)-0.7*Span/50 vc{j}(i,1)-0.7*Span/50 vc{j}(i,1) vc{j}(i,1)],[Lo+1.75*DiF/6 Hi-1.75*DiF/6 Hi-DiF/6 Lo+DiF/6],[.4 .4 .4],'EdgeColor',EdgCol,'LineWidth',1.5);
                % Mirrors
                fill([vc{j}(i,1)+1.1 vc{j}(i,1)+1.3 vc{j}(i,1)+1.5 vc{j}(i,1)+1.3],[Lo Lo-0.5 Lo-0.5 Lo],'w','EdgeColor',EdgCol,'LineWidth',1.5);
                fill([vc{j}(i,1)+1.1 vc{j}(i,1)+1.3 vc{j}(i,1)+1.5 vc{j}(i,1)+1.3],[Hi Hi+0.5 Hi+0.5 Hi],'w','EdgeColor',EdgCol,'LineWidth',1.5);
            else
                % Back Bumper
                fill([vc{j}(i,1)-0.1 vc{j}(i,1)-0.1 vc{j}(i,1)+0.1 vc{j}(i,1)+0.1],[Lo+0.5 Hi-0.5 Hi-0.5 Lo+0.5],'w','EdgeColor',EdgCol,'LineWidth',1.5);
                % Front of Truck
                fill([vc{j}(i,2)+0.7*Span/50 vc{j}(i,2)+0.7*Span/50 vc{j}(i,2) vc{j}(i,2)],[Lo+1.75*DiF/6 Hi-1.75*DiF/6 Hi-DiF/6 Lo+DiF/6],[.4 .4 .4],'EdgeColor',EdgCol,'LineWidth',1.5);
                % Mirrors
                fill([vc{j}(i,2)-1.1 vc{j}(i,2)-1.3 vc{j}(i,2)-1.5 vc{j}(i,2)-1.3],[Lo Lo-0.5 Lo-0.5 Lo],'w','EdgeColor',EdgCol,'LineWidth',1.5);
                fill([vc{j}(i,2)-1.1 vc{j}(i,2)-1.3 vc{j}(i,2)-1.5 vc{j}(i,2)-1.3],[Hi Hi+0.5 Hi+0.5 Hi],'w','EdgeColor',EdgCol,'LineWidth',1.5);
            end
            
            % Only write text if it is within the plot...
            if (vc{j}(i,1)+vc{j}(i,2))/2 > Span/15 && (vc{j}(i,1)+vc{j}(i,2))/2 < Span-Span/15
                text((vc{j}(i,1)+vc{j}(i,2))/2,5,sprintf('%i ax | %.1f t',t{j}.AX(i),t{j}.GW_TOT(i)/1000),'FontSize',11,'FontWeight','bold','VerticalAlignment','middle','HorizontalAlignment','center','Color','k')
            end
            if (vc{j}(i,1)+vc{j}(i,2))/2 > Span/8 && (vc{j}(i,1)+vc{j}(i,2))/2 < Span-Span/8 && ismember('SPEED', T.Properties.VariableNames)
                text((vc{j}(i,1)+vc{j}(i,2))/2,Hi+1,sprintf('%.0f kph',t{j}.SPEED(i)),'FontSize',11,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle','Color','w')
            end
            if (vc{j}(i,1)+vc{j}(i,2))/2 > Span/8 && (vc{j}(i,1)+vc{j}(i,2))/2 < Span-Span/8
                text((vc{j}(i,1)+vc{j}(i,2))/2,Lo-0.75,sprintf('%s',TrNames(TrTyps == t{j}.CLASS(i))),'FontSize',9,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle','Color','w')
            end
        end
        
        % Add all wheel locations (column 5 has actual wheel locations, column 1 would be approximate)
        hold on
        scatter(q{j}(:,5),(Hi-DiF/6)*ones(size(q{j},1),1),DiF*6,'filled','s','MarkerFaceColor','k')
        scatter(q{j}(:,5),(Lo+DiF/6)*ones(size(q{j},1),1),DiF*6,'filled','s','MarkerFaceColor','k')

    end
    
    % If we are NOT at the bottom lane
    if ALane(j) ~= TotalLanes
        % Add line between lanes
        %if Lane.Details.NSEW(Lane.Details.ALANE(j)) == Lane.Details.NSEW(Lane.Details.ALANE(j)+1)
        if Lane.Details.NSEW(j) == Lane.Details.NSEW(ALane == ALane(j)+1)
            % Dashed if between same direction
            hold on
            xlim([-0.5 Span+0.5]); ylim([0 10])
            %dashline([0 Span],[0.1 0.1],Span/10,Span/10,Span/10,Span/10,'w','LineWidth',2);
            dashline([0 Span],[0.1 0.1],5,5,5,5,'w','LineWidth',2);
            %dashedline([0 Span],[0.1 0.1],Span/10,'w','LineWidth',2);
        else
            % Solid if change of direction
            yline(0.1,'-w','LineWidth',2);
        end
        
        % Turn off x axis completely
        set(gca,'xticklabel',[],'xcolor','none','xtick',[])
        
        
        
    end

    % Set axis limits
    xlim([-0.5 Span+0.5]); ylim([0 10])
    % Label y axis
    if isfield(Lane,'Details')
        ylabel(['Lane ' num2str(Lane.Details.LANE(j))]); set(gca,'ytick',[]); set(gca,'yticklabel',[])
    else
        ylabel(['Lane ' num2str(j)]); set(gca,'ytick',[]); set(gca,'yticklabel',[])
    end
    % Turn y axis off
    yaxish = gca; yaxish.YRuler.Axle.Visible = 'off';

    % Set the background lane color to grey, like asphalt
    set(gca,'Color',[.6 .6 .6])
end

% Set Position
set(gcf,'Position',[100 100 900 750])

% Left Bottom Width Height
Left = sp(1).Position(1);
Width = sp(1).Position(3);
GapSize = (sp(end-1).Position(2) - sp(end).Position(4) - sp(end).Position(2))/2;
NewHeight = sp(1).Position(4) + GapSize;
NewHeightIn = sp(1).Position(4) + GapSize + GapSize;

% Remove gap between InfLine and Axle Loads
set(sp(end), 'Position',   [Left   sp(end).Position(2)                     Width    NewHeight]);
set(sp(end-1), 'Position', [Left   sp(end-1).Position(2) - GapSize         Width    NewHeight]);

% Remove gap between lanes
for i = 1:TotalLanes
    if i == 1
        set(sp(i), 'Position',   [Left   sp(i).Position(2) - GapSize       Width    NewHeight]);
    elseif i == TotalLanes
        set(sp(i), 'Position',   [Left   sp(i).Position(2)                 Width    NewHeight]);
    else
        set(sp(i), 'Position',   [Left   sp(i).Position(2) - GapSize       Width    NewHeightIn]);
    end
end

end
    



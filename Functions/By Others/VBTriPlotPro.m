function [FigNum] = VBTriPlotPro(xdata,ydata,PDets,Title,Type,FigNum,FigTitle,LaneTrDistr)
%VBTRIPLOT Makes a 3-panel plot in the style of AGB 2002-005

% Use FigNum to cascade figures
figure('Name',FigTitle,'NumberTitle','off','Position',[200+FigNum*25 200+FigNum*25 730 400]);
FigNum = FigNum + 1;

% What is xdata is different for different cases? We duplicate when it is
% the same always, to allow for the fact that it
if ~iscell(xdata)
    xdata2 = xdata; clear xdata
    for m = 1:3
        for i = 1:size(ydata,2)
            xdata{m,i} = xdata2;
        end
    end
end

%Define positions for all elements
Pos.PosGraph(1,:)=[0.08 0.12 0.18 0.62];
Pos.PosGraph(2,:)=[0.32 0.12 0.18 0.62];
Pos.PosGraph(3,:)=[0.56 0.12 0.18 0.62];
Pos.PosSysStat(1,:)=[0.08 0.79 0.18 0.18];
Pos.PosSysStat(2,:)=[0.32 0.79 0.18 0.18];
Pos.PosSysStat(3,:)=[0.56 0.79 0.18 0.18];
Pos.PosTypePont =[0.78 0.79 0.182 0.18];
Pos.PosTrafi =[0.78 0.12 0.18 0.62];


Len = 0;
for m = 1:3
    
    subplot('Position',Pos.PosGraph(m,:))
    hold on
    for i = 1:size(ydata,2)
        if strcmp(PDets.DN{m,i},'off') || strcmp(PDets.DN{m,i},'none')
            HandV = 'off';
        else, HandV = 'on'; end
        if isequal(PDets.MFC{m}(i,:),[.98 .98 .98])
            plot(xdata{m,i},ydata{m,i},'-s','Color','k','MarkerEdgeColor',PDets.MEC{m}(i,:),'MarkerFaceColor','none','MarkerSize',PDets.MS{m}(i),'DisplayName',PDets.DN{m,i},'HandleVisibility',HandV)
        else
            plot(xdata{m,i},ydata{m,i},'-s','Color','k','MarkerEdgeColor',PDets.MEC{m}(i,:),'MarkerFaceColor',PDets.MFC{m}(i,:),'MarkerSize',PDets.MS{m}(i),'DisplayName',PDets.DN{m,i},'HandleVisibility',HandV)
        end
        % Get xdataL{m} - the longest xdata
        if length(xdata{m,i}) > Len
            xdataL{m} = xdata{m,i};
            Len = length(xdata{m,i});
        end
    end
    
    % Set tick details, x-axis label, and title
    ytickformat('%.2f'); yticks(0:0.1:1); set(gca,'TickDir','out'); set(gca,'YGrid','on'); xticks(xdataL{m}); 
    xlabel('Span (m)')
    title(Title{m})
    % get handle of current, set box property to off and remove background color
    a = gca; set(a,'box','off','color','none');
    % create new, empty axes with box but without ticks
    b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
    % set original axes as active, and link axes in case of zooming
    axes(a); linkaxes([a b]);
    % Set axis limits
    ylim([0 1]); xlim([xdataL{m}(1) xdataL{m}(end)])
     
    % Y-axis only for the leftmost (first) plot
    if m == 1
        if strcmp(Type,'WIM')
            ylabel('E_{WIM}/E_{SIA}')
        else
            ylabel('E_{SIM}/E_{SIA}')
        end
        yh = get(gca,'ylabel'); % handle to the label object
        p = get(yh,'position'); % get the current position property
        p(1) = 0.9*p(1);          % double the distance,
        % negative values put the label below the axis
        set(yh,'position',p)   % set the new position
        legend('location','best')
    end
    Len = 0;
    
    %Import images systèmes statiques
    if contains(FigTitle,'Slab')
        if contains(FigTitle,'Fixed')
            if strcmp(Title(m),'M+')||strcmp(Title(m),'MxMid')
                img = imread('Graphs Drawings/Efforts/100% M+.png');
            else
                img = imread('Graphs Drawings/Efforts/100% M- V.png');
            end
        else
            if strcmp(Title(m),'M+')||strcmp(Title(m),'MxMid')
                img = imread('Graphs Drawings/Efforts/50% M+.png');
            else
                img = imread('Graphs Drawings/Efforts/50% M- V.png');
            end   
        end
    else
        if strcmp(Title(m),'M+')
                img = imread('Graphs Drawings/Efforts/M+.png');
        elseif strcmp(Title(m),'M-')
                img = imread('Graphs Drawings/Efforts/M-.png');
        elseif strcmp(Title(m),'V')
                img = imread('Graphs Drawings/Efforts/V.png');
        end
    end
    subplot('Position',Pos.PosSysStat(m,:));
    image(img);
    set(gca,'TickLength',[0,0]);
    set(gca,'XTick',[],'YTick',[]);
               
end

    %Import image type pont
    if contains(FigTitle,'Box')
    img = imread('Graphs Drawings/TypesPonts/Box.png');
    elseif contains(FigTitle,'Twin') && contains(FigTitle,'Conc')
    img = imread('Graphs Drawings/TypesPonts/TwinC.png');
    elseif contains(FigTitle,'Twin')
    img = imread('Graphs Drawings/TypesPonts/TwinA.png');
    elseif contains(FigTitle,'Multi')
    img = imread('Graphs Drawings/TypesPonts/Multi.png');
    elseif contains(FigTitle,'Slab')
    img = imread('Graphs Drawings/TypesPonts/Slab.png');
    end
    subplot('Position',Pos.PosTypePont);
    image(img);
    set(gca,'TickLength',[0,0]);
    set(gca,'XTick',[],'YTick',[]);
    
    %Trafic distr
    LaneTrDistr = strsplit(char(LaneTrDistr),',');
    [~,TrSize] = size(LaneTrDistr);
    for i=1:TrSize
    TrDist(i) = str2num(char(LaneTrDistr(i)));
    end
    clear TrSize LaneTrDistr
    
    %Import image trafic
    subplot('Position',Pos.PosTrafi);
    if contains(FigTitle,'Uni')
        if contains(FigTitle,'2L')||contains(FigTitle,'Wid12')
        img = imread('Graphs Drawings/DispoVoies/Uni11.png');
        image(img);
        hold on;
        text(150,69,append(int2str(TrDist(1)),'%'),'Color','g','FontSize',15);
        text(150,122.5,append(int2str(TrDist(2)),'%'),'Color','g','FontSize',15);
        elseif contains(FigTitle,'3L')||contains(FigTitle,'Wid15')
        img = imread('Graphs Drawings/DispoVoies/Uni111.png');
        image(img);
        hold on;
        text(150,69,append(int2str(TrDist(1)),'%'),'Color','g','FontSize',15);
        text(150,122.5,append(int2str(TrDist(2)),'%'),'Color','g','FontSize',15);
        text(150,177.5,append(int2str(TrDist(3)),'%'),'Color','g','FontSize',15);
        elseif contains(FigTitle,'4L')||contains(FigTitle,'Wid18')
        img = imread('Graphs Drawings/DispoVoies/Uni1111.png');
        image(img);
        hold on;
        text(150,69,append(int2str(TrDist(1)),'%'),'Color','g','FontSize',15);
        text(150,122.5,append(int2str(TrDist(2)),'%'),'Color','g','FontSize',15);
        text(150,177.5,append(int2str(TrDist(3)),'%'),'Color','g','FontSize',15);
        text(150,232.5,append(int2str(TrDist(4)),'%'),'Color','g','FontSize',15);
        end
    elseif contains(FigTitle,'Bi')
        if contains(FigTitle,'2L')||contains(FigTitle,'Wid12')
        img = imread('Graphs Drawings/DispoVoies/Bi12.png');
        image(img);
        hold on;
        text(150,69,append(int2str(TrDist(1)),'%'),'Color','g','FontSize',15);
        text(150,122.5,append(int2str(TrDist(2)),'%'),'Color','y','FontSize',15);
        elseif contains(FigTitle,'3L')||contains(FigTitle,'Wid15')
        img = imread('Graphs Drawings/DispoVoies/Bi112.png');
        image(img);
        hold on;
        text(150,69,append(int2str(TrDist(1)),'%'),'Color','g','FontSize',15);
        text(150,122.5,append(int2str(TrDist(2)),'%'),'Color','g','FontSize',15);
        text(150,177.5,append(int2str(TrDist(3)),'%'),'Color','y','FontSize',15);
        elseif contains(FigTitle,'4L')||contains(FigTitle,'Wid18')
        img = imread('Graphs Drawings/DispoVoies/Bi1122.png');
        image(img);
        hold on;
        text(150,69,append(int2str(TrDist(1)),'%'),'Color','g','FontSize',15);
        text(150,122.5,append(int2str(TrDist(2)),'%'),'Color','g','FontSize',15);
        text(150,177.5,append(int2str(TrDist(3)),'%'),'Color','y','FontSize',15);
        text(150,232.5,append(int2str(TrDist(4)),'%'),'Color','y','FontSize',15);
        end
    elseif contains(FigTitle,'PUN')
        if contains(FigTitle,'3L')||contains(FigTitle,'Wid12')
        img = imread('Graphs Drawings/DispoVoies/PUN111.png');
        image(img);
        hold on;
        text(150,46,append(int2str(TrDist(1)),'%'),'Color','g','FontSize',15);
        text(150,100.5,append(int2str(TrDist(2)),'%'),'Color','g','FontSize',15);
        text(150,154,append(int2str(TrDist(3)),'%'),'Color','g','FontSize',15);
        elseif contains(FigTitle,'4L')||contains(FigTitle,'Wid15')
        img = imread('Graphs Drawings/DispoVoies/PUN1111.png');
        image(img);
        hold on;
        text(150,46,append(int2str(TrDist(1)),'%'),'Color','g','FontSize',15);
        text(150,100.5,append(int2str(TrDist(2)),'%'),'Color','g','FontSize',15);
        text(150,154,append(int2str(TrDist(3)),'%'),'Color','g','FontSize',15);
        text(150,206,append(int2str(TrDist(4)),'%'),'Color','g','FontSize',15);
        elseif contains(FigTitle,'5L')||contains(FigTitle,'Wid18')
        img = imread('Graphs Drawings/DispoVoies/PUN11111.png');
        image(img);
        hold on;
        text(150,46,append(int2str(TrDist(1)),'%'),'Color','g','FontSize',15);
        text(150,100.5,append(int2str(TrDist(2)),'%'),'Color','g','FontSize',15);
        text(150,154,append(int2str(TrDist(3)),'%'),'Color','g','FontSize',15);
        text(150,206,append(int2str(TrDist(4)),'%'),'Color','g','FontSize',15);
        text(150,257,append(int2str(TrDist(5)),'%'),'Color','g','FontSize',15);
        end
    elseif contains(FigTitle,'Chan')
         if contains(FigTitle,'4L')||contains(FigTitle,'Wid12')
        img = imread('Graphs Drawings/DispoVoies/Chan1122.png');
        image(img);
        hold on;
        text(150,43,append(int2str(TrDist(1)),'%'),'Color','g','FontSize',15);
        text(150,95.5,append(int2str(TrDist(2)),'%'),'Color','g','FontSize',15);
        text(150,150.5,append(int2str(TrDist(3)),'%'),'Color','y','FontSize',15);
        text(150,205.5,append(int2str(TrDist(4)),'%'),'Color','y','FontSize',15);
        elseif contains(FigTitle,'5L')||contains(FigTitle,'Wid15')
        img = imread('Graphs Drawings/DispoVoies/Chan11122.png');
        image(img);
        hold on;
        text(150,43,append(int2str(TrDist(1)),'%'),'Color','g','FontSize',15);
        text(150,95.5,append(int2str(TrDist(2)),'%'),'Color','g','FontSize',15);
        text(150,150.5,append(int2str(TrDist(3)),'%'),'Color','g','FontSize',15);
        text(150,205.5,append(int2str(TrDist(4)),'%'),'Color','y','FontSize',15);
        text(150,260.5,append(int2str(TrDist(5)),'%'),'Color','y','FontSize',15);
        elseif contains(FigTitle,'6L')||contains(FigTitle,'Wid18')
        img = imread('Graphs Drawings/DispoVoies/Chan111222.png');
        image(img);
        hold on;
        text(150,43,append(int2str(TrDist(1)),'%'),'Color','g','FontSize',15);
        text(150,95.5,append(int2str(TrDist(2)),'%'),'Color','g','FontSize',15);
        text(150,150.5,append(int2str(TrDist(3)),'%'),'Color','g','FontSize',15);
        text(150,205.5,append(int2str(TrDist(4)),'%'),'Color','y','FontSize',15);
        text(150,260.5,append(int2str(TrDist(5)),'%'),'Color','y','FontSize',15);
        text(150,315,append(int2str(TrDist(6)),'%'),'Color','y','FontSize',15);
        end
    end
        
    if contains(FigTitle,'Uni')||contains(FigTitle,'Bi')
        if contains(FigTitle,'p1')
        plot(50,43,'*','Color','r','LineWidth',6);
        text(60,23,'p1','Color','r','FontSize',15);
        elseif contains(FigTitle,'p2')
        plot(50,69,'*','Color','r','LineWidth',6);
        text(60,49,'p2','Color','r','FontSize',15);
        elseif contains(FigTitle,'p3')
        plot(50,95,'*','Color','r','LineWidth',6);
        text(60,75,'p3','Color','r','FontSize',15);
        end
    elseif contains(FigTitle,'PUN')
        if contains(FigTitle,'p1')
        plot(50,19.5,'*','Color','r','LineWidth',6);
        text(60,0,'p1','Color','r','FontSize',15);
        elseif contains(FigTitle,'p2')
        plot(50,46,'*','Color','r','LineWidth',6);
        text(60,26,'p2','Color','r','FontSize',15);
        elseif contains(FigTitle,'p3')
        plot(50,73,'*','Color','r','LineWidth',6);
        text(60,53,'p3','Color','r','FontSize',15);
        elseif contains(FigTitle,'p4')
        plot(50,100.5,'*','Color','r','LineWidth',6);
        text(60,80.5,'p4','Color','r','FontSize',15);
        elseif contains(FigTitle,'p5')
        plot(50,128,'*','Color','r','LineWidth',6);
        text(60,108,'p5','Color','r','FontSize',15);
        end
     elseif contains(FigTitle,'Chan')
        if contains(FigTitle,'p1')
        plot(50,17.5,'*','Color','r','LineWidth',6);
        text(60,0,'p1','Color','r','FontSize',15);
        elseif contains(FigTitle,'p2')
        plot(50,43,'*','Color','r','LineWidth',6);
        text(60,23,'p2','Color','r','FontSize',15);
        elseif contains(FigTitle,'p3')
        plot(50,68,'*','Color','r','LineWidth',6);
        text(60,48,'p3','Color','r','FontSize',15);
        elseif contains(FigTitle,'p4')
        plot(50,95.5,'*','Color','r','LineWidth',6);
        text(60,75.5,'p4','Color','r','FontSize',15);
        elseif contains(FigTitle,'p5')
        plot(50,123,'*','Color','r','LineWidth',6);
        text(60,103,'p5','Color','r','FontSize',15);
        end
    end
    set(gca,'TickLength',[0,0]);
    set(gca,'XTick',[],'YTick',[]);
end


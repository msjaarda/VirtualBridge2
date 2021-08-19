function [FigNum] = VBTriPlot(xdata,ydata,PDets,Title,Type,FigNum,FigTitle)
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

Len = 0;
for m = 1:3
    
    subplot(1,3,m)
    hold on
    for i = 1:size(ydata,2)
        plot(xdata{m,i},ydata{m,i},'-s','Color','k','MarkerEdgeColor',PDets.MEC{m}(i,:),'MarkerFaceColor',PDets.MFC{m}(i,:),'MarkerSize',PDets.MS{m}(i),'DisplayName',PDets.DN{m,i})
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
end

end

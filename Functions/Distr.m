function [Distr_Excel] = Distr(PDC,TrName,TrTyps,TrAxPerGr,TrTypPri,Location,Year,plotflag)

FaceAlpa = 0.7;

%[TrDistr, P2, A1, B1, a1, b1, A2, B2, a2, b2, p, mu1, mu2, sig1, sig2] = deal(zeros(length(TrTyps),1));
[TrDistr, P2, A1, B1, a1, b1, A2, B2, a2, b2, a3, b3, p1, p2, mu1, mu2, mu3, sig1, sig2, sig3] = deal(zeros(length(TrTyps),1));


pdf_betamixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
                       p*betapdf(x,mu1,sigma1) + (1-p)*betapdf(x,mu2,sigma2);
pdf_betamixture3 = @(x,p1,p2,mu1,mu2,mu3,sigma1,sigma2,sigma3) ...
                       p1*betapdf(x,mu1,sigma1) + p2*betapdf(x,mu2,sigma2) + (1-p1-p2)*betapdf(x,mu3,sigma3);
                   
%Three = ones(length(TrTyps),1);                   
Three = zeros(length(TrTyps),1);  
if plotflag == 1
    figure('units','normalized','outerposition',[0 0 1 1])
    color = cell(length(TrTyps),1);
end

for i = 1:13%length(TrTyps)
    
    TrTyp = TrTyps(i);
    
    if plotflag == 1
        subplot(4,4,i);
        color{i} = ((13-i)/(13))*(200/255)*[1, 1, 1];
    end
    
% Reinstate in the future... but for now we need to emulate just bimodal
%     try    
%         [x, P2(i), A1(i), B1(i), a1(i), b1(i), A2(i), B2(i), a2(i), b2(i), a3(i), b3(i), p1(i), p2(i), mu1(i), mu2(i), mu3(i), sig1(i), sig2(i), sig3(i)] = PdfweightThree(PDC,TrTyp);
%     catch
%         Three(i) = 0;
%         [x, P2(i), A1(i), B1(i), a1(i), b1(i), A2(i), B2(i), a2(i), b2(i), p(i), mu1(i), mu2(i), sig1(i), sig2(i)] = Pdfweight(PDC,TrTyp);
%     end
    
    [x, P2(i), A1(i), B1(i), a1(i), b1(i), A2(i), B2(i), a2(i), b2(i), p(i), mu1(i), mu2(i), sig1(i), sig2(i)] = Pdfweight(PDC,TrTyp);
    
    if plotflag == 1
        histogram(x,'Normalization','pdf','BinWidth',5,'EdgeColor','none','FaceColor',color{i},'FaceAlpha',FaceAlpa);
    end
    
    xgrid = linspace(0.7*min(x),1.1*max(x),200);
    
    if Three(i) == 0
        pdfgridx = pdf_betamixture(xgrid/1000,1-P2(i)/100,a1(i),a2(i),b1(i),b2(i));
    else
        pdfgridx = pdf_betamixture3(xgrid/1000,p1(i),p2(i),a1(i),a2(i),a3(i),b1(i),b2(i),b3(i));
    end
    
    if plotflag == 1
        hold on
        plot(xgrid,pdfgridx/1000,'-','LineWidth',2)
    
        title(['Type ' TrName{i}])
        xlabel('Weight (kN)')
        ylabel('PDF')
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        xlim([0 700])
    end
    
    TrDistr(i) = sum(PDC.CLASS == TrTyp)/sum(PDC.CLASS > 0);

end

if plotflag == 1
    sgtitle([Location ' ' num2str(Year) ' Truck Weight Summary']);

    pdfgridr = zeros(1,300);
    xgrid = linspace(0,600,300);
    for i = 1:length(TrTyps)
        if Three(i) == 0
            pdfgrid = pdf_betamixture(xgrid/1000,1-P2(i)/100,a1(i),a2(i),b1(i),b2(i));
        else
            pdfgrid = pdf_betamixture3(xgrid/1000,p1(i),p2(i),a1(i),a2(i),a3(i),b1(i),b2(i),b3(i));
        end
        pdfgridr = pdfgridr + pdfgrid*TrDistr(i);
    end


    % Plot each one plus all weights together (x2), plus TrDistr (bar/pie)

    % Pie chart for distribution
    subplot(4,4,16)
    
    if length(TrTyps) == 12

    Labels = {sprintf('%s', TrName{1}) sprintf('%s', TrName{2})...
        sprintf('%s', TrName{3}) sprintf('')...
        sprintf('') sprintf('')...
        sprintf('%s', TrName{7}) sprintf('')...
        sprintf('') sprintf('')...
        sprintf('%s', TrName{11}) sprintf('%s', TrName{12})};
    else
    
    Labels = {sprintf('%s', TrName{1}) sprintf('%s', TrName{2})...
        sprintf('%s', TrName{3}) sprintf('')...
        sprintf('') sprintf('')...
        sprintf('%s', TrName{7}) sprintf('')...
        sprintf('') sprintf('')...
        sprintf('%s', TrName{11}) sprintf('%s', TrName{12})...
        sprintf('%s', TrName{13})};
    end

    h = pie(TrDistr(1:13), Labels);
    if length(TrTyps) == 12

    colormap([color{1}; color{2}; color{3}; color{4}; color{5}; color{6}; color{7}; color{8};...
        color{9}; color{10}; color{11}; color{12}]);
    else
        colormap([color{1}; color{2}; color{3}; color{4}; color{5}; color{6}; color{7}; color{8};...
            color{9}; color{10}; color{11}; color{12}; color{13}]);
    end

    set(findobj(h, '-property', 'FaceAlpha'), 'FaceAlpha', FaceAlpa);
    %title('Distribution of Types')
    
    subplot(4,4,[14,15])

    count = cell(length(TrTyps),1);
    scaling = 2.5;

    for i = 1:13%length(TrTyps)
        x = PDC.GW_TOT(PDC.CLASS == TrTyps(i))/102;
        count{i} = histcounts(x,'BinEdges',0:scaling:500);
    end
    % Add Unclassified
    x = PDC.GW_TOT(PDC.CLASS == 0)/102;
    %count{length(TrTyps)+1} = histcounts(x,'BinEdges',0:scaling:500);
    count{13+1} = histcounts(x,'BinEdges',0:scaling:500);
    
    y = cell2mat(count);
    Su = sum(sum(y));
    y = y/(Su*scaling);

    % figure(2)

    h = bar(y',1,'stacked','FaceAlpha', FaceAlpa,'EdgeColor','none');

    for i = 1:13%length(TrTyps)
        h(i).FaceColor = color{i};
    end
    % Add Unclassified
    h(13+1).FaceColor = 'w';
    h(13+1).EdgeColor = 'k';

    %h(length(TrTyps)+1).FaceColor = 'w';
    %h(length(TrTyps)+1).EdgeColor = 'k';
    
    xticks = get(gca,'xtick');
    scaling  = 2.5;
    newlabels = arrayfun(@(x) sprintf('%.0f', scaling * x), xticks, 'un', 0);
    set(gca,'xticklabel',newlabels);

    title('All Trucks > 6t')
    xlabel('Weight (kN)')
    ylabel('PDF')
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    xlim([0 500/scaling])
    hold on
    plot(xgrid/scaling,pdfgridr/1000,'-','LineWidth',2,'Color','r')

end

if length(TrTyps) == 12
PlatPct = 0.1*ones(12,1);
else
    PlatPct = 0.1*ones(length(TrTyps),1);
end

if plotflag == 1
    figure('units','normalized','outerposition',[0 0 1 .3])
    color = cell(length(TrTyps),1);
end

for i = 14:length(TrTyps)
    
    TrTyp = TrTyps(i);
    
    if plotflag == 1
        subplot(1,4,i-13);
         color{i} = ((i-13)/(4))*(200/255)*[1, 1, 1];
    end
    
% Reinstate in the future... but for now we need to emulate just bimodal
%     try    
%         [x, P2(i), A1(i), B1(i), a1(i), b1(i), A2(i), B2(i), a2(i), b2(i), a3(i), b3(i), p1(i), p2(i), mu1(i), mu2(i), mu3(i), sig1(i), sig2(i), sig3(i)] = PdfweightThree(PDC,TrTyp);
%     catch
%         Three(i) = 0;
%         [x, P2(i), A1(i), B1(i), a1(i), b1(i), A2(i), B2(i), a2(i), b2(i), p(i), mu1(i), mu2(i), sig1(i), sig2(i)] = Pdfweight(PDC,TrTyp);
%     end
    
    [x, P2(i), A1(i), B1(i), a1(i), b1(i), A2(i), B2(i), a2(i), b2(i), p(i), mu1(i), mu2(i), sig1(i), sig2(i)] = Pdfweight(PDC,TrTyp);
    
    if plotflag == 1
        histogram(x,'Normalization','pdf','BinWidth',5,'EdgeColor','none','FaceColor',color{i},'FaceAlpha',FaceAlpa);
    end
    
    xgrid = linspace(0.7*min(x),1.1*max(x),200);
    
    if Three(i) == 0
        pdfgridx = pdf_betamixture(xgrid/1000,1-P2(i)/100,a1(i),a2(i),b1(i),b2(i));
    else
        pdfgridx = pdf_betamixture3(xgrid/1000,p1(i),p2(i),a1(i),a2(i),a3(i),b1(i),b2(i),b3(i));
    end
    
    if plotflag == 1
        hold on
        plot(xgrid,pdfgridx/1000,'-','LineWidth',2)
    
        title(['Type ' TrName{i}])
        xlabel('Weight (kN)')
        ylabel('PDF')
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        xlim([0 1000])
    end
    
    TrDistr(i) = sum(PDC.CLASS == TrTyp)/sum(PDC.CLASS > 0);

end


% Excel-ready output
Distr_Excel = table(TrName,TrAxPerGr,TrTypPri,TrDistr,PlatPct,P2,A1,B1,a1,b1,A2,B2,a2,b2);

end



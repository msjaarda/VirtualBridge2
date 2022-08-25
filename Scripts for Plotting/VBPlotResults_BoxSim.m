clear, clc, close all

% Get VBResults using VBOutput2Struct
[VBResults] = VBOutput2Struct('BoxSIM');

% Format is:
% VBResults.AQ.Type.SubType.Width.Layout.Support.Trans.AE.Traffic
% VBResults.x.Type.SubType.Width.Layout.Support.Trans.AE.Traffic



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FName = 'Box Girder Bridge Bidirectional 2L'; FigNum = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General Info
Type = 'Box';         SubType = 'Stand';
Width = 'Wid12';     Layout = 'Bi';
Support = 'Simp'; Trans = 'p0';
% Subplot Info (m 1:3 is the subplots)
AE(1) = "Mn"; AE(2) = "Mp"; AE(3) = "V";
% Series Info (i 1:4 are the series')
Traffic(1) = "Denges2019";   Traffic(2) = "Oberburen2019";
Traffic(3) = "Ceneri2019"; Traffic(4) = "Gotthard2019";
% Other deets
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
% Set colours
D = linspecer(24);      

% Prepare plot data for VBTriPlot
for m = 1:3
    % All ClassTypes
    for i = 1:length(Traffic)
        % Get ydata
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        % xdata should always be equal... give a warning if one isn't?
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        % DisplayName
        PDets.DN{m,i} = Traffic(i);
        % Marker Size, Edge Color, Face Color
        PDets.MEC{m}(i,:) = [0 0 0]; PDets.MFC{m}(i,:) = D(i*6-3,:);
        PDets.MS{m}(i) = 5;
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum,FName);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FName = 'Box Girder Bridge PUN 3L';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changes
Layout = 'PUN';

clear ydata, clear xdata, clear PDets
for m = 1:3
    for i = 1:length(Traffic)
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        PDets.DN{m,i} = Traffic(i);
        PDets.MEC{m}(i,:) = [0 0 0]; PDets.MFC{m}(i,:) = D(i*6-3,:);
        PDets.MS{m}(i) = 5;
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum,FName);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FName = 'Box Girder Bridge Bidirectional 4L';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changes
Width = 'Wid18';     Layout = 'Bi';
Traffic(4) = [];

clear ydata, clear xdata, clear PDets
for m = 1:3
    for i = 1:length(Traffic)
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        PDets.DN{m,i} = Traffic(i);
        PDets.MEC{m}(i,:) = [0 0 0]; PDets.MFC{m}(i,:) = D(i*6-3,:);
        PDets.MS{m}(i) = 5;
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum,FName);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FName = 'Box Girder Bridge PUN 5L';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changes
Layout = 'PUN';

clear ydata, clear xdata, clear PDets
for m = 1:3
    for i = 1:length(Traffic)
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        PDets.DN{m,i} = Traffic(i);
        PDets.MEC{m}(i,:) = [0 0 0]; PDets.MFC{m}(i,:) = D(i*6-3,:);
        PDets.MS{m}(i) = 5;
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum,FName);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FName = 'Box Girder Bridge Unidirectional 2L';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changes
Layout = 'Uni'; Width = 'Wid12';
Traffic(4) = "Gotthard2019";


clear ydata, clear xdata, clear PDets
for m = 1:3
    for i = 1:length(Traffic)
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        PDets.DN{m,i} = Traffic(i);
        PDets.MEC{m}(i,:) = [0 0 0]; PDets.MFC{m}(i,:) = D(i*6-3,:);
        PDets.MS{m}(i) = 5;
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum,FName);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FName = 'Box Girder Bridge Unidirectional 2L 2008';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changes
Traffic(1) = "Schafisheim2008";
Traffic(2) = "Felsenau2008";
Traffic(3) = "Oberburen2008";
Traffic(4) = "Mattstetten2008";
Traffic(5) = "Ceneri2008";
Traffic(6) = "Trubbach2008";
Traffic(7) = "Gotthard2008";

for m = 1:3
    for i = 1:length(Traffic)
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        PDets.DN{m,i} = Traffic(i);
        PDets.MEC{m}(i,:) = [0 0 0]; PDets.MFC{m}(i,:) = D(i*3,:);
        PDets.MS{m}(i) = 5;
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum,FName);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FName = 'Box Girder Bridge Unidirectional 2L 2018';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changes
clear Traffic
Traffic(1) = "Denges2018";
Traffic(3) = "Oberburen2018";
Traffic(2) = "StMaurice2018";
Traffic(5) = "Ceneri2018";
Traffic(4) = "Trubbach2018";
Traffic(6) = "Gotthard2018";

clear ydata, clear xdata, clear PDets
for m = 1:3
    for i = 1:length(Traffic)
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        PDets.DN{m,i} = Traffic(i);
        PDets.MEC{m}(i,:) = [0 0 0]; PDets.MS{m}(i) = 5;
        if i == 6
            PDets.MFC{m}(i,:) = D(i*3+3,:);
        else
            PDets.MFC{m}(i,:) = D(i*3,:);
        end
        
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum,FName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FName = 'Box Girder Bridge Unidirectional 2L Ceneri';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Changes
clear Traffic
Traffic(1) = "Ceneri2008";
Traffic(2) = "Ceneri2009";
Traffic(3) = "Ceneri2010";
Traffic(4) = "Ceneri2011";
Traffic(5) = "Ceneri2012";
Traffic(6) = "Ceneri2013";
Traffic(7) = "Ceneri2014";
Traffic(8) = "Ceneri2015";
Traffic(9) = "Ceneri2016";
Traffic(10) = "Ceneri2017";
Traffic(11) = "Ceneri2018";
Traffic(12) = "Ceneri2019";
Traffic(13) = "Ceneri2020";

clear ydata, clear xdata, clear PDets
for m = 1:3
    for i = 1:length(Traffic)
        ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic(i));
        PDets.DN{m,i} = Traffic(i);
        PDets.MEC{m}(i,:) = [0 0 0]; PDets.MS{m}(i) = 5;
        % Shades of grey
        PDets.MFC{m}(i,:) = [1 1 1]*i/length(Traffic);
    end
end

% Plot
FigNum = VBTriPlot(xdata,ydata,PDets,Title,'WIM',FigNum,FName);
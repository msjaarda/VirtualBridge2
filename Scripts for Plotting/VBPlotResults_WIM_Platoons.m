clear, clc%, close all

% Global Definitions
FigNum = 0; D = linspecer(30);      % Set colours
SubType = 'Stand'; Support = 'Simp'; Trans = 'p0'; Layout = 'Uni';
% Subplot Info (m)
AE(1) = "Mn"; AE(2) = "Mp"; AE(3) = "V";
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
% Get VBResults using VBOutput2Struct
[VBResultsPlat] = VBOutput2Struct('Platooningno35');

%load('VBResultsP.mat') %PlatooningFull (Flowing All Stations)

%load('VBResultsPF.mat') %Platooning Just 408
%load('VBResultsPFNo35.mat') %Platooningno35 Just 408 (

load('VBResults.mat')
%load('VBResultsPFTwin.mat')
% Class
Class = 'Class';

Export = true;

%Type = 'Twin'; 
Type = 'Box';

%Width = 'Wid9';  
Width = 'Wid12';

% Series Info (i)
%Traffic = fieldnames(VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).Mn);
Traffic{1} = 'Ceneri408';
% Platooning Info
PSize = 'Size4'; PRate = 'Rate4'; 
PFolDist{1} = 'FolDist50'; PFolDist{2} = 'FolDist25'; PFolDist{3} = 'FolDist50'; PFolDist{4} = 'FolDist75'; PFolDist{5} = 'FolDist100';

FName = [Type ' Girder Uni2L ' PSize(end) ' Tr Platoons @ ' PRate(end) '0 % (' Class ')'];
for m = 1:3
    for i = 1:5
        % Get ydata and xdata
        if i == 1
            ydata{m,i} = VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{1}).(Class)(1:2:end);
            PDets.MFC{m}(i,:) = D(i*5,:);
            PDets.DN{m,i} = Traffic{1};
        elseif i == 2
            ydata{m,i} = VBResultsPlat.P.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{1}).(PSize).(PRate).(PFolDist{i}).(Class);
            PDets.MFC{m}(i,:) = D(i*5,:);
            PDets.DN{m,i} = '2.5m Fol Dist';
        elseif i == 3
            ydata{m,i} = VBResultsPlat.P.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{1}).(PSize).(PRate).(PFolDist{i}).(Class);
            PDets.MFC{m}(i,:) = D(i*5,:);
            PDets.DN{m,i} = '5.0m Fol Dist';
        elseif i == 4
            ydata{m,i} = VBResultsPlat.P.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{1}).(PSize).(PRate).(PFolDist{i}).(Class);
            PDets.MFC{m}(i,:) = D(i*5,:);
            PDets.DN{m,i} = '7.5m Fol Dist';
        else
            ydata{m,i} = VBResultsPlat.P.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{1}).(PSize).(PRate).(PFolDist{i}).(Class);
            PDets.MFC{m}(i,:) = D(i*5,:);
            PDets.DN{m,i} = '10m Fol Dist';
        end
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{1}).(Class)(1:2:end);
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5;
        PDets.MEC{m}(i,:) = [0 0 0];
    end
end

% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'WIM',FigNum,FName);

if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end

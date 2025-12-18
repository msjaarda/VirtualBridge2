clear, clc, close all

% Global Definitions
FigNum = 0; D = linspecer(30);      % Set colours
SubType = 'Stand'; Support = 'Simp'; Trans = 'p0'; Layout = 'Uni';
% Subplot Info (m)
AE(1) = "Mn"; AE(2) = "Mp"; AE(3) = "V";
Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
% Get VBResults using VBOutput2Struct
%[VBResultsPlat] = VBOutput2Struct('Platooningno35');

%load('VBResultsP.mat') %PlatooningFull (Flowing All Stations)

load('VBResultsPF.mat') %Platooning Just 408
%load('VBResultsPFNo35.mat') %Platooningno35 Just 408 (

load('VBResults.mat')
%load('VBResultsPFTwin.mat')
% Class
Class = 'Class';

Export = false;

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
            PDets.DN{m,i} = ['Ceneri Base'];
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
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'PLAT',FigNum,FName);

if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end



% Try with Jammed results

% Need to find a way to get ESIA for each of these cases...
% VBResultsPJ_V2 includes the basemax result (but we don't have ESIA!)
% shouldn't be too hard. use 

% Read Input File
FName = 'Input/VBWIMqInput_WIMPlatoonFull.xlsx';
BaseData = VBReadInputFile(FName);

% Update analysis data for current row of BaseData
[Num,Lane,ILData,~,~] = VBUpdateData(BaseData(1,:));

E = VBGetECode(ILData,1);
for i = 1:4
    ESIA.V(i) = E(i).SIA.Total;
end
for i = 1:4
    ESIA.Mp(i) = E(i+4).SIA.Total;
end
for i = 1:4
    ESIA.Mn(i) = E(i+8).SIA.Total;
end
% 
% ESIA.V = ESIA.Total(1:4);
% ESIA.Mp = ESIA.Total(5:8);
% ESIA.Mn = ESIA.Total(9:12);

% Get VBResults using VBOutput2Struct
%[VBResultsPlat] = VBOutput2Struct('Platooningno35');

%load('VBResultsP.mat') %PlatooningFull (Flowing All Stations)

load('VBResultsPJ_V2.mat') %Platooning Just 408
%load('VBResultsPFNo35.mat') %Platooningno35 Just 408 (

load('VBResults.mat')
%load('VBResultsPFTwin.mat')
% Class
Class = 'Class';


%Type = 'Twin'; 
Type = 'Box';

%Width = 'Wid9';  
Width = 'Wid12';

% Series Info (i)
%Traffic = fieldnames(VBResults.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).Mn);
%Traffic{1} = 'Ceneri2017';
Trafficx{1} = 'Ceneri2017';
% Platooning Info
PSize = 'Size4'; PRate = 'Rate4'; 
PFolDist{1} = 'FolDist50'; PFolDist{2} = 'FolDist25'; PFolDist{3} = 'FolDist50'; PFolDist{4} = 'FolDist75'; PFolDist{5} = 'FolDist100';


for m = 1:3
    for i = 1:5
        % Get ydata and xdata
        if i == 1
            ydata{m,i} = VBResultsPlat.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Trafficx{1}).(PSize).(PRate).(PFolDist{i}).(Class)./ESIA.(AE(m))';
            PDets.MFC{m}(i,:) = D(i*5,:);
            PDets.DN{m,i} = ['Ceneri Base'];
        elseif i == 2
                          % Normalize, see below... BUT we have to use the
                          % Jammed BaseData... not the flowing as here
            ydata{m,i} = ydata{m,1}.*VBResultsPlat.P.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Trafficx{1}).(PSize).(PRate).(PFolDist{i}).(Class)';
            PDets.MFC{m}(i,:) = D(i*5,:);
            PDets.DN{m,i} = '2.5m Fol Dist';
        elseif i == 3
            ydata{m,i} = ydata{m,1}.*VBResultsPlat.P.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Trafficx{1}).(PSize).(PRate).(PFolDist{i}).(Class)';
            PDets.MFC{m}(i,:) = D(i*5,:);
            PDets.DN{m,i} = '5.0m Fol Dist';
        elseif i == 4
            ydata{m,i} = ydata{m,1}.*VBResultsPlat.P.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Trafficx{1}).(PSize).(PRate).(PFolDist{i}).(Class)';
            PDets.MFC{m}(i,:) = D(i*5,:);
            PDets.DN{m,i} = '7.5m Fol Dist';
        else
            ydata{m,i} = ydata{m,1}.*VBResultsPlat.P.AQ.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Trafficx{1}).(PSize).(PRate).(PFolDist{i}).(Class)';
            PDets.MFC{m}(i,:) = D(i*5,:);
            PDets.DN{m,i} = '10m Fol Dist';
        end
        xdata{m,i} = VBResults.x.(Type).(SubType).(Width).(Layout).(Support).(Trans).(AE(m)).(Traffic{1}).(Class)(1:2:end);
        % Marker Size, Edge Color, Face Color
        PDets.MS{m}(i) = 5;
        PDets.MEC{m}(i,:) = [0 0 0];
    end
end

FName = [Type ' Girder Uni2L ' PSize(end) ' Tr Platoons @ ' PRate(end) '0 % (' Class ') v2 Jammed'];


% Plot
FigNum = VBTriPlotPro(xdata,ydata,PDets,Title,'PLAT',FigNum,FName);

if Export
    exportgraphics(gcf,FName + ".jpg",'Resolution',600);
end

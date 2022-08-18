% ------------------------------------------------------------------------
%                              VBInputGen
% ------------------------------------------------------------------------
% Summarize traffic characteristics from WIM or VirtualWIM
%       - Output is a MATSimInput spreadsheet to be used with MATSim
%       - Add classification of heavy vehicles (use ClassifyOW)
%       - Take a look at Classify to see about eliminating 11bis

% Find a way to include statistics on charts (best fit params, schematics)

% Check certain hard coded values including WitAx and Allos, as well as
% FolDist by plotting following distances less than a certain amount in the
% WIM files (consider more than just classified vehicles)

tic, clear, clc, hold off, close all, format long g, rng('shuffle')

% ----- INPUT -----

% Year, #, Station Name, string, save and plot toggles
%Year = 2017; SName = 'Ger9625'; Save = 0; PlotFits = 1;
%Year = 2016; SName = 'Ast11450'; Save = 0; PlotFits = 1;
Year = 2019; SName = ''; PlotFits = 1; Save = 0;
%Year = 2019;
% ----- ENDIN -----

% try
% load(['PrunedS1 WIM/',SName,'/',SName,'_',num2str(Year),'.mat']);
% %PD = 1; % REMOVE, just for silly dependancy analysis problem
% Station = unique(PD.ZST);
% catch
% if necessary
%load(['WIM/',SName,'.mat']);
%load('Grauholz_2017.mat')
%load('C:\Users\mjsja\Desktop\SwissTraffic2\WIMVaud\Vufflens.mat')
%load('C:\Users\mjsja\Desktop\SwissTraffic2\WIM\9625.mat')
load('C:\Users\mjsja\Desktop\SwissTraffic2\WIM\408.mat')
PD1 = PDs;
load('C:\Users\mjsja\Desktop\SwissTraffic2\WIM\409.mat')
PDs = [PD1; PDs];
%PD(PD.GW_TOT == 0,:) = [];
%PD = PD(year(PD.DTS) == 2019,:);
PDs = PDs(year(PDs.DTS) == Year,:);
%PD = 1; % REMOVE, just for silly dependancy analysis problem
Station = unique(PDs.SITE);

% end

% % Lookup from weeklycount excel file
% InputFile = 'Misc/WorkdayCounts.xlsx'; CountingData = readtable(InputFile);
% % Get Total Workday Vehicles (x2 dir) from counting station data
% try
%     TotalVeh = CountingData{Year-1999,SName};
% catch
%     TotalVeh = 100000;
% end

% Let the Classify function add the .CLASS column to PD
%PDC = Classify(PD);
PDC = PDs;

% % Initialize variables specific to each station... must go RAW for WeightFlux
% [InitialNum, NumDays, TotW] = deal(zeros(length(Station),1));
% 
% RawAvail = isfile(['Raw WIM/' SName '/' SName num2str(Station(1)) '_' num2str(Year) '.mat']);
% 
% if RawAvail
%     for i = 1:length(Station)
%         InputFile = ['Raw WIM/' SName '/' SName num2str(Station(i)) '_' num2str(Year) '.mat'];
%         % Load InputFile
%         load(InputFile)
% 
%         % Get initial number of vehicles
%         InitialNum(i) = height(RD);
%         % Get number of days station is open
%         NumDays(i) = length(unique(RD.JJJJMMTT));
%         % Get total weight of all axles
%         y = RD{:,14:23} > 30 &  RD{:,14:23} < 35000;
%         AW = RD{:,14:23};
%         TotW(i) = sum(AW(y),'all');
%     end
% 
%     % Get yearly weight flux (kN)
%     PerYearkN = ((TotW/102)./NumDays)*365;
% else
%     for i = 1:length(Station)
%         % Get number of days station is open
%         NumDays(i) = length(unique(PD.JJJJMMTT(PD.ZST == Station(i)))); 
%     end
% end

% Stage 2 Custom Pruning (two more of our own)
% 1. Disqualification by weight (under 6 tonnes)
PDC = PDC(PDC.GW_TOT > 6000,:);
% 2. Disqualification by Swiss10 Class (exclude 2,3,4,6)
% Note all years have proper Sw10 data... therefore exclude, not include
try
if sum(PDC.CS == 2 | PDC.CS == 3 | PDC.CS == 4 | PDC.CS == 6) < 0.7*height(PDC)
    PDC(PDC.CS == 2 | PDC.CS == 3 | PDC.CS == 4 | PDC.CS == 6,:) = [];
end
catch
end
%PDC = PDC(PDC.CS == 1 | PDC.CS == 5 | PDC.CS == 8 | PDC.CS == 9 | PDC.CS == 10,:);
% Get rid of overweight classes and class 11bis
PDC.CLASS(PDC.CLASS > 39 & PDC.CLASS < 50) = 0;
PDC.CLASS(PDC.CLASS == 119) = 0;

% % Get Daytype
% PDC = Daytype(PDC,Year);
% PerWorkdayTraf = sum(PDC.GW_TOT(PDC.Daytype < 6))/sum(PDC.GW_TOT);
% 
% % Find average weight of classified vehicles
% MTW = mean(PDC.GW_TOT(PDC.CLASS > 0))/102;
% if RawAvail
%     % Full weight flux
%     WF = sum(PerYearkN);
%     % Number of trucks to get weight flux
%     NV = WF/MTW;
% else
%     NV = height(PDC)*length(Station)*365/sum(NumDays);
% end
% 
% % Number of vehicles on weekdays
% NVW = NV*PerWorkdayTraf;
% % Number of vehicles per weekday
% NVpW = NVW/(365*5/7);
% % Truck Rate
% TR = NVpW/TotalVeh;
% % Classified weight flux percentage
% ClassWFp = 100*sum(PDC.GW_TOT(PDC.CLASS > 0))/sum(TotW);
% if RawAvail
%     % Percent of WIM used for veh stats
%     TotalR = 100*sum(PDC.CLASS > 0)/sum(InitialNum);
% end

% % Output basic info from analysis
% for i = 1:length(Station)
%     fprintf('Total Days Open: %i, %i\n',NumDays(i), Station(i))
% end
% if RawAvail
%     fprintf('Yearly Weight (kN): %i\n',sum(PerYearkN))
%     fprintf('Classified Weight Flux Percentage: %.2f\n',ClassWFp)
%     fprintf('Percent of WIM Logs used for Veh Stats: %.2f\n\n',TotalR)
% end
% 
% fprintf('Truck Rate (Percentage): %.2f\n\n',100*TR)

% Get Input file!

% 7 = r, 8 = a
% TrName{1} = '11'; TrName{2} = '12'; TrName{3} = '22';
% TrName{4} = '111'; TrName{5} = '1111r'; TrName{6} = '112r'; 
% TrName{7} = '1211r'; TrName{8} = '122'; TrName{9} = '1112r';    
% TrName{10} = '112a'; TrName{11} = '113a'; TrName{12} = '123a';   
% TrName = TrName';

TrName{1} = '11'; TrName{2} = '12'; TrName{3} = '22'; TrName{4} = '23';
TrName{5} = '111'; TrName{6} = '1111r'; TrName{7} = '112r'; 
TrName{8} = '1211r'; TrName{9} = '122'; TrName{10} = '1112r';    
TrName{11} = '112a'; TrName{12} = '113a'; TrName{13} = '123a';   
TrName = TrName';

TrTyps = [11; 12; 22; 23; 111; 11117; 1127; 12117; 122; 11127; 1128; 1138; 1238];
TrAxPerGr = [11; 12; 22; 23; 111; 1111; 112; 1211; 122; 1112; 112; 113; 123];
TrTypPri = [21; 21; 21; 21; 321; 2341; 231; 2341; 231; 2341; 321; 321; 321];

% TrTyps = [11; 12; 22; 111; 11117; 1127; 12117; 122; 11127; 1128; 1138; 1238];
% TrAxPerGr = [11; 12; 22; 111; 1111; 112; 1211; 122; 1112; 112; 113; 123];
% TrTypPri = [21; 21; 21; 321; 2341; 231; 2341; 231; 2341; 321; 321; 321];

BetAx_Excel = BetAx(PDC,TrName,TrTyps,TrAxPerGr,TrTypPri,SName,Year,PlotFits);
%exportgraphics(gcf,SName + "BetAx" + ".jpg",'Resolution',600);
LinFit_Excel = LinFit(PDC,TrName,TrTyps,TrAxPerGr,TrTypPri,SName,Year,PlotFits);
%exportgraphics(gcf,SName + "LinFit" + ".jpg",'Resolution',600);
Distr_Excel = Distr(PDC,TrName,TrTyps,TrAxPerGr,TrTypPri,SName,Year,PlotFits);
%exportgraphics(gcf,SName + "Distr" + ".jpg",'Resolution',600);

% Summarize Weight by Axles
[STaTr,AllAx] = AxleStats(PDC,TrAxPerGr,TrTyps,SName,Year,PlotFits);
%exportgraphics(gcf,SName + "AxStats" + ".jpg",'Resolution',600);

% Write to Input File
if Save == 1
    % Copy blank file and save as new name
    filename = strcat('Input/MATSimInput', SName, num2str(Year),'_',datestr(now,'mmmdd-yy HHMM'), '.xlsx');
    copyfile('Input\Input_blank.xlsx',filename)

    writetable(BetAx_Excel,filename,'Sheet','TrBetAx','Range','A1');
    writetable(Distr_Excel,filename,'Sheet','TrDistr','Range','A1');
    writetable(LinFit_Excel,filename,'Sheet','TrLinFit','Range','A1');
    
    C = {TR};
    writecell(C,filename,'Sheet','BaseData','Range','C2');
    
end

%LaneDistBr(PDC,TrTyps,TrAxPerGr,Station)

% If importing new TrLib to TrLib

load('TrLib.mat')
TrLib.([SName,num2str(Year)]).TrDistr = Distr_Excel;
TrLib.([SName,num2str(Year)]).TrLinFit = LinFit_Excel;
TrLib.([SName,num2str(Year)]).TrAllo = TrLib.Ceneri2018.TrAllo;
TrLib.([SName,num2str(Year)]).TrBetAx = BetAx_Excel;
TrLib.([SName,num2str(Year)]).TrWitAx = TrLib.Ceneri2018.TrWitAx;
TrLib.([SName,num2str(Year)]).TrDistr.PlatPct = [];

TrLib = orderfields(TrLib);
save('Misc/TrLib.mat','TrLib')
%ExportOpen
% MANUALLY Save



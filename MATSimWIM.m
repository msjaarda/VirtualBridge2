% ------------------------------------------------------------------------
%                            MATSim2019forWIM
% ------------------------------------------------------------------------
% Run real traffic over a bridge to find maximum load effects
%       - To be used with MATSimInput spreadsheet for Inf Lines

% Initial commands
tic, clear, clc, close all, format long g, rng('shuffle'); st = now;

% Input File Name
InputFile = 'Input/MATSimInputx.xlsx';

% Read simulation data from Input File
[BaseData,LaneData,~,~] = ReadInputFile(InputFile); OverMax = [];

% Get Influence Line Details
[InfLanes,InfNames,UniqInf,UniqInfs,UniqInfi,Num.InfCases,Infx,Infv,IntInfv,MaxInfv] = GetInfLines(LaneData,BaseData);

% WIM File
SName = 'Ceneri'; Stage2Prune = false; ClassOnly = false; SpaceSaver = 100; Year = 2015;

% Load WIM File [OR VWIM File]
%load(['PrunedS1 WIM/',SName,'/',SName,'_',num2str(Year),'.mat']);
load('WIM_Jan14 1130.mat');

% Only needs to be done to WIM

% Add row for Class, Daytime, and Daycount
PDC = Classify(PD); PDC = Daytype(PDC,Year); clear('PD')

% We treat each station separately [edit Stations(1) to Stations(2)]
Stations = unique(PDC.ZST); Station = Stations(1);

PDCx = PDC(PDC.ZST == Station,:); clear('PDC')

if Stage2Prune
    PDCx = PruneWIM2(PDCx,0);
end

if ClassOnly
    PDCx(PDCx.CLASS == 0,:) = [];
end

% Convert PDC to AllTrAx
[PDCx, AllTrAx, TrLineUp] = WIMtoAllTrAx(PDCx,SpaceSaver);

% Round TrLineUp first row, move unrounded to fifth row
TrLineUp(:,5) = TrLineUp(:,1); TrLineUp(:,1) = round(TrLineUp(:,1));
        
for t = 1:Num.InfCases
    % Subject Influence Line to Truck Axle Stream
    [MaxLE,DLF,BrStInd,AxonBr,FirstAxInd,FirstAx] = GetMaxLE(AllTrAx,Infv,InfLanes,[80 20],BaseData.RunDyn,t,UniqInfs,UniqInfi);
    % Record Maximums
    OverMax = [OverMax; [t, Year, MaxLE, DLF, BrStInd, FirstAxInd, FirstAx]];
end

% [Delete vehicles involved in maximum case, and re-analyze]

% Convert Results to Table
OverMAXT = array2table(OverMax,'VariableNames',{'InfCase','Year','MaxLE','MaxDLF','MaxBrStInd','MaxFirstAxInd','MaxFirstAx'});

% Get simulation time
Time = GetSimTime();

% Get ESIA from function
ESIA = GetESia(IntInfv,MaxInfv,InfLanes,UniqInfs);

T = Apercu(PDCx,[Year ' ' SName ' Station ' Station],Infx,Infv,BrStInd,TrLineUp,MaxLE/ESIA,DLF);



% Optional Outputs 

% figure(2)
% scatter(LowYear:HighYear,OverMAXT.MaxLE)
% hold on
% plot([LowYear HighYear],[ESIA ESIA])
% hold on
% plot([LowYear HighYear],.72*[ESIA ESIA])
% ylim([0 14000])
% xlabel('Year')
% ylabel('Moment (kNm)')
% title('Ceneri Station 408 Maximums')

%fprintf('Percent of ESIA = %.2f\n\n',MaxLE/ESIA)

% ------------------------------------------------------------------------
%                            VBDet
% ------------------------------------------------------------------------
% Run deterministic vehicles over a bridge to find maximum load effects

% Initial commands
clear, clc, format long g, rng('shuffle'), close all; 

% Input Information --------------------

% Read Input File
BaseData = VBReadInputFile('VBWIMInputOFROUDetx.xlsx');

% % Manual BaseData Input
% BaseData = table;
% % Roadway Info
% BaseData.LaneDir = {'1,2'};
% % Influence Line Info
% BaseData.ILs = {'AGBBox.Mp.S60'};  BaseData.ILRes = 0.1;
% % Analysis Info
% BaseData.RunDyn = 0; % 1.3 added manually
% % Analysis Type
% BaseData.AnalysisType = "Det";
% % Apercu
% BaseData.Apercu = true;
% % Save and Folder
% BaseData.Save = 0;
% BaseData.Folder = '/AGB2002_VB';
% % Apercu Title
% BaseData.ApercuTitle = sprintf('%s','Deterministic Analysis');
% g = 1;

% Input Complete   ---------------------

for g = 1:height(BaseData)
    
    % Obtain Influence Line Info
    [Num,Lane,ILData,~,~] = VBUpdateData(BaseData(g,:));
    
    % Initialize
    OverMax = [];
    
    FName = 'VBTrial.mat'; % 'Det60t.mat'
    load(FName)
    
    InAxs = contains(PDC.Properties.VariableNames, 'AWT');
        
    % Apply Factors according to AGB 2002/005 (1.1 and 1.5)
    PDC{PDC.GW_TOT == 60000,InAxs} = 1.1*PDC{PDC.GW_TOT == 60000,InAxs};
    PDC{PDC.GW_TOT == 40000,InAxs} = 1.5*PDC{PDC.GW_TOT == 40000,InAxs};
    PDCx = PDC;
    
    % Convert PDC to AllTrAx
    [PDCx, AllTrAx, TrLineUp] = VBWIMtoAllTrAx(PDCx,0,Lane,BaseData.ILRes(g));
    
    % Round TrLineUp first row, move unrounded to fifth row
    TrLineUp(:,5) = TrLineUp(:,1); TrLineUp(:,1) = round(TrLineUp(:,1)/BaseData.ILRes(g));
    
    % TrLineUp [ 1: AllTrAxIndex  2: AxleValue  3: Truck#  4: LaneID  5: Station(m)]
    TrLineUp = array2table(TrLineUp,'VariableNames',{'ATAIndex','AxleValue','TrNum','LaneID','mStation'});
    
    for t = 1:Num.InfCases
        
        % Subject Influence Line to Truck Axle Stream
        [MaxLE,DLF,BrStInd,R] = VBGetMaxLE(AllTrAx,ILData(t).v,0);
        % Record Maximums
        % Add AGB 1.3 DLF
        OverMax = [OverMax; [t, 1.3*MaxLE, 1.3, BrStInd]];
        
    end
    
    if BaseData.Apercu(g)
        % Display Apercu
        T = VBApercuv2(PDCx,sprintf('%s %i','Deterministic Analysis',g),ILData(t),BrStInd,table2array(TrLineUp),1.3,Lane,BaseData.ILRes(g));
    end
    
    % Convert Results to Table
    OverMaxT = array2table(OverMax,'VariableNames',{'InfCase','MaxLE','DLF','BrStInd'});
    
    
    % Optional save of OutInfo (used for deterministic AGB matching)
    OutInfo.Name = datestr(now,'mmmdd-yy HHMMSS'); OutInfo.BaseData = BaseData(g,:);
    OutInfo.OverMax = OverMax; OutInfo.OverMaxT = OverMaxT;
    OutInfo.ILData = ILData;
    
    % Create folders where there are none
    if BaseData.Save(g) == 1
        CreateFolders(BaseData.Folder{g},BaseData.VWIM(g),BaseData.Apercu(g),BaseData.Save(g))
    end
    
    if BaseData.Save(g) == 1
        save(['Output' BaseData.Folder{g} '/' OutInfo.Name], 'OutInfo','-v7.3')
    end
    
end % g
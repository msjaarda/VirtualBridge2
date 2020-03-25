% ------------------------------------------------------------------------
%                              MATSim2019
% ------------------------------------------------------------------------
% Simulate traffic over a bridge to find maximum load effects
%       - To be used with Simple or Complete MATSimInput spreadsheet
%       - Platooning only supported in Lane 1

% Initializing commands
clear, clc, close all, format long g, rng('shuffle');

% Input File or Folder Name 
InputF = 'MATSimInputAGBSimple43.xlsx'; Folder_Name = '/AGB4p3'; %'/' b4

% Read in simulation data
[BaseData,LaneData,TrData,FolDist] = ReadInputFile(['Input/' InputF]);

% Each row of BaseData represents one analysis
for g = 1:height(BaseData)

% Update analysis data for current row of BaseData
[Num.Lanes,Lane,LaneData,TrData,FolDist] = UpdateData(BaseData(g,:),LaneData,TrData,FolDist);   
    
% Get key variables from imported data
[BatchSize,Num.Batches,FixVars,PlatPct] = GetKeyVars(BaseData(g,:),TrData.TrDistr);

% Get Influence Line Details
[Inf,Num.InfCases,Inf.x,Inf.v,ESIA] = GetInfLines(LaneData,BaseData(g,:),Num.Lanes);

% Define truck type specific properties
[Num.TrTyp,TrTyp] = GetTrPpties(TrData);

% Get Per Lane Rates including TRs and TransProbs
[Lane.NumVeh, TransPrTT, TransPrCC, Surplus, Lane.DistCu] = PerLaneRates(FolDist,BaseData(g,:),Num,FixVars,TrData,TrTyp.NumAxPerGr,BatchSize,Lane);

% Initialize empty vars
[VirtualWIM, OverMax, ApercuOverMax, TableNames] = deal([]);

% Initialize parpool if necessary and initialize progress bar
if BaseData.Parallel(g) > 0, gcp; clc; end, m = StartProgBar(BaseData.NumSims(g), Num.Batches, g, height(BaseData)); tic; st = now;

parfor (v = 1:BaseData.NumSims(g), BaseData.Parallel(g)*100)
%for v = 1:BaseData.NumSims(g)  
    
    % Initialize variables outside of batch simulation loop
    LaneAxLineUp = cell(Num.Lanes,1); LaneVehLineUp = cell(Num.Lanes,1); ApercuMax = [];
    [MaxVec, Maxk, MaxDLF, MaxBrStInd, MaxFirstAxInd, MaxFirstAx] = deal(zeros(1,Num.InfCases));
            
    for k = 1:Num.Batches
        
        % ----- START OF RANDOM TRAFFIC GENERATION -----
          
        for q = 1:Num.Lanes
            
            if Lane.TrDistr(q) ~= 0 || FixVars.CarWgt > 0

                % 1) Flo.Veh/Trans | Get Flow of Cars (0s) and Trucks (#s 1 to Num.TrTyp)
                Flo = GetFloVeh(Lane.NumVeh,TransPrTT,TransPrCC,BaseData.BunchFactor(g),q,Lane.DistCu(:,q));
                
                % 2) Flo.Wgt | Get Truck / Axle Weights (kN)
                Flo = GetFloWgt(Num.TrTyp,Lane.NumVeh(q),Flo,TrData.TrDistr);
                
                % 2)bis Flo.Plat/PTrail/PLead/PPrime/Swap | Modify Flo for Platooning
                if BaseData.RunPlat(g) == 1
                    Flo = SwapforPlatoons(Flo,BaseData(g,:),q,Num.TrTyp,PlatPct,TrData.TrDistr,BatchSize,Surplus,Lane);
                end
                
                % 3) Flo.Dist | Get Intervehicle Distances (m) TC is the Car following a Truck (<<<Truck<<Car)
                Flo = GetFloDist(FolDist,FixVars,Flo,BaseData.PlatFolDist(g));
                
                % 4) Assemble Axle Loads and Axle Spacings vectors - populate Axle Weights (kN) and Inter-Axle Distances (m) within
                [Flo, LaneAxLineUp{q}, LaneVehLineUp{q}, TableNames] = GetLaneLineUp(TrTyp,Lane.Dir,q,Flo,FixVars,k,v,BaseData(g,:),TrData);
                
            end
            
        end % END OF LANE SPECIFIC TRAFFIC GENERATION
        
        % Assemble lane specific data, flip for direction, sort by axle position
        [AxLineUp, VehLineUp] = AssembleLineUps(LaneAxLineUp,LaneVehLineUp,BatchSize,Num.Lanes,BaseData(g,:),Lane.Dir,FixVars);
        
        % Log Virtual WIM if necessary... cummulative
        if BaseData.VWIM(g) == 1
            VirtualWIM = [VirtualWIM; VehLineUp];
        end
        
        % ----- END OF RANDOM TRAFFIC GENERATION -----
        % ----- START OF LOAD EFFECT CALCULATION -----
        
        if BaseData.Analysis(g) == 1
            [AllTrAx] = GetAllTrAx(AxLineUp,BaseData.ILRes(g),Lane,FixVars);
            for t = 1:Num.InfCases
                % Subject Influence Line to Axle Stream, lane specific influence line procedure included in GetMaxLE
                [MaxLE,MaxLEStatic,DLF,BrStInd,AxonBr,FirstAxInd,FirstAx] = GetMaxLE(AllTrAx,Inf,BaseData.RunDyn(g),t);
                % Update Maximums if they are exceeded
                if MaxLE > MaxVec(t)
                    [MaxVec(t),Maxk(t),MaxDLF(t),MaxBrStInd(t),MaxFirstAxInd(t),MaxFirstAx(t)] = UpMaxes(MaxLE,k,DLF,BrStInd,FirstAxInd,FirstAx);            
                    % Save results for Apercu
                    if BaseData.Apercu(g) == 1
                        % Must account for ILRes here... /ILRes
                        ApercuMax{t} = VehLineUp(VehLineUp(:,1)/BaseData.ILRes(g)>(BrStInd - 20) & VehLineUp(:,1)/BaseData.ILRes(g)<(BrStInd + Inf.x(end) + 20),:);
                    end
                end
            end 
        end
        
        % ----- END OF LOAD EFFECT CALCULATION -----
        
        % Update progress bar
        UpProgBar(m, st, v, k, BaseData.NumSims(g), Num.Batches)
       
    end % END OF TRAFFIC BATCH

    % Log overall maximum cases into OverMax and ApercuOverMax if necessary
    for i = 1:Num.InfCases
        OverMax = [OverMax; [i, v, MaxVec(i), Maxk(i), MaxDLF(i), MaxBrStInd(i), MaxFirstAxInd(i), MaxFirstAx(i)]];
        % Save VWIM to ApercuOverMax, and add column for InfCase
        if BaseData.Apercu(g) == 1
            ApercuOverMax = [ApercuOverMax; [ApercuMax{i}, repmat(i,size(ApercuMax{i},1),1)]]; 
        end
    end
       
end % END OF SIMULATION

% Get simulation time
[Time] = GetSimTime();

% In the future could add InfCaseName
OverMAXT = array2table(OverMax,'VariableNames',{'InfCase','SimNum','MaxLE','BatchNum','MaxDLF','MaxBrStInd','MaxFirstAxInd','MaxFirstAx'});

% Reshape OverMax results for output
OverMax = sortrows(OverMax); OverMax = OverMax(:,3); OverMax = reshape(OverMax,BaseData.NumSims(g),Num.InfCases);

% Duplicate if only 1 Sim to avoid statistical errors
if BaseData.NumSims(g) == 1
    OverMax = [OverMax; OverMax];
end

% Get ESIM and Ratio
ESIM = 1.1*prctile(OverMax,99); Ratio = ESIM./ESIA.Total;

% Print Summary Stats to Command Window
PrintSummary(BaseData(g,:),BatchSize,PlatPct,TrData,Num,VirtualWIM,Time,Lane.TrDistr)

% Create folders where there are none
CreateFolders(Folder_Name)

TName = datestr(now,'mmmdd-yy HHMM');

% Write results to a file (put into function)
if BaseData.Save(g) == 1
    SaveSummary(TName,strcat('Output', Folder_Name, '/MATSimOutput', InputF(12:end-5),'.xlsx'),TrData,BaseData(g,:),Time,Inf.UniqNames,FolDist,LaneData,ESIM,ESIA.Total,Ratio,OverMax);
end

% Convert VirtualWIMs to tables and save if necessary
if BaseData.VWIM(g) == 1
    PD = array2table(VirtualWIM,'VariableNames',TableNames);    
    save(['VirtualWIM' Folder_Name '/WIM_' TName], 'PD')
end

% Covert Apercus to tables and save if necessary
if BaseData.Apercu(g) == 1
    PD = array2table(ApercuOverMax,'VariableNames',[TableNames 'InfCase']);   
    save(['Apercu' Folder_Name '/AWIM_' TName], 'PD')
end

% Save structure variable with essential simulation information
OutInfo.Name = TName; OutInfo.BaseData = BaseData(g,:);
OutInfo.ESIA = ESIA; OutInfo.ESIM = ESIM;
OutInfo.Mean = mean(OverMax); OutInfo.Std = std(OverMax);
OutInfo.OverMax = OverMax; OutInfo.OverMAXT = OverMAXT;
OutInfo.InfNames = Inf.Names; OutInfo.PlatPct = max(PlatPct);
OutInfo.LaneData = LaneData;

if BaseData.Save(g) == 1
    save(['Output' Folder_Name '/' OutInfo.Name], 'OutInfo')
end

end

% Run Apercu to see critical case
if BaseData.Apercu(g) == 1
    [T, OverMx, AllTrAxx] = GetApercu(PD,OverMAXT,Num.InfCases,Inf,BaseData.RunDyn(g),ESIA.Total,Lane.Dir,BaseData.ILRes(g));
end




% ------------------------------------------------------------------------
%                             AxleStatsBasic
% ------------------------------------------------------------------------
% Assemble Simple, Tandem, and Tridem Axles, using geometry, to gain
% information on Q1 and maximum axle loads
% Differs from function AxleStats because it finds axle groups in
% non-classified vehicles. Optional variable save at end (large vars)

 % Initial commands
clear, clc, tic, format long g, rng('shuffle'), close all;

% Notes
% - Add new stations to existing variables manually (load together and cat)
% - Remember that we are limited to 25t - larger getts tossed by S1Prune
% - Make sure ClassOW is set to true inside Classify.m and no 11bis

% Input Information --------------------
                      
% Traffic Info
%Year = 2003:2019;
%Year = 2010:2011;
%SName = {'Ceneri', 'Denges', 'Gotthard', 'Oberburen'};
%SName = {'Trubbach','Mattstetten'};
SName = {'Gotthard'};

% Toggles
AxleStatsT = 0;
AxleStatsPlot = 0;
Stage2Prune = true;

ILRes = 0.2;   % Needed for WIMtoAllTrAx
    
% Input Complete   ---------------------

AxSingle = array2table(NaN(200000000,3));
AxSingle.Properties.VariableNames = {'AWT1kN','SITE','CLASS'};
AxSingle.Time(:) = NaT;
AxSingleNum = 1;

AxTandem = array2table(NaN(100000000,5));
AxTandem.Properties.VariableNames = {'AWT1kN','AWT2kN','W1_2M','SITE','CLASS'};
AxTandem.Time(:) = NaT;
AxTandemNum = 1;

AxTridem = array2table(NaN(100000000,7));
AxTridem.Properties.VariableNames = {'AWT1kN','AWT2kN','AWT3kN','W1_2M','W2_3M','SITE','CLASS'};
AxTridem.Time(:) = NaT;
AxTridemNum = 1;
    
% For each station to be analyzed
for r = 1:length(SName)
     
    if strcmp(SName{r},'Gotthard')
        LaneDir = [1 2];
    else
        LaneDir = [1 1];
    end
    
    % For each year to be analyzed
    %for i = 1:length(Year)
        
        
        % To track progress
        %fprintf('\nLocation: %s\n    Year: %i\n\n',SName{r},Year(i))
        
        % Load File
        PD = load(['Pruned WIM/',SName{r},'.mat']);
        
        % Add Datetime
        % SWITCH BELOW TO NEW
        %PDC = Classify(PD.PD);  %PD = Daytype(PD,Year(i));
        %PDC = AddDatetime(PDC,1);
        
        % Fix Oberburen pre 2006 station naming issue (SITE = 410 both dirs)
        if strcmp(SName{r},'Oberburen') && Year(i) < 2006
            PDC.SITE(PDC.FS < 3) = 415;
            PDC.SITE(PDC.SITE == 410) = 416;
        end
        if strcmp(SName{r},'Trubbach') && Year(i) < 2006
            PDC.SITE(PDC.FS < 3) = 417;
            PDC.SITE(PDC.SITE == 407) = 418;
        end
        if strcmp(SName{r},'Mattstetten') && Year(i) < 2006
            PDC.SITE(PDC.FS < 3) = 413;
            PDC.SITE(PDC.SITE == 407) = 414;
        end

        % We treat each station separately
        Stations = unique(PDC.SITE);
        
        % For each station
        for w = 1:length(Stations)
            
            Station = Stations(w);
            clc
            fprintf('\nLocation: %s\n    Station: %i\n\n',SName{r},Station)
            
            PDCx = PDC(PDC.SITE == Station,:);
            
            % Further trimming if necessary
            if Stage2Prune
                PDCx = PruneWIM2(PDCx,0);
            end
            
            % Convert PDC to AllTrAx (must be greater than 0 to actually Spacesave!
            % The 4 here is strategic... combined with min veh length
            [PDCx, AllTrAx, TrLineUp] = VBWIMtoAllTrAx(PDCx,4,LaneDir,ILRes);
                       
            % Round TrLineUp first row, move unrounded to fifth row
            TrLineUp(:,5) = TrLineUp(:,1); TrLineUp(:,1) = round(TrLineUp(:,1)/ILRes);
            
            % TrLineUp [       1             2         3        4         5       ]
            %            AllTrAxIndex    AxleValue   Truck#   LaneID  Station(m)
            % Distances refer to distance in front
                        
            % Treat the lanes separately & incorporate direction
            Lanes = unique(PDCx.FS);
            
            % Complication vehicle recognition scheme (created through trial and error)
            for p = 1:length(Lanes)
                
                % LaneB is lane boolean
                LaneB = TrLineUp(:,4) == Lanes(p);
                
                if LaneDir(p) == 1
                    
                    TrLineUp(LaneB,6) = [10; diff(TrLineUp(LaneB,5))];
                    % 7th column is Single Axles [ >= 2.4m infront and behind ]
                    % Or if the row above is a diff veh and the row below is >= 2.4m, or if the row below is a diff veh and row above is > 2.4m
                    TrLineUp(LaneB,7) = (TrLineUp(LaneB,6) >= 2.4 & circshift(TrLineUp(LaneB,6),-1) >= 2.4) | (circshift(TrLineUp(LaneB,3),1) ~= TrLineUp(LaneB,3) & circshift(TrLineUp(LaneB,6),-1) >= 2.4)...
                        | (circshift(TrLineUp(LaneB,3),-1) ~= TrLineUp(LaneB,3) & TrLineUp(LaneB,6) >= 2.4);
                    % 8th column is Tridem Axles
                    %                       cannot be single            cannot begin before the third axle                       first gap less than 2.4               second gap less than 2.4            third gap larger than 2.4
                    TrLineUp(LaneB,8) = TrLineUp(LaneB,7) == 0 & circshift(TrLineUp(LaneB,3),2) == TrLineUp(LaneB,3) & circshift(TrLineUp(LaneB,6),-1) < 2.4 & circshift(TrLineUp(LaneB,6),-2) < 2.4 & circshift(TrLineUp(LaneB,6),-3) >= 2.4;
                    
                else
                    
                    % Fix order before performing diff
                    TrLineUpx = TrLineUp(LaneB,:);
                    [TrLineUpxy, b] = sortrows(TrLineUpx,1);
                    TrLineUpxy(:,6) = [diff(TrLineUpxy(:,5)); 10];
                    
                    TrLineUp(LaneB,:) = TrLineUpxy(b,:);

                    % 7th column is Single Axles [ >= 2.4m infront and behind ]
                    % Or if the row above is a diff veh and the row below is >= 2.4m, or if the row below is a diff veh and row above is > 2.4m
                    TrLineUp(LaneB,7) = (TrLineUp(LaneB,6) >= 2.4 & circshift(TrLineUp(LaneB,6),-1) >= 2.4) | (circshift(TrLineUp(LaneB,3),1) ~= TrLineUp(LaneB,3) & circshift(TrLineUp(LaneB,6),-1) >= 2.4)...
                        | (circshift(TrLineUp(LaneB,3),-1) ~= TrLineUp(LaneB,3) & TrLineUp(LaneB,6) >= 2.4);
                    % 8th column is Tridem Axles
                    %                       cannot be single            cannot begin before the third axle                       first gap less than 2.4               second gap less than 2.4            third gap larger than 2.4
                    TrLineUp(LaneB,8) = TrLineUp(LaneB,7) == 0 & circshift(TrLineUp(LaneB,3),2) == TrLineUp(LaneB,3) & circshift(TrLineUp(LaneB,6),-1) < 2.4 & circshift(TrLineUp(LaneB,6),-2) < 2.4 & circshift(TrLineUp(LaneB,6),-3) >= 2.4;

                end
            end
                
            % Set the 7th column to 3 whenever we have a tridem
            TrLineUp(TrLineUp(:,8)==1,7) = 3;
            TrLineUp(circshift(TrLineUp(:,8)==1,1),7) = 3;
            TrLineUp(circshift(TrLineUp(:,8)==1,2),7) = 3;
            
            for p = 1:length(Lanes)
                % LaneB is lane boolean
                LaneB = TrLineUp(:,4) == Lanes(p);
                % 9th column is Tandem Axles
                %                cannot be single or tridem  neither can the following one                      gap less than 2.4         gap infront larger than 2         second gap larger than 2
                TrLineUp(LaneB,9) = TrLineUp(LaneB,7) == 0 & circshift(TrLineUp(LaneB,7),-1) == 0 & circshift(TrLineUp(LaneB,6),-1) < 2.4 & TrLineUp(LaneB,6) >= 2;% & circshift(TrLineUp(LaneB,6),-2) >= 1.2;
            end
                
            % Clean up 3 in a rows
            TrLineUp(TrLineUp(:,9) == 1 & circshift(TrLineUp(:,9),-1) == 1 & circshift(TrLineUp(:,9),1) == 1,9) = 0;
            
            % Set the 7th column to 2 whenever we have a tandem
            TrLineUp(TrLineUp(:,9)==1,7) = 2;
            TrLineUp(circshift(TrLineUp(:,9)==1,1),7) = 2;
            
            % TrLineUp [       1             2         3        4         5            6            7           8         9       ]
            %            AllTrAxIndex    AxleValue   Truck#   LaneID  Station(m) StationDiff(m) SingleFlag TridemFlag TandemFlag          
            % Distances refer to distance in front
            
            % Optional: Verify with AxleStats
            if AxleStatsT
                TrTyps = [11; 12; 22; 23; 111; 11117; 1127; 12117; 122; 11127; 1128; 1138; 1238];
                TrAxPerGr = [11; 12; 22; 23; 111; 1111; 112; 1211; 122; 1112; 112; 113; 123];
                Vec = [0 1 2 1 0 0 1 1 2 1 1 0 1];
                
                [STaTr,AllAx] = AxleStats(PDCx,TrAxPerGr,TrTyps,[SName{r} ' ' num2str(Station)],Year(i),AxleStatsPlot);
            end
            
            % PDCx Index is found in TrLineUp(:,3)
            
            % Legacy (kept because it is informative)
            % All [Q Class]
            %Axles.([SName{r} num2str(Station)])(Year(i)-2000).All = [TrLineUp(:,2) PDCx.CLASS(TrLineUp(:,3))];
            % Single [Q Class]
            %Axles.([SName{r} num2str(Station)])(Year(i)-2000).Single = [TrLineUp(TrLineUp(:,7)==1,2) PDCx.CLASS(TrLineUp(TrLineUp(:,7)==1,3))];
            % Tandem [Total Q1 Q2 Space Class]                                                            Total                                                       Axle 1                          Axle 2                                                          Space btwn                                             Class
            %Axles.([SName{r} num2str(Station)])(Year(i)-2000).Tandem = [(TrLineUp(TrLineUp(:,9)==1,2) + TrLineUp(circshift(TrLineUp(:,9)==1,1),2)) TrLineUp(TrLineUp(:,9)==1,2) TrLineUp(circshift(TrLineUp(:,9)==1,1),2) abs(TrLineUp(circshift(TrLineUp(:,9)==1,1),5)-TrLineUp(TrLineUp(:,9)==1,5)) PDCx.CLASS(TrLineUp(TrLineUp(:,9)==1,3))]; 
            % Tridem [Total Q1 Q2 Q3 Space1 Space2 Class]
            %Axles.([SName{r} num2str(Station)])(Year(i)-2000).Tridem = [(TrLineUp(TrLineUp(:,8)==1,2) + TrLineUp(circshift(TrLineUp(:,8)==1,1),2) + TrLineUp(circshift(TrLineUp(:,8)==1,2),2))  TrLineUp(TrLineUp(:,8)==1,2)  TrLineUp(circshift(TrLineUp(:,8)==1,1),2)  TrLineUp(circshift(TrLineUp(:,8)==1,2),2) abs(TrLineUp(circshift(TrLineUp(:,8)==1,1),5)-TrLineUp(TrLineUp(:,8)==1,5)) abs(TrLineUp(circshift(TrLineUp(:,8)==1,2),5)-TrLineUp(circshift(TrLineUp(:,8)==1,1),5)) PDCx.CLASS(TrLineUp(TrLineUp(:,8)==1,3))];
            
            Ind = TrLineUp(TrLineUp(:,7)==1,3);
            Len = length(Ind);
            AxSingle.AWT1kN(AxSingleNum:AxSingleNum+Len-1) = TrLineUp(TrLineUp(:,7)==1,2);
            AxSingle.SITE(AxSingleNum:AxSingleNum+Len-1) = PDCx.SITE(Ind);
            AxSingle.CLASS(AxSingleNum:AxSingleNum+Len-1) = PDCx.CLASS(Ind);
            AxSingle.Time(AxSingleNum:AxSingleNum+Len-1) = PDCx.Time(Ind);
            AxSingleNum = AxSingleNum+Len;
            
            Ind = TrLineUp(TrLineUp(:,9)==1,3);
            Len = length(Ind);
            AxTandem.AWT1kN(AxTandemNum:AxTandemNum+Len-1) = TrLineUp(TrLineUp(:,9)==1,2);
            AxTandem.AWT2kN(AxTandemNum:AxTandemNum+Len-1) = TrLineUp(circshift(TrLineUp(:,9)==1,1),2);
            AxTandem.W1_2M(AxTandemNum:AxTandemNum+Len-1) = abs(TrLineUp(circshift(TrLineUp(:,9)==1,1),5)-TrLineUp(TrLineUp(:,9)==1,5));
            AxTandem.SITE(AxTandemNum:AxTandemNum+Len-1) = PDCx.SITE(Ind);
            AxTandem.CLASS(AxTandemNum:AxTandemNum+Len-1) = PDCx.CLASS(Ind);
            AxTandem.Time(AxTandemNum:AxTandemNum+Len-1) = PDCx.Time(Ind);
            AxTandemNum = AxTandemNum+Len;
            
            Ind = TrLineUp(TrLineUp(:,8)==1,3);
            Len = length(Ind);
            AxTridem.AWT1kN(AxTridemNum:AxTridemNum+Len-1) = TrLineUp(TrLineUp(:,8)==1,2);
            AxTridem.AWT2kN(AxTridemNum:AxTridemNum+Len-1) = TrLineUp(circshift(TrLineUp(:,8)==1,1),2);
            AxTridem.AWT3kN(AxTridemNum:AxTridemNum+Len-1) = TrLineUp(circshift(TrLineUp(:,8)==1,2),2);
            AxTridem.W1_2M(AxTridemNum:AxTridemNum+Len-1) = abs(TrLineUp(circshift(TrLineUp(:,8)==1,1),5)-TrLineUp(TrLineUp(:,8)==1,5));
            AxTridem.W2_3M(AxTridemNum:AxTridemNum+Len-1) = abs(TrLineUp(circshift(TrLineUp(:,8)==1,2),5)-TrLineUp(circshift(TrLineUp(:,8)==1,1),5));
            AxTridem.SITE(AxTridemNum:AxTridemNum+Len-1) = PDCx.SITE(Ind);
            AxTridem.CLASS(AxTridemNum:AxTridemNum+Len-1) = PDCx.CLASS(Ind);
            AxTridem.Time(AxTridemNum:AxTridemNum+Len-1) = PDCx.Time(Ind);
            AxTridemNum = AxTridemNum+Len;
            
        end
    %end
end

AxSingle(isnan(AxSingle.SITE),:) = [];
AxTandem(isnan(AxTandem.SITE),:) = [];
AxTridem(isnan(AxTridem.SITE),:) = [];

% Saving is manual
toc

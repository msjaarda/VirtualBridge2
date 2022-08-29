function E = VBGetECode(ILData,ILRes)
% VBGetECode takes an influence line (or set) in the form of 
% To be decided if we want ILData to be simply the line (ILData.v) or the
% structure complete with values (v) and Name
% Warning!!!! ILData must have the same resolution as ILRes otherwise
% results will be wrong!!

for i = 1:length(ILData)

    CSA.Axle(i) = 0; % initialise the number of axle on the bridge for CSA code
    
    % Interpolate around influence lines to figure out next biggest max
    for j = 1:size(ILData(i).v,2)
        yILv = ILData(i).v(:,j);
        
        if j == 1
            PL = 300; % SIA
            PLBR(1) = 210; %210; %140; % Bruhwiler values from 15 october 2019 document, updated with the one found with blockmax axles types
            PLBR(2) = 175; %175; %110;
            PLBR(3) = 135; %135; %95;
            PL12 = 120; % 12 tonnes per axle (load model for 41, 48 and 72 truck)
            PLEU = 300; % EUROCODE
            PLCSA(1) = 50; PLCSA(2) = 125; PLCSA(3) = 125; PLCSA(4) = 175; PLCSA(5) = 150; % CSA (Canadian code)
            PLAAS(1) = 35.6; PLAAS(2) = 142.5; PLAAS(3) = 111.5; % AAASHTO (American code)
        elseif j == 2
            PL = 200; % SIA
            PLBR(1) = 0;
            PLBR(2) = 150; %150; %85;
            PLBR(3) = 110; %110; %70;
            PLEU = 200; % EUROCODE
        elseif j == 3
            PL = 0; % SIA
            PLBR(2) = 0; %0;
            PLBR(3) = 0; %0;
            PLEU = 100; % EUROCODE
        else
            PLEU = 0; % EUROCODE
        end
        
        % Concentrated loads
        Conc.SIA = zeros(round(1.2/ILRes)+1,1); % SIA Load Model
        Conc.BR1 = PLBR(1); % Brühwiler Load Model, 1 axle
        Conc.BR2 = zeros(round(1.2/ILRes)+1,1); % Brühwiler Load Model, 2 axles
        Conc.BR3 = zeros(round(2.4/ILRes)+1,1); % Brühwiler Load Model, 3 axles
        Conc.K41 = zeros(round(7.2/ILRes)+1,1); % KUBA-ST Load Model, 41 truck type
        Conc.K48 = zeros(round(9.1/ILRes)+1,1); % KUBA-ST Load Model, 48 truck type
        Conc.K72 = zeros(round(14/ILRes)+1,1); % KUBA-ST Load Model, 72 truck type
        Conc.KSPTR = zeros(round(12.5/ILRes)+1,1); % KUBA-ST Load Model, combination of special transport with 15 meters gap
        Conc.EURO = zeros(round(1.2/ILRes)+1,1); % EUROCODE Load Model 1
        Conc.CSA = zeros(round(18/ILRes)+1,1); % CSA truck CL-625 (Canadian Code)
        Conc.AAS1 = zeros(round(13.5/ILRes)+1,1); % AASHTO truck HS20-44 (American Code)
        Conc.AAS2 = zeros(round(1.2/ILRes)+1,1); % AASHTO Tandem (American Code)
                
           % SIA Load model (Swiss Code)  
           Conc.SIA(1) = PL;
           Conc.SIA(end) = Conc.SIA(end) + PL;
           
           % Brühwiler Load Model, 2 axles
           Conc.BR2(1) = PLBR(2);
           Conc.BR2(end) = Conc.BR2(end) + PLBR(2);
           
           % EUROCODE Load Model 1
           Conc.EURO(1) = PLEU;
           Conc.EURO(end) = Conc.EURO(end) + PLEU;
                        
        if max(size(Conc.BR3)) == 1
            Conc.BR3(1) = PLBR(3)*3;
        elseif max(size(Conc.BR3)) == 2
            Conc.BR3(1) = PLBR(3)*1.5; Conc.BR3(end) = PLBR(3)*1.5;
        else
            Conc.BR3(1) = PLBR(3); Conc.BR3(end) = PLBR(3);
            try
                Conc.BR3((end+1)/2) = PLBR(3);
            catch
                Conc.BR3((end)/2) = PLBR(3)*.5; Conc.BR3((end)/2+1) = PLBR(3)*.5;
            end
        end
        
        if j == 1 % Load model for the 41 crane, due to resolution must add previous load if applied in the same spot (correction for low resolution)
           Conc.K41(1) = PL12;
           Conc.K41(round(2.4/ILRes)+1) = Conc.K41(round(2.4/ILRes)+1) + PL12;
           Conc.K41(round(4/ILRes)+1) = Conc.K41(round(4/ILRes)+1) + PL12;
           Conc.K41(round(5.6/ILRes)+1) = Conc.K41(round(5.6/ILRes)+1) + PL12;
           Conc.K41(end) = Conc.K41(end) + PL12;
           
           % Load model for the 48 crane (correction for low resolution)
           Conc.K48(1) = PL12;
           Conc.K48(round(2.6/ILRes)+1) = Conc.K48(round(2.6/ILRes)+1) + PL12;
           Conc.K48(round(4.2/ILRes)+1) = Conc.K48(round(4.2/ILRes)+1) + PL12;
           Conc.K48(round(5.9/ILRes)+1) = Conc.K48(round(5.9/ILRes)+1) + PL12;
           Conc.K48(round(7.5/ILRes)+1) = Conc.K48(round(7.5/ILRes)+1) + PL12;
           Conc.K48(end) = Conc.K48(end) + PL12;
           
           % Load model for the 48 crane (correction for low resolution)
           Conc.K72(1) = PL12;
           Conc.K72(round(3.2/ILRes)+1) = Conc.K72(round(3.2/ILRes)+1) + PL12;
           Conc.K72(round(4.6/ILRes)+1) = Conc.K72(round(4.6/ILRes)+1) + PL12;
           Conc.K72(round(9.8/ILRes)+1) = Conc.K72(round(9.8/ILRes)+1) + PL12;
           Conc.K72(round(11.2/ILRes)+1) = Conc.K72(round(11.2/ILRes)+1) + PL12;
           Conc.K72(round(12.6/ILRes)+1) = Conc.K72(round(12.6/ILRes)+1) + PL12;
           Conc.K72(end) = Conc.K72(end) + PL12;
           
           % Load model for combination of all special transports (correction for low resolution)
           Conc.KSPTR(1) = PL12;
           Conc.KSPTR(round(3.2/ILRes)+1) = Conc.KSPTR(round(3.2/ILRes)+1) + PL12;
           Conc.KSPTR(round(4.6/ILRes)+1) = Conc.KSPTR(round(4.6/ILRes)+1) + PL12;
           Conc.KSPTR(round(9.7/ILRes)+1) = Conc.KSPTR(round(9.7/ILRes)+1) + PL12;
           Conc.KSPTR(round(11.1/ILRes)+1) = Conc.KSPTR(round(11.1/ILRes)+1) + PL12;
           Conc.KSPTR(end) = Conc.KSPTR(end) + PL12;
           
        else
           Conc.K41 = Conc.SIA;
           Conc.K48 = Conc.SIA;
           Conc.K72 = Conc.SIA;
           Conc.KSPTR = Conc.SIA;
        end
        
           % CSA Load model (Canadian Code) (correction for low resolution)
           Conc.CSA(1) = PLCSA(1);
           Conc.CSA(round(3.6/ILRes)+1) = Conc.CSA(round(3.6/ILRes)+1) + PLCSA(2);
           Conc.CSA(round(4.8/ILRes)+1) = Conc.CSA(round(4.8/ILRes)+1) + PLCSA(3);
           Conc.CSA(round(11.4/ILRes)+1) = Conc.CSA(round(11.4/ILRes)+1) + PLCSA(4);
           Conc.CSA(end) = Conc.CSA(end) + PLCSA(5);
           
           % AASHTO Load model (American Code) (correction for low resolution)
           Conc.AAS1(1) = PLAAS(1);
           Conc.AAS1(round(4.3/ILRes)+1) = Conc.AAS1(round(4.3/ILRes)+1) + PLAAS(2);
           Conc.AAS1(end) = Conc.AAS1(end) + PLAAS(2);
           Conc.AAS2(1) = PLAAS(3);
           Conc.AAS2(end) = Conc.AAS2(end) + PLAAS(3);
                   
        % Save the biggest contribution, save also the position for load models with q gaps
        MaxInf.vCONV(j,i) = max(conv(Conc.SIA,yILv));
        [MaxInf.vCONVBR1(j,i),MaxInf.vCONVBR1Posi(j,i)] = max(conv(Conc.BR1,yILv));
        [MaxInf.vCONVBR2(j,i),MaxInf.vCONVBR2Posi(j,i)] = max(conv(Conc.BR2,yILv));
        [MaxInf.vCONVBR3(j,i),MaxInf.vCONVBR3Posi(j,i)] = max(conv(Conc.BR3,yILv));
        [MaxInf.vCONV41(j,i),MaxInf.vCONV41Posi(j,i)] = max(conv(Conc.K41,yILv));
        [MaxInf.vCONV48(j,i),MaxInf.vCONV48Posi(j,i)] = max(conv(Conc.K48,yILv));
        [MaxInf.vCONV72(j,i),MaxInf.vCONV72Posi(j,i)] = max(conv(Conc.K72,yILv));
        [MaxInf.vCONVSPTR(j,i),MaxInf.vCONVSPTRPosi(j,i)] = max(conv(Conc.KSPTR,yILv));
        MaxInf.vCONVEURO(j,i) = max(conv(Conc.EURO,yILv));
        [MaxInf.vCONVCSA(j,i),MaxInf.vCONVCSAPosi(j,i)] = max(conv(Conc.CSA,yILv));
        LongeurPontRes = (size(ILData(i).v,1)-1); % Length of the bridge with resolution factor
        % Find wich axles are determinant for CSA dynamic factor
        CSA.OnBridge = [(MaxInf.vCONVCSAPosi(j,i)<=LongeurPontRes+1)&&(MaxInf.vCONVCSAPosi(j,i)>1);(MaxInf.vCONVCSAPosi(j,i)-3.6/ILRes<=LongeurPontRes+1)&&(MaxInf.vCONVCSAPosi(j,i)-3.6/ILRes>0);(MaxInf.vCONVCSAPosi(j,i)-4.8/ILRes<=LongeurPontRes+1)&&(MaxInf.vCONVCSAPosi(j,i)-4.8/ILRes>0);(MaxInf.vCONVCSAPosi(j,i)-11.4/ILRes<=LongeurPontRes+1)&&(MaxInf.vCONVCSAPosi(j,i)-11.4/ILRes>0);(MaxInf.vCONVCSAPosi(j,i)-18/ILRes<=LongeurPontRes+1)&&(MaxInf.vCONVCSAPosi(j,i)-18/ILRes>0)];
        CSA.OnBridge = CSA.OnBridge.*[0.8;0.8;0.8;1;1]; % weight of the axle for the combination CSA dynamic factor
        CSA.Axle(i) = CSA.Axle(i)+round(sum(CSA.OnBridge)); % Number of CSA axle on the bridge
        % Loop for AASHTO, to modify the rear axle position
        MaxInf.vCONVAAS1(j,i) = max(conv(Conc.AAS1,yILv));
        for k=1:round(4.9/ILRes)
            Conc.AAS1(end-k) = Conc.AAS1(end-k+1); Conc.AAS1(end-k+1) = 0;
            MaxTempAAS1 = max(conv(Conc.AAS1,yILv));
            MaxInf.vCONVAAS1(j,i) = max(MaxInf.vCONVAAS1(j,i),MaxTempAAS1);
        end
        MaxInf.vCONVAAS2(j,i) = max(conv(Conc.AAS2,yILv)); % AASHTO tandem (American Code), not sure if you must add other tandem at 8 to 12 meters
        Conc.AAS3 = 0.9*Conc.AAS1; % AASHTO 2 trucks HS20-44 90% for Mn (American Code)
        MaxInf.vCONVAAS3(j,i) = 0;
        [MaxInf.vCONVAAS31(j,i),MaxInf.vCONVAAS3Posi(j,i)] = max(conv(Conc.AAS3,yILv));
        yILvtemp = yILv; yILvtemp(max(MaxInf.vCONVAAS3Posi(j,i)+round(-8.6/ILRes-15.2/ILRes),1):min(MaxInf.vCONVAAS3Posi(j,i)+round(15.2/ILRes),end)) = 0;
        MaxInf.vCONVAAS32(j,i) = max(conv(Conc.AAS3,yILvtemp));
        MaxInf.vCONVAAS3(j,i) = MaxInf.vCONVAAS31(j,i)+MaxInf.vCONVAAS32(j,i);
    end
    CSA.AllDLA = [1.4;1.3;1.25]; % CSA all dynamic factors
    CSA.DLA(i) = CSA.AllDLA(min(CSA.Axle(i),3)); % Current dynamic factor for the infl studied
 end

% Assign integral values into IntInfv (each InfCase)
for i = 1:length(ILData)
    x = 0:ILRes:(length(ILData(i).v)-1)*ILRes;
    A = ILData(i).v;
    A(A<0) = 0;
    IntInf.v(:,i) = trapz(x,A);
    % Create gap for load models who need it
    for j = 1:size(A,2)
        IntInf.vBR1(j,i) = IntInf.v(j,i)-(j==1)*trapz(x(max(MaxInf.vCONVBR1Posi(j,i)-round(2/ILRes),1):min(MaxInf.vCONVBR1Posi(j,i)+round(2/ILRes),end)),A(max(MaxInf.vCONVBR1Posi(j,i)-round(2/ILRes),1):min(MaxInf.vCONVBR1Posi(j,i)+round(2/ILRes),end),j));
        IntInf.vBR2(j,i) = IntInf.v(j,i)-(j==1||j==2)*trapz(x(max(MaxInf.vCONVBR2Posi(j,i)-round(3.2/ILRes),1):min(MaxInf.vCONVBR2Posi(j,i)+round(2/ILRes),end)),A(max(MaxInf.vCONVBR2Posi(j,i)-round(3.2/ILRes),1):min(MaxInf.vCONVBR2Posi(j,i)+round(2/ILRes),end),j));
        IntInf.vBR3(j,i) = IntInf.v(j,i)-(j==1||j==2)*trapz(x(max(MaxInf.vCONVBR3Posi(j,i)-round(4.4/ILRes),1):min(MaxInf.vCONVBR3Posi(j,i)+round(2/ILRes),end)),A(max(MaxInf.vCONVBR3Posi(j,i)-round(4.4/ILRes),1):min(MaxInf.vCONVBR3Posi(j,i)+round(2/ILRes),end),j));
        IntInf.v41(j,i) = IntInf.v(j,i)-(j==1)*trapz(x(max(MaxInf.vCONV41Posi(j,i)-round(8.9/ILRes),1):min(MaxInf.vCONV41Posi(j,i)+round(4.3/ILRes),end)),A(max(MaxInf.vCONV41Posi(j,i)-round(8.9/ILRes),1):min(MaxInf.vCONV41Posi(j,i)+round(4.3/ILRes),end),j));
        IntInf.v48(j,i) = IntInf.v(j,i)-(j==1)*trapz(x(max(MaxInf.vCONV48Posi(j,i)-round(11/ILRes),1):min(MaxInf.vCONV48Posi(j,i)+round(4.9/ILRes),end)),A(max(MaxInf.vCONV48Posi(j,i)-round(11/ILRes),1):min(MaxInf.vCONV48Posi(j,i)+round(4.9/ILRes),end),j));
        IntInf.v72(j,i) = IntInf.v(j,i)-(j==1)*trapz(x(max(MaxInf.vCONV72Posi(j,i)-round(15/ILRes),1):min(MaxInf.vCONV72Posi(j,i)+round(1.5/ILRes),end)),A(max(MaxInf.vCONV72Posi(j,i)-round(15/ILRes),1):min(MaxInf.vCONV72Posi(j,i)+round(1.5/ILRes),end),j));
        IntInf.vSPTR(j,i) = IntInf.v(j,i)-(j==1)*trapz(x(max(MaxInf.vCONVSPTRPosi(j,i)-round(13.5/ILRes),1):min(MaxInf.vCONVSPTRPosi(j,i)+round(1.5/ILRes),end)),A(max(MaxInf.vCONVSPTRPosi(j,i)-round(13.5/ILRes),1):min(MaxInf.vCONVSPTRPosi(j,i)+round(1.5/ILRes),end),j));
    end
end

% Define ESIA details
LaneWidth = 3; % meters, hard coded
% Distributed loads CHANGE TO 2.5 NORMALLY (3 FOR EXPERIMENT)
qk = 2.5*ones(width(ILData(1).v),1);
qkBR = qk;
qk05 = qk; % kN/m2 with alpha = 1 (needed alpha = 0.5)
qk(1) = 9; % kN/m2
qkBR(1) = 3.6; % kN/m2
qk05(1) = 9; % kN/m2 with alpha = 1 (needed alpha = 0.5)
qkCSA = 3*ones(width(ILData(1).v),1); % kN/m2 Canadian code
qkAAS = 3.11*ones(width(ILData(1).v),1); % kN/m2 American code
% Alpha is 1 to make ratios easier (note that it is 0.9 in the code)
Alpha = 1;

% On 25.03.2021 Matt and Lucas used LucasInfluenceLine to show that this
% method underpredicts ESIA for twin girder bridges because in Lucas' code he
% shifts the point loads Q1 and Q2 to the edge, and I do not. TM did the
% same as Lucas.

% Calculate ESIA for each InfCase
for i = 1:length(ILData)
    Maxv.CONV = MaxInf.vCONV(:,i); % SIA CODE
    Maxv.CONVBR1 = MaxInf.vCONVBR1(:,i); % Brühwiler load model 1
    Maxv.CONVBR2 = MaxInf.vCONVBR2(:,i); % Brühwiler load model 2
    Maxv.CONVBR3 = MaxInf.vCONVBR3(:,i); % Brühwiler load model 3
    Maxv.CONV41 = MaxInf.vCONV41(:,i); % KUBA-ST model for 41 truck
    Maxv.CONV48 = MaxInf.vCONV48(:,i); % KUBA-ST model for 48 truck
    Maxv.CONV72 = MaxInf.vCONV72(:,i); % KUBA-ST model for 72 truck
    Maxv.CONVSPTR = MaxInf.vCONVSPTR(:,i); % KUBA-ST model for special transport truck (generalised truck)
    Maxv.CONVSPTR_RmvL1 = Maxv.CONVSPTR; 
    Maxv.CONVSPTR_RmvL1(1) = 0; % KUBA-ST model for special transport truck (generalised truck), removed truck from 1st Lane (to analyse accompaniment)
    Maxv.CONVEURO = MaxInf.vCONVEURO(:,i); % EUROCODE
    Maxv.CONVCSA = MaxInf.vCONVCSA(:,i); % CSA CODE
    Maxv.CONVAAS1 = MaxInf.vCONVAAS1(:,i); % AASHTO CODE, truck HS20-44
    Maxv.CONVAAS2 = MaxInf.vCONVAAS2(:,i); % AASHTO CODE, tandem
    Maxv.CONVAAS3 = MaxInf.vCONVAAS3(:,i); % AASHTO CODE, 0.9*2truck HS20-44 for Mn
    
    Int.v = IntInf.v(:,i); % SIA CODE
    Int.vBR1 = IntInf.vBR1(:,i); % Brühwiler load model 1
    Int.vBR2 = IntInf.vBR2(:,i); % Brühwiler load model 2
    Int.vBR3 = IntInf.vBR3(:,i); % Brühwiler load model 3
    Int.v41 = IntInf.v41(:,i); % KUBA-ST model for 41 truck
    Int.v48 = IntInf.v48(:,i); % KUBA-ST model for 48 truck
    Int.v72 = IntInf.v72(:,i); % KUBA-ST model for 72 truck
    Int.vSPTR = IntInf.vSPTR(:,i); % KUBA-ST model for special transport truck (generalised truck)
    Int.vEURO = IntInf.v(:,i); % EUROCODE
    Int.vCSA = IntInf.v(:,i); % CSA CODE
    Int.vAAS = IntInf.v(:,i); % AASHTO CODE

    % SIA (Swiss code)
    E(i).SIA.Total = 1.5*Alpha*(sum(Maxv.CONV)+Int.v'*qk*LaneWidth);
    E(i).SIA.EQ = Maxv.CONV;
    E(i).SIA.Eq = Int.v.*qk*LaneWidth;
    
    % EUROCODE (European code) - load model n°1
    E(i).EURO.Total = 1.35*Alpha*(sum(Maxv.CONVEURO)+Int.vEURO'*qk*LaneWidth);
    E(i).EURO.EQ = Maxv.CONVEURO;
    E(i).EURO.Eq = Int.vEURO.*qk*LaneWidth;
    
    % CSA (Canadian code)
    MultiLaneFactorCSA = [1;0.9;0.8;0.7;0.6;0.55];
    LiveLoadFactorCSA = 1.7;
    E(i).CSA.Total = 0;
    E(i).CSA.EQm1 = MultiLaneFactorCSA(min(size(ILData(i).v,2),6))*Maxv.CONVCSA*CSA.DLA(i);
    E(i).CSA.EQm2 = MultiLaneFactorCSA(min(size(ILData(i).v,2),6))*(Int.vCSA.*qkCSA*LaneWidth+0.8*Maxv.CONVCSA);
    E(i).CSA.Total = max([sum(E(i).CSA.EQm1),sum(E(i).CSA.EQm2)])*LiveLoadFactorCSA;
        
    % AASHTO (American code)
    AASHTODLA = 1.33; % AASHTO Dynamic Load
    MultiLaneFactorAASHTO = [1.2;1.0;0.85;0.65];
    LiveLoadFactorAASHTO = 1.75;
    E(i).AASHTO.Total = 0;
    E(i).AASHTO.EQm1 = MultiLaneFactorAASHTO(min(size(ILData(i).v,2),4))*Maxv.CONVAAS1*AASHTODLA+Int.vAAS.*qkAAS*LaneWidth;
    E(i).AASHTO.EQm2 = MultiLaneFactorAASHTO(min(size(ILData(i).v,2),4))*Maxv.CONVAAS2*AASHTODLA+Int.vAAS.*qkAAS*LaneWidth;
    E(i).AASHTO.EQm3 = MultiLaneFactorAASHTO(min(size(ILData(i).v,2),4))*Maxv.CONVAAS3*AASHTODLA+Int.vAAS.*qkAAS*LaneWidth*0.9;
    E(i).AASHTO.Total = max([sum(E(i).AASHTO.EQm1),sum(E(i).AASHTO.EQm2),sum(E(i).AASHTO.EQm3)])*LiveLoadFactorAASHTO;
    
    % KUBA-ST (ASTRA load model) 
    E(i).KUBA_T41.Total= 1.5*Alpha*(sum(Maxv.CONV41)+Int.v41'*qk05*LaneWidth);
    E(i).KUBA_T41.EQ = Maxv.CONV41;
    E(i).KUBA_T41.Eq = Int.v41.*qk05*LaneWidth;
    E(i).KUBA_T48.Total= 1.5*Alpha*(sum(Maxv.CONV48)+Int.v48'*qk05*LaneWidth);
    E(i).KUBA_T48.EQ = Maxv.CONV48;
    E(i).KUBA_T48.Eq = Int.v48.*qk05*LaneWidth;
    E(i).KUBA_T72.Total= 1.5*Alpha*(sum(Maxv.CONV72)+Int.v72'*qk05*LaneWidth);
    E(i).KUBA_T72.EQ = Maxv.CONV72;
    E(i).KUBA_T72.Eq = Int.v72.*qk05*LaneWidth;
    E(i).KUBA_SPTR.Total= 1.5*Alpha*(sum(Maxv.CONVSPTR)+Int.vSPTR'*qk05*LaneWidth);
    E(i).KUBA_SPTR.EQ = Maxv.CONVSPTR;
    E(i).KUBA_SPTR.Eq = Int.vSPTR.*qk05*LaneWidth;
    % EKUBA but the load from the KUBA truck in 1st lane is removed
    E(i).KUBA_SPTR_RmvL1.Total= 1.5*Alpha*(sum(Maxv.CONVSPTR_RmvL1)+Int.vSPTR'*qk05*LaneWidth);
    E(i).KUBA_SPTR_RmvL1.EQ = Maxv.CONVSPTR_RmvL1;
    E(i).KUBA_SPTR_RmvL1.Eq = Int.vSPTR.*qk05*LaneWidth;
    
    % EBRU (Professor Eugen Brühwiler load model (EPFL))
    [E(i).BRU.Total,posi] = max([1.5*(sum(Maxv.CONVBR1)+Int.vBR1'*qkBR*LaneWidth);1.5*(sum(Maxv.CONVBR2)+Int.vBR2'*qkBR*LaneWidth);1.5*(sum(Maxv.CONVBR3)+Int.vBR3'*qkBR*LaneWidth)]);
    E(i).BRU.LoadMod = append('Model n°',int2str(posi));
    if posi == 1
        E(i).BRU.EQ = Maxv.CONVBR1;
        E(i).BRU.Eq = Int.vBR1.*qkBR*LaneWidth;
    elseif posi == 2
        E(i).BRU.EQ = Maxv.CONVBR2;
        E(i).BRU.Eq = Int.vBR2.*qkBR*LaneWidth;
    elseif posi == 3
        E(i).BRU.EQ = Maxv.CONVBR3;
        E(i).BRU.Eq = Int.vBR3.*qkBR*LaneWidth;
    else
        error('error with EBRU');
    end
        
end

end


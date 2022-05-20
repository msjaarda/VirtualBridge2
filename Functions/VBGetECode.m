function E = VBGetECode(ILData,ILRes)
% VBGetECode takes an influence line (or set) in the form of 
% To be decided if we want ILData to be simply the line (ILData.v) or the
% structure complete with values (v) and Name

for i = 1:length(ILData)

    % Interpolate around influence lines to figure out next biggest max
    for j = 1:size(ILData(i).v,2)
        yILv = ILData(i).v(:,j);
        
        if j == 1
            PL = 300; % SIA
            PLBR(1) = 210; %140; % Bruhwiler values from 15 october 2019 document, updated with the one found with blockmax axles types
            PLBR(2) = 170; %110;
            PLBR(3) = 170; %95;
            PL12 = 120; % 12 tonnes per axle (load model for 41, 48 and 72)
            PLEU = 300; % EUROCODE
        elseif j == 2
            PL = 200; % SIA
            PLBR(1) = 0;
            PLBR(2) = 150; %85;
            PLBR(3) = 150; %70;
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
        Conc.BR2 = zeros(round(1.2/ILRes)+1,1); % Brühwiler Load Model, 2 axle
        Conc.BR3 = zeros(round(2.4/ILRes)+1,1); % Brühwiler Load Model, 3 axle
        Conc.K41 = zeros(round(7.2/ILRes)+1,1); % KUBA-ST Load Model, 41 truck type
        Conc.K48 = zeros(round(9.1/ILRes)+1,1); % KUBA-ST Load Model, 48 truck type
        Conc.K72 = zeros(round(14/ILRes)+1,1); % KUBA-ST Load Model, 72 truck type
        Conc.KSPTR = zeros(round(12.5/ILRes)+1,1); % KUBA-ST Load Model, combination of special transport with 15 meters gap
        Conc.EURO = zeros(round(1.2/ILRes)+1,1); % EUROCODE Load Model 1
        
        if max(size(Conc.SIA)) == 1
            Conc.SIA(1) = PL*2;
            Conc.BR2(1) = PLBR(2)*2;
            Conc.EURO(1) = PLEU*2;
        else
            Conc.SIA(1) = PL; Conc.SIA(end) = PL;
            Conc.BR2(1) = PLBR(2); Conc.BR2(end) = PLBR(2);
            Conc.EURO(1) = PLEU; Conc.EURO(end) = PLEU;
        end
        
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
        
    end
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
% Distributed loads
qk = 2.5*ones(width(ILData(1).v),1);
qkBR = qk;
qk05 = qk; % kN/m2 with alpha = 1 (needed alpha = 0.5)
qk(1) = 9; % kN/m2
qkBR(1) = 3.6; % kN/m2
qk05(1) = 9; % kN/m2 with alpha = 1 (needed alpha = 0.5)
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
    
    Int.v = IntInf.v(:,i); % SIA CODE
    Int.vBR1 = IntInf.vBR1(:,i); % Brühwiler load model 1
    Int.vBR2 = IntInf.vBR2(:,i); % Brühwiler load model 2
    Int.vBR3 = IntInf.vBR3(:,i); % Brühwiler load model 3
    Int.v41 = IntInf.v41(:,i); % KUBA-ST model for 41 truck
    Int.v48 = IntInf.v48(:,i); % KUBA-ST model for 48 truck
    Int.v72 = IntInf.v72(:,i); % KUBA-ST model for 72 truck
    Int.vSPTR = IntInf.vSPTR(:,i); % KUBA-ST model for special transport truck (generalised truck)
    Int.vEURO = IntInf.v(:,i); % EUROCODE

    % SIA (Swiss code)
    E(i).SIA.Total = 1.5*Alpha*(sum(Maxv.CONV)+Int.v'*qk*LaneWidth);
    E(i).SIA.EQ = Maxv.CONV;
    E(i).SIA.Eq = Int.v.*qk*LaneWidth;
    
    % EUROCODE (European code) - load model n°1
    E(i).EURO.Total = 1.35*Alpha*(sum(Maxv.CONVEURO)+Int.vEURO'*qk*LaneWidth);
    E(i).EURO.EQ = Maxv.CONVEURO;
    E(i).EURO.Eq = Int.vEURO.*qk*LaneWidth;
    
    % CSA (Canadian code)
    
    % AASHTO (American code)
    
    % KUBA-ST (ASTRA load model) 
    E(i).EKUBA_T41.Total = 1.5*Alpha*(sum(Maxv.CONV41)+Int.v41'*qk05*LaneWidth);
    E(i).EKUBA_T41.EQ = Maxv.CONV41;
    E(i).EKUBA_T41.Eq = Int.v41.*qk05*LaneWidth;
    E(i).EKUBA_T48.Total = 1.5*Alpha*(sum(Maxv.CONV48)+Int.v48'*qk05*LaneWidth);
    E(i).EKUBA_T48.EQ = Maxv.CONV48;
    E(i).EKUBA_T48.Eq = Int.v48.*qk05*LaneWidth;
    E(i).EKUBA_T72.Total = 1.5*Alpha*(sum(Maxv.CONV72)+Int.v72'*qk05*LaneWidth);
    E(i).EKUBA_T72.EQ = Maxv.CONV72;
    E(i).EKUBA_T72.Eq = Int.v72.*qk05*LaneWidth;
    E(i).EKUBA_SPTR.Total = 1.5*Alpha*(sum(Maxv.CONVSPTR)+Int.vSPTR'*qk05*LaneWidth);
    E(i).EKUBA_SPTR.EQ = Maxv.CONVSPTR;
    E(i).EKUBA_SPTR.Eq = Int.vSPTR.*qk05*LaneWidth;
    % EKUBA but the load from the KUBA truck in 1st lane is removed
    E(i).EKUBA_SPTR_RmvL1.Total = 1.5*Alpha*(sum(Maxv.CONVSPTR_RmvL1)+Int.vSPTR'*qk05*LaneWidth);
    E(i).EKUBA_SPTR_RmvL1.EQ = Maxv.CONVSPTR_RmvL1;
    E(i).EKUBA_SPTR_RmvL1.Eq = Int.vSPTR.*qk05*LaneWidth;
    
    % EBRU (Professor Eugen Brühwiler load model (EPFL))
    [E(i).EBRU.Total,posi] = max([1.5*(sum(Maxv.CONVBR1)+Int.vBR1'*qkBR*LaneWidth);1.5*(sum(Maxv.CONVBR2)+Int.vBR2'*qkBR*LaneWidth);1.5*(sum(Maxv.CONVBR3)+Int.vBR3'*qkBR*LaneWidth)]);
    E(i).EBRU.LoadMod = append('Model n°',int2str(posi));
    if posi == 1
        E(i).EBRU.EQ = Maxv.CONVBR1;
        E(i).EBRU.Eq = Int.vBR1.*qkBR*LaneWidth;
    elseif posi == 2
        E(i).EBRU.EQ = Maxv.CONVBR2;
        E(i).EBRU.Eq(:,i) = Int.vBR2.*qkBR*LaneWidth;
    elseif posi == 3
        E(i).EBRU.EQ = Maxv.CONVBR3;
        E(i).EBRU.Eq = Int.vBR3.*qkBR*LaneWidth;
    else
        error('error with EBRU');
    end
        
end

end


function PDs = SwapforPlatoonsWIM(PDs,PlatSize,FolDist,TrTyps,PlatPct)
% See if we can perform a swapping routine to form platoons out of WIM traffic

% Initialize
[PDs.Plat, PDs.Trail, PDs.Lead, PDs.Prime, PDs.Swap] = deal(false(height(PDs),1));

% Extract Lane 1... we can put back in later and sort based on DTS
PD = PDs(PDs.LANE == 1,:);

% Get number of each TrTyp in PDs
N = histcounts(categorical(PDs.CLASS),categorical(TrTyps))';

% Get number of each truck type
NumTrTyp = length(TrTyps);

% Get total number of platoons that will result
for i = 1:NumTrTyp
    sumN(i) = N(i)*PlatPct(i)/PlatSize;
end
sumN = round(cumsum(sumN))';

for i = 1:NumTrTyp
    
    % Get all Platoon Candidates
    PrimeCand = PD.CLASS == TrTyps(i) & PD.GW_TOT > prctile(PD.GW_TOT(PD.CLASS == TrTyps(i)),35) & PD.Plat == false & PD.Trail == false & PD.Lead == false;
    % TRY NO 35!! Change after 15/12/21
    %PrimeCand = PD.CLASS == TrTyps(i) & PD.Plat == false & PD.Trail == false & PD.Lead == false;
    PrimeCand(1:4) = false;  PrimeCand(end-PlatSize-1:end) = false;
    
    % New Indexed Swapping Procedure
    
    % Initialize P1s
    P1s = false(height(PD),1);
    
    % Select all lead vehicles at once - here we must make sure that we
    % respect: must be at least 5 vehicles infront before next one
    % PlatCand is the size of PD...
    % Out of PrimeCand == 1 we choose [ round(N(i)*PlatPct(i)/PlatSize) ]
    NumPlats = min([round(N(i)*PlatPct(i)/PlatSize) round(sum(PrimeCand)*.9./PlatSize)]);
    
    if NumPlats == round(sum(PrimeCand)*.9./PlatSize)
        fprintf('\nWarning: Not enough vehicles to form platoon for Truck Type %i, %i formed instead of %i \n',i,round(sum(PrimeCand)*.9./PlatSize),round(N(i)*PlatPct(i)/PlatSize))
        fprintf('(90%% of candidates in Lane 1) \n')
    end
    
    % Actual Index #s
    P1Inds = randsample(find(PrimeCand),NumPlats);
    
    for k = 1:100
    
        % Place into P1s
        P1s(P1Inds) = true;
        
        % Now lets try to figure out how many violations we have
        AA = find(P1s);
        % Get difference between indexes
        BB = [5; diff(AA)];
        
        % Find where < 5 - Vio = Violations (# of violations is the length)
        Vio = find(BB < 5);
        
        if isempty(Vio), break, end
        
        % Set violators to not be lead vehicles
        P1s(AA(Vio)) = false;
        
        % Set anything within 5 to not be candidates (0 is actual P1s)
        for j = 0:PlatSize+2
            PrimeCand(circshift(P1s,j)) = false;
        end
        
        % Redo Vios
        P1Inds = randsample(find(PrimeCand),length(Vio));
    
        if k == 100
            fprintf('Took more than 100 tries for Truck Type %i... could have illegal platoons',i)
        end
    end
    
    % Now P1s are set...
    PD.Plat(P1s) = true;
    PD.Prime(P1s) = true;
    PD.Lead(circshift(P1s,-1)) = true;
    
    % Complete the rest of the plat
    for u = 1:PlatSize-1
        
        % Next we must set P2s
        PlatCand = PD.CLASS == TrTyps(i) & PD.GW_TOT > prctile(PD.GW_TOT(PD.CLASS == TrTyps(i)),35) & PD.Plat == false & PD.Trail == false & PD.Lead == false;
        % Actual Index #s
        PInds = randsample(find(PlatCand),NumPlats);
        
        % We try to place them with their closest match
        PInds = sort(PInds);
        
        % Get Indexes of vehicles following Prime (change circshift in loop)
        Swap = find(circshift(P1s == true,u));
        
        % Label as a swap (reverse of names)
        PD.Swap(PInds) = true;
        PD.Plat(Swap) = true;
        if u == PlatSize-1
            PD.Trail(Swap+1) = true;
        end
        
        % Swap Axle:Class
        PD([PInds, Swap],find(strcmpi(PD.Properties.VariableNames,'AX')):find(strcmpi(PD.Properties.VariableNames,'CLASS'))) = PD([Swap, PInds],find(strcmpi(PD.Properties.VariableNames,'AX')):find(strcmpi(PD.Properties.VariableNames,'CLASS')));
        
        % Redefine PD.DTS of swapped vehicle
        PD.DTS(Swap) = PD.DTS(Swap-1) + seconds((PD.LENTH(Swap-1)/100 + FolDist)./(PD.SPEED(Swap)*0.2777777777778));
        if sum(PD.DTS > datetime(2010,1,1)) > 0
            g = 7;
        end

    end
    %fprintf('Percent Done: %.2f%%\n',100*(sumN(i)./sumN(end)))
end

PDs(PDs.LANE == 1,:) = PD;
PDs = sortrows(PDs,2);

end
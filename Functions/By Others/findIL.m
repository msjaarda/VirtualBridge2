% findIL gets influence lines at the end of all branches

function [NumInfCases, ILData] = findIL(ILs,ILRes,NumLanes,RType) % 05/10/22 Lucas PUN added, PUN = 0 (no PUN infl lanes), PUN = 1 (only PUN infl lanes)

% Initialize
NumInfCases = 0; ILData = struct();
load Misc/ILLib.mat

for i = 1:length(ILs)
     
    TName = ['ILLib.' ILs{i}];
    [NumInfCases, ILData] = recurser(NumInfCases,ILData,TName,ILRes,NumLanes,ILLib,RType);
    
end
end

function [NumInfCases, ILData] = recurser(NumInfCases,ILData,TName,ILRes,NumLanes,ILLib,RType)

% Detect if we are at the end
if isnumeric(eval([TName '(:,1)']))
    ILt = eval([TName '(:,2:end)']);
    if size(ILt,2) == NumLanes || contains(TName,'Box') % Added a check to keep only the right Infl Lane in accordance with the lane configuration Lucas 27/09/22
        NumInfCases = NumInfCases + 1;
        ILx = eval([TName '(:,1)']);
        ILv = eval([TName '(:,2:end)']);
        % Round to refinement level (ILRes)
        RoundedILx = ILx(1):ILRes:ILx(end); RoundedILx = RoundedILx';
        % Now we interpolate the influence lines and populate IL.v
        ILData(NumInfCases).v = interp1(ILx,ILv,RoundedILx);
        % Add to IL.Name
        ILData(NumInfCases).Name = TName;
%         if size(ILData(NumInfCases).v,2) < NumLanes
%             % Can choose to duplicate, or add zeros. The rule will be that
%             % if size(ILData.v{Num.InfCases},2) == 1, we duplicate, but if
%             % it is greater, we add zeros with a notification
%             if size(ILData(NumInfCases).v,2) == 1
%                 ILData(NumInfCases).v = repmat(ILData(NumInfCases).v,1,NumLanes);
%             else
%                 for t = size(ILData(NumInfCases).v,2) + 1:NumLanes % No longer necessary Lucas 27/09/22 
%                     ILData(NumInfCases).v(:,t) = 0;
%                     Turn back on later... was annoying!
%                     fprintf('\nWARNING Lane mismatch for IL: %s, ILs with zeros added',TName)
%                 end
%             end
        %end
    end
else
    FNames = fieldnames(eval(TName));
    if any(strcmp(FNames,'PUN')+strcmp(FNames,'Uni')+strcmp(FNames,'Bi')) %Check if we are in traffic dir, if yes, check if only PUN or no PUN infl lanes have to be loaded + if PUN~=0 or 1, load everything (if needed)
        if strcmp(RType,'PUN')
            FNames = FNames(strcmp(FNames,'PUN'));
        elseif strcmp(RType,'Bi') 
            FNames = FNames(strcmp(FNames,'Bi'));
        elseif strcmp(RType,'Uni')
            FNames = FNames(strcmp(FNames,'Uni'));
        end    
    end
    for j = 1:length(FNames)
        T1Name = [TName '.' FNames{j}];
        [NumInfCases, ILData] = recurser(NumInfCases,ILData,T1Name,ILRes,NumLanes,ILLib,RType);
    end
end
end


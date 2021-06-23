function [FName,Section,Config,Dist,AE,Title,PLoc,X,XTIX] = VBAGBFigDets(Fig)
%AGBFigDets gives figure details according to AGB 2002/005 Report

if Fig == 1
    
    % Assign Figure Name
    FName = 'Figure 4.2 Box Girder, Bidirectional';
    
    % Set Plot Parameters
    Section = 'Box'; % Box, Twin, TwinRed, TwinExp, TwinCon
    Config = 'Bi';   % Bi, Mo
    Dist = 'Split'; % Split, Stand, ExFast, ExSlow
    
    % Set Action Effects
    AE{1} = 'Mn'; AE{2} = 'Mp'; AE{3} = 'V';
    % Set Titles
    Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
    PLoc = 0;
    
elseif Fig == 2
    
    % Assign Figure Name
    FName = 'Figure 4.3 Box Girder, Motorway';
    
    % Set Plot Parameters
    Section = 'Box'; % Box, Twin, TwinRed, TwinExp, TwinCon
    Config = 'Mo';   % Bi, Mo
    Dist{1} = 'ExFast';  Dist{2} = 'Stand'; Dist{3} = 'ExSlow'; % Split, Stand, ExFast, ExSlow
    
    % Set Action Effects
    AE{1} = 'Mp'; AE{2} = 'Mp'; AE{3} = 'Mp';
    % Set Titles
    Title{1} = 'M+ 85%-15%'; Title{2} = 'M+ 96%-4%'; Title{3} = 'M+ 100%-0%';
    PLoc = 0;
    
elseif Fig == 3
    
    % Assign Figure Name
    FName = 'Figure 4.4 Twin Girder, Bidirectional';
    
    % Set Plot Parameters
    Section = 'Twin'; % Box, Twin, TwinRed, TwinExp, TwinCon
    Config = 'Bi';   % Bi, Mo
    Dist = 'Split'; % Split, Stand, ExFast, ExSlow
    
    % Set Action Effects
    AE{1} = 'Mn'; AE{2} = 'Mp'; AE{3} = 'V';
    % Set Titles
    Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
    PLoc = 0;
    
elseif Fig == 4
    
    % Assign Figure Name
    FName = 'Figure 4.5 Twin Girder, Bidirectional';
    
    % Set Plot Parameters
    Section{1} = 'TwinRed'; Section{2} = 'TwinExp'; Section{3} = 'TwinCon'; % Box, Twin, TwinRed, TwinExp, TwinCon
    Config = 'Bi';   % Bi, Mo
    Dist = 'Split'; % Split, Stand, ExFast, ExSlow
    
    % Set Action Effects
    AE = 'Mp';
    % Set Titles
    Title{1} = 'M+ Reduced'; Title{2} = 'M+ Expanded'; Title{3} = 'M+ Concrete';
    PLoc = 0;
    
elseif Fig == 5
    
    % Assign Figure Name
    FName = 'Figure 4.6 Twin Girder, Motorway';
    
    % Set Plot Parameters
    Section = 'Twin'; % Box, Twin, TwinRed, TwinExp, TwinCon
    Config = 'Mo';   % Bi, Mo
    Dist{1} = 'ExFast';  Dist{2} = 'Stand'; Dist{3} = 'ExSlow'; % Split, Stand, ExFast, ExSlow
    
    % Set Action Effects
    AE = 'Mp';
    % Set Titles
    Title{1} = 'M+ 85%-15%'; Title{2} = 'M+ 96%-4%'; Title{3} = 'M+ 100%-0%';
    PLoc = 0;
    
elseif Fig == 6
    
    % Assign Figure Name
    FName = 'Figure 4.7 Multi Girder, Bidirectional';
    
    % Set Plot Parameters
    Section = 'Multi'; % Box, Twin, TwinRed, TwinExp, TwinCon
    Config = 'Bi';   % Bi, Mo
    Dist = 'Split'; % Split, Stand, ExFast, ExSlow
    PLoc{1} = 'P1'; PLoc{2} = 'P2'; PLoc{3} = 'P3';
    
    % Set Action Effects
    AE = 'Mp';
    % Set Titles
    Title{1} = 'Girder 1'; Title{2} = 'Girder 2'; Title{3} = 'Girder 3';
    
elseif Fig == 7
    
    % Assign Figure Name
    FName = 'Figure 4.8 Multi Girder, Bidirectional';
    
    % Set Plot Parameters
    Section = 'Multi'; % Box, Twin, TwinRed, TwinExp, TwinCon
    Config = 'Bi';   % Bi, Mo
    Dist = 'Split'; % Split, Stand, ExFast, ExSlow
    PLoc{1} = 'P1'; PLoc{2} = 'P2'; PLoc{3} = 'P3';
    
    % Set Action Effects
    AE = 'Mn';
    % Set Titles
    Title{1} = 'Girder 1'; Title{2} = 'Girder 2'; Title{3} = 'Girder 3';
    
 elseif Fig == 8
    
    % Assign Figure Name
    FName = 'Figure 4.9 Slab, Bidirectional';
    
    % Set Plot Parameters
    Section = 'SlabSemi'; % Box, Twin, TwinRed, TwinExp, TwinCon
    Config = 'Bi';   % Bi, Mo
    Dist = 'Split'; % Split, Stand, ExFast, ExSlow
    PLoc{1} = 'p1'; PLoc{2} = 'p2'; PLoc{3} = 'p3';
    
    % Set Action Effects
    AE = 'Mp';
    % Set Titles
    Title{1} = 'Lane Edge'; Title{2} = 'Lane Centre'; Title{3} = 'Bridge Centre';
    
elseif Fig == 9
    
    % Assign Figure Name
    FName = 'Figure 4.10 Slab, Bidirectional';
    
    % Set Plot Parameters
    Section = 'SlabSemi'; % Box, Twin, TwinRed, TwinExp, TwinCon
    Config = 'Bi';   % Bi, Mo
    Dist = 'Split'; % Split, Stand, ExFast, ExSlow
    PLoc{1} = 'p1'; PLoc{2} = 'p2'; PLoc{3} = 'p3';
    
    % Set Action Effects
    AE = 'Mn';
    % Set Titles
    Title{1} = 'Lane Edge'; Title{2} = 'Lane Centre'; Title{3} = 'Bridge Centre';
    
elseif Fig == 10
    
    % Assign Figure Name
    FName = 'Figure 4.11 Slab, Bidirectional';
    
    % Set Plot Parameters
    Section{1} = 'SlabPin'; Section{2} = 'SlabSemi'; Section{3} = 'SlabFix';% Box, Twin, TwinRed, TwinExp, TwinCon
    Config = 'Bi';   % Bi, Mo
    Dist = 'Split'; % Split, Stand, ExFast, ExSlow
    PLoc{1} = 'p3'; PLoc{2} = 'p3'; PLoc{3} = 'p3';
    
    % Set Action Effects
    AE = 'Mp';
    % Set Titles
    Title{1} = 'Pinned p3'; Title{2} = 'Semi Rigid p3'; Title{3} = 'Fixed p3';
    
end

% Convert to cellular if necessary
if ~iscell(Section)
    temp = Section; clear Section; [Section{1:3}] = deal(temp);
end
if ~iscell(AE)
    temp = AE; clear AE; [AE{1:3}] = deal(temp);
end
if ~iscell(Dist)
    temp = Dist; clear Dist; [Dist{1:3}] = deal(temp);
end

if Fig < 6
    % X Data
    X = [10:10:80]'; XTIX = 10:10:80;
elseif Fig > 7
    X = [5:5:30]'; XTIX = 5:5:30;
else
    X = [20:10:30]'; XTIX = 15:5:35;
end    


% Storing Custom Input here... just in case we wanna bring it back

% % CUSTOM INPUT ----------
% 
% % Assign Figure Name
% FName = 'Figure 4.2 Box Girder, Bidirectional';
% 
% % Set Plot Parameters
% % Set Section
% Section{1} = 'Box';
% % Box, Twin, TwinRed, TwinExp, TwinCon
% % Set Configuration
% Config{1} = 'Bi';  
% % Bi, Mo
% % Set Distribution
% Dist{1} = 'Split'; % Split, Stand, ExFast, ExSlow
% % Set Action Effects
% AE{1} = 'Mn'; AE{2} = 'Mp'; AE{3} = 'V';
% % Set Titles
% Title{1} = 'M-'; Title{2} = 'M+'; Title{3} = 'V';
% % Set Locations
% Loc{1} = 'GD'; Loc{2} = 'MD'; Loc{3} = 'DD'; Loc{4} = 'DetD';
% Loc{5} = 'CD'; Loc{6} = 'xTCD';
% 
% % CUSTOM INPUT OVER ----------


end


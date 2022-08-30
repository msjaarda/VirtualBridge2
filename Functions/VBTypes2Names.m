function [TrTyps,TrNames] = VBTypes2Names
%VBTypes2Names So that we only have to change them here...

% UPDATE FOR CLASS++

% TrTyps = [11; 12; 22; 23; 111; 11117; 1127; 12117; 122; 11127; 1128; 1138; 1238; 41; 42; 43; 44; 45; 46; 48; 49; 0; 99];
% TrNames = ["11" "12" "22" "23" "111" "1111r" "112r" "1211r" "122" "11112r" "112a" "113a" "123a"...
%     "60t Crane" "6ax 60t" "7ax 72t" "8ax 84t" "9ax 96t" "96t Crane" "72t Crane" "72t Crane" "84t Crane" "No Class" "Empty/Light"]';
% 
% VType = table();
% VType.CLASS = TrTyps;
% VType.Name = TrNames;
% VType.Light = zeros(height(VType),1); VType.Light(end) = 1;
% VType.Class = ones(height(VType),1); VType.Class(14:end) = 0;
% VType.ClassOW = zeros(height(VType),1); VType.ClassOW(14:end-2) = 1;

TrTyps = [11; 12; 22; 23; 111; 11117; 1127; 12117; 122; 11127; 1128; 1138; 1238; 41; 48; 49; 46; 0; 99];
TrNames = ["11" "12" "22" "23" "111" "1111r" "112r" "1211r" "122" "11112r" "112a" "113a" "123a"...
    "60t Crane" "72t Crane" "84t Crane" "96t Crane" "No Class" "Empty/Light"]';

VType = table();
VType.CLASS = TrTyps;
VType.Name = TrNames;
VType.Light = zeros(height(VType),1); VType.Light(end) = 1;
VType.Class = ones(height(VType),1); VType.Class(14:end) = 0;
VType.ClassOW = zeros(height(VType),1); VType.ClassOW(14:end-2) = 1;

VType = sortrows(VType);

end


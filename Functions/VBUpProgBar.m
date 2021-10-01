function [LenPrint] = VBUpProgBar(m, st, v, k, NumSims, NumBatches, Ram, r, LenPrint)
%GETKEYVARS Grabs key variables
t = v;
Num = NumSims;


no = (now-st);
if no > 1
    marker = 'd';
elseif no*24 > 1
    no = no*24;  marker = 'h';
elseif no*24*60 > 1
    no = no*24*60; marker = 'm';
else
    no = no*3600*24; marker = 's';
end

if r == 1

    fprintf('\b\b\b\b\b\b\b\b|%s%.2f%s%s%.1f%s\n','  ',no,marker,', Used Ram : ',Ram,'%');
    LenPrint = 23 + numel(num2str(floor(no))) + numel(num2str(floor(Ram)));
    
else
  
    Remo = repmat('\b',1,LenPrint);
    fprintf(Remo);
    fprintf('|%s%.2f%s%s%.1f%s\n','  ',no,marker,', Used Ram : ',Ram,'%');
    LenPrint = 23 + numel(num2str(floor(no))) + numel(num2str(floor(Ram)));
    
end

end

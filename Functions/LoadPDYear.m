function [PDsy] = LoadPDYear(Name,Year)
%LOADPDYEAR Loads 1 year of PD (function is an effective way of clearing
%what you don't want)

load(Name);

PDsy = PDs(year(PDs.DTS) == Year,:);
            
end


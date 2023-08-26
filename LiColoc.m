function [ICQ LiSum] = LiColoc(BGdropGreen,BGdropRed,AvgGreen,AvgRed,CircSelect);


if CircSelect ~= "Not enough info"
    
    DiffG = (BGdropGreen - AvgGreen(1));
    DiffR = (BGdropRed - AvgRed(1));
    DiffProd = DiffG .* DiffR;
    Pos = DiffProd > 0;
    
    ICQ = (nnz(Pos) / numel(DiffProd)) - 0.5;
    LiSum = nansum(DiffR, 'all') * nansum(DiffG, 'all');
    
else
    
    ICQ = [];
    LiSum = [];
    
end
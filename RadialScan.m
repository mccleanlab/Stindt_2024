function [RadialMetrics] = RadialScan(Center, Radius, CircSelect, BGdropGreen, BGdropRed)


if CircSelect ~= "Not enough info"
        theta = linspace(0, 360, 4*pi*(Radius)); % More than needed to avoid gaps.
        xA = Center(1) + (Radius) * cosd(theta); %"Radius + 10" used to give a little extra buffer region
        yA = Center(2) + (Radius) * sind(theta);
        xy = round([xA', yA']);
        xy = unique(xy,'rows');
        xA = xy(:, 1);
        yA = xy(:, 2);

        
        % Take intensity profiles from center to circumference
        
        I = zeros(length(xy),round(Radius));
        for coord = 1 : length(xy)
            IG0 = improfile(BGdropGreen,...
                [Center(1) xA(coord)],[Center(2) yA(coord)])';
            IG(coord,1:length(IG0)) = IG0;
            IR0 = improfile(BGdropRed,...
                [Center(1) xA(coord)],[Center(2) yA(coord)])';
            IR(coord,1:length(IR0)) = IR0;
%             IRat0 = improfile(Red2Green,...
%                 [Center(1) xA(coord)],[Center(2) yA(coord)])';
%             IRat(coord,1:length(IRat0)) = IRat0;
        end


        %Here's your intensity profile for each. Rows are "r" (radius),
        %columns are theta (angle). As a check, notice that the top row of
        %each (r = 0) has the same value for every column.
        
        IGt = IG'; % Green
        IRt = IR'; % Red
        IRatt = IGt./IRt; % Ratio
        
        %But the closer to the center this is, the more duplicate reads
        %there are (the center intensity gets logged again for every
        %circumference read, e.g.). The interval of "correct" reads to keep
        %is closer to C/(2*pi*r), where C is the outer circumference, whose
        %size is based on the number of xy coordinates.
        
        R = size(IGt,1);
        IGtM = IGt;
        IRtM = IRt;
        
        IGtM(1,2:size(IGtM,2)) = nan; %First row manually, just take 1 reading
        IRtM(1,2:size(IGtM,2)) = nan; %First row manually, just take 1 reading

        for r = 2:(size(IGt,1)-1)
            idx = zeros(1,size(IGt,2)); %Set an empty idx at first

            if r <= R/2
                int = round(R/r); 
                idx(1:int:end) = 1; %Idx every ith to keep
                IGtM(r,~idx) = nan; %Set all else in row to nan
                IRtM(r,~idx) = nan;
            elseif r > R/2
                int = round(1/(1 - (r/R))); %Get the mirror fraction around 1/2
                idx(2:int:end) = 1; %Idx every ith to keep
                IGtM(r,idx==1) = nan; %Now ONLY set ith to nan
                IRtM(r,idx==1) = nan;
            end
            
        end
            
        IRattM = IGtM ./ IRtM;
        %Now do calculations on those profiles (mean, variance, SD) per
        %radius
        
        BasicMetics = [];
        
        for n = 1 : size(IGt,1)
            clear DiffG DiffR DiffProd Pos

            AvgIG(n) = nanmean(IGtM(n,:));
            AvgIR(n) = nanmean(IRtM(n,:));
            AvgIRat(n) = nanmean(IRattM(n,:));

            %Variance per radius
            VarIG(n) = var(IGtM(n,:));
            VarIR(n) = var(IRtM(n,:));
            VarIRat(n) = var(IRattM(n,:));

            %Std Dev per radius
            SDIG(n) = std(IGtM(n,:));
            SDIR(n) = std(IRtM(n,:));
            SDIRat(n) = std(IRattM(n,:));

            %Li Coloc per radius
            DiffG = (IGtM(n,:) - AvgIG(n));
            DiffR = (IGtM(n,:) - AvgIR(n));
            DiffProd = DiffG .* DiffR;
            Pos = DiffProd > 0;
            ICQ(n) = (nnz(Pos) / numel(DiffProd)) - 0.5;
            LiSum(n) = nansum(DiffR, 'all') * nansum(DiffG, 'all');
            
        end
        RadialMetrics.RadialAvgG{1} = AvgIG;
        RadialMetrics.RadialAvgR{1} = AvgIR;
        RadialMetrics.RadialAvgRatio{1} = AvgIRat;

        RadialMetrics.RadialVarG{1} = VarIG;
        RadialMetrics.RadialVarR{1} = VarIR;
        RadialMetrics.RadialVarRatio{1} = VarIRat;

        RadialMetrics.RadialSDG{1} = SDIG;
        RadialMetrics.RadialSDR{1} = SDIR;
        RadialMetrics.RadialSDRatio{1} = SDIRat;
        
        RadialMetrics.RadialICQ{1} = ICQ;
        RadialMetrics.RadialLiSum{1} = LiSum;
        
        %This stuff got put back into the main script
        
%         for row = 1:size(BGdropGreen,1)
%             for col = 1:size(BGdropRed,2)
%                 if round(sqrt((row-Center(1))^2+(col-Center(2))^2)) > Radius
%                     BGdropGreen(row,col) = "NaN";
%                     BGdropRed(row,col) = "NaN";
%                 end
%             end
%         end

%         RadialMetrics.TotGreen(1) = nansum(BGdropGreen,'all');     %Some colony-wide metrics
%         RadialMetrics.TotRed(1) = nansum(BGdropRed,'all');
%         RadialMetrics.AvgGreen(1) = mean(BGdropGreen,'all','omitnan');
%         RadialMetrics.AvgRed(1) = mean(BGdropRed,'all','omitnan');
        
        RadialMetrics = struct2table(RadialMetrics);
        
else 
    
    RadialMetrics = [];
end
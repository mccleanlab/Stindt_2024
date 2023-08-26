function [Center Radius CircSelect] = ColID(fillG, fillR)

    %Combine both signals for maximum binarized object
    fillCombo = fillG + fillR; 

    
    PropsC = regionprops('table',fillCombo,'Area','Centroid','MajorAxisLength');
    %Looks for range of diameters, selects the largest w/in range
    if ~isempty(PropsC.Area(PropsC.MajorAxisLength > 300 & PropsC.MajorAxisLength < 1500)) 
        rowC = max(PropsC.Area(PropsC.MajorAxisLength > 300 & PropsC.MajorAxisLength < 1500));
        Ccenter = PropsC.Centroid(PropsC.Area==rowC,:); %Assigns green circle center
        Cradius = PropsC.MajorAxisLength(PropsC.Area==rowC)/2; %Assigns green radius
    else
        Ccenter = [0,0]; %If no circle diameters w/in range, center & radius set to zero
        Cradius = 0;
    end

    PropsG = regionprops('table',fillG,'Area','Centroid','MajorAxisLength');
    %Looks for range of diameters, selects the largest w/in range
    if ~isempty(PropsG.Area(PropsG.MajorAxisLength > 300 & PropsG.MajorAxisLength < 1500)) 
        rowG = max(PropsG.Area(PropsG.MajorAxisLength > 300 & PropsG.MajorAxisLength < 1500));
        Gcenter = PropsG.Centroid(PropsG.Area==rowG,:); %Assigns green circle center
        Gradius = PropsG.MajorAxisLength(PropsG.Area==rowG)/2; %Assigns green radius
    else
        Gcenter = [0,0]; %If no circle diameters w/in range, center & radius set to zero
        Gradius = 0;
    end
    
    PropsR = regionprops('table',fillR,'Area','Centroid','MajorAxisLength');
    if ~isempty(PropsR.Area(PropsR.MajorAxisLength > 300 & PropsR.MajorAxisLength < 1500))
        rowR = max(PropsR.Area(PropsR.MajorAxisLength > 300 & PropsR.MajorAxisLength < 1500));
        Rcenter = PropsR.Centroid(PropsR.Area==rowR,:);
        Rradius = PropsR.MajorAxisLength(PropsR.Area==rowR)/2;
    else
        Rcenter = [0,0];
        Rradius = 0;
    end

    Rad = vertcat(Gradius,Rradius,Cradius);
    Cen = vertcat(Gcenter,Rcenter,Ccenter);
    
    %Decide which circle to use. Larger preferenced bc it allows wider radius plot

    MaxRad = max(Rad);
    if MaxRad ~= 0
        Chan = find(Rad==MaxRad);
        Pick = max(Chan);
        
        Center = Cen(Pick,:);
        Radius = Rad(Pick);
    else 
        Center = [0,0];
        Radius = 0;
        Pick = 0;
        CircSelect = categorical("Not enough info");
    end
    
    if Pick == 1
        CircSelect = categorical("Green");
    elseif Pick == 2
        CircSelect = categorical("Red");
    elseif Pick == 3
        CircSelect = categorical("Combo");
    end

    %Test if circles are very different from each other & flag it
    
%     if (sqrt((Rcenter(1) - Gcenter(1))^2 + (Rcenter(2) - Gcenter(2))^2) > 60 &&... 
%         sqrt((Rcenter(1) - Gcenter(1))^2 + (Rcenter(2) - Gcenter(2))^2) < 500)
%     
%         CircSelect = categorical("Check circles");
%         
%     end
    
    if Radius==0
        
        CircSelect = categorical("Not enough info");
        
    end
end
function [ newInR, newRaster ] = coarsenInput(raster,InRR,OutRR,planet)
%coarsen input raster if resolution of output is significantly coarser
%Input
%   raster - original raster
%   InRR - input raster reference
%   OutRR - output raster reference
%   planet - designation of which planet
%Output
%   newInR - raster reference of newRaster
%   newRaster - coarsened input raster

thresholdRatio = 3;
inProj = ~contains(class(InRR),'geographic','IgnoreCase',true);
outProj = ~contains(class(OutRR),'geographic','IgnoreCase',true);

if xor(inProj,outProj) % one projected, one geographic
    S = referenceEllipsoid(planet);
    R = S.MeanRadius;
    % convert the geographic pixelsize to meters
    if inProj
        hody = deg2rad(OutRR.CellExtentInLatitude)*R;
        hodx = deg2rad(OutRR.CellExtentInLongitude)*R*(cosd(mean(OutRR.LatitudeLimits)));
        hdy = InRR.CellExtentInWorldY;
        hdx = InRR.CellExtentInWorldX;
    else
        hdy = deg2rad(InRR.CellExtentInLatitude)*R;
        hdx = deg2rad(InRR.CellExtentInLongitude)*R*(cosd(mean(InRR.LatitudeLimits)));
        hody = OutRR.CellExtentInWorldY;
        hodx = OutRR.CellExtentInWorldX;
    end
else
    % either both geographic or both projected
    if inProj
        hody = OutRR.CellExtentInWorldY;
        hodx = OutRR.CellExtentInWorldX;
        hdy = InRR.CellExtentInWorldY;
        hdx = InRR.CellExtentInWorldX;
    else
        hody = OutRR.CellExtentInLatitude;
        hodx = OutRR.CellExtentInLongitude;
        hdy = InRR.CellExtentInLatitude;
        hdx = InRR.CellExtentInLongitude;
    end
end

cellRatio = mean([hodx/hdx hody/hdy]);

%use mapresize or georesize, depending on whether input is projected or
%geographic
if cellRatio>thresholdRatio
    if inProj
        [newRaster,newInR] = mapresize(raster,InRR,2/cellRatio);
    else
        [newRaster,newInR] = georesize(raster,InRR,2/cellRatio);
    end
else
    newInR = InRR;
    newRaster = raster;
end
end
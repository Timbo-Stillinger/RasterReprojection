function [newRR] = transformToAffine(RR,rotationAngle)
%transform raster reference to same with rotation, an affine transform
if contains(class(RR),'geographic','IgnoreCase',true)
    [lat11,lon11] = intrinsicToGeographic(RR,1,1);
    [lat12,lon12] = intrinsicToGeographic(RR,1,2);
    [lat21,lon21] = intrinsicToGeographic(RR,2,1);
    polyin = polyshape([lon11 lon21 lon12 lon11],[lat11 lat21 lat12 lat11]);
    polyout = rotate(polyin,rotationAngle,[lon11 lat11]);
    dx = [polyout.Vertices(3,1)-polyout.Vertices(1,1)...
        polyout.Vertices(2,1)-polyout.Vertices(1,1)];
    dy = [polyout.Vertices(3,2)-polyout.Vertices(1,2)...
        polyout.Vertices(2,2)-polyout.Vertices(1,2)];
    RM = makerefmat(lon11,lat11,dx,dy);
    % new image size
    [latend,lonend] = intrinsicToGeographic(RR,RR.RasterSize(2),...
        RR.RasterSize(1));
    [row,col] = map2pix(RM,lonend,latend);
    newSize = ceil([row col]);
    newRR = refmatToGeoRasterReference(RM,newSize,...
        RR.RasterInterpretation);
elseif contains(class(RR),'mapcells','IgnoreCase',true) ||...
        contains(class(RR),'mappostings','IgnoreCase',true)
    cx = [1 2 2 1];
    cy = [1 1 2 2];
    [x,y] = intrinsicToWorld(RR,cx,cy);
    quadin = polyshape(x,y);
    quadout = rotate(quadin,rotationAngle,[mean(x) mean(y)]);
    dydy = quadout.Vertices(4,2)-quadout.Vertices(1,2);
    dxdy = quadout.Vertices(4,1)-quadout.Vertices(1,1);
    dydx = quadout.Vertices(2,2)-quadout.Vertices(1,2);
    dxdx = quadout.Vertices(2,1)-quadout.Vertices(1,1);
    dx = [dxdy dxdx];
    dy = [dydy dydx];
    RM = makerefmat(x(1),y(1),dx,dy);
    newRR = refmatToMapRasterReference(RM,RR.RasterSize,...
        RR.RasterInterpretation);
else
    error('class input raster reference %s not recognized', class(RR))
end


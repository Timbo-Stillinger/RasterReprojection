function [ B, RRB, varargout] = rasterReprojection(A,InR,InProj,OutProj,varargin )
% [ B, RRB [,fillvalue]] = rasterReprojection(A,InR,InProj,OutProj [,Prop/Value pairs] )
%Reprojects raster from a projected, geographic, or geolocated raster (2D or 3D)
%to a different projection or geographic raster, or to the same projection with
%a different cell size
%(there is no option for output geolocated raster)
%
%INPUT
%   A - input raster (2D or 3D), any numeric type, or logical
%   InR - raster reference (geographic or mapping) for A. InR must be empty
%       if input data are geolocated, in which case lat-lon grids are
%       specified below.
%   InProj - input projection structure, [] if geographic or geolocated
%   OutProj - output projection structure, [] if geographic
% OPTIONAL INPUT
%   name-value pairs, case-insensitive, abbreviations of 3 or more letters
%       generally work, in any order specifying:
%       'planet' - planet name as a character string, case insensitive,
%           defaults to 'earth'
%       'hemisphere' - specify if input image is entire northern or southern
%           hemisphere or whole planet, choices are 'north', 'south', or
%           'both' or 'whole'
%       'method' - interpolation method, options when input data are gridded
%           (i.e. when InR is specified) are 'linear' (default), 'nearest'
%           (fastest), 'cubic', 'spline', or 'makima'; options when input
%           data are geolocated are 'linear' (default), 'nearest', or
%           'natural'
%           (if input data are logical or categorical, interpolation is 'nearest')
%       'rasterref' - output raster reference object, mapping or geographic
%           (this option allows output to exactly match another known
%           image, for example if fitting an elevation model to a satellite
%           image frame)
%       'latitude' and 'longitude' - needed when input A data are geolocated
%           (irregularly spaced cells, for example swath satellite data) -
%           matrices of latitude and longitude, same size as first 2
%           dimensions of A
%           (if input irregularly spaced data are in projected coordinates,
%           convert to lat-lon)
%       'fillvalue' - output value for the cells that do not fall within
%           the projection boundaries, or that are of unknown value in the
%           input raster (defaults are NaN for floating point, minimum for
%           signed integers, maximum for unsigned integers, or you can
%           specify a value)
%
%       The following arguments are ignored if 'rasterref' option is used
%       'pixelsize' - either a 2-element vector specifying height and width
%           of output cells, or a scalar if height=width (default is to
%           approximately match the cell size of input raster)
%       'XLimit' and 'YLimit' - each vectors of length 2: minimum & maximum
%           of output x- and y-coordinates (default is to cover extent of A)
%       'Origin' - 'ul' (upper left, default unless 'rasterref' specified),
%           other options are 'll', 'ur', or 'lr'
%       'adjust' - true or false to adjust x- and y-limits to be a multiple
%           of the pixelsize (default true unless 'rasterref' specified, in
%           which case default is false)
%       'cells' - true or false to specifiy whether inputs are cells or
%           postings (default true, also ignored if InR is a raster
%           reference)
%       'rotate' - in degrees, +ccw, if the output projection is rotated so
%           an affine transformation is needed
%
%OUTPUT
%   B output reprojected raster, same class as input A
%   RRB raster reference object for B
%       (if you want a referencing matrix also, use the optional output)
%Optional OUTPUT, in order
%   fillvalue - especially useful if input data are not floating point and
%       you want to convert them to floating point
%   RefMatrix - referencing matrix

assert(ismatrix(A) || ndims(A)==3,...
    'input array must have 2 or 3 dimensions')

% parse inputs
optargin=length(varargin);
assert (mod(optargin,2)==0,'must be even number of optional arguments')
[InRasterRef,OutRasterRef,planet,method,inLat,inLon,fillvalue] =...
    parseInput(A,InR,InProj,OutProj,varargin{:});

% coarsen input image if output is at significantly coarser
% resolution, so that the output is averaged over multiple input pixels
if ~isempty(InRasterRef) &&...
        contains(InRasterRef.RasterInterpretation,'cells','IgnoreCase',true)
    [InRasterRef,A] = coarsenInput(A,InRasterRef,OutRasterRef,planet);
end

% world coordinates in output image
    [XIntrinsic,YIntrinsic] =...
        meshgrid(1:OutRasterRef.RasterSize(2),1:OutRasterRef.RasterSize(1));
if contains(class(OutRasterRef),'map.rasterref.Map','IgnoreCase',true)
    [XWorld, YWorld] = intrinsicToWorld(OutRasterRef,XIntrinsic,YIntrinsic);
    try % minvtran fails on some projections
        [lat,lon] = minvtran(OutProj,XWorld,YWorld);
    catch % in those cases, use projinv instead, which also fails on some
        [lat,lon] = projinv(OutProj,XWorld,YWorld);
    end
elseif contains(class(OutRasterRef),'Geographic')
    [lat,lon] = intrinsicToGeographic(OutRasterRef,XIntrinsic,YIntrinsic);
else
    error('OutRasterRef class %s unrecognized',class(OutRasterRef))
end

% input coordinates that correspond to all output coordinates
if isempty(InProj) % input grid is lat-lon
    Xq = lon;
    Yq = lat;
else
    try
        [Xq,Yq] = mfwdtran(InProj,lat,lon);
    catch
        [Xq,Yq] = projfwd(InProj,lat,lon);
    end
end

geolocated = ~isempty(inLat); % otherwise geographic or projected

%set fillvalues if not already specified
origType = class(A);
if isempty(fillvalue) % fill value depending on original type
    switch origType
        case {'single','double'}
            fillvalue = NaN;
        case {'uint8','uint16','uint32','uint64'}
            fillvalue = intmax(origType);
        case {'int8','int16','int32','int64'}
            fillvalue = intmin(origType);
        case 'logical'
            fillvalue = false;
        case 'categorical'
            fillvalue = categorical(0);
        otherwise
            error('class ''%s'' not recognized',origType)
    end
else % make sure fill value is within range
    switch origType
        case {'uint8','uint16','uint32','uint64','int8','int16','int32','int64'}
            assert(fillvalue>=intmin(origType) && fillvalue<=intmax(origType),...
                'fillvalue must be >=%d and <=%d for class ''%s''',...
                intmin(origType),intmax(origType),origType)
        case 'logical'
            assert(islogical(fillvalue) || fillvalue==0 || fillvalue==1,...
                'fillvalue must be true or false')
            if ~islogical(fillvalue)
                fillvalue = fillvalue==1;
            end
        case 'categorical'
            if ~iscategorical(fillvalue)
                fillvalue = categorical(fillvalue);
            end
    end
end

% interpolate to output points specified in terms of input coordinates
% set input values to double, setting their fillvalues to NaNs
dblA = double(A);
if ~(strcmpi(origType,'single') || strcmpi(origType,'double'))
    dblA(A==fillvalue) = NaN;
end
if geolocated
    % B will be double
    B = interpolateGeolocatedRaster(inLon,inLat,dblA,Xq,Yq,method);
    t = isnan(B);
    % to keep NaNs from propagating, fill in with values from nearest-neighbor
    if any(t(:)) && ~strcmpi(method,'nearest')
        newB = interpolateGeolocatedRaster(inLon,inLat,dblA,Xq,Yq,'nearest');
        B(t) = newB(t);
    end
else
    B = interpolateRaster(dblA,InRasterRef,Xq,Yq,method);
    t = isnan(B);
    % to keep NaNs from propagating, fill in with values from nearest-neighbor
    if any(t(:)) && ~strcmpi(method,'nearest')
        newB = interpolateRaster(dblA,InRasterRef,Xq,Yq,'nearest');
        B(t) = newB(t);
    end
end

% reset to original type
t = isnan(B);
if contains(origType,'single')
    B = single(B);
    if ~isnan(fillvalue)
        B(t) = fillvalue;
    end
elseif contains(origType,'double')
    if ~isnan(fillvalue)
        B(t) = fillvalue;
    end
elseif contains(origType,'logical')
    x = B>=0.5;
    B = x;
    B(t) = fillvalue;
elseif contains(origType,'categorical')
    B = categorical(round(B));
    B(t) = fillvalue;
else % just signed or unsigned ints are the possibilities
    B = cast(round(B),origType);
    B(t) = fillvalue;
end

% set the optional output fillvalue and referencing matrix
RRB = OutRasterRef;
if nargout>2
    varargout{1} = fillvalue;
    if nargout>3
        varargout{2} = RasterRef2Refmat(RRB);
    end
end
end
function [InRR,OutRR,planet,method,inLat,inLon,fillvalue] =...
    parseReprojectionInput(A,InRR,InS,OutS,varargin)
%parse input values to produce output raster reference object

p = inputParser;
defaultPixelSize = NaN;
defaultXLimit = NaN;
defaultYLimit = defaultXLimit;
defaultAdjust = true;
defaultOrigin = 'ul';
defaultPlanet = 'earth';
defaultHemisphere = '';
defaultMethod = 'linear';
defaultRotation = 0;
defaultRR = [];
defaultCells = true;
inLat = [];
inLon = [];

addRequired(p,'A',@(x) (isnumeric(x) || islogical(x) || iscategorical(x)) &&...
    (ismatrix(x) || ndims(x)==3))
addRequired(p,'InRR',@(x) contains(class(x),'rasterref') || isempty(x))
addRequired(p,'InS',@(x) isstruct(x) || isempty(x))
addRequired(p,'OutS',@(x) isstruct(x) || isempty(x))
addParameter(p,validatestring('pixelsize',{'pix','pixel','pixelsize'}),...
    defaultPixelSize,...
    @(x) isnumeric(x) && (isscalar(x) || length(x)==2))
addParameter(p,validatestring('XLimit',{'xli','xlim','xlimit'}),...
    defaultXLimit,@(x) isnumeric(x) && length(x)==2)
addParameter(p,validatestring('YLimit',{'yli','ylim','ylimit'}),...
    defaultYLimit,@(x) isnumeric(x) && length(x)==2)
addParameter(p,validatestring('adjust',{'adj','adjust'}),...
    defaultAdjust,@islogical)
addParameter(p,validatestring('cells',{'cel','cell','cells'}),...
    defaultCells,@islogical)
addParameter(p,validatestring('origin',{'ori','org','orig','origin'}),...
    defaultOrigin,@(x) ischar(x) &&...
    (strcmpi(x,'ul') || strcmpi(x,'ll') || strcmpi(x,'ur') || strcmpi(x,'lr')))
addParameter(p,validatestring('planet',{'pla','plan','planet'}),...
    defaultPlanet,@ischar)
addParameter(p,validatestring('hemisphere',{'hem','hemi','hemisphere'}),...
    defaultHemisphere,@ischar)
addParameter(p,validatestring('method',{'met','meth','method'}),...
    defaultMethod,@ischar)
addParameter(p,validatestring('rasterref',{'ras','rast','raster','rasterref'}),...
    defaultRR,@(x) contains(class(x),'rasterref'))
addParameter(p,validatestring('latitude',{'lat','lats','latitude'}),...
    [],@(x) isnumeric(x) && ismatrix(x))
addParameter(p,validatestring('longitude',{'lon','lons','longitude'}),...
    [],@(x) isnumeric(x) && ismatrix(x))
addParameter(p,validatestring('fillvalue',{'fil','fill','fillv','fillvalue'}),...
    [],@(x) (isnumeric(x) || iscategorical(x) || islogical(x)) && isscalar(x))
addParameter(p,validatestring('rotate',{'rot','rotate','rotation'}),...
    defaultRotation,@(x) isnumeric(x) && isscalar(x))


parse(p,A,InRR,InS,OutS,varargin{:})

if isempty(p.Results.InS)
    InS = struct([]);
else
    InS = p.Results.InS;
end
if isempty(p.Results.OutS)
    OutS = struct([]);
else
    OutS = p.Results.OutS;
end

% which planet
planet = lower(p.Results.planet);

% interpolation method
if islogical(A) || iscategorical(A)
    method = 'nearest';
else
if isempty(InRR)
    expectedMethods = {'linear','nearest','natural'};
else
    expectedMethods = {'linear','nearest','cubic','spline','makima'};
end
method = validatestring(lower(p.Results.method),expectedMethods);
end

fillvalue = p.Results.fillvalue;

% bounding input coordinates
if strncmpi(p.Results.hemisphere,'nor',3)
    lonlim = [-180 180];
    latlim = [0 90];
elseif strncmpi(p.Results.hemisphere,'sou',3)
    lonlim = [-180 180];
    latlim = [-90 0];
elseif strncmpi(p.Results.hemisphere,'who',3) ||...
        strncmpi(p.Results.hemisphere,'bot',3)
    lonlim = [-180 180];
    latlim = [-90 90];
elseif isempty(p.Results.hemisphere)
    latlim = NaN;
    lonlim = NaN;
else
    error('''hemisphere'' ''%s'' not recognized',p.Results.hemisphere)
end

% check if size matches raster reference
sizeA = [size(p.Results.A,1) size(p.Results.A,2)];
if isempty(p.Results.InRR)
    inLat = double(p.Results.latitude);
    inLon = double(p.Results.longitude);
    assert(isequal(size(inLat),sizeA) &&...
        isequal(size(inLon),sizeA),...
        'if input referencing matrix or raster reference is empty, latitude/longitude of same size as input raster must be specified')
    InRR = [];
    latlim = [min(inLat(:)) max(inLat(:))];
    lonlim = [min(inLon(:)) max(inLon(:))];
else
    assert(isempty(p.Results.latitude) && isempty(p.Results.longitude),...
        'if input raster reference is specified, then input ''latitude'' and ''longitude'' must NOT be specified')
    if isempty(InS) % input is geographic, not projected
        if contains(class(p.Results.InRR),'geographic','IgnoreCase',true)
            InRR = p.Results.InRR;
        elseif contains(class(p.Results.InRR),'mapcells','IgnoreCase',true) ||...
                contains(class(p.Results.inRR),'mappostings','IgnoreCase',true)
            error('if input A is geographic, input raster reference must be Geographic, not MapCells or MapPostings')
        end
        if isnan(latlim)
            latlim = InRR.LatitudeLimits;
            lonlim = InRR.LongitudeLimits;
        end
    else % input is projected
        if contains(class(p.Results.InRR),'mapcells','IgnoreCase',true) ||...
                contains(class(p.Results.InRR),'mappostings','IgnoreCase',true)
            InRR = p.Results.InRR;
        elseif contains(class(p.Results.InRR),'geographic','IgnoreCase',true)
            error('if input A is in a projection, input raster reference must be Map, not Geographic')
        end
        if isnan(latlim)
            [XIntrinsic,YIntrinsic] =...
                meshgrid([1 InRR.RasterSize(2)],[1 InRR.RasterSize(1)]);
            [xWorld,yWorld] = intrinsicToWorld(InRR,XIntrinsic, YIntrinsic);
            try % minvtran fails on some projections if projection structure incorrect
                [latlim,lonlim] = minvtran(InS,xWorld,yWorld);
            catch % in those cases, use projinv instead, which also fails on some, like UTM
                [latlim,lonlim] = projinv(InS,xWorld,yWorld);
            end
        end
    end
end

% if output raster reference specified, other values default to it
% otherwise, parse other inputs
if isempty(p.Results.rasterref)
    % start of rows and cols
    if strncmpi(p.Results.origin,'l',1)
        startcol = 'south';
    else
        startcol = 'north';
    end
    if strcmpi(p.Results.origin(2),'l')
        startrow = 'west';
    else
        startrow = 'east';
    end
    
    % x- and y-limits
    if any(isnan(p.Results.xlimit(:))) || any(isnan(p.Results.ylimit(:)))
        if isempty(OutS)
            xlimit = [min(lonlim(:)) max(lonlim(:))];
            ylimit = [min(latlim(:)) max(latlim(:))];
        else
            try
                [x,y] = mfwdtran(OutS,latlim,lonlim);
            catch
                [x,y] = projfwd(OutS,latlim,lonlim);
            end
            xlimit = [min(x(:)) max(x(:))];
            ylimit = [min(y(:)) max(y(:))];
        end
    end
    if isnan(p.Results.xlimit)
        XLimit = xlimit;
    else
        XLimit = p.Results.xlimit;
    end
    if isnan(p.Results.ylimit)
        YLimit = ylimit;
    else
        YLimit = p.Results.ylimit;
    end
    
    % pixel size
    if isnan(p.Results.pixelsize)
        % pixel size is [height width], default is preserve original pixel
        % size (approximately)
        if isempty(InS) && isempty(OutS)
            if strcmpi(class(InRR),'map.rasterref.GeographicPostingsReference')
                pixelsize = [InRR.SampleSpacingInLatitude InRR.SampleSpacingInLongitude];
            elseif strcmpi(class(InRR),'map.rasterref.GeographicCellsReference')
                pixelsize = [InRR.CellExtentInLatitude InRR.CellExtentInLongitude];
            else
                error('class input raster reference %s not recognized', class(inRR))
            end
        elseif ~isempty(InS) && ~isempty(OutS)
            if isequal(class(InRR),'map.rasterref.MapPostingsReference')
                pixelsize = [InRR.SampleSpacingInWorldY InRR.SampleSpacingInWorldX];
            elseif isequal(class(InRR),'map.rasterref.MapCellsReference')
                pixelsize = [InRR.CellExtentInWorldY InRR.CellExtentInWorldX];
            else
                error('class input raster reference %s not recognized', class(inRR))
            end
        elseif xor(isempty(InS),isempty(OutS))
            pixelsize = [(max(YLimit)-min(YLimit))/size(A,1)...
                (max(XLimit)-min(XLimit))/size(A,2)];
        end
    else
        pixelsize = p.Results.pixelsize;
        if isscalar(pixelsize)
            pixelsize(2) = pixelsize(1);
        end
    end
    
    % adjust limits to multiple of pixel size
    if p.Results.adjust
        [XLimit,YLimit] = adjustReprojectionLimits(XLimit,YLimit,pixelsize);
    end
    nRows = (YLimit(2)-YLimit(1))/pixelsize(1);
    nCols = (XLimit(2)-XLimit(1))/pixelsize(2);
    if nRows~=round(nRows)
        YLimit(2) = YLimit(1)+round(nRows)*pixelsize(1);
        nRows = round((YLimit(2)-YLimit(1))/pixelsize(1));
    end
    if nCols~=round(nCols)
        XLimit(2) = XLimit(1)+round(nCols)*pixelsize(2);
        nCols = round((XLimit(2)-XLimit(1))/pixelsize(2));
    end
    
    % geographic or projected
    if isempty(OutS)
        if p.Results.cells
            OutRR = georefcells(YLimit,XLimit,[nRows nCols],...
                'ColumnsStartFrom',startcol,'RowsStartFrom',startrow);
        else
            OutRR = georefpostings(XLimit,YLimit,[nRows nCols],...
                'ColumnsStartFrom',startcol,'RowsStartFrom',startrow);
        end
    else
        if p.Results.cells
            OutRR = maprefcells(XLimit,YLimit,[nRows nCols],...
                'ColumnsStartFrom',startcol,'RowsStartFrom',startrow);
        else
            OutRR = maprefpostings(XLimit,YLimit,[nRows nCols],...
                'ColumnsStartFrom',startcol,'RowsStartFrom',startrow);
        end
    end
    
    % affine transformation in output raster reference?
    doAffine = p.Results.rotate~=0;
    if doAffine
        OutRR = transformToAffine(OutRR,p.Results.rotate);
    end
    
else
    OutRR = p.Results.rasterref;
end

end
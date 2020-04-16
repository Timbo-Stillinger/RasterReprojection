function RefMat = RasterRef2RefMat(R)
%R = input raster reference object
%RefMat = output referencing matrix
W = worldFileMatrix(R);
RefMat = worldFileMatrixToRefmat(W);
end
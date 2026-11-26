function [slope,relief,aspect] = terrainMetrics(z,cellSize,windowCells)
%TERRAINMETRICS Slope in degrees, square-window relief in m, downslope azimuth.
arguments
    z (:,:) {mustBeNumeric}
    cellSize (1,1) double {mustBePositive}
    windowCells (1,1) double {mustBePositive,mustBeInteger}
end
% Raster rows run southward; gradient's second output follows row direction.
[dzdx,dzdsouth]=gradient(double(z),cellSize);
slope=single(atand(hypot(dzdx,dzdsouth)));
aspect=single(mod(atan2d(-dzdx,dzdsouth),360));
aspect(slope==0)=NaN;
footprint=true(windowCells);
relief=imdilate(single(z),footprint)-imerode(single(z),footprint);
end

function prepare_dem()
%PREPARE_DEM Mosaic and project Copernicus GLO-30 to a 30 m metric grid.
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root,'vendor','TecGeoTools','Topotoolbox3','toolbox')));
files = dir(fullfile(root,'data','raw','*.tif'));
assert(numel(files)==2,'Expected the N37E101 and N37E102 source tiles.');
[z1,r1] = readgeoraster(fullfile(files(1).folder,files(1).name),'OutputType','single');
[z2,r2] = readgeoraster(fullfile(files(2).folder,files(2).name),'OutputType','single');
[lat,~] = intrinsicToGeographic(r1,ones(r1.RasterSize(1),1),(1:r1.RasterSize(1))');
[~,lon1] = intrinsicToGeographic(r1,1:r1.RasterSize(2),ones(1,r1.RasterSize(2)));
[~,lon2] = intrinsicToGeographic(r2,1:r2.RasterSize(2),ones(1,r2.RasterSize(2)));
lon = [lon1 lon2]; z = [z1 z2];
rows = lat>=37.2 & lat<=37.98;
cols = lon>=101.4 & lon<=102.5;
lat = lat(rows); lon=lon(cols); z=z(rows,cols);
assert(all(isfinite(z(:))) && all(z(:)>0),'Missing or invalid source elevations.');
crs = projcrs(32648);
[xc,yc] = projfwd(crs,[min(lat) min(lat) max(lat) max(lat)], ...
    [min(lon) max(lon) min(lon) max(lon)]);
% An interior rectangle avoids artificial nodata wedges at projected edges.
x = (ceil((max(xc([1 3]))+60)/30)*30):30:(floor((min(xc([2 4]))-60)/30)*30);
y = (floor((min(yc([3 4]))-60)/30)*30):-30:(ceil((max(yc([1 2]))+60)/30)*30);
[xx,yy]=meshgrid(x,y);
[qlat,qlon]=projinv(crs,xx,yy);
zi=interp2(lon,flipud(lat),flipud(z),qlon,qlat,'linear');
assert(all(isfinite(zi(:))),'Projection rectangle extends beyond source mosaic.');
DEM=GRIDobj(x,y,single(zi));
R=maprefcells([x(1)-15 x(end)+15],[y(end)-15 y(1)+15],size(zi),'ColumnsStartFrom','north');
R.ProjectedCRS=crs; DEM.georef=R;
DEM.name='Copernicus GLO-30'; DEM.xyunit='m'; DEM.zunit='m';
save(fullfile(root,'data','dem30.mat'),'DEM','-v7.3');
geotiffwrite(fullfile(root,'data','dem30_utm48.tif'),single(zi),R,'CoordRefSysCode',32648);
f=figure('Visible','off','Position',[100 100 1200 850]);
imagesc(lon,lat,z); axis xy equal tight; colorbar; colormap(turbo);
xlabel('Longitude (E)'); ylabel('Latitude (N)'); hold on;
plot(102.05,37.8667,'wp','MarkerFaceColor','k','MarkerSize',12);
text(102.055,37.875,'Published Jiutiaoling location','Color','w');
exportgraphics(f,fullfile(root,'..','.work','dem_overview.png'),'Resolution',150);
close(f);
fprintf('Projected DEM: %d x %d, elevation %.1f--%.1f m\n',size(zi),min(zi(:)),max(zi(:)));
end

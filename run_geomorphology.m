function run_geomorphology()
%RUN_GEOMORPHOLOGY Compute Xiying catchment and Lenglongling segment metrics.
% Run download_data.py once; then run_geomorphology from MATLAB.
root=fileparts(mfilename('fullpath')); addpath(root);
addpath(genpath(fullfile(root,'vendor','TecGeoTools','Topotoolbox3','toolbox')));
cfg=geomorph_config(); out=fullfile(root,'results');
if ~isfolder(out), mkdir(out); end
if ~isfolder(fullfile(root,'figures')), mkdir(fullfile(root,'figures')); end
if ~isfile(fullfile(root,'data','dem30.mat')), prepare_dem(); end
load(fullfile(root,'data','dem30.mat'),'DEM');
if isfile(fullfile(root,'data','hydrology.mat'))
    load(fullfile(root,'data','hydrology.mat'),'FD','A');
else
    FD=FLOWobj(DEM,'preprocess','carve','uselibtt',false);
    A=flowacc(FD);
    save(fullfile(root,'data','hydrology.mat'),'FD','A','-v7.3');
end
Sall=STREAMobj(FD,'minarea',cfg.threshold_km2*1e6/DEM.cellsize^2);
crs=projcrs(cfg.epsg); [x,y]=getcoordinates(DEM);
[xx,yy]=meshgrid(x,y);
[ox,oy]=projfwd(crs,cfg.outlet_lonlat(2),cfg.outlet_lonlat(1));
% Restrict snapping to a main channel, to avoid a nearby minor tributary.
main=A.Z(Sall.IXgrid)*900>500e6;
dd=hypot(Sall.x-ox,Sall.y-oy); dd(~main)=inf;
[snap_m,ii]=min(dd); outlet=Sall.IXgrid(ii);
assert(snap_m<1500,'Outlet main-channel snap exceeds 1.5 km.');
B=drainagebasins(FD,outlet); mask=B.Z>0;
assert(~any(mask(1,:)) && ~any(mask(end,:)) && ~any(mask(:,1)) && ~any(mask(:,end)), ...
    'Catchment intersects DEM edge; enlarge DEM first.');
S=modify(Sall,'upstreamto',outlet); T=trunk(S,A);
z=double(getnal(S,DEM)); zmin=imposemin(S,z,0.00001);
rawKsn=ksn(S,zmin,A,cfg.theta);
k=aggregate(S,rawKsn,'seglength',cfg.reachLength_m);
chi=chitransform(S,A,'a0',cfg.referenceArea_m2,'mn',cfg.theta);
valid=S.distance>=cfg.excludeOutletReach_m;
k(~valid)=NaN;
[slope,relief,aspect]=terrainMetrics(DEM.Z,DEM.cellsize,cfg.reliefWindow_cells);
% Southern boundary of the Xiying watershed is the local main divide.
ridgeY=nan(1,numel(x));
for j=1:numel(x)
    rr=find(mask(:,j),1,'last');
    if ~isempty(rr), ridgeY(j)=y(rr); end
end
boundaryPath=fullfile(root,cfg.mountainBoundaryFile);
assert(isfile(boundaryPath),'Run select_mountain_boundary first to digitize the mountain segment.');
boundary=readtable(boundaryPath);
assert(all(ismember({'X_m','Y_m'},boundary.Properties.VariableNames)), ...
    'Mountain boundary CSV must contain X_m and Y_m columns.');
mountain=inpolygon(xx,yy,boundary.X_m,boundary.Y_m);
ry=interp1(x(isfinite(ridgeY)),ridgeY(isfinite(ridgeY)),x,'linear',NaN);
north=mountain & yy>=ry;
south=mountain & yy<ry;
assert(nnz(mountain & ~(north|south))==0,'Mountain extends beyond the extracted divide.');
assert(nnz(north)>1000 && nnz(south)>1000,'Both mountain flanks must be present.');
ridgeInside=inpolygon(x,ry,boundary.X_m,boundary.Y_m);
rx=x(ridgeInside); ry=ry(ridgeInside);
assert(numel(rx)>100,'Drawn mountain polygon must enclose a substantial divide segment.');

% Independent tributary basins attached to the main trunk.
St=modify(S,'tributaryto',T);
trib=streampoi(St,'outlets','ix');
areas=double(A.Z(trib))*900/1e6;
[areas,order]=sort(areas,'descend'); trib=trib(order);
keep=find(areas>=15,8,'first'); trib=trib(keep);
L=drainagebasins(FD,trib);
sub=table();
for j=1:numel(trib)
    bm=L.Z==j; km=k(bm(S.IXgrid));
    [tl,tn]=projinv(crs,xx(trib(j)),yy(trib(j)));
    row=table(j,nnz(bm)*900/1e6,mean(double(DEM.Z(bm))),mean(double(slope(bm))), ...
        hypsometricIntegral(DEM.Z(bm)),median(km,'omitnan'),tn,tl, ...
        'VariableNames',{'ID','Area_km2','MeanElevation_m','MeanSlope_deg','HI','MedianKsn','OutletLon','OutletLat'});
    sub=[sub;row]; %#ok<AGROW>
end
writetable(sub,fullfile(out,'subbasins.csv'),'Encoding','UTF-8');
regions={mask,mountain,north,south}; names=["西营河山区流域";"冷龙岭山段";"冷龙岭北坡";"冷龙岭南坡"];
summary=table();
for j=1:4
    m=regions{j}; zz=double(DEM.Z(m));
    row=table(names(j),nnz(m)*900/1e6,min(zz),mean(zz),max(zz), ...
        mean(double(slope(m))),median(double(relief(m))), ...
        'VariableNames',{'Region','Area_km2','MinElevation_m','MeanElevation_m','MaxElevation_m','MeanSlope_deg','MedianRelief3030m_m'});
    summary=[summary;row]; %#ok<AGROW>
end
writetable(summary,fullfile(out,'terrain_summary.csv'),'Encoding','UTF-8');
[outlat,outlon]=projinv(crs,xx(outlet),yy(outlet));
% Station snapping is a check, not a substitute for the reservoir outlet.
[sx,sy]=projfwd(crs,cfg.station_lonlat(2),cfg.station_lonlat(1));
[stationSnap,si]=min(hypot(Sall.x-sx,Sall.y-sy));
stationArea=double(A.Z(Sall.IXgrid(si)))*900/1e6;
metadata=struct('epsg',cfg.epsg,'cell_size_m',30,'outlet_lon',outlon,'outlet_lat',outlat, ...
    'outlet_snap_m',snap_m,'catchment_area_km2',nnz(mask)*900/1e6, ...
    'catchment_HI',hypsometricIntegral(DEM.Z(mask)), ...
    'jiutiaoling_area_km2',stationArea,'jiutiaoling_snap_m',stationSnap, ...
    'river_length_km',sum(hypot(S.x(S.ix)-S.x(S.ixc),S.y(S.ix)-S.y(S.ixc)))/1000, ...
    'main_trunk_length_km',max(T.distance)/1000,'matlab_version',version, ...
    'mountain_boundary_file',cfg.mountainBoundaryFile,'config',cfg);
fid=fopen(fullfile(out,'parameters.json'),'w','n','UTF-8');
fprintf(fid,'%s',jsonencode(metadata,PrettyPrint=true)); fclose(fid);

% Channel attributes retain raw and processed elevations and flag excluded reaches.
[rlat,rlon]=projinv(crs,S.x,S.y);
nodes=table(S.IXgrid,rlon,rlat,S.distance,double(A.Z(S.IXgrid))*900,z,zmin,chi,rawKsn,k,valid, ...
    'VariableNames',{'GridIndex','Lon','Lat','DistanceToOutlet_m','Area_m2','ElevationRaw_m', ...
    'ElevationMonotone_m','Chi_m','KsnRaw','KsnReach1000m','ValidForKsn'});
writetable(nodes,fullfile(out,'river_nodes.csv'));

% Candidate convexities on the main trunk; no tectonic attribution is implied.
tz=double(getnal(T,DEM)); tzm=imposemin(T,tz,0.00001);
[~,kp]=knickpointfinder(T,tzm,'tol',cfg.knickpointTolerance_m,'split',false, ...
    'plot',false,'verbose',false,'uselibtt',false);
kpKeep=kp.distance>=cfg.excludeOutletReach_m;
[klat,klon]=projinv(crs,kp.x(kpKeep),kp.y(kpKeep));
knicks=table(klon,klat,kp.distance(kpKeep),kp.z(kpKeep),kp.dz(kpKeep), ...
    'VariableNames',{'Lon','Lat','DistanceToOutlet_m','Elevation_m','ConvexityOffset_m'});
writetable(knicks,fullfile(out,'candidate_knickpoints.csv'));

% Limited sensitivity study; standardized slope uses Aref=1 km2 for all theta.
sensitivity=table();
for threshold=cfg.thresholdSensitivity_km2
    Sn=STREAMobj(FD,'minarea',threshold*1e6/900);
    Sn=modify(Sn,'upstreamto',outlet);
    zn=imposemin(Sn,double(getnal(Sn,DEM)),0.00001);
    for theta=cfg.thetaSensitivity
        kn=aggregate(Sn,ksn(Sn,zn,A,theta),'seglength',cfg.reachLength_m);
        kn(Sn.distance<cfg.excludeOutletReach_m)=NaN;
        medians=nan(height(sub),1);
        labels=L.Z(Sn.IXgrid);
        for b=1:height(sub), medians(b)=median(kn(labels==b),'omitnan'); end
        rho=corr(medians,sub.MedianKsn,'Type','Spearman','Rows','complete');
        row=table(threshold,theta,numel(Sn.x),median(kn,'omitnan'), ...
            median(kn,'omitnan')/cfg.referenceArea_m2^theta,rho, ...
            'VariableNames',{'Threshold_km2','Theta','NodeCount','MedianKsn','MedianSlopeAt1km2','SubbasinRankSpearman'});
        sensitivity=[sensitivity;row]; %#ok<AGROW>
    end
end
writetable(sensitivity,fullfile(out,'sensitivity.csv'));

% Persist GIS rasters, geographic boundaries, and complete calculation state.
layers={uint8(mask),uint8(mountain)+uint8(north),single(slope),single(relief),uint16(L.Z)};
fn={'xiying_basin','lenglongling_flanks_1south_2north','slope_deg','relief_3030m','subbasins'};
for j=1:numel(layers)
    geotiffwrite(fullfile(out,[fn{j} '.tif']),layers{j},DEM.georef,'CoordRefSysCode',cfg.epsg);
end
writeBoundary(mask,x,y,crs,fullfile(out,'xiying_boundary.csv'));
writeBoundary(mountain,x,y,crs,fullfile(out,'lenglongling_boundary.csv'));
[rla,rlo]=projinv(crs,rx(:),ry(:));
writetable(table(rlo,rla,'VariableNames',{'Lon','Lat'}),fullfile(out,'divide_trace.csv'));
save(fullfile(out,'analysis.mat'),'DEM','S','T','B','L','A','k','chi','slope','relief','aspect', ...
    'mask','mountain','north','south','rx','ry','sub','summary','knicks','sensitivity','metadata','cfg','-v7.3');
disp(summary); disp(metadata);
plot_geomorphology();
end

function writeBoundary(mask,x,y,crs,file)
boundaries=bwboundaries(mask,'noholes');
sizes=cellfun(@(b)size(b,1),boundaries); [~,q]=max(sizes); ij=boundaries{q};
bx=x(ij(:,2)); by=y(ij(:,1)); [lat,lon]=projinv(crs,bx(:),by(:));
writetable(table(lon,lat,'VariableNames',{'Lon','Lat'}),file);
end

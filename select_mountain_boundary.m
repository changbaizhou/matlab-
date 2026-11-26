function select_mountain_boundary()
%SELECT_MOUNTAIN_BOUNDARY Digitize the selected mountain segment on the DEM.
root=fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root,'vendor','TecGeoTools','Topotoolbox3','toolbox')));
demFile=fullfile(root,'data','dem30.mat');
if ~isfile(demFile), prepare_dem(); end
d=load(demFile,'DEM'); DEM=d.DEM;
[x,y]=getcoordinates(DEM); z=double(DEM.Z);
f=figure('Color','w','Name','选择冷龙岭山段范围','NumberTitle','off');
imagesc(x/1000,y/1000,z); axis xy equal tight;
colormap(turbo(256)); cb=colorbar; cb.Label.String='高程 / m';
hold on; contour(x/1000,y/1000,z,3000:250:5000,'Color',[.15 .15 .15],'LineWidth',.35);
ridgeFile=fullfile(root,'results','analysis.mat');
if isfile(ridgeFile)
    r=load(ridgeFile,'rx','ry');
    plot(r.rx/1000,r.ry/1000,'w--','LineWidth',1.5);
end
title({'沿山段外围点击，双击结束并闭合多边形', ...
    '范围应包含所选山脊及其南北两侧坡面'},'FontSize',12);
xlabel('UTM 48N东坐标 / km'); ylabel('UTM 48N北坐标 / km');
xlim([205 240]); ylim([4128 4175]); box on;
h=drawpolygon('Color',[.85 .18 .08],'LineWidth',1.5);
wait(h); points=h.Position;
assert(size(points,1)>=4,'Draw at least three boundary vertices.');
boundary=table(points(:,1)*1000,points(:,2)*1000, ...
    'VariableNames',{'X_m','Y_m'});
outDir=fullfile(root,'results');
if ~isfolder(outDir), mkdir(outDir); end
cfg=geomorph_config();
writetable(boundary,fullfile(root,cfg.mountainBoundaryFile));
fprintf('Saved %d boundary vertices to %s\n',height(boundary),cfg.mountainBoundaryFile);
close(f);
end

function verify_outputs()
%VERIFY_OUTPUTS Independent unit, topology, area and output-contract checks.
root=fileparts(mfilename('fullpath')); addpath(root);
addpath(genpath(fullfile(root,'vendor','TecGeoTools','Topotoolbox3','toolbox')));
r=runtests(fullfile(root,'tests')); assertSuccess(r);
d=load(fullfile(root,'results','analysis.mat'));
cs=d.DEM.cellsize; S=d.S;
area=double(d.A.Z(S.IXgrid))*cs^2;
z=imposemin(S,double(getnal(S,d.DEM)),0.00001);
ds=hypot(S.x(S.ix)-S.x(S.ixc),S.y(S.ix)-S.y(S.ixc));
expected=ds.*((d.cfg.referenceArea_m2./area(S.ix)).^d.cfg.theta+ ...
    (d.cfg.referenceArea_m2./area(S.ixc)).^d.cfg.theta)/2;
actual=d.chi(S.ix)-d.chi(S.ixc);
assert(max(abs(actual-expected))<1e-7,'Chi edge integral or area units inconsistent.');
assert(all(actual>0),'Chi must increase upstream.');
assert(all(z(S.ix)>=z(S.ixc)),'Conditioned channel elevations increase downstream.');
assert(nnz(d.north & d.south)==0 && isequal(d.north|d.south,d.mountain), ...
    'Mountain flank masks are not a disjoint partition.');
assert(abs(d.metadata.catchment_area_km2-nnz(d.mask)*cs^2/1e6)<1e-9);
outlets=streampoi(S,'outlets','ix');
assert(numel(outlets)==1,'Selected basin river network must have one outlet.');
assert(abs(double(d.A.Z(outlets))*cs^2-nnz(d.mask)*cs^2)<=cs^2,'Outlet area disagrees with basin mask.');
zz=sort(double(d.DEM.Z(d.mask)),'descend'); n=numel(zz);
h=(zz-min(zz))/(max(zz)-min(zz));
hiIntegral=trapz([0;((1:n)'-.5)/n;1],[1;h;0]);
assert(abs(hiIntegral-d.metadata.catchment_HI)<1e-10,'HI curve disagrees with mean-elevation HI.');
assert(d.metadata.jiutiaoling_area_km2>1000 && d.metadata.jiutiaoling_area_km2<1150, ...
    'Jiutiaoling independent area check outside expected vicinity.');
assert(all(d.sub.HI>=0 & d.sub.HI<=1));
assert(height(d.sensitivity)==9,'Expected 3 x 3 parameter sensitivity cases.');
png=dir(fullfile(root,'figures','*.png')); pdf=dir(fullfile(root,'figures','*.pdf'));
assert(numel(png)==9 && numel(pdf)==9,'Expected exactly nine PNG/PDF figure pairs.');
assert(all([png.bytes]>30000) && all([pdf.bytes]>15000),'Some figures are unexpectedly small.');
for j=1:numel(png)
    im=imfinfo(fullfile(png(j).folder,png(j).name));
    assert(min(im.Width,im.Height)>=1500 && max(im.Width,im.Height)>=2200,'Figure export resolution too low.');
end
files={'run_geomorphology.m','prepare_dem.m','geomorph_config.m','terrainMetrics.m', ...
    'hypsometricIntegral.m','plot_geomorphology.m','verify_outputs.m'};
for j=1:numel(files)
    issues=checkcode(fullfile(root,files{j}),'-struct');
    fprintf('%s: %d static-analysis messages\n',files{j},numel(issues));
    for q=1:numel(issues)
        fprintf('  %d: %s\n',issues(q).line,issues(q).message);
    end
end
fprintf('VERIFIED: %d numerical tests; basin %.3f km2; HI %.6f; chi integration; flank masks; 9 figure pairs.\n', ...
    numel(r),d.metadata.catchment_area_km2,d.metadata.catchment_HI);
end

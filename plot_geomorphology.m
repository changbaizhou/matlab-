function plot_geomorphology()
%PLOT_GEOMORPHOLOGY Export nine reproducible PNG/PDF scientific figures.
root=fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root,'vendor','TecGeoTools','Topotoolbox3','toolbox')));
d=load(fullfile(root,'results','analysis.mat'));
DEM=d.DEM; S=d.S; T=d.T; L=d.L; cfg=d.cfg;
mask=d.mask; mountain=d.mountain; north=d.north; south=d.south;
slope=d.slope; relief=d.relief; rx=d.rx; ry=d.ry;
k=d.k; chi=d.chi; sub=d.sub; knicks=d.knicks;
sensitivity=d.sensitivity; metadata=d.metadata;
[x,y]=getcoordinates(DEM); x=x(:)'/1000; y=y(:)/1000;
figdir=fullfile(root,'figures');
set(groot,'defaultAxesFontName','Microsoft YaHei','defaultTextFontName','Microsoft YaHei', ...
    'defaultAxesFontSize',10,'defaultTextInterpreter','none','defaultLegendInterpreter','none');
blue=[.10 .36 .57]; orange=[.82 .34 .16];
terrain=interp1(linspace(0,1,6),[.24 .40 .34;.43 .57 .37;.70 .72 .48;.80 .74 .59;.66 .60 .54;.97 .96 .91],linspace(0,1,256));

f=newfig(); ax=axes(f); m=mask|mountain;
mapfield(ax,double(DEM.Z),m,[],terrain,'高程 / m',false);
hold(ax,'on'); outlines(ax,mask,blue,1.8); outlines(ax,mountain,orange,1.8);
[sx,sy]=STREAMobj2XY(S); plot(ax,sx/1000,sy/1000,'Color',blue,'LineWidth',.65);
plot(ax,rx/1000,ry/1000,'Color',[.35 .08 .38],'LineWidth',1.6);
crs=projcrs(cfg.epsg); [ox,oy]=projfwd(crs,metadata.outlet_lat,metadata.outlet_lon);
plot(ax,ox/1000,oy/1000,'kv','MarkerFaceColor','k','MarkerSize',7);
text(ax,ox/1000-1,oy/1000-1,'西营水库附近出口','FontSize',9,'VerticalAlignment','top','HorizontalAlignment','right');
text(ax,mean(S.x)/1000,mean(S.y)/1000+6,'西营河山区流域','Color',blue,'FontWeight','bold');
text(ax,mean(rx)/1000,mean(ry)/1000-7,'冷龙岭山段','Color',orange,'FontWeight','bold');
title(ax,{'研究范围：西营河山区流域与冷龙岭南北坡山段','Copernicus GLO-30（2021）｜30 m｜WGS84 / UTM 48N'},'FontSize',11);
decorate(ax); finish(f,'01_研究范围与地形');

f=newfig(); tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); mapfield(ax,slope,mask,[0 55],parula(256),'坡度 / °',true);
title(ax,'西营河山区流域：坡度');
ax=nexttile(tl); mapfield(ax,relief,mask,[0 1800],turbo(256),'3.03 km窗口起伏度 / m',true);
title(ax,'西营河山区流域：局部地形起伏');
finish(f,'02_流域坡度与起伏度');

f=newfig(); tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); labels=double(L.Z); labels(labels==0)=NaN;
mapfield(ax,labels,mask,[.5 height(sub)+.5],turbo(height(sub)),'子流域编号',true);
colorbar(ax,'off'); % Polygon IDs and the matching curve legend encode the categories.
hold(ax,'on'); plot(ax,sx/1000,sy/1000,'Color',[.2 .2 .2],'LineWidth',.4);
for b=1:height(sub)
    [r,c]=find(L.Z==b);
    text(ax,mean(x(c)),mean(y(r)),string(b),'FontSize',12,'FontWeight','bold','BackgroundColor','w','Margin',2);
end
title(ax,'主干流两侧独立支流流域');
ax=nexttile(tl); hold(ax,'on'); colors=turbo(height(sub));
for b=1:height(sub)
    z=sort(double(DEM.Z(L.Z==b)),'descend');
    a=linspace(0,1,numel(z)); h=(z-min(z))/(max(z)-min(z));
    take=unique(round(linspace(1,numel(z),300)));
    plot(ax,a(take),h(take),'Color',colors(b,:),'LineWidth',1.5, ...
        'DisplayName',sprintf('%d：HI=%.3f',b,sub.HI(b)));
end
xlabel(ax,'高于给定高程的面积比例 a / A'); ylabel(ax,'归一化高程 h / H');
axis(ax,[0 1 0 1]); grid(ax,'on'); legend(ax,'Location','northeast','FontSize',9);
title(ax,'面积—高程曲线'); finish(f,'03_子流域与面积高程曲线');

f=newfig(); tl=tiledlayout(f,2,1,'TileSpacing','compact','Padding','compact');
[~,loc]=ismember(T.IXgrid,S.IXgrid); tz=double(getnal(T,DEM));
ax=nexttile(tl); hold(ax,'on');
[td,order]=sort(T.distance);
plot(ax,td/1000,tz(order),'Color',blue,'LineWidth',1.5,'DisplayName','主干流原始高程');
scatter(ax,knicks.DistanceToOutlet_m/1000,knicks.Elevation_m,40,orange,'filled','DisplayName','候选凸起点（容差30 m）');
xline(ax,cfg.excludeOutletReach_m/1000,':','水库影响剔除界线','LabelVerticalAlignment','middle','HandleVisibility','off');
xlabel(ax,'沿河距出口距离 / km'); ylabel(ax,'高程 / m'); grid(ax,'on'); legend(ax,'Location','northwest');
title(ax,'西营河主干流纵剖面（按最大汇水面积追踪）');
ax=nexttile(tl); hold(ax,'on');
scatter(ax,chi/1000,double(getnal(S,DEM)),3,[.73 .77 .79],'filled','DisplayName','支流与主流节点');
[tc,order]=sort(chi(loc));
plot(ax,tc/1000,tz(order),'Color',blue,'LineWidth',1.4,'DisplayName','主干流');
xlabel(ax,'χ / km（A₀=1 km²，θ=0.45）'); ylabel(ax,'高程 / m'); grid(ax,'on'); legend(ax,'Location','northwest');
title(ax,'χ—高程剖面'); finish(f,'04_河流纵剖面与chi剖面');

f=newfig(); tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); background(ax,mask); hold(ax,'on');
scatter(ax,S.x/1000,S.y/1000,5,k,'filled'); colormap(ax,turbo); clim(ax,[0 prctile(k(isfinite(k)),98)]);
cb=colorbar(ax); cb.Label.String='k_{sn} / m^{0.9}（色标上限P98）'; cb.Label.Interpreter='tex';
title(ax,'河道陡峭指数（约1 km河段平均）');
ax=nexttile(tl); background(ax,mask); hold(ax,'on');
scatter(ax,S.x/1000,S.y/1000,5,chi/1000,'filled'); colormap(ax,parula); cb=colorbar(ax); cb.Label.String='χ / km';
title(ax,'χ空间分布（统一出口与参考面积）');
finish(f,'05_河网ksn与chi分布');

f=newfig(); tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl); mapfield(ax,double(DEM.Z),mountain,[2800 5200],terrain,'高程 / m',true);
hold(ax,'on'); plot(ax,rx/1000,ry/1000,'Color',[.35 .08 .38],'LineWidth',1.6);
for q=1:numel(cfg.swath_longitudes)
    lon=cfg.swath_longitudes(q); [cx,~]=projfwd(crs,37.5,lon); cx=cx/1000;
    row=mountain(:,max(1,min(numel(x),round(interp1(x,1:numel(x),cx)))));
    yy=y(row); plot(ax,[cx cx],[min(yy) max(yy)],'k--','LineWidth',.8);
    text(ax,cx,min(yy)-.7,sprintf('P%d',q),'HorizontalAlignment','center','FontWeight','bold');
end
text(ax,mean(x(any(north,1))),mean(y(any(north,2))),'北坡','FontWeight','bold','BackgroundColor','w');
text(ax,mean(x(any(south,1))),mean(y(any(south,2))),'南坡','FontWeight','bold','BackgroundColor','w');
title(ax,'冷龙岭山段、主分水岭与横断条带');
ax=nexttile(tl); mapfield(ax,slope,mountain,[0 55],parula(256),'坡度 / °',true);
hold(ax,'on'); plot(ax,rx/1000,ry/1000,'Color',[.35 .08 .38],'LineWidth',1.6);
title(ax,'山段坡度及分水岭两侧'); finish(f,'06_冷龙岭南北坡地形');

f=newfig(); tl=tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
vals={double(DEM.Z),double(slope),double(relief)};
edges={2500:100:5300,0:2:70,0:75:2400}; labs={'高程 / m','坡度 / °','3.03 km窗口起伏度 / m'};
for q=1:3
    ax=nexttile(tl); hold(ax,'on');
    histogram(ax,vals{q}(north),edges{q},'Normalization','probability','DisplayStyle','stairs','EdgeColor',blue,'LineWidth',1.6,'DisplayName','北坡');
    histogram(ax,vals{q}(south),edges{q},'Normalization','probability','DisplayStyle','stairs','EdgeColor',orange,'LineWidth',1.6,'DisplayName','南坡');
    xlabel(ax,labs{q}); ylabel(ax,'像元比例'); grid(ax,'on');
    title(ax,sprintf('北坡均值 %.1f；南坡均值 %.1f',mean(vals{q}(north)),mean(vals{q}(south))),'FontSize',10);
    legend(ax,'Location','best');
end
title(tl,'冷龙岭南北坡地形分布（按主分水岭分区）'); finish(f,'07_冷龙岭南北坡参数对比');

f=newfig(); tl=tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
swathRows=table();
for q=1:3
    [cx,~]=projfwd(crs,37.5,cfg.swath_longitudes(q));
    ridge=interp1(rx,ry,cx,'linear');
    cols=abs(x*1000-cx)<=cfg.swathWidth_m/2;
    zz=double(DEM.Z(:,cols)); mm=mountain(:,cols); zz(~mm)=NaN;
    rr=any(mm,2); dist=(y(rr)*1000-ridge)/1000;
    zz=zz(rr,:); lo=min(zz,[],2,'omitnan'); hi=max(zz,[],2,'omitnan'); av=mean(zz,2,'omitnan');
    ax=nexttile(tl); hold(ax,'on');
    fill(ax,[dist;flipud(dist)],[lo;flipud(hi)],[.82 .85 .81],'EdgeColor','none','DisplayName','条带最小—最大高程');
    plot(ax,dist,av,'Color',blue,'LineWidth',1.5,'DisplayName','条带平均高程');
    xline(ax,0,'--','主分水岭','HandleVisibility','off');
    xlabel(ax,'距主分水岭的北向距离 / km（南负、北正）'); ylabel(ax,'高程 / m'); grid(ax,'on');
    title(ax,sprintf('P%d：3 km宽横向条带',q));
    swathRows=[swathRows;table(repmat(q,numel(dist),1),dist,lo,av,hi, ...
        'VariableNames',{'ProfileID','NorthwardDistance_km','Min_m','Mean_m','Max_m'})]; %#ok<AGROW>
end
ax=nexttile(tl); hold(ax,'on');
% Along-strike strips are parallel to the endpoint-to-endpoint divide trend.
p0=[rx(1) ry(1)]; direction=[rx(end)-rx(1),ry(end)-ry(1)]; direction=direction/norm(direction);
[XX,YY]=meshgrid(x*1000,y*1000);
along=((XX-p0(1))*direction(1)+(YY-p0(2))*direction(2))/1000;
across=abs(-(XX-p0(1))*direction(2)+(YY-p0(2))*direction(1));
sel=mountain & across<=cfg.swathWidth_m/2;
bins=0:cfg.swathBin_m/1000:max(along(sel)); centers=(bins(1:end-1)+bins(2:end))/2;
mn=nan(size(centers)); mx=mn; av=mn;
for q=1:numel(centers)
    z=double(DEM.Z(sel & along>=bins(q) & along<bins(q+1)));
    if ~isempty(z), mn(q)=min(z); mx(q)=max(z); av(q)=mean(z); end
end
good=isfinite(av); centers=centers(good); mn=mn(good); mx=mx(good); av=av(good);
fill(ax,[centers fliplr(centers)],[mn fliplr(mx)],[.82 .85 .81],'EdgeColor','none');
plot(ax,centers,av,'Color',orange,'LineWidth',1.5); xlabel(ax,'沿山段走向距西端距离 / km'); ylabel(ax,'高程 / m'); grid(ax,'on');
title(ax,'沿走向3 km宽条带'); title(tl,{'冷龙岭山段条带地形剖面','阴影：最小—最大高程；实线：平均高程'},'FontSize',12);
finish(f,'08_冷龙岭条带地形剖面');
writetable(swathRows,fullfile(root,'results','cross_swath_profiles.csv'));
writetable(table(centers(:),mn(:),av(:),mx(:),'VariableNames',{'Distance_km','Min_m','Mean_m','Max_m'}),fullfile(root,'results','along_swath_profile.csv'));

f=newfig(); tl=tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact'); colors=lines(3);
ax=nexttile(tl); hold(ax,'on');
for q=1:3
    t=sensitivity(sensitivity.Theta==cfg.thetaSensitivity(q),:);
    plot(ax,t.Threshold_km2,t.MedianSlopeAt1km2,'-o','Color',colors(q,:),'LineWidth',1.5,'DisplayName',sprintf('θ=%.2f',cfg.thetaSensitivity(q)));
end
xlabel(ax,'河网起始汇水面积 / km²'); ylabel(ax,'参考面积1 km²处的标准化坡降'); grid(ax,'on');legend(ax,'Location','best');
title(ax,'不同阈值与参考凹度的指标变化');
ax=nexttile(tl); hold(ax,'on');
for q=1:3
    t=sensitivity(sensitivity.Theta==cfg.thetaSensitivity(q),:);
    plot(ax,t.Threshold_km2,t.SubbasinRankSpearman,'-o','Color',colors(q,:),'LineWidth',1.5,'DisplayName',sprintf('θ=%.2f',cfg.thetaSensitivity(q)));
end
ylim(ax,[-1.05 1.05]); xlabel(ax,'河网起始汇水面积 / km²'); ylabel(ax,'子流域排序与基准结果的Spearman相关');
title(ax,'子流域陡峭程度排序稳定性'); grid(ax,'on');legend(ax,'Location','best'); finish(f,'09_参数敏感性');
fprintf('Exported 9 PNG + 9 PDF figures.\n');

    function f=newfig()
        f=figure('Visible','off','Color','w','Position',[80 80 1320 820]);
    end
    function finish(f,name)
        drawnow;
        exportgraphics(f,fullfile(figdir,[name '.png']),'Resolution',300);
        exportgraphics(f,fullfile(figdir,[name '.pdf']),'ContentType','auto','Resolution',300);
        close(f);
    end
    function mapfield(ax,z,m,limits,cm,label,doOutline)
        ii=1:3:numel(y); jj=1:3:numel(x); zz=double(z(ii,jj));
        h=imagesc(ax,x(jj),y(ii),zz); axis(ax,'xy'); axis(ax,'equal');
        set(h,'AlphaData',double(m(ii,jj)&isfinite(zz))); set(ax,'Color',[.95 .95 .93]);
        colormap(ax,cm); if ~isempty(limits), clim(ax,limits); end
        cb=colorbar(ax); cb.Label.String=label;
        bounds(ax,m); xlabel(ax,'UTM 48N东坐标 / km'); ylabel(ax,'UTM 48N北坐标 / km'); box(ax,'on');
        if doOutline, hold(ax,'on'); outlines(ax,m,[.2 .2 .2],.7); end
    end
    function background(ax,m)
        ii=1:3:numel(y); jj=1:3:numel(x);
        z=double(DEM.Z(ii,jj)); grayLevel=.98-.18*min(max((z-2000)/3300,0),1);
        h=imagesc(ax,x(jj),y(ii),repmat(grayLevel,1,1,3));
        set(h,'AlphaData',double(m(ii,jj))); axis(ax,'xy'); axis(ax,'equal');
        set(ax,'Color',[.97 .97 .96]); bounds(ax,m); hold(ax,'on');
        outlines(ax,m,[.25 .25 .25],.7); box(ax,'on');
        xlabel(ax,'UTM 48N东坐标 / km'); ylabel(ax,'UTM 48N北坐标 / km');
    end
    function bounds(ax,m)
        [r,c]=find(m); xlim(ax,[min(x(c))-1.5 max(x(c))+1.5]); ylim(ax,[min(y(r))-1.5 max(y(r))+1.5]);
    end
    function outlines(ax,m,color,width)
        bb=bwboundaries(m,'noholes');
        for n=1:numel(bb)
            ij=bb{n}; plot(ax,x(ij(:,2)),y(ij(:,1)),'Color',color,'LineWidth',width);
        end
    end
    function decorate(ax)
        xl=xlim(ax); yl=ylim(ax); xx=xl(1)+.08*diff(xl); yy=yl(1)+.08*diff(yl);
        plot(ax,[xx xx+10],[yy yy],'k','LineWidth',2);
        text(ax,xx+5,yy+1,'10 km','HorizontalAlignment','center');
        xx=xl(1)+.10*diff(xl); yy=yl(2)-.13*diff(yl);
        quiver(ax,xx,yy,0,4,0,'k','LineWidth',1.5,'MaxHeadSize',.8);
        text(ax,xx,yy+4.7,'网格北','HorizontalAlignment','center','FontSize',9);
    end
end

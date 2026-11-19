function [figs,graphics] = plotFIMResults(input,outputDir)
%PLOTFIMRESULTS 四种阵列的 CDF、参数表与真实部署样本；支持 MAT 直接重绘。
% plotFIMResults('results/.../study.mat') 不调用信道、不重新抽样、不重新优化。
if ischar(input)
    path=input; data=load(path);
    if isfield(data,'study'), input=data.study;
    elseif isfield(data,'result'), input=data.result;
    else, error('MAT must contain study or result.'); end
    if nargin<2 || isempty(outputDir), outputDir=fullfile(fileparts(path),'replot'); end
end
if isfield(input,'cases'), cases=input.cases; else, cases={input}; end
assert(iscell(cases) && ~isempty(cases),'No simulation cases.');
cfg=cases{1}.cfg;
if ~exist('outputDir','var') || isempty(outputDir), outputDir=cfg.output_dir; end
if exist(outputDir,'dir')~=7, mkdir(outputDir); end
figs={}; names={}; graphics.files={}; graphics.failures={};
try
    if strcmp(cases{1}.experiment,'cdf')
        fields={'throughput_bps','snr','sinr','sumRate'};
        captions={'UE THP (Mbit/s), CP included',...
            'Rate-equivalent SNR (dB), thermal noise only',...
            'Rate-equivalent SINR (dB), interference included','Sum SE (bit/s/Hz), CP included'};
        stems={'THP_CDF','SNR_CDF','SINR_CDF','sum_SE_CDF'};
        for j=1:numel(fields)
            figs{end+1}=metricCDF(cases,fields{j},captions{j}); names{end+1}=stems{j};
        end
        figs{end+1}=gainCDF(cases,false); names{end+1}='gain_over_fixed_CDF';
        figs{end+1}=gainCDF(cases,true); names{end+1}='subarray_minus_PCA_CDF';
    else
        figs{end+1}=snrSweep(cases); names{end+1}='SNR_sweep_SE_THP';
    end
catch err
    graphics.failures{end+1}=['Performance: ' err.message];
    warning('FIM:PerformancePlot','Performance plots failed: %s',err.message);
end
for c=1:numel(cases)
    r=cases{c}; K=r.cfg.K; stem=fullfile(outputDir,sprintf('K%d',K));
    try
        [~,tables]=reportFIMParameters(r,stem,r.cfg.plot.show_parameter_table,false);
        graphics.files=[graphics.files {[stem '_parameters.txt'],[stem '_parameters.csv']}];
        for j=1:numel(tables)
            figs{end+1}=tables{j}; names{end+1}=sprintf('K%d_parameters_%02d',K,j);
        end
    catch err
        graphics.failures{end+1}=sprintf('K%d parameters: %s',K,err.message);
    end
    if r.cfg.plot.show_geometry
        try
            figs{end+1}=arrayGeometry(r); names{end+1}=sprintf('K%d_four_array_geometry',K);
        catch err
            graphics.failures{end+1}=sprintf('K%d geometry: %s',K,err.message);
        end
    end
    if r.cfg.plot.show_scene
        indices=unique(r.cfg.plot.scene_indices(:).');
        indices=indices(indices<=numel(r.scene_cell));
        for mc=indices
            try
                figs{end+1}=deploymentScene(r,mc); names{end+1}=sprintf('K%d_BS_UE_scene_MC%04d',K,mc);
            catch err
                graphics.failures{end+1}=sprintf('K%d scene %d: %s',K,mc,err.message);
            end
        end
    end
end
for j=1:numel(figs)
    set(figs{j},'Visible',cfg.plot.figure_visible);
    if ~cfg.save_figures, continue; end
    path=fullfile(outputDir,[names{j} '.png']);
    try
        drawnow;
        device='-dpng';
        if exist('OCTAVE_VERSION','builtin')~=0
            device='-dpngcairo'; % 避免 gnuplot 的 PostScript 路径丢失中文字符。
        end
        print(figs{j},path,device,sprintf('-r%d',cfg.plot.export_dpi));
        assert(exist(path,'file')==2,'PNG was not created.');
        graphics.files{end+1}=path;
    catch firstError
        try
            saveas(figs{j},path);
            assert(exist(path,'file')==2,'PNG fallback did not create a file.');
            graphics.files{end+1}=path;
        catch
            graphics.failures{end+1}=[names{j} ': ' firstError.message];
        end
    end
end
graphics.figure_names=names; graphics.requested_figures=numel(figs);
if ~isempty(graphics.failures)
    for j=1:numel(graphics.failures), warning('FIM:GraphicsIncomplete','%s',graphics.failures{j}); end
end
if cfg.plot.close_after_save
    for j=1:numel(figs), close(figs{j}); end
end
end

function fig=metricCDF(cases,field,caption)
[labels,colors,styles,fields]=schemes(); C=numel(cases); values=cell(C,4); allValues=[];
for c=1:C
    r=cases{c}; assert(isscalar(r.snr_dB),'CDF requires one SNR point per case.');
    for j=1:4
        d=r.(fields{j}); value=d.(field)(:,:,1); value=value(:);
        assert(isreal(value) && ~any(isnan(value)) && all(value>=0),'Invalid CDF samples.');
        if strcmp(field,'throughput_bps'), value=value/1e6; end
        values{c,j}=value; allValues=[allValues;value];
    end
end
if any(strcmp(field,{'snr','sinr'}))
    finitePositive=allValues(allValues>0 & isfinite(allValues));
    if isempty(finitePositive), display=[-100 100];
    else, db=10*log10(finitePositive); display=[floor(min(db))-10 ceil(max(db))+10]; end
    for i=1:numel(values)
        v=10*log10(values{i});
        countZero=sum(v==-Inf); countInf=sum(v==Inf);
        v(v==-Inf)=display(1); v(v==Inf)=display(2); values{i}=v;
        if countZero || countInf
            fprintf('%s display only: %d zero -> %.1f dB, %d Inf -> %.1f dB.\n',field,countZero,display(1),countInf,display(2));
        end
    end
end
bounds=commonBounds(values); fig=newFigure([field ' empirical CDF'],C);
for c=1:C
    subplot(ceil(C/min(C,3)),min(C,3),c);
    for j=1:4, drawCDF(values{c,j},colors(j,:),styles{j}); hold on; end
    r=cases{c};
    if strcmp(field,'sumRate')
        ttl=sprintf('K=%d | %d MC samples',r.cfg.K,r.cfg.MC);
    else
        ttl=sprintf('K=%d | MC=%d | %d UE samples',r.cfg.K,r.cfg.MC,r.cfg.MC*r.cfg.K);
    end
    finishAxes(caption,'Empirical CDF',ttl);
    legend(labels{:},'Location','best'); ylim([0 1]); xlim(bounds);
end
end

function fig=gainCDF(cases,againstPCA)
[labels,colors,styles,fields]=schemes(); C=numel(cases); values=cell(C,3); excluded=zeros(C,1);
for c=1:C
    r=cases{c};
    if againstPCA
        values{c,1}=r.fim.sumRate(:,1)-r.pca.sumRate(:,1);
    else
        baseline=r.fixed.sumRate(:,1); valid=baseline>0; excluded(c)=sum(~valid);
        for j=2:4
            d=r.(fields{j}); values{c,j-1}=100*(d.sumRate(valid,1)-baseline(valid))./baseline(valid);
        end
    end
end
bounds=commonBounds(values); fig=newFigure('Paired gain empirical CDF',C);
for c=1:C
    subplot(ceil(C/min(C,3)),min(C,3),c);
    if againstPCA
        drawCDF(values{c,1},colors(4,:),'-');
        caption='Subarray minus PCA sum SE (bit/s/Hz)';
        detail='positive favors subarray';
    else
        for j=1:3, drawCDF(values{c,j},colors(j+1,:),styles{j+1}); hold on; end
        if any(cellfun(@(v) ~isempty(v),values(c,:)))
            legend(labels{2:4},'Location','best');
        else
            text(0.5,0.5,'No positive fixed baseline: relative gain undefined','Units','normalized','HorizontalAlignment','center');
        end
        caption='Paired sum-SE gain over fixed (%)';
        detail=sprintf('%d undefined zero baselines omitted only here',excluded(c));
    end
    finishAxes(caption,'Empirical CDF',sprintf('K=%d | %d paired MC | %s',cases{c}.cfg.K,cases{c}.cfg.MC,detail));
    ylim([0 1]); xlim(bounds);
end
end

function fig=snrSweep(cases)
[labels,colors,styles,fields]=schemes(); C=numel(cases); fig=newFigure('SNR sweep',C);
for c=1:C
    r=cases{c}; rateField='sumRate'; thpField='sumThroughput_bps';
    if strcmp(r.experiment,'worstuser'), rateField='worstRate'; thpField='worstThroughput_bps'; end
    for p=1:2
        subplot(2,C,(p-1)*C+c);
        for j=1:4
            d=r.(fields{j});
            if p==1, v=mean(d.(rateField),1); else, v=mean(d.(thpField),1)/1e6; end
            plot(r.snr_dB,v,'Color',colors(j,:),'LineStyle',styles{j},'LineWidth',1.7); hold on;
        end
        yl='SE (bit/s/Hz)'; if p==2, yl='THP (Mbit/s)'; end
        finishAxes('Per-UE input transmit SNR (dB)',yl,sprintf('K=%d | %s | %s',r.cfg.K,r.experiment,r.cfg.power.snr_reference));
        legend(labels{:},'Location','best');
    end
end
end

function fig=arrayGeometry(r)
cfg=r.cfg; mc=min(cfg.plot.geometry_mc_index,cfg.MC); s=min(cfg.plot.frequency_snr_index,numel(r.snr_dB));
g=cell(1,4); g{1}=r.geometry.fixed;
pcfg=cfg; pcfg.deform.family='pca'; pcfg.comparison.include_pca=false;
[g{2},valid]=fimGeometry('build',r.pca_model,reshape(r.pca_parameters(mc,:,s),1,[]),pcfg);
assert(valid,'Saved PCA parameters are infeasible.');
selection=reshape(r.orientation_selection(mc,:,s),[],1);
g{3}=fimGeometry('orientation',g{1},r.orientation_states,selection,cfg);
[g{4},valid]=fimGeometry('build',r.deformation_model,reshape(r.deformation_parameters(mc,:,s),1,[]),cfg);
assert(valid,'Saved subarray parameters are infeasible.');
[labels,colors]=schemes(); allPos=[];
for j=1:4, allPos=[allPos;g{j}.positions_m*100]; end
lo=min(allPos,[],1); hi=max(allPos,[],1); span=max(max(hi-lo),1);
lo=lo-0.08*span; hi=hi+0.08*span;
for d=1:3
    if hi(d)-lo(d)<0.25*span, center=(hi(d)+lo(d))/2; lo(d)=center-0.125*span; hi(d)=center+0.125*span; end
end
fig=figure('Name',sprintf('Four array forms | K%d MC%d',cfg.K,mc),'NumberTitle','off','Color','w','Position',[40 40 1300 950]);
for j=1:4
    subplot(2,2,j); p=g{j}.positions_m*100;
    plot3(p(:,1),p(:,2),p(:,3),'o','Color',colors(j,:),'MarkerFaceColor',colors(j,:),'MarkerSize',3); hold on;
    ids=unique(round(linspace(1,cfg.N,min(cfg.N,32)))); arrowLength=0.5*cfg.c0/cfg.fc_Hz*100;
    quiver3(p(ids,1),p(ids,2),p(ids,3),arrowLength*g{j}.normals(ids,1),arrowLength*g{j}.normals(ids,2),arrowLength*g{j}.normals(ids,3),0,'Color',colors(j,:));
    axis equal; xlim([lo(1) hi(1)]); ylim([lo(2) hi(2)]); zlim([lo(3) hi(3)]); view(38,24); grid on; box on;
    xlabel('Local x (cm)'); ylabel('Local y (cm)'); zlabel('Local z (cm)'); title(labels{j},'Interpreter','none');
end
end

function fig=deploymentScene(r,mc)
scene=r.scene_cell{mc}; cfg=r.cfg; p=scene.positions_m; bs=scene.bs.position_m;
assert(size(p,1)==cfg.K && size(p,2)==3,'Invalid saved UE deployment.');
fig=figure('Name',sprintf('BS UE deployment | K%d MC%d',cfg.K,mc),'NumberTitle','off','Color','w','Position',[40 50 1400 680]);
blue=[0.05 0.35 0.72]; red=[0.72 0.08 0.12];
subplot(1,2,1);
ue=plot(p(:,1),p(:,2),'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',7); hold on;
bh=plot(bs(1),bs(2),'s','Color',red,'MarkerFaceColor',red,'MarkerSize',10);
if strcmp(cfg.users.mode,'random')
    az=cfg.bs.azimuth_deg+linspace(cfg.users.azimuth_offsets_deg(1),cfg.users.azimuth_offsets_deg(2),120);
    for radius=cfg.users.radius_range_m
        plot(bs(1)+radius*sind(az),bs(2)+radius*cosd(az),':','Color',[0.6 0.65 0.7]);
    end
end
for k=1:cfg.K
    plot([bs(1) p(k,1)],[bs(2) p(k,2)],':','Color',[0.7 0.73 0.77]);
    text(p(k,1),p(k,2),sprintf('  UE%d',k),'Interpreter','none');
end
axis equal; grid on; box on; xlabel('Global x (m)'); ylabel('Global y (m)');
title(sprintf('Saved deployment | K=%d | MC=%d',cfg.K,mc)); legend([bh ue],{'BS','UE'},'Location','best');
subplot(1,2,2);
ue=plot3(p(:,1),p(:,2),p(:,3),'o','Color',blue,'MarkerFaceColor',blue,'MarkerSize',7); hold on;
bh=plot3(bs(1),bs(2),bs(3),'s','Color',red,'MarkerFaceColor',red,'MarkerSize',10);
plot3([bs(1) bs(1)],[bs(2) bs(2)],[cfg.bs.ground_height_m bs(3)],'-','Color',red,'LineWidth',1.8);
for k=1:cfg.K
    plot3([bs(1) p(k,1)],[bs(2) p(k,2)],[bs(3) p(k,3)],':','Color',[0.7 0.73 0.77]);
    text(p(k,1),p(k,2),p(k,3),sprintf('  UE%d',k),'Interpreter','none');
end
normal=scene.bs.R_global_from_array*[0;1;0]; scale=0.15*max(1,max(scene.distance_3d_m));
quiver3(bs(1),bs(2),bs(3),scale*normal(1),scale*normal(2),scale*normal(3),0,'Color',red,'LineWidth',2);
axis equal; grid on; box on; view(38,24); xlabel('Global x (m)'); ylabel('Global y (m)'); zlabel('Global z (m)');
title(sprintf('BS height %.4g m | downtilt %.4g deg | dotted links are schematic',cfg.bs.height_m,cfg.bs.downtilt_deg));
legend([bh ue],{'BS','UE'},'Location','best');
end

function [labels,colors,styles,fields]=schemes()
labels={'Fixed planar','PCA','Pattern reconfigurable','Subarray pose reconfigurable'};
colors=[0.30 0.32 0.35;0.05 0.35 0.72;0.78 0.48 0.05;0.72 0.08 0.12];
styles={'-','--',':','-.'}; fields={'fixed','pca','orientation','fim'};
end

function drawCDF(values,color,style)
if isempty(values), return; end
values=sort(values(:)); n=numel(values); x=[values(1);values]; F=[0;(1:n)'/n];
stairs(x,F,'Color',color,'LineStyle',style,'LineWidth',1.8);
end

function bounds=commonBounds(values)
v=[]; for i=1:numel(values), v=[v;values{i}(:)]; end
v=v(isfinite(v));
if isempty(v), bounds=[-1 1]; return; end
lo=min(v); hi=max(v); margin=max(1e-3,0.05*(hi-lo));
if hi==lo, margin=max(0.1,0.05*abs(lo)); end
bounds=[lo-margin hi+margin];
end

function fig=newFigure(name,C)
fig=figure('Name',name,'NumberTitle','off','Color','w','Position',[50 60 min(1800,650*min(C,3)) min(1000,600*ceil(C/3))]);
end

function finishAxes(xl,yl,ttl)
grid on; box on; xlabel(xl,'Interpreter','none'); ylabel(yl,'Interpreter','none'); title(ttl,'Interpreter','none'); set(gca,'FontSize',10);
end

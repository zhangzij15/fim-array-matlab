function [rows,figs] = reportFIMParameters(result,outputStem,makeFigure,printConsole)
%REPORTFIMPARAMETERS 实际仿真参数表；长表分页，另保存完整 TXT/CSV。
if nargin<2, outputStem = ''; end
if nargin<3, makeFigure = true; end
if nargin<4, printConsole = true; end
cfg = result.cfg; figs = {}; rows = {};
M = size(result.fixed.rate,1); K = cfg.K;
lambda = cfg.c0/cfg.fc_Hz;
s = min(cfg.plot.frequency_snr_index,numel(result.snr_dB));
mc = min(cfg.plot.geometry_mc_index,M);
rows(end+1,:) = {'比较方案','固定平面 / PCA / 方向图重构 / 子阵位姿重构'};
rows(end+1,:) = {'统计规模',sprintf('K=%d；MC=%d；全用户样本数=%d；每条 CDF 按同一 SNR 单独统计',K,M,K*M)};
rows(end+1,:) = {'基站阵列',sprintf('%d x %d = %d 阵元；载频 %.5g GHz；波长 %.4g mm',cfg.Nh,cfg.Nv,cfg.N,cfg.fc_Hz/1e9,lambda*1e3)};
rows(end+1,:) = {'阵元间距 / 初始口径',sprintf('水平 %.4g lambda；纵向 %.4g lambda；中心跨度 %.4g x %.4g m',cfg.dh_lambda,cfg.dv_lambda,(cfg.Nh-1)*cfg.dh_lambda*lambda,(cfg.Nv-1)*cfg.dv_lambda*lambda)};
rows(end+1,:) = {'基站部署',sprintf('BS=%s m；方位 %.4g deg；下倾 %.4g deg',mat2str(cfg.bs.position_m,6),cfg.bs.azimuth_deg,cfg.bs.downtilt_deg)};
rows(end+1,:) = {'用户位置分布',sprintf('%s；水平距离 %s m；相对扇区角 %s deg；UE 高度 %.4g m',cfg.users.mode,mat2str(cfg.users.radius_range_m),mat2str(cfg.users.azimuth_offsets_deg),cfg.users.height_m)};
rows(end+1,:) = {'OFDM 资源',sprintf('%d RB x %d 子载波/RB；SCS %.5g kHz；实际带宽 %.6g MHz',cfg.ofdm.num_rb,cfg.ofdm.subcarriers_per_rb,cfg.ofdm.scs_Hz/1e3,cfg.ofdm.active_bandwidth_Hz/1e6)};
rows(end+1,:) = {'频率计算',sprintf('%s；计算 %d / %d 实际活动频点；FFT=%d',cfg.ofdm.frequency_mode,cfg.ofdm.num_tones,cfg.ofdm.num_active,cfg.ofdm.nfft)};
rows(end+1,:) = {'CP 与速率系数',sprintf('CP/Tu=%.6g；CP=%.6g us；计入速率的系数=%.6g',cfg.ofdm.cp_ratio,cfg.ofdm.cp_duration_s*1e6,cfg.ofdm.cp_efficiency)};
rows(end+1,:) = {'接收机 / 优化目标',[upper(cfg.receiver) ' / ' cfg.opt_mode]};
rows(end+1,:) = {'输入 SNR / 功率口径',[mat2str(result.snr_dB,7) ' dB；' cfg.power.snr_reference]};
if strcmp(cfg.power.snr_reference,'total_over_tone_noise')
    tx = cfg.noiseVar*10.^(result.snr_dB/10);
    rows(end+1,:) = {'每 UE 总功率',sprintf('%s W；各 UE 独立等功率分配到自身实际占用子载波',mat2str(tx,6))};
end
rows(end+1,:) = {'单子载波噪声',sprintf('%.7g；与发射功率、外部协方差使用相同功率单位',cfg.noiseVar)};
rows(end+1,:) = {'路径损耗',sprintf('启用=%d；指数 %.5g；参考距离 %.4g m；参考损耗 %s dB',cfg.pathloss.enabled,cfg.pathloss.exponent,cfg.pathloss.reference_distance_m,mat2str(cfg.pathloss.reference_loss_dB,7))};
rows(end+1,:) = {'信道与方向图','CDL 簇级标量模型；完美 CSI；所有形态复用同一传播路径；方向图只归一化一次'};
if result.pattern.is_test_pattern
    rows(end+1,:) = {'输入标记','合成验证输入 TEST ONLY；本次结果不代表实测天线增益'};
end
rows(end+1,:) = {'方向图来源',cfg.pattern.file};
rows(end+1,:) = {'方向图主瓣诊断',sprintf('处理后 HPBW：水平 %.5g deg，俯仰 %.5g deg',result.pattern.diagnostics.hpbw_az_deg,result.pattern.diagnostics.hpbw_el_deg)};
rows(end+1,:) = {'子阵划分 / 参数数',sprintf('%d x %d 块；每块 %d x %d 阵元；%d 个 d/yaw/tilt 参数',cfg.subarray.m,cfg.subarray.n,cfg.Nh/cfg.subarray.m,cfg.Nv/cfg.subarray.n,3*cfg.subarray.m*cfg.subarray.n)};
rows(end+1,:) = {'子阵位移网格',[mat2str(cfg.subarray.displacement_grid_lambda,5) ' lambda；沿原阵面法线']};
rows(end+1,:) = {'子阵角度网格',['yaw=' mat2str(cfg.subarray.yaw_grid_deg) ' deg；tilt=' mat2str(cfg.subarray.tilt_grid_deg) ' deg']};
rows(end+1,:) = {'子阵搜索与间距',sprintf('坐标搜索 %d 轮；最小中心间距 %.4g lambda',cfg.subarray.sweeps,cfg.subarray.min_spacing_lambda)};
rows(end+1,:) = {'方向图重构',sprintf('位置固定；2N=%d 个独立角度；yaw=%s deg；tilt=%s deg；搜索 %d 轮',2*cfg.N,mat2str(cfg.orientation.yaw_grid_deg),mat2str(cfg.orientation.tilt_grid_deg),cfg.orientation.sweeps)};
rows(end+1,:) = {'SNR 与 SINR','SNR 只计热噪声；SINR 计多用户与外部干扰；二者均在同一实际合并器后计算，并先频率汇总'};
rows(end+1,:) = {'PCA 控制与区域',sprintf('2N=%d 个 x/z 参数；%s；方向图朝向固定',2*cfg.N,cfg.pca.region_mode)};
rows(end+1,:) = {'PCA 位移网格',['x=' mat2str(cfg.pca.x_grid_lambda,5) '；z=' mat2str(cfg.pca.z_grid_lambda,5) ' lambda']};
rows(end+1,:) = {'PCA 搜索与间距',sprintf('坐标搜索 %d 轮；最小中心间距 %.4g lambda',cfg.pca.sweeps,cfg.pca.min_spacing_lambda)};
rows(end+1,:) = {'公平比较条件','同一 MC 的 UE / 多径 / 衰落 / RB / 功率 / 接收机；每种机构一套参数供全部用户和频点共用'};
rows(end+1,:) = {'搜索结论范围','保留所选目标的固定基线；有限网格局部搜索；不保证每个 UE 或某种机构始终占优'};
rows(end+1,:) = {'机械模型范围','理想刚体或移动相位中心及点间距；未计算实体碰撞、互耦、机械能耗'};
rows(end+1,:) = {'速率统计','先对各频点计算 log2(1+SINR)，再频率加权；吞吐量计入 CP；等效 SINR 不含 CP'};
rows(end+1,:) = {'随机种子',sprintf('首批 seed=%g；分批研究需同时使用保存的 batch_seeds / batch_MC 复现',cfg.seed)};
if isfield(result,'batch_seeds')
    rows(end+1,:) = {'分批统计',sprintf('%d 个批次；每批最多 %d MC；均值按实际批次 MC 数加权',numel(result.batch_seeds),result.batch_MC)};
end
rows(end+1,:) = {'所需 CP 最大值',sprintf('固定 / 子阵 / PCA = %.6g / %.6g / %.6g us',max(result.requiredCP_s(:))*1e6,max(result.fimRequiredCP_s(:))*1e6,max(result.pcaRequiredCP_s(:))*1e6)};
rows(end+1,:) = {'所选实例',sprintf('下列位置与机构参数来自 MC=%d、SNR=%.6g dB',mc,result.snr_dB(s))};
triples = reshape(result.deformation_parameters(mc,:,s),3,[]).';
for b=1:size(triples,1)
    rows(end+1,:) = {sprintf('子阵 B%d',b),sprintf('d=%.6g lambda；yaw=%.6g deg；tilt=%.6g deg',triples(b,1),triples(b,2),triples(b,3))};
end
scene = result.scene_cell{mc};
for k=1:K
    rows(end+1,:) = {sprintf('UE%d 实际位置',k),sprintf('%s m；3D 距离 %.5g m；实际施加路损 %.5g dB',mat2str(scene.positions_m(k,:),6),scene.distance_3d_m(k),scene.pathloss_dB(k))};
end
if printConsole
    for i=1:size(rows,1), fprintf('%s: %s\n',rows{i,1},rows{i,2}); end
end
if ~isempty(outputStem)
    parent=fileparts(outputStem); if ~isempty(parent) && exist(parent,'dir')~=7, mkdir(parent); end
    fid=fopen([outputStem '_parameters.txt'],'w','n','UTF-8');
    assert(fid~=-1,'Cannot write parameter TXT.'); cleanup=onCleanup(@() fclose(fid));
    for i=1:size(rows,1), fprintf(fid,'%s: %s\n',rows{i,1},rows{i,2}); end
    fprintf(fid,'\n完整执行配置：\n'); writeValue(fid,cfg,'cfg');
    clear cleanup;
    fid=fopen([outputStem '_parameters.csv'],'w','n','UTF-8');
    assert(fid~=-1,'Cannot write parameter CSV.'); cleanup=onCleanup(@() fclose(fid));
    fwrite(fid,uint8([239 187 191]),'uint8'); fprintf(fid,'"参数","值"\r\n');
    for i=1:size(rows,1)
        fprintf(fid,'"%s","%s"\r\n',strrep(rows{i,1},'"','""'),strrep(rows{i,2},'"','""'));
    end
    clear cleanup;
end
if ~makeFigure, return; end
perPage=22; pages=ceil(size(rows,1)/perPage); font=chooseFont(cfg.plot.font_name);
for page=1:pages
    fig=figure('Name',sprintf('K%d actual simulation parameters %d/%d',K,page,pages),...
        'NumberTitle','off','Color','w','Position',[40 40 1800 1100],'Visible',cfg.plot.figure_visible);
    ax=axes('Parent',fig,'Position',[0 0 1 1]); axis(ax,[0 1 0 1]); axis(ax,'off'); hold(ax,'on');
    rectangle('Position',[0 0.91 1 0.09],'FaceColor',[0.58 0.04 0.08],'EdgeColor','none');
    text(ax,0.025,0.955,sprintf('固定 / PCA / 方向图 / 子阵位姿 | 实际参数 | K=%d',K),...
        'FontName',font,'FontSize',23,'FontWeight','bold','Color','w','Interpreter','none');
    text(ax,0.025,0.875,sprintf('MC=%d | %s | 第 %d / %d 页',M,upper(cfg.receiver),page,pages),...
        'FontName',font,'FontSize',13,'Color',[0.3 0.3 0.3],'Interpreter','none');
    indices=(page-1)*perPage+1:min(page*perPage,size(rows,1)); height=0.8/perPage;
    for j=1:numel(indices)
        y=0.85-j*height; fill=[0.96 0.966 0.973]; if mod(j,2)==0, fill=[1 1 1]; end
        rectangle('Position',[0.02 y 0.96 height],'FaceColor',fill,'EdgeColor',[0.85 0.86 0.88]);
        i=indices(j);
        text(ax,0.03,y+height/2,wrapText(rows{i,1},26,2),'FontName',font,'FontSize',12,...
            'FontWeight','bold','VerticalAlignment','middle','Interpreter','none');
        text(ax,0.235,y+height/2,wrapText(rows{i,2},112,2),'FontName',font,'FontSize',11,...
            'VerticalAlignment','middle','Interpreter','none');
    end
    text(ax,0.025,0.025,'图中超长字段以 ... 省略；完整配置、位置和参数见配套 TXT / CSV。',...
        'FontName',font,'FontSize',11,'Color',[0.4 0.4 0.4],'Interpreter','none');
    figs{end+1}=fig;
end
end

function out=wrapText(in,maxUnits,maxLines)
chars=regexp(in,'[\s\S]','match'); parts={}; piece=''; units=0;
for i=1:numel(chars)
    width=1+any(double(chars{i})>127);
    if units+width>maxUnits && ~isempty(piece), parts{end+1}=piece; piece=''; units=0; end
    piece=[piece chars{i}]; units=units+width;
end
parts{end+1}=piece;
if numel(parts)>maxLines
    parts=parts(1:maxLines); tail=regexp(parts{end},'[\s\S]','match');
    while ~isempty(tail) && numel([tail{:}])>maxUnits-3, tail(end)=[]; end
    parts{end}=[tail{:} '...'];
end
out=strjoin(parts,sprintf('\n'));
end

function font=chooseFont(preferred)
font=preferred;
try
    names=listfonts;
    candidates={preferred,'Noto Sans CJK SC','Noto Sans CJK JP','WenQuanYi Zen Hei','Arial Unicode MS','DejaVu Sans'};
    for i=1:numel(candidates)
        if any(strcmpi(names,candidates{i})), font=candidates{i}; return; end
    end
catch
end
end

function writeValue(fid,v,prefix)
if isstruct(v)
    names=fieldnames(v);
    for i=1:numel(names), writeValue(fid,v.(names{i}),[prefix '.' names{i}]); end
elseif iscell(v)
    for i=1:numel(v), writeValue(fid,v{i},sprintf('%s{%d}',prefix,i)); end
elseif isnumeric(v) || islogical(v)
    fprintf(fid,'%s = %s\n',prefix,mat2str(v,17));
elseif isa(v,'function_handle')
    fprintf(fid,'%s = @%s\n',prefix,func2str(v));
elseif ischar(v)
    fprintf(fid,'%s = %s\n',prefix,v);
end
end

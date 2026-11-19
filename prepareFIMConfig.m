function cfg = prepareFIMConfig(cfg)
%PREPAREFIMCONFIG 验证输入并派生频率网格；不对各阵列/频点单独归一化信道。
% 此函数可重复调用；全频与分组模式始终按 active_indices 重新生成网格。
% 使用 configFIMSimulation() 创建的新 cfg 已显式带 rb/ratio/物理部署选项。
cfg = fillMissingDefaults(cfg,configFIMSimulation());
cfg = fillMissingDefaults(cfg,struct('irc',struct('external_covariance',[])));
pcaDefaults = struct('x_grid_lambda',[-0.5 0 0.5],'z_grid_lambda',[-0.5 0 0.5],...
    'min_spacing_lambda',0.4,'sweeps',1,'region_mode','per_element',...
    'x_bounds_lambda',[],'z_bounds_lambda',[]);
cfg = fillMissingDefaults(cfg,struct('pca',pcaDefaults,'comparison',struct('include_pca',false)));
checkBool(cfg.comparison.include_pca,'comparison.include_pca');
checkGrid(cfg.orientation.yaw_grid_deg,'orientation.yaw_grid_deg');
checkGrid(cfg.orientation.tilt_grid_deg,'orientation.tilt_grid_deg');
assert(all(abs(cfg.orientation.yaw_grid_deg(:))<=180) && all(abs(cfg.orientation.tilt_grid_deg(:))<=90),'Invalid orientation angle ranges.');
checkInteger(cfg.orientation.sweeps,0,Inf,'orientation.sweeps');
ints = {'Nh','Nv','K','MC','progress_every','fast_MC'};
for i = 1:numel(ints)
    checkInteger(cfg.(ints{i}),1,Inf,ints{i});
end
checkInteger(cfg.fast_sweeps,0,Inf,'fast_sweeps');
checkInteger(cfg.deform.cache_limit,1,Inf,'deform.cache_limit');
checkInteger(cfg.seed,0,2^32-1,'seed');
positive = {'dh_lambda','dv_lambda','fc_Hz','c0','noiseVar'};
for i = 1:numel(positive)
    checkScalar(cfg.(positive{i}),0,Inf,positive{i},true);
end
checkScalar(cfg.snr_dB,-Inf,Inf,'snr_dB',false);
assert(isnumeric(cfg.snr_dB_list) && isreal(cfg.snr_dB_list) && isvector(cfg.snr_dB_list) && ~isempty(cfg.snr_dB_list) && all(isfinite(cfg.snr_dB_list(:))),'Invalid snr_dB_list.');
checkRange(cfg.user_az_range,-Inf,Inf,'user_az_range');
checkRange(cfg.user_el_range,-90,90,'user_el_range');
checkScalar(cfg.objective_tol,0,Inf,'objective_tol',false);
checkScalar(cfg.check.norm_tol,0,1,'check.norm_tol',true);
checkScalar(cfg.check.objective_tol,0,Inf,'check.objective_tol',false);
checkScalar(cfg.check.zf_rcond_switch,0,1,'check.zf_rcond_switch',true);
checkScalar(cfg.check.svd_rel_tol,0,1,'check.svd_rel_tol',false);
checkScalar(cfg.opt.gram_rcond_min,0,1,'opt.gram_rcond_min',true);
assert(cfg.K<=cfg.Nh*cfg.Nv,'This simulator requires K <= N.');
assert(any(strcmp(cfg.opt_mode,{'sum_rate','max_min_rate'})),'Invalid opt_mode.');
assert(strcmp(cfg.fim.mode,'deformation'),'Only physical subarray/PCA geometry is supported.');
assert(any(strcmp(cfg.deform.family,{'subarray','pca'})),'Invalid deform.family.');
assert(~cfg.comparison.include_pca || (strcmp(cfg.fim.mode,'deformation') && strcmp(cfg.deform.family,'subarray')),...
    'Three-array comparison requires deformation/subarray mode; use simulateFIM or run_fim_comparison.');
assert(any(strcmp(cfg.output.mat_format,{'default','-v7','-v7.3'})),'Invalid output.mat_format.');
if strcmp(cfg.fim.mode,'deformation') && strcmp(cfg.deform.family,'subarray')
    checkInteger(cfg.subarray.m,1,cfg.Nh,'subarray.m');
    checkInteger(cfg.subarray.n,1,cfg.Nv,'subarray.n');
    assert(mod(cfg.Nh,cfg.subarray.m)==0 && mod(cfg.Nv,cfg.subarray.n)==0,...
        'Subarray counts m/n must divide Nh/Nv exactly.');
    checkGrid(cfg.subarray.displacement_grid_lambda,'subarray.displacement_grid_lambda');
    checkGrid(cfg.subarray.yaw_grid_deg,'subarray.yaw_grid_deg');
    checkGrid(cfg.subarray.tilt_grid_deg,'subarray.tilt_grid_deg');
    assert(all(abs(cfg.subarray.yaw_grid_deg(:))<=180),'Subarray yaw angles must lie in [-180,180] degrees.');
    assert(all(abs(cfg.subarray.tilt_grid_deg(:))<=90),'Subarray tilt angles must lie in [-90,90] degrees.');
    checkInteger(cfg.subarray.sweeps,0,Inf,'subarray.sweeps');
    checkScalar(cfg.subarray.min_spacing_lambda,0,Inf,'subarray.min_spacing_lambda',false);
end
if cfg.comparison.include_pca || (strcmp(cfg.fim.mode,'deformation') && strcmp(cfg.deform.family,'pca'))
    checkGrid(cfg.pca.x_grid_lambda,'pca.x_grid_lambda');
    checkGrid(cfg.pca.z_grid_lambda,'pca.z_grid_lambda');
    checkInteger(cfg.pca.sweeps,0,Inf,'pca.sweeps');
    checkScalar(cfg.pca.min_spacing_lambda,0,Inf,'pca.min_spacing_lambda',false);
    assert(any(strcmp(cfg.pca.region_mode,{'per_element','shared'})),...
        'pca.region_mode must be per_element or shared.');
    if ~isempty(cfg.pca.x_bounds_lambda)
        checkRange(cfg.pca.x_bounds_lambda,-Inf,Inf,'pca.x_bounds_lambda');
    end
    if ~isempty(cfg.pca.z_bounds_lambda)
        checkRange(cfg.pca.z_bounds_lambda,-Inf,Inf,'pca.z_bounds_lambda');
    end
end
assert(any(strcmp(cfg.opt.score_mode,{'gram_verified','exact'})),'Invalid opt.score_mode.');
if strcmp(cfg.receiver,'mmse_irc')
    cfg.receiver = 'mmse-irc';
end
assert(ischar(cfg.receiver) && any(strcmp(cfg.receiver,{'zf','mmse-irc'})),...
    'receiver must be zf or mmse-irc.');
cfg.irc.external_covariance = validateExternalCovariance(cfg.irc.external_covariance,cfg.Nh*cfg.Nv);
assert(strcmp(cfg.throughput_mode,'shannon'),'This release implements Shannon spectral efficiency only.');
assert(any(strcmp(cfg.power.snr_reference,{'per_subcarrier','total_over_tone_noise'})),'Invalid power.snr_reference.');
checkBool(cfg.fast_mode,'fast_mode');
checkBool(cfg.front_only,'front_only');
checkBool(cfg.save_results,'save_results');
checkBool(cfg.save_figures,'save_figures');
checkBool(cfg.pattern.use_phase,'pattern.use_phase');
assert(~cfg.pattern.use_phase,'Pattern phase needs calibration; this release uses power gain only.');
checkBool(cfg.output.save_per_tone,'output.save_per_tone');
checkBool(cfg.output.save_realizations,'output.save_realizations');
checkBool(cfg.output.save_geometry,'output.save_geometry');
checkBool(cfg.plot.show_geometry,'plot.show_geometry');
checkBool(cfg.plot.show_scene,'plot.show_scene');
checkBool(cfg.plot.enabled,'plot.enabled');
checkBool(cfg.plot.show_parameter_table,'plot.show_parameter_table');
checkBool(cfg.plot.show_per_tone_cdf,'plot.show_per_tone_cdf');
checkBool(cfg.plot.show_frequency_profile,'plot.show_frequency_profile');
checkInteger(cfg.plot.frequency_snr_index,1,Inf,'plot.frequency_snr_index');
checkInteger(cfg.plot.geometry_mc_index,1,Inf,'plot.geometry_mc_index');
checkInteger(cfg.plot.export_dpi,50,600,'plot.export_dpi');
checkBool(cfg.plot.close_after_save,'plot.close_after_save');
assert(any(strcmp(cfg.plot.figure_visible,{'on','off'})),'plot.figure_visible must be on or off.');
assert(ischar(cfg.cdl.generator) && ~isempty(cfg.cdl.generator),'cdl.generator must be a named function.');
checkGrid(cfg.plot.scene_indices,'plot.scene_indices');
assert(all(cfg.plot.scene_indices(:)>=1 & cfg.plot.scene_indices(:)==fix(cfg.plot.scene_indices(:))),'Scene indices must be positive integers.');
assert(ischar(cfg.pattern.file) && ischar(cfg.output_dir) && ischar(cfg.cdl.profile_mat),'File paths must be character vectors.');

assert(any(strcmp(cfg.pattern.alignment,{'manual','peak'})),'Invalid pattern.alignment.');
checkScalar(cfg.pattern.expected_az_hpbw_deg,0,360,'pattern.expected_az_hpbw_deg',true);
checkScalar(cfg.pattern.expected_el_hpbw_deg,0,180,'pattern.expected_el_hpbw_deg',true);

% 基站全局部署和 UE 模式验证；指定坐标必须与 K 一一对应，不自动删补用户。
checkGrid(cfg.bs.position_xy_m,'bs.position_xy_m');
assert(numel(cfg.bs.position_xy_m)==2,'bs.position_xy_m must be [x y].');
cfg.bs.position_xy_m = cfg.bs.position_xy_m(:).';
checkScalar(cfg.bs.ground_height_m,-Inf,Inf,'bs.ground_height_m',false);
checkScalar(cfg.bs.height_m,0,Inf,'bs.height_m',false);
checkScalar(cfg.bs.azimuth_deg,-Inf,Inf,'bs.azimuth_deg',false);
checkScalar(cfg.bs.downtilt_deg,-90,90,'bs.downtilt_deg',false);
assert(any(strcmp(cfg.users.mode,{'specified','random','angles'})),'Invalid users.mode.');
checkRange(cfg.users.radius_range_m,0,Inf,'users.radius_range_m');
assert(cfg.users.radius_range_m(2)>0,'Random-user maximum radius must be positive.');
checkRange(cfg.users.azimuth_offsets_deg,-Inf,Inf,'users.azimuth_offsets_deg');
checkScalar(cfg.users.height_m,0,Inf,'users.height_m',false);
checkScalar(cfg.users.angular_distance_m,0,Inf,'users.angular_distance_m',true);
assert(any(strcmp(cfg.users.sampling,{'area','radius'})),'Invalid users.sampling.');
if strcmp(cfg.users.mode,'specified')
    pos = cfg.users.positions_m;
    assert(isnumeric(pos) && isreal(pos) && isequal(size(pos),[cfg.K 3]) && all(isfinite(pos(:))),'users.positions_m must have exactly K rows and 3 columns; update positions when changing K.');
end
checkBool(cfg.pathloss.enabled,'pathloss.enabled');
checkScalar(cfg.pathloss.exponent,0,Inf,'pathloss.exponent',false);
checkScalar(cfg.pathloss.reference_distance_m,0,Inf,'pathloss.reference_distance_m',true);
checkScalar(cfg.pathloss.distance_floor_m,0,Inf,'pathloss.distance_floor_m',true);
if ~isempty(cfg.pathloss.reference_loss_dB)
    checkScalar(cfg.pathloss.reference_loss_dB,-Inf,Inf,'pathloss.reference_loss_dB',false);
end
scene = fimScene('prepare',cfg);
cfg.bs.position_m = scene.bs.position_m;
cfg.bs.R_global_from_array = scene.bs.R_global_from_array;
cfg.bs.boresight_global = scene.bs.boresight_global;
if strcmp(cfg.users.mode,'specified')
    delta = bsxfun(@minus,cfg.users.positions_m,cfg.bs.position_m);
    assert(all(sum(delta.^2,2)>0),'A specified UE cannot coincide with the BS position.');
end


% CDL 返回的角度定义、侧别需由用户确认；这里仅验证已选枚举。
cfg.cdl = validateCDLConfig(cfg.cdl,cfg.fc_Hz);
assert(any(strcmp(cfg.cdl.eoa_convention,{'elevation','zenith'})),'Invalid EoA convention.');
assert(any(strcmp(cfg.cdl.bs_angle_side,{'arrival','departure'})),'Invalid BS angle side.');
assert(any(strcmp(cfg.cdl.center_mode,{'power_mean','explicit'})),'Invalid CDL center_mode.');
if strcmp(cfg.cdl.center_mode,'explicit')
    v = cfg.cdl.reference_az_el_deg;
    assert(isnumeric(v) && isreal(v) && numel(v)==2 && all(isfinite(v(:))) && abs(v(2))<=90,'Set reference [standard az, elevation].');
end
checkScalar(cfg.cdl.frequency_offset_Hz,-Inf,Inf,'cdl.frequency_offset_Hz',false);
assert(iscell(cfg.cdl.extra_args) && mod(numel(cfg.cdl.extra_args),2)==0,'CDL arguments must be name-value pairs.');
for i = 1:2:numel(cfg.cdl.extra_args)
    assert(ischar(cfg.cdl.extra_args{i}),'CDL parameter names must be character vectors.');
    assert(~strcmpi(cfg.cdl.extra_args{i},'CenterFrequency'),'Set CenterFrequency only via cfg.fc_Hz.');
end
R = cfg.pattern.native_to_body;
assert(isnumeric(R) && isreal(R) && isequal(size(R),[3 3]) && all(isfinite(R(:))) && norm(R'*R-eye(3),'fro')<1e-10 && abs(det(R)-1)<1e-10,'native_to_body must be a proper rotation.');

% 调试参数不降低频率分辨率：分组近似必须显式设置 frequency_mode。
if cfg.fast_mode
    cfg.MC = min(cfg.MC,cfg.fast_MC);
    if strcmp(cfg.fim.mode,'deformation') && strcmp(cfg.deform.family,'subarray')
        cfg.subarray.sweeps = min(cfg.subarray.sweeps,cfg.fast_sweeps);
    end
    cfg.pca.sweeps = min(cfg.pca.sweeps,cfg.fast_sweeps);
end
cfg.N = cfg.Nh*cfg.Nv;
cfg.snr_dB_list = cfg.snr_dB_list(:).';
checkBool(cfg.ofdm.enabled,'ofdm.enabled');
checkBool(cfg.ofdm.include_cp_overhead,'ofdm.include_cp_overhead');
checkBool(cfg.ofdm.frequency_dependent_steering,'ofdm.frequency_dependent_steering');
checkScalar(cfg.ofdm.scs_Hz,0,Inf,'ofdm.scs_Hz',true);
checkInteger(cfg.ofdm.nfft,2,Inf,'ofdm.nfft');
assert(any(strcmp(cfg.ofdm.frequency_mode,{'full','grouped'})),'Invalid ofdm.frequency_mode.');
assert(any(strcmp(cfg.ofdm.resource_mode,{'rb','indices'})),'Invalid ofdm.resource_mode.');
assert(any(strcmp(cfg.ofdm.cp_mode,{'ratio','samples'})),'Invalid ofdm.cp_mode.');
cfg.ofdm.sample_rate_Hz = cfg.ofdm.nfft*cfg.ofdm.scs_Hz;
cfg.ofdm.useful_duration_s = 1/cfg.ofdm.scs_Hz;
if strcmp(cfg.ofdm.cp_mode,'ratio')
    checkScalar(cfg.ofdm.cp_ratio,0,1,'ofdm.cp_ratio',false);
    % CP 比例是精确时间定义；样点数只作为外部时域实现的取整建议。
    cfg.ofdm.cp_samples_equivalent = cfg.ofdm.nfft*cfg.ofdm.cp_ratio;
    cfg.ofdm.cp_samples = round(cfg.ofdm.cp_samples_equivalent);
    cfg.ofdm.cp_duration_s = cfg.ofdm.cp_ratio/cfg.ofdm.scs_Hz;
else
    checkInteger(cfg.ofdm.cp_samples,0,cfg.ofdm.nfft,'ofdm.cp_samples');
    cfg.ofdm.cp_samples_equivalent = cfg.ofdm.cp_samples;
    cfg.ofdm.cp_duration_s = cfg.ofdm.cp_samples/cfg.ofdm.sample_rate_Hz;
    cfg.ofdm.cp_ratio = cfg.ofdm.cp_samples/cfg.ofdm.nfft;
end
cfg.ofdm.symbol_duration_s = cfg.ofdm.useful_duration_s+cfg.ofdm.cp_duration_s;

if ~cfg.ofdm.enabled
    % 单频点模式保持所有用户共频，不计 CP，不使用 RB 调度网格。
    assert(cfg.fc_Hz+cfg.cdl.frequency_offset_Hz>0,'Narrowband physical frequency must be positive.');
    cfg.ofdm.freq_offsets_Hz = cfg.cdl.frequency_offset_Hz;
    cfg.ofdm.freq_weights = 1;
    cfg.ofdm.frequency_counts = 1;
    cfg.ofdm.selected_indices = cfg.cdl.frequency_offset_Hz/cfg.ofdm.scs_Hz;
    cfg.ofdm.num_tones = 1;
    cfg.ofdm.num_active = 1;
    cfg.ofdm.active_bandwidth_Hz = cfg.ofdm.scs_Hz;
    cfg.ofdm.frequency_span_Hz = cfg.ofdm.scs_Hz;
    cfg.ofdm.cp_efficiency = 1;
    cfg.ofdm.is_frequency_approximation = false;
    cfg.ofdm.user_active_mask = true(cfg.K,1);
    cfg.ofdm.user_sample_mask = true(cfg.K,1);
    cfg.ofdm.user_num_active = ones(cfg.K,1);
    cfg.ofdm.rb_of_active = 0;
    cfg.ofdm.rb_of_sample = 0;
    cfg.ofdm.user_rb_mask = [];
    return
end

% RB 模式显式定义每个用户的分配；indices 模式保持旧版所有用户共享网格。
if strcmp(cfg.ofdm.resource_mode,'rb')
    checkInteger(cfg.ofdm.num_rb,1,Inf,'ofdm.num_rb');
    checkInteger(cfg.ofdm.subcarriers_per_rb,1,Inf,'ofdm.subcarriers_per_rb');
    A = cfg.ofdm.num_rb*cfg.ofdm.subcarriers_per_rb;
    if isempty(cfg.ofdm.rb_start_index)
        firstIndex = -floor(A/2);
    else
        checkInteger(cfg.ofdm.rb_start_index,-Inf,Inf,'ofdm.rb_start_index');
        firstIndex = cfg.ofdm.rb_start_index;
    end
    idx = firstIndex+(0:A-1);
    rbOfTone = floor((0:A-1)/cfg.ofdm.subcarriers_per_rb)+1;
    rbMask = buildRBAllocation(cfg);
    userMask = rbMask(:,rbOfTone);
else
    idx = cfg.ofdm.active_indices;
    assert(isnumeric(idx) && isreal(idx) && isvector(idx) && ~isempty(idx) && all(isfinite(idx(:))) && all(idx(:)==fix(idx(:))),'active_indices must be nonempty integer FFT indices.');
    idx = sort(double(idx(:).'));
    assert(numel(unique(idx))==numel(idx),'active_indices must not contain duplicates.');
    assert(isempty(cfg.ofdm.rb_allocation),'rb_allocation is only supported when resource_mode is rb.');
    A = numel(idx);
    rbOfTone = zeros(1,A);
    rbMask = [];
    userMask = true(cfg.K,A);
end
assert(all(idx>=-floor(cfg.ofdm.nfft/2) & idx<=ceil(cfg.ofdm.nfft/2)-1),'Active resource indices fall outside signed FFT bins; increase NFFT or revise RB start/count.');
assert(all(cfg.fc_Hz+idx*cfg.ofdm.scs_Hz>0),'Every physical subcarrier frequency must be positive.');
cfg.ofdm.active_indices = idx;
cfg.ofdm.num_active = A;
cfg.ofdm.active_bandwidth_Hz = A*cfg.ofdm.scs_Hz;
cfg.ofdm.frequency_span_Hz = (idx(end)-idx(1)+1)*cfg.ofdm.scs_Hz;
if strcmp(cfg.cdl.model,'rays')
    assert(cfg.ofdm.frequency_span_Hz<=2e9,'3GPP CDL rays support at most a 2 GHz frequency span.');
end
cfg.ofdm.user_active_mask = userMask;
cfg.ofdm.user_num_active = sum(userMask,2);
cfg.ofdm.user_rb_mask = rbMask;
cfg.ofdm.rb_of_active = rbOfTone;
if strcmp(cfg.ofdm.frequency_mode,'full')
    ordinals = 1:A;
    counts = ones(1,A);
else
    checkInteger(cfg.ofdm.num_frequency_samples,1,A,'ofdm.num_frequency_samples');
    F = cfg.ofdm.num_frequency_samples;
    if strcmp(cfg.ofdm.resource_mode,'rb')
        assert(F>=cfg.ofdm.num_rb,'RB grouped sampling requires num_frequency_samples >= num_rb, so each RB has at least one sample.');
        [ordinals,counts] = groupWithinRBs(cfg.ofdm.num_rb,cfg.ofdm.subcarriers_per_rb,F);
    else
        % 兼容旧版：按显式活动频点的排序序号等分，再选每组实际中间频点。
        [ordinals,counts] = middleOrdinals(A,F,0);
    end
end
cfg.ofdm.selected_indices = idx(ordinals);
cfg.ofdm.frequency_counts = counts;
cfg.ofdm.freq_offsets_Hz = cfg.ofdm.selected_indices*cfg.ofdm.scs_Hz;
cfg.ofdm.freq_weights = counts/A;
cfg.ofdm.num_tones = numel(ordinals);
cfg.ofdm.user_sample_mask = userMask(:,ordinals);
cfg.ofdm.rb_of_sample = rbOfTone(ordinals);
cfg.ofdm.is_frequency_approximation = numel(ordinals)<A;
cfg.ofdm.cp_efficiency = 1;
if cfg.ofdm.include_cp_overhead
    cfg.ofdm.cp_efficiency = 1/(1+cfg.ofdm.cp_ratio);
end
assert(sum(counts)==A && abs(sum(cfg.ofdm.freq_weights)-1)<1e-12,'Frequency quadrature weights are inconsistent.');
represented = double(cfg.ofdm.user_sample_mask)*counts(:);
assert(isequal(represented,cfg.ofdm.user_num_active),'Grouped sampling must preserve every user allocation exactly.');
end

function cdl = validateCDLConfig(cdl,fc_Hz)
profile = cdl.profile;
assert((ischar(profile) && isrow(profile)) || (isstring(profile) && isscalar(profile)),...
    'cdl.profile must be CDL-A through CDL-E.');
profile = upper(char(profile));
if startsWith(profile,'NR-')
    profile = profile(4:end);
end
assert(any(strcmp(profile,{'CDL-A','CDL-B','CDL-C','CDL-D','CDL-E'})),...
    'cdl.profile must be CDL-A through CDL-E.');
cdl.profile = profile;
assert(ischar(cdl.model) && any(strcmp(cdl.model,{'cluster','rays'})),...
    'cdl.model must be cluster or rays.');
checkBool(cdl.align_to_ue,'cdl.align_to_ue');
assert(ischar(cdl.ray_coupling) && any(strcmp(cdl.ray_coupling,{'random','aligned'})),...
    'cdl.ray_coupling must be random or aligned.');
checkScalar(cdl.rx_slant_deg,-Inf,Inf,'cdl.rx_slant_deg',false);
checkScalar(cdl.tx_slant_deg,-Inf,Inf,'cdl.tx_slant_deg',false);
checkGrid(cdl.velocity_mps,'cdl.velocity_mps');
assert(numel(cdl.velocity_mps)==3,'cdl.velocity_mps must contain [vx vy vz] in global axes.');
cdl.velocity_mps = double(cdl.velocity_mps(:).');
checkScalar(cdl.time_s,0,Inf,'cdl.time_s',false);
if ~isempty(cdl.k_factor_dB)
    checkScalar(cdl.k_factor_dB,-Inf,Inf,'cdl.k_factor_dB',false);
    assert(strcmp(cdl.model,'rays'),'cdl.k_factor_dB requires the CDL rays model.');
    assert(~isempty(cdl.profile_mat) || any(strcmp(cdl.profile,{'CDL-D','CDL-E'})),...
        'cdl.k_factor_dB requires a CDL-D or CDL-E profile.');
end
scaling = cdl.angle_scaling_deg;
assert(isnumeric(scaling) && isreal(scaling),'cdl.angle_scaling_deg must be real numeric.');
if ~isempty(scaling)
    assert(isequal(size(scaling),[2 4]),...
        'cdl.angle_scaling_deg must be 2 x 4 for [AoD AoA ZoD ZoA].');
    selected = scaling(:,~all(isnan(scaling),1));
    assert(all(isfinite(selected(:))) && all(selected(2,:)>=0),...
        'Each angle-scaling column needs a finite mean and a nonnegative RMS spread, or two NaNs.');
    zenith = scaling(1,3:4);
    assert(all(isnan(zenith) | (zenith>=0 & zenith<=180)),...
        'Mean zenith angles must be in [0, 180] degrees.');
    assert(strcmp(cdl.model,'rays'),'cdl.angle_scaling_deg requires the CDL rays model.');
end
if strcmp(cdl.model,'rays')
    assert(strcmp(cdl.bs_angle_side,'arrival'),'The uplink CDL rays model requires cdl.bs_angle_side=arrival.');
    checkScalar(fc_Hz,0.5e9,100e9,'fc_Hz for the 3GPP CDL rays model',false);
end
end

function mask = buildRBAllocation(cfg)
B = cfg.ofdm.num_rb;
if isempty(cfg.ofdm.rb_allocation)
    mask = true(cfg.K,B);
    return
end
allocation = cfg.ofdm.rb_allocation;
assert(iscell(allocation) && numel(allocation)==cfg.K,'rb_allocation must be empty or have exactly K cells.');
mask = false(cfg.K,B);
for k = 1:cfg.K
    list = allocation{k};
    assert(isnumeric(list) && isreal(list) && isvector(list) && ~isempty(list) && all(isfinite(list(:))) && all(list(:)==fix(list(:))) && all(list(:)>=1 & list(:)<=B),'Each UE must receive at least one valid 1-based RB index.');
    assert(numel(unique(list))==numel(list),'A UE RB list must not contain duplicates.');
    mask(k,list(:)) = true;
end
end

function [ordinals,counts] = groupWithinRBs(B,L,F)
% 每个 RB 至少一组，余下样点均匀分配；最终组数严格等于请求 F。
% 一组永不跨 RB，因此不会把不同用户活动集合混在同一个频率代表点中。
perRB = floor(F/B)*ones(1,B);
remainder = mod(F,B);
% MLang 对空数组广播支持有限；F 为 B 的整数倍时不做空切片算术。
if remainder>0
    perRB(1:remainder) = perRB(1:remainder)+1;
end
ordinals = zeros(1,F);
counts = zeros(1,F);
next = 1;
for b = 1:B
    [localOrdinals,localCounts] = middleOrdinals(L,perRB(b),(b-1)*L);
    target = next:next+perRB(b)-1;
    ordinals(target) = localOrdinals;
    counts(target) = localCounts;
    next = next+perRB(b);
end
end

function [ordinals,counts] = middleOrdinals(A,F,offset)
edges = floor((0:F)*A/F);
counts = diff(edges);
ordinals = floor((edges(1:end-1)+1+edges(2:end))/2)+offset;
end

function checkInteger(x,lo,hi,name)
assert(isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x==fix(x) && x>=lo && x<=hi,['Invalid ' name '.']);
end

function checkScalar(x,lo,hi,name,strictLower)
assert(isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x>=lo && x<=hi,['Invalid ' name '.']);
if strictLower
    assert(x>lo,['Invalid ' name '.']);
end
end

function checkRange(x,lo,hi,name)
assert(isnumeric(x) && isreal(x) && numel(x)==2 && all(isfinite(x(:))) && x(2)>=x(1) && all(x(:)>=lo & x(:)<=hi),['Invalid ' name '.']);
end

function checkBool(x,name)
assert((isnumeric(x) || islogical(x)) && isreal(x) && isscalar(x) && (x==0 || x==1),['Invalid logical flag ' name '.']);
end

function cfg = fillMissingDefaults(cfg,defaults)
% 不覆盖已有字段，仅为旧版 cfg 补充缺失的结构层和叶节点。
assert(isstruct(cfg) && isscalar(cfg),'cfg must be a scalar structure.');
names = fieldnames(defaults);
for i = 1:numel(names)
    name = names{i};
    if ~isfield(cfg,name)
        cfg.(name) = defaults.(name);
    elseif isstruct(defaults.(name))
        assert(isstruct(cfg.(name)) && isscalar(cfg.(name)),['Invalid configuration structure: ' name]);
        cfg.(name) = fillMissingDefaults(cfg.(name),defaults.(name));
    end
end
end

function C = validateExternalCovariance(C,N)
% 固定天线端口协方差；不重复加入热噪声或已建模用户的协方差。
if isempty(C)
    assert(isnumeric(C),'irc.external_covariance must be numeric.');
    C = [];
    return
end
assert(isnumeric(C) && isequal(size(C),[N N]) && all(isfinite(C(:))),...
    'irc.external_covariance must be empty or a finite N x N matrix.');
scale = max(abs(C(:)));
if scale==0
    C = [];
    return
end
% 在归一化尺度验证，避免瓦特量级输入被绝对容差误判。
A = C/scale;
assert(norm(A-A','fro')<=1e-10*norm(A,'fro'),...
    'irc.external_covariance must be Hermitian.');
A = (A+A')/2;
[V,D] = eig(A);
d = real(diag(D));
assert(min(d)>=-1e-10*max(abs(d)),...
    'irc.external_covariance must be positive semidefinite.');
% 只移除容差内的浮点负特征值，不接受真实负功率。
if any(d<0)
    d(d<0) = 0;
    A = V*diag(d)*V';
    A = (A+A')/2;
end
C = scale*A;
end

function checkGrid(v,name)
assert(isnumeric(v) && isreal(v) && isvector(v) && ~isempty(v) && all(isfinite(v(:))),['Invalid ' name '.']);
end

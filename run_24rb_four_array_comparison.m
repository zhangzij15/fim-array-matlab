%RUN_24RB_FOUR_ARRAY_COMPARISON 24 RB、四种阵列形态的轻量配对仿真。
% 比较：固定阵列 / PCA / 逐阵元方向图重构 / 刚性子阵位姿重构(RSA/FIM)。
% 本入口按128阵元(16×8)配置；每个子阵为4×4阵元，共4×2=8个子阵。
% 轻量化只减少MC、频率代表点和搜索网格，不改变24 RB物理带宽及每UE总功率。

%% 常用选项：通常只修改这一段
UE_CASE = [1 10];
MC_PER_CASE = 20;              % 首次运行建议设为1；正式统计可增大。
SAMPLES_PER_RB = 1;             % 1=24代表频点；2=48；12=完整288频点。
RANDOM_SEED = 20260928;
SHOW_FIGURES = false;           % false：隐藏图窗但仍保存PNG。
RESUME_DIRECTORY = '';          % 空=新任务；续跑时填原任务目录。
PATTERN_TXT = '';                % 空值使用通用配置的合成方向图，可填写真实 TXT 路径。
CDL_PROFILE_MAT = '';            % 填入已有 rawCDL MAT 可跳过生成器。
CDL_PROFILE = 'CDL-C';           % 默认生成器支持 CDL-A、CDL-B、CDL-C、CDL-D、CDL-E。
CDL_MODEL = 'cluster';           % 'rays' 用标准子径模式；'cluster' 用于复现已有研究。

assert(isnumeric(UE_CASE) && isvector(UE_CASE) && ~isempty(UE_CASE) && ...
    all(ismember(UE_CASE,[1 10])) && numel(unique(UE_CASE))==numel(UE_CASE), ...
    'UE_CASE must be 1, 10, or [1 10].');
assert(isnumeric(MC_PER_CASE) && isscalar(MC_PER_CASE) && isfinite(MC_PER_CASE) && ...
    MC_PER_CASE>=1 && MC_PER_CASE==fix(MC_PER_CASE), ...
    'MC_PER_CASE must be a positive integer.');
assert(isnumeric(SAMPLES_PER_RB) && isscalar(SAMPLES_PER_RB) && ...
    any(SAMPLES_PER_RB==[1 2 3 4 6 12]), ...
    'Use 1, 2, 3, 4, 6 or 12 samples per RB.');
assert(exist('runFIMStudy','file')==2 && exist('configFIMSimulation','file')==2, ...
    'Run this script from the compact four-array project directory.');

%% 载入基础配置并覆盖本次实验条件
cfg = configFIMSimulation();
cfg.cdl.profile = CDL_PROFILE;
cfg.cdl.model = CDL_MODEL;
if ~isempty(PATTERN_TXT)
    cfg.pattern.file = PATTERN_TXT;
end
if ~isempty(CDL_PROFILE_MAT)
    cfg.cdl.profile_mat = CDL_PROFILE_MAT;
end

% 阵列：16个水平阵元×8个竖直阵元=128阵元。
cfg.Nh = 16;
cfg.Nv = 8;
cfg.fc_Hz = 6.7e9;
cfg.c0 = 299792458;
cfg.dh_lambda = 0.5;
cfg.dv_lambda = 1.5;

% BS/UE部署。
cfg.bs.position_xy_m = [0 0];
cfg.bs.ground_height_m = 0;
cfg.bs.height_m = 20;
cfg.bs.azimuth_deg = 0;
% 机械下倾角沿用通用配置，为四种方案共用的安装角。
cfg.users.mode = 'random';
cfg.users.sampling = 'area';
cfg.users.radius_range_m = [20 200];
cfg.users.azimuth_offsets_deg = [-60 60];
cfg.users.height_m = 1.5;

%% OFDM：24 RB×12子载波=288活动子载波，8.64 MHz
cfg.ofdm.enabled = true;
cfg.ofdm.resource_mode = 'rb';
cfg.ofdm.num_rb = 24;
cfg.ofdm.subcarriers_per_rb = 12;
cfg.ofdm.scs_Hz = 30e3;
cfg.ofdm.nfft = 512;             % 512点FFT容纳288个活动子载波；不能写成51。
cfg.ofdm.rb_start_index = [];
cfg.ofdm.rb_allocation = [];
cfg.ofdm.num_frequency_samples = 24*SAMPLES_PER_RB;
if SAMPLES_PER_RB==12
    cfg.ofdm.frequency_mode = 'full';
else
    cfg.ofdm.frequency_mode = 'grouped';
end
cfg.ofdm.frequency_dependent_steering = true;
cfg.ofdm.cp_mode = 'ratio';
cfg.ofdm.cp_ratio = 0.12;
cfg.ofdm.include_cp_overhead = true;

%% 功率、噪声、路径损耗和接收机
Ptx_dBm = 23;
cfg.noiseVar = 1.380649e-23*290*10^(7/10)*cfg.ofdm.scs_Hz;
cfg.power.snr_reference = 'total_over_tone_noise';
cfg.snr_dB = 10*log10(10^((Ptx_dBm-30)/10)/cfg.noiseVar);
cfg.snr_dB_list = cfg.snr_dB;
cfg.pathloss.enabled = true;
cfg.pathloss.exponent = 2.2;
cfg.pathloss.reference_distance_m = 1;
cfg.pathloss.distance_floor_m = 1;
cfg.pathloss.reference_loss_dB = 20*log10(4*pi*cfg.fc_Hz/cfg.c0)+35;
cfg.receiver = 'mmse-irc';
cfg.irc.external_covariance = [];
cfg.throughput_mode = 'shannon';
cfg.opt_mode = 'sum_rate';

%% 粗搜索设置
% 4×2个子阵；每个子阵4×4阵元。法向位移固定为0，只搜索yaw/tilt。
cfg.subarray.m = 4;
cfg.subarray.n = 2;
cfg.subarray.displacement_grid_lambda = 0;
cfg.subarray.yaw_grid_deg = [-30 0 30];
cfg.subarray.tilt_grid_deg = [-20 0 20];
cfg.subarray.min_spacing_lambda = 0;
cfg.subarray.sweeps = 1;

% PCA：每个阵元独立在原阵面内移动，阵元朝向固定；不施加最小间距约束。
cfg.pca.region_mode = 'per_element';
cfg.pca.x_grid_lambda = [-0.25 0 0.25];
cfg.pca.z_grid_lambda = [-1.5 0 1.5];
cfg.pca.x_bounds_lambda = [];
cfg.pca.z_bounds_lambda = [];
cfg.pca.min_spacing_lambda = 0;
cfg.pca.sweeps = 1;

% 逐阵元方向图重构：阵元中心不动，只选择方向图朝向。
cfg.orientation.yaw_grid_deg = [-30 0 30];
cfg.orientation.tilt_grid_deg = [-20 0 20];
cfg.orientation.sweeps = 1;

cfg.fim.mode = 'deformation';
cfg.deform.family = 'subarray';
cfg.comparison.include_pca = true;
cfg.deform.cache_limit = 4;
cfg.fast_mode = false;
cfg.seed = RANDOM_SEED;
cfg.progress_every = 1;

%% 批量MC、检查点与输出
cfg.study.K_values = UE_CASE;
cfg.study.MC_per_case = MC_PER_CASE;
cfg.study.batch_MC = 1;
cfg.study.resume_directory = RESUME_DIRECTORY;
cfg.study.make_plots = true;
cfg.output_dir = fullfile(cfg.output_dir,'four_array_24rb_lowcost');
cfg.save_results = true;
cfg.save_figures = true;
cfg.output.save_per_tone = false;
cfg.output.save_realizations = false;
cfg.output.save_geometry = false;

%% 图形输出
cfg.plot.enabled = true;
cfg.plot.show_parameter_table = true;
cfg.plot.show_geometry = true;
cfg.plot.geometry_mc_index = 1;
cfg.plot.show_scene = true;
cfg.plot.scene_indices = 1:min(3,MC_PER_CASE);
cfg.plot.show_per_tone_cdf = false;
cfg.plot.show_frequency_profile = false;
cfg.plot.frequency_snr_index = 1;
cfg.plot.export_dpi = 150;
if SHOW_FIGURES
    cfg.plot.figure_visible = 'on';
else
    cfg.plot.figure_visible = 'off';
end
cfg.plot.close_after_save = ~SHOW_FIGURES;

%% 正式输入预检查：避免跑到信道初始化阶段才发现文件路径错误
assert(exist(cfg.pattern.file,'file')==2, ...
    'Antenna pattern TXT not found. Please edit PATTERN_TXT or cfg.pattern.file: %s',cfg.pattern.file);
if ~isempty(CDL_PROFILE_MAT)
    assert(exist(CDL_PROFILE_MAT,'file')==2, ...
        'CDL_PROFILE_MAT not found: %s',CDL_PROFILE_MAT);
else
    fprintf('CDL profile %s (%s mode) will be generated by %s.\n', ...
        cfg.cdl.profile,cfg.cdl.model,cfg.cdl.generator);
end

fprintf('Four arrays: K=%s; MC/case=%d; N=%d (%d×%d); computed tones=%d/288.\n', ...
    mat2str(UE_CASE),MC_PER_CASE,cfg.Nh*cfg.Nv,cfg.Nh,cfg.Nv,cfg.ofdm.num_frequency_samples);
fprintf('Physical bandwidth=8.64 MHz; each UE total power=23 dBm; PCA min spacing=0.\n');
[study,figs,files] = runFIMStudy(cfg);

% 数值结果先由runFIMStudy保存；随后检查主要PNG是否成功导出。
requiredGraphics = {'THP_CDF.png','SNR_CDF.png'};
missingGraphics = {};
for imageIndex = 1:numel(requiredGraphics)
    if exist(fullfile(study.directory,requiredGraphics{imageIndex}),'file')~=2
        missingGraphics{end+1} = requiredGraphics{imageIndex};
    end
end
for caseIndex = 1:numel(study.K_values)
    K = study.K_values(caseIndex);
    if isempty(dir(fullfile(study.directory,sprintf('K%d_parameters_*.png',K))))
        missingGraphics{end+1} = sprintf('K%d parameter-table PNG',K);
    end
    for sceneIndex = cfg.plot.scene_indices
        sceneName = sprintf('K%d_BS_UE_scene_MC%04d.png',K,sceneIndex);
        if exist(fullfile(study.directory,sceneName),'file')~=2
            missingGraphics{end+1} = sceneName;
        end
    end
end
if ~isempty(missingGraphics)
    error('FIM:RequiredGraphicsMissing', ...
        'Numeric results saved in %s, but required images are missing: %s. Replot with plotFIMResults.', ...
        study.directory,strjoin(missingGraphics,', '));
end
fprintf('\nRequired CDF, parameter-table and BS/UE PNGs saved in:\n%s\n',study.directory);
% 单独重绘：plotFIMResults(fullfile(study.directory,'study.mat'));
% 四方案使用同一MC路径、功率和接收机；不预设任何性能增益排序。

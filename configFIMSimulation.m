function cfg = configFIMSimulation()
%CONFIGFIMSIMULATION 四方对比唯一配置：固定 / PCA / 方向图重构 / 子阵位姿重构。
% 修改本文件后运行 run_fim_comparison。
root = fileparts(mfilename('fullpath'));

% 阵列与部署；阵面本地 x 水平、y 法向、z 竖直。
cfg.Nh = 16; cfg.Nv = 8;
cfg.dh_lambda = 0.5; cfg.dv_lambda = 3;
cfg.fc_Hz = 6.7e9; cfg.c0 = 299792458;
cfg.K = 1; % 直接调用 simulateFIM 时使用；批量研究由 study.K_values 覆盖。
cfg.bs.position_xy_m = [0 0]; cfg.bs.ground_height_m = 0;
cfg.bs.height_m = 20; cfg.bs.azimuth_deg = 0;
% 24° 由 288 RB 的固定阵列校准选出，作为四种方案共用的安装角。
cfg.bs.downtilt_deg = 24;
cfg.users.mode = 'random';
cfg.users.positions_m = []; % specified 模式需提供 K x 3 全局绝对坐标（m）。
cfg.users.radius_range_m = [10 150];
cfg.users.azimuth_offsets_deg = [-60 60];
cfg.users.height_m = 1.5; cfg.users.sampling = 'area';
cfg.users.angular_distance_m = 100;
cfg.user_az_range = [-60 60]; cfg.user_el_range = [-10 10];

% OFDM：24 RB、288 个实际子载波、8.64 MHz；full 为全部频点。
cfg.ofdm.enabled = true;
cfg.ofdm.resource_mode = 'rb'; cfg.ofdm.num_rb = 24;
cfg.ofdm.subcarriers_per_rb = 12; cfg.ofdm.scs_Hz = 30e3;
cfg.ofdm.nfft = 512; cfg.ofdm.rb_start_index = [];
cfg.ofdm.rb_allocation = []; % 所有 UE 共享全部 RB；指定用户分配可用 K 个 cell。
cfg.ofdm.active_indices = -144:143;
cfg.ofdm.frequency_mode = 'full'; cfg.ofdm.num_frequency_samples = 24;
cfg.ofdm.frequency_dependent_steering = true;
cfg.ofdm.cp_mode = 'ratio'; cfg.ofdm.cp_ratio = 0.12;
cfg.ofdm.cp_samples = 61; cfg.ofdm.include_cp_overhead = true;

% 每 UE 总功率 23 dBm；改变 SCS 时也应重算单子载波噪声。
Ptx_dBm = 23;
cfg.noiseVar = 1.380649e-23*290*10^(7/10)*cfg.ofdm.scs_Hz;
cfg.power.snr_reference = 'total_over_tone_noise';
cfg.snr_dB = 10*log10(10^((Ptx_dBm-30)/10)/cfg.noiseVar);
cfg.snr_dB_list = cfg.snr_dB+(-15:5:15);
cfg.pathloss.enabled = true; cfg.pathloss.exponent = 2.2;
cfg.pathloss.reference_distance_m = 1; cfg.pathloss.distance_floor_m = 1;
cfg.pathloss.reference_loss_dB = 20*log10(4*pi*cfg.fc_Hz/cfg.c0)+35;
cfg.receiver = 'mmse-irc'; % 可改为 zf；四种阵列始终使用相同接收机。
cfg.irc.external_covariance = []; % 接收端口固定的外部干扰协方差（可选）。
cfg.throughput_mode = 'shannon'; cfg.opt_mode = 'sum_rate';

% 子阵：4 x 2 块，每块沿原法向移动并绕自身中心转动。
cfg.subarray.m = 4; cfg.subarray.n = 2;
cfg.subarray.displacement_grid_lambda = -0.5:0.25:0.5;
cfg.subarray.yaw_grid_deg = -30:10:30;
cfg.subarray.tilt_grid_deg = -20:10:20;
cfg.subarray.min_spacing_lambda = 0.4; cfg.subarray.sweeps = 1;
cfg.orientation.yaw_grid_deg = -30:10:30;
cfg.orientation.tilt_grid_deg = -20:10:20;
cfg.orientation.sweeps = 1;
% 方向图重构：固定阵元位置，逐阵元独立选择上述 yaw/tilt。
% PCA：逐阵元在原阵面 x-z 内移动，朝向固定。
cfg.pca.region_mode = 'per_element';
cfg.pca.x_grid_lambda = [-0.25 0 0.25];
cfg.pca.z_grid_lambda = -1.5:0.5:1.5;
cfg.pca.x_bounds_lambda = []; cfg.pca.z_bounds_lambda = [];
cfg.pca.min_spacing_lambda = 0.4; cfg.pca.sweeps = 1;

% 当前显式使用合成方向图；可指定真实 TXT 覆盖，路径错误时直接报错。
cfg.pattern.file = fullfile(root,'patterns','Synthetic_6p7GHz_H120_V24.txt');
cfg.pattern.native_to_body = eye(3); cfg.pattern.alignment = 'peak';
cfg.pattern.expected_az_hpbw_deg = 120; cfg.pattern.expected_el_hpbw_deg = 24;
cfg.pattern.use_phase = false; cfg.front_only = true;
cfg.cdl.profile_mat = ''; % 填含 rawCDL 或 result.rawCDL 的 MAT 可跳过生成器。
cfg.cdl.generator = 'generateCDLBase'; % 基础 MATLAB 读取随工程提供的标准 CDL 簇参数。
cfg.cdl.extra_args = {'DelaySpread',100};
cfg.cdl.eoa_convention = 'elevation'; cfg.cdl.bs_angle_side = 'arrival';
cfg.cdl.center_mode = 'power_mean'; cfg.cdl.reference_az_el_deg = [];
cfg.cdl.frequency_offset_Hz = 0;
cfg.cdl.profile = 'CDL-C'; % 生成模板时可选 CDL-A 至 CDL-E。
cfg.cdl.model = 'cluster'; % 'rays' 使用 3GPP 子径；'cluster' 用于复现簇级结果。
cfg.cdl.align_to_ue = true; % false 按模板的固定全局角度生成。
cfg.cdl.ray_coupling = 'random'; % 'aligned' 用于固定子径配对的核对。
cfg.cdl.rx_slant_deg = 0; cfg.cdl.tx_slant_deg = 0;
cfg.cdl.velocity_mps = [0 0 0]; % 所有 UE 共用的全局 [vx vy vz]。
cfg.cdl.time_s = 0; % 相对于随机初相的快照时刻。
cfg.cdl.k_factor_dB = []; % 仅 CDL-D/E 子径模式，全模型 K 因子。
cfg.cdl.angle_scaling_deg = []; % [AoD AoA ZoD ZoA] 的目标均值与 RMS 角扩展。

% 默认做 K=1 / K=10 配对研究；每组20 MC、每批1 MC，支持断点续算。
cfg.study.K_values = [1 10]; cfg.study.MC_per_case = 20;
cfg.study.batch_MC = 1; cfg.study.resume_directory = '';
cfg.study.make_plots = true;
cfg.MC = 3; cfg.seed = 20260928; cfg.progress_every = 1;
cfg.fast_mode = false; cfg.fast_MC = 2; cfg.fast_sweeps = 1;

% 图形与输出。scene_indices 指向真实保存的 MC 样本，不重新抽取场景。
cfg.output_dir = fullfile(root,'results');
cfg.save_results = true; cfg.save_figures = true;
cfg.output.save_per_tone = false;
cfg.output.save_realizations = false; cfg.output.save_geometry = false;
cfg.output.mat_format = 'default';
cfg.plot.enabled = true; cfg.plot.show_geometry = true;
cfg.plot.show_scene = true; cfg.plot.scene_indices = [1 2 3];
cfg.plot.show_parameter_table = true;
cfg.plot.show_per_tone_cdf = false; cfg.plot.show_frequency_profile = false;
cfg.plot.frequency_snr_index = 1; cfg.plot.geometry_mc_index = 1;
cfg.plot.font_name = 'Microsoft YaHei'; cfg.plot.figure_visible = 'on';
cfg.plot.export_dpi = 180; cfg.plot.close_after_save = false;

% 内部标签与数值选项；主流程只比较这四种阵列。
cfg.fim.mode = 'deformation'; cfg.deform.family = 'subarray';
cfg.comparison.include_pca = true; cfg.deform.cache_limit = 4;
cfg.objective_tol = 1e-10;
cfg.opt.score_mode = 'gram_verified'; cfg.opt.gram_rcond_min = 1e-8;
cfg.check.norm_tol = 1e-10; cfg.check.objective_tol = 1e-8;
cfg.check.zf_rcond_switch = 1e-8; cfg.check.svd_rel_tol = 0;
end

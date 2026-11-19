function [study,figs,files] = run_24rb_pca_no_spacing(cfg,resumeDirectory,scale)
%RUN_24RB_PCA_NO_SPACING PCA无间距约束的四方案入口，默认轻量计算。
% 直接运行：run_24rb_pca_no_spacing                 % light：每组最多5 MC
% 先跑通：  run_24rb_pca_no_spacing([],[],'smoke')  % 每组最多1 MC
% 稍大规模：run_24rb_pca_no_spacing([],[],'balanced') % 每组最多20 MC
% 完整配置：run_24rb_pca_no_spacing([],[],'full')   % 沿用本地MC/频点/网格
% 指定单个K：cfg=configFIMSimulation(); cfg.study.K_values=10;
%           run_24rb_pca_no_spacing(cfg);
% 续算需保持scale及本地配置一致；第二参数显式给出本任务目录。
% light/smoke使用每RB一个代表频点、粗搜索；balanced每RB两个代表频点。
% 各模式均保持阵元数/初始孔径、实际RB/带宽、每UE功率及正式输入接口。
% PCA最小间距强制为0，允许理想相位中心重合；仍保留位移边界和CP检查。
if nargin<1 || isempty(cfg)
    cfg=configFIMSimulation();
end
if nargin<2 || isempty(resumeDirectory)
    resumeDirectory='';
end
if nargin<3 || isempty(scale)
    scale='light';
end
assert(isstruct(cfg) && isscalar(cfg) && isfield(cfg,'pca') && ...
    isfield(cfg,'study') && isfield(cfg,'output_dir'),...
    'Use configFIMSimulation from the compact four-array project.');
if isa(resumeDirectory,'string') && isscalar(resumeDirectory)
    resumeDirectory=char(resumeDirectory);
end
assert(ischar(resumeDirectory),'resumeDirectory must be a character-vector path.');
if isa(scale,'string') && isscalar(scale), scale=char(scale); end
assert(ischar(scale),'scale must be light, smoke, balanced or full.');
scale=lower(strtrim(scale));
switch scale
    case 'light'
        maxMC=5; samplesPerRB=1;
    case 'smoke'
        maxMC=1; samplesPerRB=1;
    case 'balanced'
        maxMC=20; samplesPerRB=2;
    case 'full'
        maxMC=[]; samplesPerRB=[];
    otherwise
        error('FIM:Scale','scale must be light, smoke, balanced or full.');
end

cfg.pca.min_spacing_lambda=0;
% 显式的第二参数才启用续算，防止继承cfg中指向旧0.4lambda任务的续跑目录。
cfg.study.resume_directory=resumeDirectory;
if ~strcmp(scale,'full')
    assert(cfg.ofdm.enabled && strcmp(cfg.ofdm.resource_mode,'rb'),...
        'Reduced scales require enabled OFDM with resource_mode=rb.');
    cfg.fast_mode=false; % MC上限由本入口明确设置，不使用隐式fast模式。
    cfg.study.MC_per_case=min(cfg.study.MC_per_case,maxMC);
    % 保留原batch_MC，使同一模板和种子下的小MC组沿用原批次序列前缀。
    cfg.ofdm.frequency_mode='grouped';
    cfg.ofdm.num_frequency_samples=cfg.ofdm.num_rb*min(samplesPerRB,cfg.ofdm.subcarriers_per_rb);
    % 只保留每个原网格的两端和零点，保留原最大位移/转角范围。
    cfg.pca.x_grid_lambda=coarseGrid(cfg.pca.x_grid_lambda);
    cfg.pca.z_grid_lambda=coarseGrid(cfg.pca.z_grid_lambda);
    cfg.orientation.yaw_grid_deg=coarseGrid(cfg.orientation.yaw_grid_deg);
    cfg.orientation.tilt_grid_deg=coarseGrid(cfg.orientation.tilt_grid_deg);
    cfg.subarray.displacement_grid_lambda=coarseGrid(cfg.subarray.displacement_grid_lambda);
    cfg.subarray.yaw_grid_deg=coarseGrid(cfg.subarray.yaw_grid_deg);
    cfg.subarray.tilt_grid_deg=coarseGrid(cfg.subarray.tilt_grid_deg);
    cfg.pca.sweeps=min(cfg.pca.sweeps,1);
    cfg.orientation.sweeps=min(cfg.orientation.sweeps,1);
    cfg.subarray.sweeps=min(cfg.subarray.sweeps,1);
    outputTag=['pca_no_spacing_' scale];
else
    outputTag='pca_no_spacing'; % 与上一版完整无间距入口的目录和续算计划兼容。
end
cfg.output_dir=fullfile(cfg.output_dir,outputTag);
fprintf('PCA rerun: minimum element spacing DISABLED (0 lambda).\n');
fprintf('Computation scale: %s; configured MC per case: %s.\n',scale,mat2str(cfg.study.MC_per_case));
if ~strcmp(scale,'full')
    fprintf('Frequency approximation: %d representative tones over %d physical active tones.\n',...
        cfg.ofdm.num_frequency_samples,cfg.ofdm.num_rb*cfg.ofdm.subcarriers_per_rb);
    fprintf('Coarse search grids; original element count, aperture, bandwidth and per-UE power retained.\n');
end
if isempty(resumeDirectory)
    fprintf('Starting a new study under: %s\n',cfg.output_dir);
else
    fprintf('Resuming the explicitly selected study; saved configuration must match.\n');
end
[study,figs,files]=runFIMStudy(cfg);
end

function grid=coarseGrid(input)
assert(isnumeric(input) && isreal(input) && ~isempty(input) && ...
    isvector(input) && all(isfinite(input(:))),'Candidate grids must contain finite real values.');
grid=unique([min(input(:)) 0 max(input(:))]);
end

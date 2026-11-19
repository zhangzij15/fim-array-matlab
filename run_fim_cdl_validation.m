function report = run_fim_cdl_validation()
%RUN_FIM_CDL_VALIDATION Analytic and four-array checks for terrestrial ray CDL.
savedRng = rng;
cleanup = onCleanup(@() rng(savedRng));
rng(20261120,'twister');

profiles = {'CDL-A','CDL-B','CDL-C','CDL-D','CDL-E'};
expectedRays = [460 460 480 261 281];
expectedXPR_dB = [10 8 7 11 8];
expectedMaxNormalizedDelay = [9.6586 4.7834 8.6523 12.525 20.6419];
for i = 1:numel(profiles)
    raw = generateCDLBase(profiles{i},'DelaySpread',100);
    rays = expandCDLRays(raw,struct('ray_coupling','aligned'));
    assert(numel(rays.power)==expectedRays(i) && nnz(rays.isLOS)==(i>=4));
    assert(abs(sum(rays.power)-1)<1e-14);
    assert(max(abs(accumarray(rays.cluster,rays.power)-raw.pow(:)))<1e-14);
    assert(abs(max(rays.tau)-expectedMaxNormalizedDelay(i)*100e-9)<1e-19);
    assert(all(abs(rays.xpr-10^(expectedXPR_dB(i)/10))<1e-14));
end

raw = generateCDLBase('CDL-E','DelaySpread',60);
options = struct('k_factor_dB',3,'angle_scaling_deg',[0 NaN NaN 90;0 NaN NaN 0]);
rays = expandCDLRays(raw,options);
assert(abs(10*log10(sum(rays.power(rays.isLOS))/sum(rays.power(~rays.isLOS)))-3)<1e-12);
delaySpread = sqrt(sum(rays.power.*(rays.tau-sum(rays.power.*rays.tau)).^2));
assert(abs(delaySpread-60e-9)<1e-19);
assert(max(abs(rays.Udeparture(2,:)))<1e-14 && all(rays.Udeparture(1,:)>=0));
assert(max(abs(rays.Uarrival(3,:)))<1e-14);

cfg = integrationConfig();
checkAnalyticLOS(cfg);
for i = 1:numel(profiles)
    cfg.cdl.profile = profiles{i};
    for k = [1 10]
        cfg.K = k;
        [r,~,~] = simulateFIM(cfg,'cdf');
        assert(strcmp(r.rawCDL.name,['NR-' profiles{i}]));
        paths = r.realizations{1}.paths;
        assert(numel(paths)==k && all(cellfun(@(p) numel(p.tau)==expectedRays(i),paths)));
        for scheme = {'fixed','pca','orientation','fim'}
            stats = r.(scheme{1});
            assert(all(isfinite(stats.sinr(:))) && all(stats.sinr(:)>=0));
            assert(all(isfinite(stats.sumRate(:))) && ...
                all(stats.sumRate(:)>=r.fixed.sumRate(:)-cfg.check.objective_tol));
        end
    end
    fprintf('PASS: %s, %d rays, four arrays at K=1 and K=10.\n',profiles{i},expectedRays(i));
end
report = struct('passed',true,'profiles',{profiles},'ray_counts',expectedRays, ...
    'K_values',[1 10],'num_four_array_cases',2*numel(profiles));
fprintf('PASS: CDL ray counts, powers, delay scaling, LOS, polarization, Doppler and four-array integration.\n');
fprintf('PASS: CDL time series preserves its time axis, Doppler, RNG state and CP checks.\n');
end

function cfg = integrationConfig()
cfg = configFIMSimulation();
cfg.cdl.model = 'rays';
cfg.K = 1;
cfg.MC = 1;
cfg.Nh = 4;
cfg.Nv = 4;
cfg.dv_lambda = 0.75;
cfg.ofdm.num_rb = 1;
cfg.ofdm.frequency_mode = 'grouped';
cfg.ofdm.num_frequency_samples = 3;
cfg.subarray.m = 2;
cfg.subarray.n = 2;
cfg.subarray.yaw_grid_deg = [-10 0 10];
cfg.subarray.displacement_grid_lambda = [0 0.15];
cfg.subarray.tilt_grid_deg = 0;
cfg.orientation.yaw_grid_deg = [-10 0 10];
cfg.orientation.tilt_grid_deg = 0;
cfg.pca.x_grid_lambda = [-0.1 0 0.1];
cfg.pca.z_grid_lambda = 0;
cfg.save_results = false;
cfg.save_figures = false;
cfg.plot.enabled = false;
cfg.output.save_realizations = true;
end

function checkAnalyticLOS(cfg)
cfg.bs.downtilt_deg = 0;
cfg.users.mode = 'specified';
cfg.users.positions_m = [0 100 cfg.bs.height_m];
cfg.cdl.profile = 'CDL-D';
cfg.cdl.velocity_mps = [0 -2 0];
cfg.pathloss.enabled = false;
cfg = prepareFIMConfig(cfg);
raw = generateCDLBase('CDL-D','DelaySpread',100);
raw.pow(:) = 0;
raw.pow(1) = 1;
state = fimChannel('realization',cfg,struct('ray_template',raw));
model = fimGeometry('prepare',cfg);
[geom,valid] = fimGeometry('build',model,model.zero,cfg);
assert(valid);
geom.positions_m(:,2) = linspace(0,0.15,cfg.N);
pattern = struct('phi',0:360,'theta',0:180,'G',ones(181,361), ...
    'native_to_body',eye(3),'front_only',false);
t = 0.0125;
H0 = fimChannel('geometry',state,pattern,geom,cfg);
H1 = fimChannel('geometry',fimChannel('time',state,t),pattern,geom,cfg);
excessPosition = geom.positions_m(:,2);
nu = reshape(cfg.ofdm.freq_offsets_Hz,1,1,[]);
expected = exp(-1i*2*pi*cfg.fc_Hz*100/cfg.c0).* ...
    exp(1i*2*pi*(cfg.fc_Hz+nu).*excessPosition/cfg.c0);
assert(max(abs(H0(:)-expected(:)))<1e-10,'LOS distance or array phase mismatch.');
doppler = 2*cfg.fc_Hz/cfg.c0;
assert(max(abs(H1(:)-H0(:)*exp(1i*2*pi*doppler*t)))<1e-10,'LOS Doppler phase mismatch.');
times = [t;0;2*t;t];
checkDopplerSeries(state,pattern,geom,cfg,times,H0,doppler);

stationaryCfg = cfg;
stationaryCfg.cdl.velocity_mps = [0 0 0];
stationary = fimChannel('realization',stationaryCfg,struct('ray_template',raw));
stationaryTrace = fimChannel('series',stationary,pattern,geom,stationaryCfg,times);
stationaryError = stationaryTrace-H0;
assert(max(abs(stationaryError(:)))<1e-12, ...
    'A stationary CDL realization must produce a constant trace.');

crossPolarized = state;
crossPolarized.paths{1}.rx_slant_deg = 90;
Hcross = fimChannel('geometry',crossPolarized,pattern,geom,cfg);
assert(max(abs(Hcross(:)))<1e-14,'Orthogonal polarizations must reject the single specular ray.');

diffuseCfg = cfg;
diffuseCfg.cdl.profile = 'CDL-C';
diffuseCfg.K = 2;
diffuseCfg.users.positions_m = [cfg.users.positions_m; 10 80 cfg.bs.height_m];
diffuseCfg = prepareFIMConfig(diffuseCfg);
diffuse = fimChannel('realization',diffuseCfg,struct('ray_template',generateCDLBase('CDL-C')));
checkDopplerSeries(diffuse,pattern,geom,diffuseCfg,times,[],[]);
end

function checkDopplerSeries(state,pattern,geom,cfg,times,H0,doppler)
% Starting with an advanced snapshot must not offset absolute series times.
state = fimChannel('time',state,3*max(times));
previous = rng;
[trace,info] = fimChannel('series',state,pattern,geom,cfg,times);
assert(isequal(rng,previous),'A time series must reuse the existing realization without drawing randomness.');
assert(info.cp_feasible && isequal(info.time_s,times.') && ...
    isequal(size(trace),[cfg.N cfg.K cfg.ofdm.num_tones numel(times)]));
repeatedError = trace(:,:,:,1)-trace(:,:,:,end);
assert(max(abs(repeatedError(:)))<1e-12,'Repeated times must reproduce the same channel.');
for i = 1:numel(times)
    actual = trace(:,:,:,i);
    snapshot = fimChannel('time',state,times(i));
    expected = fimChannel('geometry',snapshot,pattern,geom,cfg);
    assert(max(abs(actual(:)-expected(:)))<1e-10,'Series and geometry snapshots disagree.');
    if ~isempty(H0)
        assert(max(abs(actual(:)-H0(:)*exp(1i*2*pi*doppler*times(i))))<1e-10, ...
            'LOS time series has an incorrect Doppler phase.');
    end
end
invalidGeom = geom;
invalidGeom.positions_m(end,2) = cfg.c0*(cfg.ofdm.cp_duration_s + 1e-3);
[invalid,invalidInfo] = fimChannel('series',state,pattern,invalidGeom,cfg,times);
assert(isempty(invalid) && ~invalidInfo.cp_feasible && isequal(invalidInfo.time_s,times.'), ...
    'A time series must report infeasible CP and preserve its time axis.');
end

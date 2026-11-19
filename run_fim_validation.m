function report = run_fim_validation()
%RUN_FIM_VALIDATION 无实测输入的统一数值自检；所有临时输入自动清理。
report.subarray = sa_test_run_subarray_sanity_checks();
report.receiver = rx_test_run_receiver_sanity_checks();
report.pca = pc_test_run_pca_sanity_checks();
report.fourway = fourWayValidation();
report.passed = report.subarray.passed && report.receiver.passed && report.pca.passed && report.fourway.passed;
report.num_groups = numel(report.subarray.checks)+numel(report.receiver.checks)+numel(report.pca.checks)+numel(report.fourway.checks);
fprintf('PASS: %d numerical validation groups. Synthetic fixtures only.\n',report.num_groups);
end

function report=fourWayValidation()
testDir=tempname; mkdir(testDir); cleanup=onCleanup(@() removeFourWayDirectory(testDir));
patternPath=fullfile(testDir,'TEST_four_array_pattern.txt'); writeTestPattern(patternPath);
rawCDL=struct('name','TEST_ONLY_four_array_three_path_fixture','pow',[0.55 0.3 0.15],...
    'taus',[0 40 120]*1e-9,'AoA',[70 95 110]*pi/180,'EoA',[-3 0 5]*pi/180,...
    'numClusters',3,'test_fixture_only',true);
rawCDL.AoD=rawCDL.AoA; rawCDL.EoD=rawCDL.EoA;
profilePath=fullfile(testDir,'TEST_rawCDL.mat'); save(profilePath,'rawCDL');
cfg=configFIMSimulation(); cfg.Nh=4; cfg.Nv=4; cfg.K=2; cfg.MC=2;
cfg.dh_lambda=0.5; cfg.dv_lambda=0.5; cfg.seed=20261102;
cfg.users.mode='random'; cfg.pathloss.enabled=false; cfg.noiseVar=1;
cfg.power.snr_reference='per_subcarrier'; cfg.snr_dB=0;
cfg.ofdm.num_rb=2; cfg.ofdm.frequency_mode='full';
cfg.pattern.file=patternPath; cfg.cdl.profile_mat=profilePath;
cfg.subarray.m=2; cfg.subarray.n=2; cfg.subarray.sweeps=1;
cfg.subarray.displacement_grid_lambda=0; cfg.subarray.yaw_grid_deg=0; cfg.subarray.tilt_grid_deg=0;
cfg.pca.x_grid_lambda=0; cfg.pca.z_grid_lambda=0;
cfg.orientation.yaw_grid_deg=0; cfg.orientation.tilt_grid_deg=0; cfg.orientation.sweeps=1;
cfg.save_results=true; cfg.save_figures=false; cfg.plot.enabled=false;
cfg.output.save_realizations=true; cfg.output.save_geometry=true;
cfg.output_dir=fullfile(testDir,'zero');
[zero,~,~]=simulateFIM(cfg,'cdf');
for name={'pca','orientation','fim'}
    d=zero.(name{1}); fourClose(d.rate,zero.fixed.rate,1e-10,'Four-way zero-grid identity');
end
checks={'Four-way zero-grid planar baseline identity'};
cfg.subarray.displacement_grid_lambda=[-0.2 0 0.2]; cfg.subarray.yaw_grid_deg=[-10 0 10];
cfg.pca.x_grid_lambda=[-0.2 0 0.2]; cfg.pca.z_grid_lambda=[-0.2 0 0.2];
cfg.orientation.yaw_grid_deg=[-10 0 10]; cfg.output_dir=fullfile(testDir,'moving');
[moving,~,~]=simulateFIM(cfg,'cdf');
fourClose(moving.fixed.rate,zero.fixed.rate,0,'Searches must not change fading or the fixed baseline');
checks{end+1}='Four searches share unchanged UE/fading/reference channel';
for mc=1:cfg.MC
    pcfg=moving.cfg; pcfg.deform.family='pca'; pcfg.comparison.include_pca=false;
    [pg,valid]=fimGeometry('build',moving.pca_model,moving.pca_parameters(mc,:,1),pcfg);
    assert(valid,'Saved PCA geometry invalid.');
    og=fimGeometry('orientation',moving.geometry.fixed,moving.orientation_states,moving.orientation_selection(mc,:,1),moving.cfg);
    geometries={moving.geometry.fixed,pg,og,moving.geometry.fim{mc,1}};
    names={'fixed','pca','orientation','fim'};
    fourClose(og.positions_m,moving.geometry.fixed.positions_m,0,'Pattern reconfiguration keeps physical positions fixed');
    for j=1:4
        [H,info]=fimChannel('geometry',moving.realizations{mc},moving.pattern,geometries{j},moving.cfg);
        assert(info.cp_feasible,'Rebuilt four-way geometry violates CP.');
        metric=uplinkOFDMRate(H,moving.snr_dB,moving.cfg); d=moving.(names{j});
        fourClose(metric.rate.',d.rate(mc,:,1),1e-10,'Saved four-way metrics rebuild from selected geometry');
        assert(all(metric.snrPerTone(:)>=metric.sinrPerTone(:)-1e-10),'Thermal SNR must be at least actual SINR.');
    end
end
checks{end+1}='Four selected geometries reconstruct receiver scores and preserve physical CP';
checks{end+1}='Pattern reconfiguration changes orientation only; SNR and SINR are distinct';
one=cfg; one.K=1; one.MC=1; one.save_results=false;
[singleResult,~,~]=simulateFIM(one,'cdf');
for name={'fixed','pca','orientation','fim'}
    d=singleResult.(name{1}); fourClose(d.snr,d.sinr,1e-12,'Single-user white-noise SNR/SINR identity');
end
checks{end+1}='Single-user white-noise SNR equals SINR';

% Independently known pattern-rotation gain, without relying on the main search model.
known=prepareFIMConfig(one); known.power.snr_reference='per_subcarrier';
[P,T]=meshgrid(0:360,0:180);
pat=struct('phi',0:360,'theta',0:180,'G',1+0.5*sind(T).*cosd(P),'native_to_body',eye(3),'front_only',false);
realization.paths={struct('U',[0;1;0],'alpha',1,'tau',0)};
st=struct('yaw',[0;-90],'tilt',[0;0],'count',2,'zero',1);
candidate=fimChannel('candidates',realization,pat,st,known);
for receiver={'zf','mmse-irc'}
    known.receiver=receiver{1};
    [~,selection,metric]=fimDeformation('orientation',candidate,st,0,known);
    assert(all(selection==2),'Pattern search missed a known best direction.');
    expected=known.ofdm.cp_efficiency*log2(1+1.5*known.N);
    fourClose(metric.rate,expected,1e-10,'Analytical pattern-rotation gain');
end
checks{end+1}='ZF/MMSE pattern search recovers an independent analytical directional gain';

% Unequal last batch, saved statistics, and one missing checkpoint exercise actual resume.
batch=cfg; batch.output_dir=fullfile(testDir,'batch');
batch.subarray.displacement_grid_lambda=0; batch.subarray.yaw_grid_deg=0;
batch.pca.x_grid_lambda=0; batch.pca.z_grid_lambda=0; batch.orientation.yaw_grid_deg=0;
batch.study.K_values=[1 2]; batch.study.MC_per_case=[3 3]; batch.study.batch_MC=2;
batch.study.make_plots=false; batch.study.resume_directory='';
[first,~,~]=runFIMStudy(batch);
for c=1:2
    r=first.cases{c}; total=zeros(size(r.fixed.meanRatePerToneNoCP)); M=0;
    for b=1:numel(first.batch_files{c})
        loaded=load(first.batch_files{c}{b},'checkpoint'); part=loaded.checkpoint.result;
        total=total+part.cfg.MC*part.fixed.meanRatePerToneNoCP; M=M+part.cfg.MC;
    end
    fourClose(total/M,r.fixed.meanRatePerToneNoCP,1e-12,'Last batch must have its actual MC weight');
end
delete(first.batch_files{1}{2});
batch.study.resume_directory=first.directory;
[resumed,~,~]=runFIMStudy(batch);
for c=1:2
    for name={'fixed','pca','orientation','fim'}
        a=first.cases{c}.(name{1}); b=resumed.cases{c}.(name{1});
        fourClose(a.rate,b.rate,0,'Resume rates'); fourClose(a.snr,b.snr,0,'Resume SNR');
        fourClose(a.throughput_bps,b.throughput_bps,0,'Resume throughput');
    end
end
checks{end+1}='Four-way batch means use actual MC counts, including the smaller final batch';
checks{end+1}='Missing-checkpoint resume reproduces all four rates/SNR/throughput exactly';
bad=cfg; bad.pattern.file=fullfile(testDir,'missing_measured_pattern.txt'); rejected=false;
try, fimChannel('prepare',prepareFIMConfig(bad)); catch, rejected=true; end
assert(rejected,'Missing formal input must not silently fall back to synthetic data.');
checks{end+1}='Missing formal input rejects clearly; temporary fixtures are removed on return';
report=struct('passed',true,'checks',{checks},'test_fixture_only',true,'temporary_files_removed_on_return',true);
end

function removeFourWayDirectory(path)
if exist(path,'dir')==7, rmdir(path,'s'); end
end

function fourClose(a,b,tolerance,label)
assert(isequal(size(a),size(b)),['Size mismatch: ' label]);
assert(all(isfinite(a(:))) && all(isfinite(b(:))),['Nonfinite comparison: ' label]);
assert(max(abs(a(:)-b(:)))<=tolerance*max(1,max(abs(b(:)))),['Failed: ' label]);
end
function report = sa_test_run_subarray_sanity_checks()
%RUN_SUBARRAY_SANITY_CHECKS 刚性子阵平移/转动的独立数值夹具检查。
% 不读取实测 TXT，不调用 Julia，不创建测试数据或性能结论文件。
% 检查公开几何/信道/接收/优化接口；人工输入不代表真实天线性能。
previousRng = rng;
cleanup = onCleanup(@() rng(previousRng));
rng(20261030,'twister');
checks = {};
cfg = sa_test_fixtureConfig();
model = fimGeometry('prepare',cfg);
[flat,valid] = fimGeometry('build',model,model.zero,cfg);
assert(valid,'Zero subarray state must be feasible.');
assert(model.parameter_count==3*cfg.subarray.m*cfg.subarray.n,'Three variables per subarray are required.');
sa_test_assertClose(flat.positions_m,model.flat_positions_m,0,'Exact planar zero state');
sa_test_assertClose(flat.normals,repmat([0 1 0],cfg.N,1),0,'Exact planar normals');
for k = 1:cfg.N
    sa_test_assertClose(flat.R(:,:,k),eye(3),0,'Exact planar local frames');
end
% Nh=4,Nv=4,m=2,n=2: horizontal element and horizontal block indices vary first.
expectedIDs = [1;1;2;2;1;1;2;2;3;3;4;4;3;3;4;4];
sa_test_assertClose(double(model.subarray_ids(:)),expectedIDs,0,'Horizontal-first rectangular partition');
for b = 1:4
    assert(isequal(sort(model.element_groups{b}(:)),find(expectedIDs==b)),'Element group differs from block ID map.');
end
checks{end+1} = 'Exact flat limit, 3mn variables, complete horizontal-first rectangular partition';

% Mixed signs and simultaneous rotations: centers move only along ORIGINAL +y.
blockParameters = [0.2 20 10;-0.2 -20 -10;0.1 -10 15;-0.1 10 -15];
parameters = reshape(blockParameters.',1,[]);
[moved,valid] = fimGeometry('build',model,parameters,cfg);
assert(valid,'Moderate rigid-subarray fixture should be feasible.');
for b = 1:4
    idx = model.element_groups{b};
    originalCenter = mean(flat.positions_m(idx,:),1);
    movedCenter = mean(moved.positions_m(idx,:),1);
    expectedCenter = originalCenter+[0 blockParameters(b,1)*model.lambda_m 0];
    sa_test_assertClose(movedCenter,expectedCenter,1e-12,'Signed displacement along original normal');
    sa_test_assertClose(moved.subarray_centers_m(b,:),expectedCenter,1e-12,'Saved moved block center');
    sa_test_assertClose(sa_test_pairDistances(moved.positions_m(idx,:)),sa_test_pairDistances(flat.positions_m(idx,:)),1e-12,'Rigid block pair distances');
    yaw = blockParameters(b,2);
    tilt = blockParameters(b,3);
    expectedNormal = [sind(yaw)*cosd(tilt),cosd(yaw)*cosd(tilt),sind(tilt)];
    sa_test_assertClose(moved.normals(idx,:),repmat(expectedNormal,numel(idx),1),1e-12,'Shared block yaw and tilt normal');
    for k = reshape(idx,1,[])
        R = moved.R(:,:,k);
        sa_test_assertClose(R'*R,eye(3),1e-12,'Block rotation orthogonality');
        sa_test_assertClose(det(R),1,1e-12,'Right-handed rigid block frame');
        sa_test_assertClose(R,moved.R(:,:,idx(1)),0,'One rigid orientation per block');
    end
end
sa_test_assertClose(moved.subarray_parameters,blockParameters,0,'Stored parameter ordering and units');
checks{end+1} = 'Both translation signs, original-normal center motion, rigid distances and common block frames';

% Rigid blocks must not inherit continuous-membrane strain/curvature/element excursion bounds.
strict = cfg;
strict.deform.max_displacement_lambda = 0;
strict.deform.max_normal_angle_deg = 0;
strict.deform.max_strain = 0;
strict.deform.max_curvature_per_lambda = 0;
[~,valid] = fimGeometry('build',model,parameters,strict);
assert(valid,'Continuous-surface limits must not reject a valid rigid-subarray state.');
bad = parameters;
bad(1) = 2;
sa_test_assertRejectedGeometry(model,bad,cfg,'Center translation range');
bad = parameters;
bad(2) = 80;
sa_test_assertRejectedGeometry(model,bad,cfg,'Yaw range');
bad = parameters;
bad(3) = 80;
sa_test_assertRejectedGeometry(model,bad,cfg,'Tilt range');
spacing = cfg;
spacing.subarray.min_spacing_lambda = 0.75;
sa_test_assertRejectedGeometry(model,model.zero,spacing,'Minimum element separation');
invalidPartition = cfg;
invalidPartition.subarray.m = 3;
sa_test_assertRejectedConfig(invalidPartition,'Nondivisible horizontal partition');
invalidPartition = cfg;
invalidPartition.subarray.n = cfg.Nv+1;
sa_test_assertRejectedConfig(invalidPartition,'Too many vertical blocks');
checks{end+1} = 'Separate subarray limits, invalid ranges/partition and actual element-spacing rejection';

% Whole-array block: cardinal rotations give simple independent expected positions.
one = cfg;
one.subarray.m = 1;
one.subarray.n = 1;
one.subarray.yaw_grid_deg = [-90 0 90];
one.subarray.tilt_grid_deg = [-90 0 90];
one = prepareFIMConfig(one);
oneModel = fimGeometry('prepare',one);
[yaw90,valid] = fimGeometry('build',oneModel,[0.2 90 0],one);
assert(valid,'Whole-array 90-degree yaw remains rigid and separated.');
P = oneModel.flat_positions_m;
expected = [zeros(one.N,1),0.2*oneModel.lambda_m-P(:,1),P(:,3)];
sa_test_assertClose(yaw90.positions_m,expected,1e-12,'Yaw rotates positions about own translated center');
[tilt90,valid] = fimGeometry('build',oneModel,[0.2 0 90],one);
assert(valid,'Whole-array 90-degree tilt remains rigid and separated.');
expected = [P(:,1),0.2*oneModel.lambda_m-P(:,3),zeros(one.N,1)];
sa_test_assertClose(tilt90.positions_m,expected,1e-12,'Tilt rotates positions about own translated center');
elements = cfg;
elements.subarray.m = elements.Nh;
elements.subarray.n = elements.Nv;
elements = prepareFIMConfig(elements);
elementsModel = fimGeometry('prepare',elements);
p = elementsModel.zero;
p(2:3:end) = 20;
p(3:3:end) = -10;
[pointBlocks,valid] = fimGeometry('build',elementsModel,p,elements);
assert(valid,'Single-element blocks should be valid.');
sa_test_assertClose(pointBlocks.positions_m,elementsModel.flat_positions_m,0,'Single-element center rotation leaves positions fixed');
checks{end+1} = 'Whole-array and oneUE-element block limits; cardinal rigid rotations change physical positions';

% Deployment is composed outside block motion, and never rotates the local lift axis twice.
deployed = cfg;
deployed.bs.position_xy_m = [13 -7];
deployed.bs.height_m = 19;
deployed.bs.azimuth_deg = 23;
deployed.bs.downtilt_deg = 11;
deployed = prepareFIMConfig(deployed);
deployedModel = fimGeometry('prepare',deployed);
[deployedGeom,valid] = fimGeometry('build',deployedModel,parameters,deployed);
assert(valid,'Deployed block fixture should be feasible.');
scene = fimScene('prepare',deployed);
expected = bsxfun(@plus,moved.positions_m*scene.bs.R_global_from_array.',scene.bs.position_m);
sa_test_assertClose(deployedGeom.positions_global_m,expected,1e-12,'Composed BS position and orientation');
for k = 1:cfg.N
    sa_test_assertClose(deployedGeom.R_global(:,:,k),scene.bs.R_global_from_array*moved.R(:,:,k),1e-12,'Composed BS/block frames');
end
checks{end+1} = 'Global deployment composes with local block motion exactly once';

% Zero blocks reproduce the already-tested planar channel interface.
pattern = sa_test_fixturePattern(false);
realization = sa_test_fixtureRealization(cfg.K);
states = struct('yaw',0,'tilt',0,'zero',1,'count',1);
legacy = planarReference(realization,pattern,cfg);
legacy = reshape(legacy,cfg.N,cfg.K,cfg.ofdm.num_tones);
[Hflat,info] = fimChannel('geometry',realization,pattern,flat,cfg);
assert(info.cp_feasible,'Numerical fixture needs sufficient CP.');
sa_test_assertClose(Hflat,legacy,1e-11,'Subarray planar channel versus fixed array');

% Independent oneUE-path formula checks actual rigid positions at every frequency.
oneUE = cfg;
oneUE.K = 1;
oneUE = prepareFIMConfig(oneUE);
path = struct('U',[0.6;0.8;0],'alpha',0.7+0.2i,'tau',50e-9);
singleRealization = struct('paths',{{path}});
[Hmove,info] = fimChannel('geometry',singleRealization,sa_test_fixturePattern(true),moved,oneUE);
assert(info.cp_feasible,'Physical oneUE-path fixture should satisfy CP.');
for f = 1:oneUE.ofdm.num_tones
    nu = oneUE.ofdm.freq_offsets_Hz(f);
    expected = path.alpha*exp(1i*2*pi*(oneUE.fc_Hz+nu)*(moved.positions_m*path.U/oneUE.c0))*exp(-1i*2*pi*nu*path.tau);
    sa_test_assertClose(Hmove(:,1,f),expected,1e-11,'Actual rigid-block spatial and excess-delay phases');
end
[HsingleFlat,~] = fimChannel('geometry',singleRealization,sa_test_fixturePattern(true),flat,oneUE);
a = uplinkOFDMRate(HsingleFlat,3,oneUE);
b = uplinkOFDMRate(Hmove,3,oneUE);
sa_test_assertClose(a.rate,b.rate,1e-11,'Constant-gain oneUE-path phase changes cannot manufacture gain');
checks{end+1} = 'Planar channel identity, physical wideband phases and no artificial phase-only oneUE-path gain';

% Known directional case proves the optimizer actually changes block orientation.
directional = one;
directional.K = 1;
directional.subarray.displacement_grid_lambda = 0;
directional.subarray.tilt_grid_deg = 0;
directional.subarray.sweeps = 1;
directional.opt.score_mode = 'exact';
directional = prepareFIMConfig(directional);
directionalModel = fimGeometry('prepare',directional);
frontPath = struct('U',[0;1;0],'alpha',1,'tau',0);
frontRealization = struct('paths',{{frontPath}});
work = fimDeformation('prepare',frontRealization,pattern,directionalModel,directional);
baseline = uplinkOFDMRate(work.Hfixed,0,directional);
[~,chosen,best] = fimDeformation('optimize',work,0,directional);
sa_test_assertClose(chosen.parameters(2),-90,1e-11,'Known directional best yaw');
sa_test_assertClose(best.rate,log2(1+1.3*directional.N),1e-10,'Analytical oneUE-user directional array gain');
assert(best.sumRate>baseline.sumRate+0.01,'Subarray search failed to find a known positive gain.');
checks{end+1} = 'Optimizer recovers a known directional gain and physical best block yaw';

% Full receiver and Gram-verified coordinate searches must agree for BOTH objectives.
for mode = {'sum_rate','max_min_rate'}
    exact = cfg;
    exact.opt_mode = mode{1};
    exact.opt.score_mode = 'exact';
    work = fimDeformation('prepare',realization,pattern,model,exact);
    reference = uplinkOFDMRate(work.Hfixed,2,exact);
    [Hselected,selected,metric,trace] = fimDeformation('optimize',work,2,exact);
    [Hrebuilt,selectedInfo] = fimChannel('geometry',realization,pattern,selected,exact);
    assert(selected.feasible && selectedInfo.cp_feasible,'Selected blocks must be feasible.');
    sa_test_assertClose(Hselected,Hrebuilt,1e-11,'One selected block configuration constructs all frequency channels');
    assert(sa_test_objective(metric,mode{1})>=sa_test_objective(reference,mode{1})-1e-9,'Selected sa_test_objective fell below fixed baseline.');
    assert(all(diff(trace.values)>=-1e-9),'Accepted block-coordinate steps must be monotonic.');
    accelerated = exact;
    accelerated.opt.score_mode = 'gram_verified';
    gramWork = fimDeformation('prepare',realization,pattern,model,accelerated);
    [~,~,gramMetric] = fimDeformation('optimize',gramWork,2,accelerated);
    sa_test_assertClose(sa_test_objective(metric,mode{1}),sa_test_objective(gramMetric,mode{1}),1e-9,'Exact versus Gram-verified block sa_test_objective');
end
checks{end+1} = 'Both objectives preserve baseline, agree across scoring modes, and share one geometry across tones';

% Unequal RB scheduling must use each UE own tone count and full receiver scoring.
rb = cfg;
rb.ofdm.resource_mode = 'rb';
rb.ofdm.num_rb = 3;
rb.ofdm.rb_allocation = {1,2:3};
rb.power.snr_reference = 'total_over_tone_noise';
rb.subarray.yaw_grid_deg = 0;
rb.subarray.tilt_grid_deg = 0;
rb.subarray.sweeps = 1;
rb = prepareFIMConfig(rb);
rbModel = fimGeometry('prepare',rb);
rbWork = fimDeformation('prepare',realization,pattern,rbModel,rb);
reference = uplinkOFDMRate(rbWork.Hfixed,10,rb);
[~,~,metric,trace] = fimDeformation('optimize',rbWork,10,rb);
sa_test_assertClose(metric.txPowerTotal,10*rb.noiseVar*ones(rb.K,1),1e-11,'Per-UE fixed total transmit power');
sa_test_assertClose(metric.userNumActive,[12;24],0,'Unequal per-user resource counts');
assert(all(metric.sinrPerTone(~rb.ofdm.user_sample_mask)==0),'Unscheduled tones must not contribute SINR.');
assert(trace.gram_scores==0,'Nonshared RB allocation must use scheduling-aware receiver evaluation.');
assert(metric.sumRate>=reference.sumRate-1e-9,'RB-aware block search lost baseline.');
checks{end+1} = 'Unequal RBs retain per-UE total power, inactive zeros and allocation-aware scoring';

report = struct('passed',true,'checks',{checks},'test_fixture_only',true,'measured_pattern_tested',false,'actual_julia_cdl_tested',false);
fprintf('PASS: %d subarray numerical check groups; synthetic fixtures only.\n',numel(checks));
for i = 1:numel(checks)
    fprintf('  %d. %s\n',i,checks{i});
end
end

function cfg = sa_test_fixtureConfig()
cfg = configFIMSimulation();
cfg.receiver='zf'; cfg.noiseVar=1; cfg.power.snr_reference='per_subcarrier';
cfg.comparison.include_pca=false; cfg.subarray.m=1; cfg.subarray.n=1;
cfg.orientation.yaw_grid_deg=0; cfg.orientation.tilt_grid_deg=0;
cfg.Nh = 4;
cfg.Nv = 4;
cfg.dh_lambda = 0.5;
cfg.dv_lambda = 0.5;
cfg.K = 2;
cfg.users.mode = 'angles';
cfg.bs.position_xy_m = [0 0];
cfg.bs.ground_height_m = 0;
cfg.bs.height_m = 0;
cfg.bs.azimuth_deg = 0;
cfg.bs.downtilt_deg = 0;
cfg.pathloss.enabled = false;
cfg.MC = 1;
cfg.fast_mode = false;
cfg.fim.mode = 'deformation';
cfg.deform.family = 'subarray';
cfg.subarray.m = 2;
cfg.subarray.n = 2;
cfg.subarray.displacement_grid_lambda = [-0.3 0 0.3];
cfg.subarray.yaw_grid_deg = [-30 0 30];
cfg.subarray.tilt_grid_deg = [-20 0 20];
cfg.subarray.min_spacing_lambda = 0.05;
cfg.subarray.sweeps = 1;
cfg.ofdm.resource_mode = 'indices';
cfg.ofdm.active_indices = [-1 0 1];
cfg.ofdm.rb_allocation = [];
cfg.ofdm.frequency_mode = 'full';
cfg.ofdm.include_cp_overhead = false;
cfg.save_results = false;
cfg.save_figures = false;
cfg.plot.enabled = false;
cfg = prepareFIMConfig(cfg);
end

function pattern = sa_test_fixturePattern(constant)
[P,T] = meshgrid(0:360,0:180);
G = ones(size(P));
if ~constant
    G = 1+0.3*sind(T).*cosd(P)+0.2*cosd(T).^2;
end
pattern = struct('phi',0:360,'theta',0:180,'G',G,'native_to_body',eye(3),'front_only',false);
end

function realization = sa_test_fixtureRealization(K)
realization.paths = cell(K,1);
for k = 1:K
    U = [0.15+0.2*k,-0.4;0.9,0.8;0.1,0.15*k];
    U = bsxfun(@rdivide,U,sqrt(sum(U.^2,1)));
    realization.paths{k} = struct('U',U,'alpha',[0.8+0.1i*k;-0.2+0.4i],'tau',[0;40e-9]);
end
end

function distances = sa_test_pairDistances(P)
distances = [];
for a = 1:size(P,1)
    for b = a+1:size(P,1)
        distances(end+1,1) = norm(P(a,:)-P(b,:));
    end
end
end

function sa_test_assertRejectedGeometry(model,parameters,cfg,label)
rejected = false;
try
    [~,valid] = fimGeometry('build',model,parameters,cfg);
    rejected = ~valid;
catch
    rejected = true;
end
assert(rejected,['Expected geometry rejection: ' label]);
end

function sa_test_assertRejectedConfig(cfg,label)
rejected = false;
try
    prepareFIMConfig(cfg);
catch
    rejected = true;
end
assert(rejected,['Expected config rejection: ' label]);
end

function value = sa_test_objective(metric,mode)
if strcmp(mode,'sum_rate')
    value = metric.sumRate;
else
    value = metric.worstRate;
end
end

function sa_test_assertClose(actual,expected,tolerance,label)
assert(isequal(size(actual),size(expected)),['Shape mismatch: ' label]);
if isempty(actual)
    return
end
errorSize = max(abs(actual(:)-expected(:)));
scale = max(1,max(abs(expected(:))));
assert(errorSize<=tolerance*scale,['Failed: ' label]);
end

function report = rx_test_run_receiver_sanity_checks()
%RUN_RECEIVER_SANITY_CHECKS MMSE-IRC 与 ZF 的独立数值和优化一致性检查。
% 人工信道/方向图夹具；不读取实测 TXT，不调用 Julia，不写测试数据。
% 参考解在天线域逐用户构造干扰协方差，独立于生产代码的小矩阵实现。
previousRng = rng;
cleanup = onCleanup(@() rng(previousRng));
rng(20261030,'twister');
checks = {};

% 1. 默认兼容、别名和非法输入。
defaults = configFIMSimulation();
legacy = defaults;
legacy.receiver = 'zf';
if isfield(legacy,'irc')
    legacy = rmfield(legacy,'irc');
end
legacy = prepareFIMConfig(legacy);
assert(strcmp(legacy.receiver,'zf') && isempty(legacy.irc.external_covariance),'Old configurations must gain an empty IRC covariance.');
cfg = rx_test_fixtureConfig(6,3,5);
alias = cfg;
alias.receiver = 'mmse_irc';
alias = prepareFIMConfig(alias);
assert(strcmp(alias.receiver,'mmse-irc'),'MMSE-IRC alias must normalize before use.');
bad = cfg; bad.receiver = 'unsupported-receiver'; rx_test_assertRejected(bad,'Unknown receiver');
bad = cfg; bad.irc.external_covariance = eye(5); rx_test_assertRejected(bad,'Wrong covariance dimension');
bad = cfg; bad.irc.external_covariance = eye(6); bad.irc.external_covariance(1,2) = 1;
rx_test_assertRejected(bad,'Non-Hermitian covariance');
bad = cfg; bad.irc.external_covariance = diag([-0.1 ones(1,5)]);
rx_test_assertRejected(bad,'Indefinite external covariance');
bad = cfg; bad.noiseVar = 2e-15;
bad.irc.external_covariance = bad.noiseVar*diag([-0.1 ones(1,5)]);
rx_test_assertRejected(bad,'Indefinite covariance cannot pass because its absolute power is tiny');
bad = cfg; bad.irc.external_covariance = eye(6); bad.irc.external_covariance(1,1) = NaN;
rx_test_assertRejected(bad,'Nonfinite covariance');
bad = cfg; bad.irc.external_covariance = repmat(eye(6),[1 1 2]);
rx_test_assertRejected(bad,'Unsupported frequency-dependent covariance');
checks{end+1} = 'Legacy ZF cfg and missing IRC defaults, MMSE alias, invalid receiver/covariance rejection';

% 2. 白噪声：每个频点/用户与天线域的独立干扰抑制参考解一致。
H = (randn(cfg.N,cfg.K,cfg.ofdm.num_tones)+1i*randn(cfg.N,cfg.K,cfg.ofdm.num_tones))/sqrt(2);
for snr_dB = [-60 0 40]
    metric = uplinkOFDMRate(H,snr_dB,cfg);
    reference = rx_test_antennaReference(H,snr_dB,cfg);
    rx_test_assertClose(metric.sinrPerTone,reference,2e-8,'White-noise antenna-space reference');
    assert(all(isfinite(metric.sinrPerTone(:))) && all(metric.sinrPerTone(:)>=0),'White-noise SINR must stay finite/nonnegative.');
end
checks{end+1} = 'Low/medium/high SNR match independent per-user antenna-space MMSE reference';

% 极小信道不能在增广 QR 的上半块提取中丢失有效数字。
tiny = rx_test_fixtureConfig(4,1,1);
tiny.noiseVar = 2e-15;
v = [1;1i;0.3;0.2i];
tiny.irc.external_covariance = tiny.noiseVar*(v*v');
tiny = prepareFIMConfig(tiny);
for amplitude = [1e-8 1e-16]
    tinyH = amplitude*[1+2i;-3+0.5i;0.2-1i;2];
    tinyMetric = uplinkOFDMRate(tinyH,-40,tiny);
    tinyReference = rx_test_antennaReference(tinyH,-40,tiny);
    assert(tinyReference>0 && tinyMetric.sinrPerTone>0,'Tiny nonzero channel must retain a positive SINR.');
    relativeError = abs(tinyMetric.sinrPerTone-tinyReference)/tinyReference;
    assert(relativeError<1e-10,'Tiny-channel SINR lost relative precision.');
end
checks{end+1} = 'Tiny complex channels retain nonzero SINR with an independent relative-error check';

% 3. 空间有色外部干扰，包含实际物理功率单位和复数互相关。
colored = cfg;
colored.noiseVar = 2e-15;
B = randn(cfg.N,2)+1i*randn(cfg.N,2);
colored.irc.external_covariance = colored.noiseVar*(2*(B*B')+diag(0:cfg.N-1));
colored = prepareFIMConfig(colored);
for snr_dB = [-15 5 35]
    metric = uplinkOFDMRate(H,snr_dB,colored);
    reference = rx_test_antennaReference(H,snr_dB,colored);
    rx_test_assertClose(metric.sinrPerTone,reference,2e-8,'Colored interference antenna-space reference');
    assert(all(metric.externalInterference(:)>=0),'External interference power must be nonnegative.');
    denominator = metric.interference+metric.noise+metric.externalInterference;
    rx_test_assertClose(metric.desired./denominator,metric.sinrPerTone,1e-11,'Actual three-part denominator');
end
checks{end+1} = 'Complex colored covariance in physical units and actual desired/interference/noise powers';

% 4. 单用户退化为广义匹配滤波，方向与协方差相关而非仅 H。
oneUE = rx_test_fixtureConfig(4,1,2);
v = [1;1i;0.3;0.2i];
oneUE.irc.external_covariance = 7*(v*v');
oneUE = prepareFIMConfig(oneUE);
h = [0.5+0.2i;1;0.3i;-0.7];
singleH = repmat(h,[1 1 oneUE.ofdm.num_tones]);
singleMetric = uplinkOFDMRate(singleH,4,oneUE);
expected = 10^(4/10)*real(h'*((eye(4)+oneUE.irc.external_covariance)\h));
rx_test_assertClose(singleMetric.sinrPerTone,expected*ones(1,2),1e-11,'Single-user generalized matched filter');
assert(all(singleMetric.interference(:)==0),'Single UE has no served-user interference.');
checks{end+1} = 'Single-user MMSE-IRC equals the analytical generalized matched filter';

% 5. 完全相关、零列和全零信道不被正则化伪装为满秩信道。
singular = rx_test_fixtureConfig(4,2,2);
singularH = ones(4,2,2);
metric = uplinkOFDMRate(singularH,0,singular);
rx_test_assertClose(metric.sinrPerTone,0.8*ones(2,2),1e-11,'Collinear-user residual interference');
assert(all(metric.receiver.rank==1) && all(metric.receiver.rank_deficient),'Rank diagnostics must describe H, not its regularized Gram.');
singularH(:,2,:) = 0;
metric = uplinkOFDMRate(singularH,0,singular);
rx_test_assertClose(metric.sinrPerTone,[4 4;0 0],1e-11,'One zero user column');
zeroMetric = uplinkOFDMRate(zeros(4,2,2),0,singular);
assert(all(zeroMetric.sinrPerTone(:)==0) && all(zeroMetric.rate==0),'Zero channel must have zero SINR and rate.');
assert(all(zeroMetric.receiver.rank==0),'All-zero channel rank must remain zero.');
checks{end+1} = 'Collinear, zero-user and zero-channel cases retain true rank and residual interference';

% 6. 对相同信道/功率/协方差，最大 SINR 线性合并不劣于原 ZF。
for useColored = 0:1
    comparison = cfg;
    if useColored
        comparison.irc.external_covariance = colored.irc.external_covariance/colored.noiseVar;
    end
    comparison = prepareFIMConfig(comparison);
    zfCfg = comparison;
    zfCfg.receiver = 'zf';
    for snr_dB = [-20 0 25]
        mmseMetric = uplinkOFDMRate(H,snr_dB,comparison);
        zfMetric = uplinkOFDMRate(H,snr_dB,zfCfg);
        scale = max(1,max(zfMetric.sinrPerTone(:)));
        assert(all(mmseMetric.sinrPerTone(:)>=zfMetric.sinrPerTone(:)-1e-8*scale),'MMSE-IRC underperformed ZF on the same H and covariance.');
    end
end
checks{end+1} = 'For identical H/power/covariance, each MMSE-IRC stream SINR is no lower than ZF';

% 7. 资源不均分：真实占用频点数分功率，未调度 RB 严格为零。
scheduled = cfg;
scheduled.ofdm.resource_mode = 'rb';
scheduled.ofdm.num_rb = 4;
scheduled.ofdm.rb_allocation = {[1 2],[2 3],3};
scheduled.power.snr_reference = 'total_over_tone_noise';
scheduled.noiseVar = colored.noiseVar;
scheduled.irc.external_covariance = colored.irc.external_covariance;
scheduled.ofdm.include_cp_overhead = true;
scheduled = prepareFIMConfig(scheduled);
flatH = H(:,:,1);
scheduledH = repmat(flatH,[1 1 scheduled.ofdm.num_tones]);
metric = uplinkOFDMRate(scheduledH,13,scheduled);
rx_test_assertClose(metric.sinrPerTone,rx_test_antennaReference(scheduledH,13,scheduled),1e-10,'Scheduled users and unequal per-tone power');
rx_test_assertClose(metric.txPowerTotal/scheduled.noiseVar,10^(13/10)*ones(3,1),1e-11,'Each UE keeps its own total transmit power');
rx_test_assertClose(metric.userNumActive,[24;24;12],0,'Unequal allocated tone counts');
assert(all(metric.sinrPerTone(~scheduled.ofdm.user_sample_mask)==0),'Unscheduled UE/tone SINR must remain zero.');
assert(all(metric.receiver.num_scheduled(37:48)==0),'Last RB must remain inactive.');
assert(all(metric.noise(:,37:48)==0,'all') && all(metric.externalInterference(:,37:48)==0,'all'),'Inactive RBs must not be evaluated as scheduled streams.');
manualRate = scheduled.ofdm.cp_efficiency*sum(log1p(metric.sinrPerTone)/log(2),2)/48;
rx_test_assertClose(metric.rate,manualRate,1e-11,'CP-adjusted scheduled full-band spectral efficiency');
rx_test_assertClose(metric.sumThroughput_bps,scheduled.ofdm.scs_Hz*48*sum(manualRate),1e-11,'Active-band throughput units');
checks{end+1} = 'Unequal shared/disjoint/inactive RBs preserve per-UE total power, CP and bandwidth units';

% 8. 常频响条件下，分组是精确等价；分功率仍按48个真实频点。
grouped = scheduled;
grouped.ofdm.frequency_mode = 'grouped';
grouped.ofdm.num_frequency_samples = 7;
grouped = prepareFIMConfig(grouped);
groupH = repmat(flatH,[1 1 grouped.ofdm.num_tones]);
groupMetric = uplinkOFDMRate(groupH,13,grouped);
rx_test_assertClose(groupMetric.rate,metric.rate,1e-11,'Full versus weighted grouped receiver rates');
rx_test_assertClose(groupMetric.txPowerPerTone,metric.txPowerPerTone,1e-11,'Grouped samples do not change power distribution');
rx_test_assertClose(groupMetric.sumThroughput_bps,metric.sumThroughput_bps,1e-11,'Full versus grouped throughput');
checks{end+1} = 'Grouped quadrature matches full-frequency results for a flat channel with unequal RB allocations';

% 9. 改变统一功率单位不能改变接收结果，也不能被绝对容差吞掉协方差。
unitCfg = cfg;
unitCfg.irc.external_covariance = colored.irc.external_covariance/colored.noiseVar;
unitCfg = prepareFIMConfig(unitCfg);
scaledCfg = unitCfg;
scaledCfg.noiseVar = 1e-18;
scaledCfg.irc.external_covariance = unitCfg.irc.external_covariance*1e-18;
scaledCfg = prepareFIMConfig(scaledCfg);
a = uplinkOFDMRate(H,7,unitCfg);
b = uplinkOFDMRate(H,7,scaledCfg);
rx_test_assertClose(a.sinrPerTone,b.sinrPerTone,1e-10,'Covariance/thermal/transmit unit scaling invariance');
checks{end+1} = 'Joint covariance, thermal-noise and transmit-power unit scaling leaves SINR invariant';

% 11. 新子阵形态：两目标、两评分设置共享同一接收口径和物理形态。
subarray = rx_test_fixtureConfig(4,2,3);
subarray.fim.mode = 'deformation';
subarray.deform.family = 'subarray';
subarray.subarray.m = 2;
subarray.subarray.n = 1;
subarray.subarray.displacement_grid_lambda = 0;
subarray.subarray.yaw_grid_deg = [-30 0 30];
subarray.subarray.tilt_grid_deg = 0;
subarray.subarray.sweeps = 1;
subarray.subarray.min_spacing_lambda = 0.05;
subarray.irc.external_covariance = diag([0.1 0.3 0.2 0.4]);
subarray = prepareFIMConfig(subarray);
model = fimGeometry('prepare',subarray);
pattern = rx_test_directionalPattern();
realization.paths = cell(2,1);
realization.paths{1} = struct('U',[0.2;sqrt(0.96);0],'alpha',1,'tau',0);
realization.paths{2} = struct('U',[-0.3;sqrt(0.91);0],'alpha',0.8+0.2i,'tau',20e-9);
for objectiveMode = {'sum_rate','max_min_rate'}
    exact = subarray;
    exact.opt_mode = objectiveMode{1};
    exact.opt.score_mode = 'exact';
    work = fimDeformation('prepare',realization,pattern,model,exact);
    baseline = uplinkOFDMRate(work.Hfixed,0,exact);
    [Hselected,selected,metricExact,traceExact] = fimDeformation('optimize',work,0,exact);
    gated = exact;
    gated.opt.score_mode = 'gram_verified';
    workGated = fimDeformation('prepare',realization,pattern,model,gated);
    [Hgated,selectedGated,metricGated,traceGated] = fimDeformation('optimize',workGated,0,gated);
    rx_test_assertClose(selected.parameters,selectedGated.parameters,0,'Selected block parameters agree');
    rx_test_assertClose(Hselected,Hgated,0,'Selected subarray channel agrees across score settings');
    rx_test_assertClose(metricExact.rate,metricGated.rate,1e-11,'Subarray exact and gated scores agree');
    assert(traceExact.gram_scores==0 && traceGated.gram_scores==0,'MMSE subarray must not use ZF-only Gram scoring.');
    assert(rx_test_objective(metricExact,exact.opt_mode)>=rx_test_objective(baseline,exact.opt_mode)-1e-10,'MMSE subarray must preserve planar baseline.');
    assert(all(diff(traceGated.values)>=-1e-10),'MMSE subarray rx_test_objective must be monotone.');
    [rebuilt,info] = fimChannel('geometry',realization,pattern,selected,exact);
    assert(info.cp_feasible && selected.feasible,'MMSE selected block geometry must remain feasible.');
    rx_test_assertClose(rebuilt,Hselected,1e-12,'Saved one geometry rebuilds all selected frequencies');
end
checks{end+1} = 'Rigid-subarray optimizer shares MMSE scoring, preserves baseline and rebuilds one broadband shape';

report = struct('passed',true,'checks',{checks},'test_fixture_only',true,...
    'measured_pattern_tested',false,'actual_julia_cdl_tested',false);
fprintf('PASS: %d receiver numerical check groups; synthetic fixtures only.\n',numel(checks));
for i = 1:numel(checks)
    fprintf('  %d. %s\n',i,checks{i});
end
end

function cfg = rx_test_fixtureConfig(N,K,F)
cfg = configFIMSimulation();
cfg.receiver='zf'; cfg.noiseVar=1; cfg.power.snr_reference='per_subcarrier';
cfg.comparison.include_pca=false; cfg.subarray.m=1; cfg.subarray.n=1;
cfg.orientation.yaw_grid_deg=0; cfg.orientation.tilt_grid_deg=0;
cfg.Nh = N;
cfg.Nv = 1;
cfg.subarray.m=1; cfg.subarray.n=1;
cfg.K = K;
cfg.receiver = 'mmse-irc';
cfg.noiseVar = 1;
cfg.irc.external_covariance = [];
cfg.power.snr_reference = 'per_subcarrier';
cfg.users.mode = 'angles';
cfg.pathloss.enabled = false;
cfg.bs.height_m = 0;
cfg.bs.downtilt_deg = 0;
cfg.MC = 1;
cfg.ofdm.enabled = true;
cfg.ofdm.resource_mode = 'indices';
cfg.ofdm.active_indices = -floor(F/2)+(0:F-1);
cfg.ofdm.rb_allocation = [];
cfg.ofdm.frequency_mode = 'full';
cfg.ofdm.include_cp_overhead = false;
cfg.save_results = false;
cfg.save_figures = false;
cfg.plot.enabled = false;
cfg = prepareFIMConfig(cfg);
end

function expected = rx_test_antennaReference(H,snr_dB,cfg)
% 独立参考：逐用户排除自身后在 N x N 天线域求解，不使用 Gram 逆。
expected = zeros(cfg.K,cfg.ofdm.num_tones);
q = 10^(snr_dB/10)*ones(cfg.K,1);
if strcmp(cfg.power.snr_reference,'total_over_tone_noise')
    q = q./cfg.ofdm.user_num_active(:);
end
R0 = eye(cfg.N);
if ~isempty(cfg.irc.external_covariance)
    R0 = R0+cfg.irc.external_covariance/cfg.noiseVar;
end
for f = 1:cfg.ofdm.num_tones
    active = find(cfg.ofdm.user_sample_mask(:,f));
    for k = reshape(active,1,[])
        R = R0;
        for j = reshape(active,1,[])
            if j~=k
                R = R+q(j)*(H(:,j,f)*H(:,j,f)');
            end
        end
        h = H(:,k,f);
        expected(k,f) = max(0,q(k)*real(h'*(R\h)));
    end
end
end

function pattern = rx_test_directionalPattern()
[P,T] = meshgrid(0:360,0:180);
pattern = struct('phi',0:360,'theta',0:180,...
    'G',1+0.5*sind(T).*cosd(P),'native_to_body',eye(3),'front_only',false);
end

function rx_test_assertRejected(cfg,label)
rejected = false;
try
    prepareFIMConfig(cfg);
catch
    rejected = true;
end
assert(rejected,['Expected config rejection: ' label]);
end

function value = rx_test_objective(metric,mode)
if strcmp(mode,'sum_rate')
    value = metric.sumRate;
else
    value = metric.worstRate;
end
end

function rx_test_assertClose(actual,expected,tolerance,label)
assert(isequal(size(actual),size(expected)),['Shape mismatch: ' label]);
assert(all(isfinite(actual(:))) && all(isfinite(expected(:))),['Nonfinite comparison: ' label]);
if isempty(actual)
    return
end
errorSize = max(abs(actual(:)-expected(:)));
scale = max(1,max(abs(expected(:))));
assert(errorSize<=tolerance*scale,['Failed: ' label]);
end

function report = pc_test_run_pca_sanity_checks()
%RUN_PCA_SANITY_CHECKS PCA平面位移、物理信道和三阵列配对仿真的数值检查。
% 使用人工夹具，不读取实测TXT或调用Julia，不代表真实天线的性能结论。
% 完整流程的临时TXT/MAT及结果在函数结束（包括报错）时自动删除。
oldRng = rng;
rngCleanup = onCleanup(@() rng(oldRng));
rng(20261031,'twister');
checks = {};
cfg = pc_test_fixtureConfig();
model = fimGeometry('prepare',cfg);
[flat,valid] = fimGeometry('build',model,model.zero,cfg);
assert(valid,'PCA zero state must be feasible.');
assert(model.parameter_count==2*cfg.N,'PCA must have two displacement variables per element.');
pc_test_assertClose(flat.positions_m,model.flat_positions_m,0,'Exact flat PCA positions');
pc_test_assertClose(flat.normals,repmat([0 1 0],cfg.N,1),0,'Exact unchanged PCA normals');
for k = 1:cfg.N
    pc_test_assertClose(flat.R(:,:,k),eye(3),0,'Exact unchanged PCA local frames');
end
checks{end+1} = 'Exact zero displacement, 2N controls, identity element frames';

p = model.zero;
p(1:2:end) = linspace(-0.1,0.1,cfg.N);
p(2:2:end) = linspace(0.1,-0.1,cfg.N);
[moved,valid] = fimGeometry('build',model,p,cfg);
assert(valid,'Moderate PCA displacement should be feasible.');
expected = flat.positions_m;
expected(:,1) = expected(:,1)+model.lambda_m*p(1:2:end).';
expected(:,3) = expected(:,3)+model.lambda_m*p(2:2:end).';
pc_test_assertClose(moved.positions_m,expected,1e-13,'Independent x/z displacement in metres');
pc_test_assertClose(moved.positions_m(:,2),zeros(cfg.N,1),0,'No normal displacement');
pc_test_assertClose(moved.R,flat.R,0,'Movement never rotates the pattern');
pc_test_assertClose(moved.normals,flat.normals,0,'Movement never rotates normals');
one = model.zero;
one(3:4) = [0.1 -0.1];
[oneMoved,valid] = fimGeometry('build',model,one,cfg);
assert(valid,'Single-element move should be feasible.');
others = [1 3:cfg.N];
pc_test_assertClose(oneMoved.positions_m(others,:),flat.positions_m(others,:),0,'Other elements stay fixed');
strict = cfg;
strict.deform.max_displacement_lambda = 0;
strict.deform.max_normal_angle_deg = 0;
strict.deform.max_strain = 0;
strict.deform.max_curvature_per_lambda = 0;
[~,valid] = fimGeometry('build',model,p,strict);
assert(valid,'PCA independent motion must not inherit continuous-surface mechanical limits.');
checks{end+1} = 'Independent signed x/z motion, no y motion or rotation, no membrane constraints';

bad = model.zero;
bad(1) = 0.31;
pc_test_assertRejectedGeometry(model,bad,cfg,'x displacement bound');
bad = model.zero;
bad(2) = -0.31;
pc_test_assertRejectedGeometry(model,bad,cfg,'z displacement bound');
collision = model.zero;
collision(1) = 0.25;
collision(3) = -0.25;
pc_test_assertRejectedGeometry(model,collision,cfg,'actual element collision');
tooClose = cfg;
tooClose.pca.min_spacing_lambda = 0.51;
pc_test_assertRejectedGeometry(model,model.zero,tooClose,'minimum centre distance');
shared = cfg;
shared.pca.region_mode = 'shared';
shared.pca.x_bounds_lambda = [];
shared.pca.z_bounds_lambda = [];
shared = prepareFIMConfig(shared);
sharedModel = fimGeometry('prepare',shared);
bad = sharedModel.zero;
bad(1) = -0.1;
pc_test_assertRejectedGeometry(sharedModel,bad,shared,'shared aperture clips outward left-edge move');
inside = sharedModel.zero;
inside(1) = 0.1;
[~,valid] = fimGeometry('build',sharedModel,inside,shared);
assert(valid,'Shared aperture must admit an interior element position.');
badCfg = cfg;
badCfg.pca.x_grid_lambda = [0 NaN];
pc_test_assertRejectedConfig(badCfg,'nonfinite PCA grid');
badCfg = cfg;
badCfg.pca.sweeps = -1;
pc_test_assertRejectedConfig(badCfg,'negative PCA sweeps');
badCfg = cfg;
badCfg.pca.region_mode = 'unrecognized';
pc_test_assertRejectedConfig(badCfg,'unrecognized PCA region mode');
checks{end+1} = 'Local bounds, shared aperture, all-pair spacing and invalid configuration rejection';

deployed = cfg;
deployed.bs.position_xy_m = [13 -7];
deployed.bs.height_m = 19;
deployed.bs.azimuth_deg = 23;
deployed.bs.downtilt_deg = 11;
deployed = prepareFIMConfig(deployed);
deployedModel = fimGeometry('prepare',deployed);
[deployedGeom,valid] = fimGeometry('build',deployedModel,p,deployed);
assert(valid,'Deployed PCA fixture must remain feasible.');
scene = fimScene('prepare',deployed);
expected = bsxfun(@plus,moved.positions_m*scene.bs.R_global_from_array.',scene.bs.position_m);
pc_test_assertClose(deployedGeom.positions_global_m,expected,1e-12,'Deployment applied exactly once');
delta = deployedGeom.positions_global_m-deployedGeom.flat_positions_global_m;
pc_test_assertClose(delta*scene.bs.R_global_from_array(:,2),zeros(cfg.N,1),1e-12,'Global move remains parallel to installed array plane');
for k = 1:cfg.N
    pc_test_assertClose(deployedGeom.R_global(:,:,k),scene.bs.R_global_from_array,1e-12,'Global element orientation remains fixed');
end
checks{end+1} = 'Installed-array plane and fixed element orientation under BS rotation and translation';

pattern = pc_test_fixturePattern();
realization = pc_test_fixtureRealization(cfg.K);
states = struct('yaw',0,'tilt',0,'zero',1,'count',1);
legacy = planarReference(realization,pattern,cfg);
legacy = reshape(legacy,cfg.N,cfg.K,cfg.ofdm.num_tones);
[Hflat,info] = fimChannel('geometry',realization,pattern,flat,cfg);
assert(info.cp_feasible,'Fixture must have sufficient CP.');
pc_test_assertClose(Hflat,legacy,1e-11,'PCA flat channel equals old fixed UPA');
oneUE = cfg;
oneUE.K = 1;
oneUE = prepareFIMConfig(oneUE);
u = [0.3;sqrt(0.75);0.4];
path = struct('U',u,'alpha',0.7+0.2i,'tau',50e-9);
singleRealization = struct('paths',{{path}});
[Hs,~] = fimChannel('geometry',singleRealization,pattern,flat,oneUE);
[Hm,~] = fimChannel('geometry',singleRealization,pattern,moved,oneUE);
for f = 1:oneUE.ofdm.num_tones
    frequency = oneUE.fc_Hz+oneUE.ofdm.freq_offsets_Hz(f);
    phase = exp(1i*2*pi*frequency*((moved.positions_m-flat.positions_m)*u)/oneUE.c0);
    pc_test_assertClose(Hm(:,1,f),Hs(:,1,f).*phase,1e-11,'Single path wideband phase uses actual displacement once');
end
pc_test_assertClose(abs(Hm),abs(Hs),1e-12,'Fixed element pattern preserves oneUE-path magnitudes');
for receiver = {'zf','mmse-irc'}
    oneUE.receiver = receiver{1};
    a = uplinkOFDMRate(Hs,3,oneUE);
    b = uplinkOFDMRate(Hm,3,oneUE);
    pc_test_assertClose(a.rate,b.rate,1e-10,'Single-user oneUE-path phase motion cannot manufacture array gain');
end
checks{end+1} = 'Old fixed-channel identity, exact wideband displacement phase and no artificial oneUE-path gain';

% A scalar two-path channel has constructive positions inside the grid.
known = cfg;
known.Nh = 2;
known.Nv = 2;
known.dv_lambda = 1;
known.K = 1;
known.pca.x_grid_lambda = [-0.25 0 0.25];
known.pca.z_grid_lambda = 0;
known.receiver = 'mmse-irc';
known.opt.score_mode = 'exact';
known = prepareFIMConfig(known);
knownModel = fimGeometry('prepare',known);
path = struct('U',[0.5 -0.5;sqrt(0.75) sqrt(0.75);0 0],...
    'alpha',[1;1],'tau',[0;0]);
knownRealization = struct('paths',{{path}});
work = fimDeformation('prepare',knownRealization,pattern,knownModel,known);
reference = uplinkOFDMRate(work.Hfixed,0,known);
[~,chosen,best] = fimDeformation('optimize',work,0,known);
assert(any(chosen.parameters~=0),'Known multipath fixture requires a nonzero selected displacement.');
assert(best.sumRate>reference.sumRate+0.01,'PCA optimizer must recover known constructive multipath gain.');
checks{end+1} = 'Optimization moves phase centres and recovers a known multipath gain';

for receiver = {'zf','mmse-irc'}
    for mode = {'sum_rate','max_min_rate'}
        exact = cfg;
        exact.receiver = receiver{1};
        exact.opt_mode = mode{1};
        exact.opt.score_mode = 'exact';
        work = fimDeformation('prepare',realization,pattern,model,exact);
        reference = uplinkOFDMRate(work.Hfixed,2,exact);
        [Hselected,selected,metric,trace] = fimDeformation('optimize',work,2,exact);
        [rebuilt,selectedInfo] = fimChannel('geometry',realization,pattern,selected,exact);
        assert(selected.feasible && selectedInfo.cp_feasible,'Selected PCA geometry must be feasible.');
        pc_test_assertClose(Hselected,rebuilt,1e-11,'One PCA geometry constructs every UE and tone channel');
        pc_test_assertClose(metric.rate,pc_test_uplinkRate(rebuilt,2,exact),1e-11,'Selected metric uses actual configured receiver');
        [referenceParameters,referenceMetric] = pc_test_fullChannelCoordinateReference(realization,pattern,model,2,exact);
        pc_test_assertClose(selected.parameters,referenceParameters,0,'Row-update search matches independent full-channel coordinate enumeration');
        pc_test_assertClose(metric.rate,referenceMetric.rate,1e-11,'Row-update search rate equals independent full-channel enumeration');
        assert(pc_test_objective(metric,mode{1})>=pc_test_objective(reference,mode{1})-1e-9,'PCA search lost fixed baseline.');
        assert(all(diff(trace.values)>=-1e-9),'Accepted PCA updates must be monotonic.');
        accelerated = exact;
        accelerated.opt.score_mode = 'gram_verified';
        work = fimDeformation('prepare',realization,pattern,model,accelerated);
        [~,~,other] = fimDeformation('optimize',work,2,accelerated);
        pc_test_assertClose(pc_test_objective(metric,mode{1}),pc_test_objective(other,mode{1}),1e-9,'Exact and verified PCA scoring agree');
    end
end
checks{end+1} = 'ZF/MMSE-IRC and both objectives: actual-channel rebuild, monotone search and baseline retention';

colored = cfg;
colored.receiver = 'mmse-irc';
colored.irc.external_covariance = 0.2*ones(cfg.N)+diag((1:cfg.N)/20);
colored.ofdm.resource_mode = 'rb';
colored.ofdm.num_rb = 3;
colored.ofdm.rb_allocation = {1,2:3};
colored.power.snr_reference = 'total_over_tone_noise';
colored = prepareFIMConfig(colored);
coloredModel = fimGeometry('prepare',colored);
work = fimDeformation('prepare',realization,pattern,coloredModel,colored);
[Hselected,~,metric] = fimDeformation('optimize',work,3,colored);
pc_test_assertClose(metric.rate,pc_test_uplinkRate(Hselected,3,colored),1e-11,'Correlated interference and unequal RB full evaluation');
pc_test_assertClose(metric.userNumActive,[12;24],0,'Unequal active resources retained');
assert(all(metric.sinrPerTone(~colored.ofdm.user_sample_mask)==0),'Unscheduled PCA tones must be zero.');
pc_test_assertClose(metric.txPowerTotal,10^(3/10)*colored.noiseVar*ones(cfg.K,1),1e-11,'PCA motion does not change UE transmit power');
checks{end+1} = 'Correlated external interference, unequal RB allocation and per-UE total power';

% Each individual row has zero delay spread, but the whole array can exceed CP.
cp = cfg;
cp.K = 1;
cp.receiver = 'mmse-irc';
cp.pca.z_grid_lambda = 0;
cp = prepareFIMConfig(cp);
cpModel = fimGeometry('prepare',cp);
[cpFlat,~] = fimGeometry('build',cpModel,cpModel.zero,cp);
cpPath = struct('U',[1;0;0],'alpha',1,'tau',0);
cpRealization = struct('paths',{{cpPath}});
[~,cpInfo] = fimChannel('geometry',cpRealization,pattern,cpFlat,cp);
cp.ofdm.cp_duration_s = cpInfo.delay_span_s;
outward = cpModel.zero;
outward(1) = -0.3;
[outGeom,valid] = fimGeometry('build',cpModel,outward,cp);
assert(valid,'Outward CP fixture must satisfy the motion and spacing limits.');
[~,outInfo] = fimChannel('geometry',cpRealization,pattern,outGeom,cp);
assert(~outInfo.cp_feasible,'Common array delay span must reject an outward edge move.');
work = fimDeformation('prepare',cpRealization,pattern,cpModel,cp);
[~,~,~,trace] = fimDeformation('optimize',work,0,cp);
assert(trace.cp_rejections>0,'Row-based PCA search must check the whole-array common CP, not one row only.');
checks{end+1} = 'Common array CP span rejects moves even when every individual row has zero spread';

integration = pc_test_runPairedIntegration();
checks{end+1} = 'Three-way zero-grid identity, common realizations, direct channel rebuild and MAT save/reload';
report = struct('passed',true,'checks',{checks},'integration',integration,...
    'test_fixture_only',true,'measured_pattern_tested',false,...
    'actual_julia_cdl_tested',false,'graphics_tested',false);
fprintf('PASS: %d PCA numerical check groups; temporary synthetic inputs removed on return.\n',numel(checks));
for i = 1:numel(checks)
    fprintf('  %d. %s\n',i,checks{i});
end
end

function cfg = pc_test_fixtureConfig()
cfg = configFIMSimulation();
cfg.receiver='zf'; cfg.noiseVar=1; cfg.power.snr_reference='per_subcarrier';
cfg.comparison.include_pca=false; cfg.subarray.m=1; cfg.subarray.n=1;
cfg.orientation.yaw_grid_deg=0; cfg.orientation.tilt_grid_deg=0;
cfg.Nh = 4; cfg.Nv = 2; cfg.K = 2;
cfg.dh_lambda = 0.5; cfg.dv_lambda = 0.5;
cfg.users.mode = 'angles';
cfg.bs.position_xy_m = [0 0];
cfg.bs.ground_height_m = 0; cfg.bs.height_m = 0;
cfg.bs.azimuth_deg = 0; cfg.bs.downtilt_deg = 0;
cfg.pathloss.enabled = false;
cfg.noiseVar = 1;
cfg.MC = 1;
cfg.fast_mode = false;
cfg.fim.mode = 'deformation';
cfg.deform.family = 'pca'; cfg.comparison.include_pca = false;
cfg.pca.x_grid_lambda = [-0.3 0 0.3];
cfg.pca.z_grid_lambda = [-0.3 0 0.3];
cfg.pca.min_spacing_lambda = 0.05;
cfg.pca.sweeps = 1;
cfg.pca.region_mode = 'per_element';
cfg.pca.x_bounds_lambda = [];
cfg.pca.z_bounds_lambda = [];
cfg.ofdm.resource_mode = 'indices';
cfg.ofdm.active_indices = [-1 0 1];
cfg.ofdm.rb_allocation = [];
cfg.ofdm.frequency_mode = 'full';
cfg.ofdm.include_cp_overhead = false;
cfg.save_results = false;
cfg.save_figures = false;
cfg.plot.enabled = false;
cfg = prepareFIMConfig(cfg);
end

function pattern = pc_test_fixturePattern()
[P,T] = meshgrid(0:360,0:180);
G = 1+0.3*sind(T).*cosd(P)+0.2*cosd(T).^2;
pattern = struct('phi',0:360,'theta',0:180,'G',G,...
    'native_to_body',eye(3),'front_only',false);
end

function realization = pc_test_fixtureRealization(K)
realization.paths = cell(K,1);
for k = 1:K
    U = [0.15+0.2*k,-0.4;0.9,0.8;0.1,0.15*k];
    U = bsxfun(@rdivide,U,sqrt(sum(U.^2,1)));
    realization.paths{k} = struct('U',U,'alpha',[0.8+0.1i*k;-0.2+0.4i],'tau',[0;40e-9]);
end
end

function report = pc_test_runPairedIntegration()
testDir = tempname;
mkdir(testDir);
assert(exist(testDir,'dir')==7,'Cannot create PCA fixture directory.');
cleanup = onCleanup(@() pc_test_removeTestDirectory(testDir));
patternPath = fullfile(testDir,'TEST_PCA_pattern.txt');
writeTestPattern(patternPath);
rawCDL.name = 'TEST_ONLY_PCA_fixture_NOT_standard_CDL_C';
rawCDL.pow = [0.55 0.30 0.15];
rawCDL.taus = [0 40 120]*1e-9;
rawCDL.AoA = [70 95 110]*pi/180;
rawCDL.EoA = [-3 0 5]*pi/180;
rawCDL.AoD = rawCDL.AoA;
rawCDL.EoD = rawCDL.EoA;
rawCDL.numClusters = 3;
rawCDL.DS = sqrt(sum(rawCDL.pow.*(rawCDL.taus-sum(rawCDL.pow.*rawCDL.taus)).^2));
rawCDL.KF = -Inf; rawCDL.XPR = 0; rawCDL.test_fixture_only = true;
profilePath = fullfile(testDir,'TEST_PCA_rawCDL.mat');
save(profilePath,'rawCDL');
cfg = pc_test_fixtureConfig();
cfg.pattern.file = patternPath;
cfg.cdl.profile_mat = profilePath;
cfg.receiver = 'mmse-irc';
cfg.MC = 2;
cfg.seed = 20261031;
cfg.snr_dB_list = [-2 4];
cfg.output.save_realizations = true;
cfg.output.save_geometry = true;
cfg.save_results = true;
cfg.output_dir = fullfile(testDir,'zero');
cfg.subarray.m = 2; cfg.subarray.n = 1;
cfg.subarray.displacement_grid_lambda = 0;
cfg.subarray.yaw_grid_deg = 0; cfg.subarray.tilt_grid_deg = 0;
cfg.subarray.min_spacing_lambda = 0.05;
cfg.subarray.sweeps = 1;
cfg.pca.x_grid_lambda = 0; cfg.pca.z_grid_lambda = 0;
[zero,~,files] = simulateFIM(cfg,'sumrate');
pc_test_assertClose(zero.fixed.rate,zero.fim.rate,0,'Subarray zero grid exact rate identity');
pc_test_assertClose(zero.fixed.rate,zero.pca.rate,0,'PCA zero grid exact rate identity');
assert(all(zero.pca_parameters(:)==0),'Zero PCA grid must retain zero positions.');
cfg.output_dir = fullfile(testDir,'moving');
cfg.subarray.displacement_grid_lambda = [-0.2 0 0.2];
cfg.subarray.yaw_grid_deg = [-10 0 10];
cfg.pca.x_grid_lambda = [-0.2 0 0.2];
cfg.pca.z_grid_lambda = [-0.2 0 0.2];
[moving,~,movingFiles] = simulateFIM(cfg,'sumrate');
pc_test_assertClose(moving.fixed.rate,zero.fixed.rate,0,'Search-space changes cannot alter paired fixed baseline');
for mc = 1:cfg.MC
    for k = 1:cfg.K
        a = moving.realizations{mc}.paths{k};
        b = zero.realizations{mc}.paths{k};
        pc_test_assertClose(a.U,b.U,0,'Shared path directions independent of search space');
        pc_test_assertClose(a.alpha,b.alpha,0,'Shared fading independent of search space');
        pc_test_assertClose(a.tau,b.tau,0,'Shared excess delays independent of search space');
    end
    for s = 1:numel(cfg.snr_dB_list)
        geometries = {moving.geometry.fixed,moving.geometry.fim{mc,s},moving.geometry.pca{mc,s}};
        names = {'fixed','fim','pca'};
        for j = 1:3
            [H,info] = fimChannel('geometry',moving.realizations{mc},moving.pattern,geometries{j},moving.cfg);
            assert(info.cp_feasible,'Every compared geometry must satisfy physical CP.');
            metric = uplinkOFDMRate(H,moving.snr_dB(s),moving.cfg);
            stats = moving.(names{j});
            pc_test_assertClose(stats.rate(mc,:,s),metric.rate.',1e-10,'All three rates reconstruct from the same saved realization');
        end
        assert(moving.pca.objective(mc,s)>=moving.fixed.objective(mc,s)-1e-9,'PCA result must retain baseline pc_test_objective.');
        assert(moving.fim.objective(mc,s)>=moving.fixed.objective(mc,s)-1e-9,'Subarray result must retain baseline pc_test_objective.');
    end
end
stored = pc_test_loadResultMat(movingFiles);
pc_test_assertClose(stored.pca.rate,moving.pca.rate,0,'Saved PCA statistics reload exactly');
pc_test_assertClose(stored.pca_parameters,moving.pca_parameters,0,'Saved PCA controls reload exactly');
assert(isstruct(pc_test_loadResultMat(files)),'Zero-grid result must also save correctly.');
assert(moving.pattern.is_test_pattern && moving.rawCDL.test_fixture_only,'Synthetic source tags must survive.');
report = struct('passed',true,'MC',cfg.MC,'SNR_points',numel(cfg.snr_dB_list),...
    'temporary_files_removed_on_return',true);
end

function result = pc_test_loadResultMat(files)
result = [];
for i = 1:numel(files)
    [~,~,extension] = fileparts(files{i});
    if strcmpi(extension,'.mat')
        load(files{i},'result');
        break
    end
end
assert(isstruct(result),'Expected saved result MAT.');
end

function pc_test_removeTestDirectory(path)
if exist(path,'dir')==7
    rmdir(path,'s');
end
end

function rate = pc_test_uplinkRate(H,snr,cfg)
metric = uplinkOFDMRate(H,snr,cfg);
rate = metric.rate;
end

function [p,metric] = pc_test_fullChannelCoordinateReference(realization,pattern,model,snr,cfg)
% Deliberately rebuild every candidate in full: independent of PCA row caching.
p = model.zero;
[geom,valid] = fimGeometry('build',model,p,cfg);
assert(valid,'Reference zero state must be feasible.');
[H,info] = fimChannel('geometry',realization,pattern,geom,cfg);
assert(info.cp_feasible,'Reference zero state must satisfy CP.');
metric = uplinkOFDMRate(H,snr,cfg);
value = pc_test_objective(metric,cfg.opt_mode);
for sweep = 1:cfg.pca.sweeps
    changed = false;
    for d = 1:model.parameter_count
        best = value;
        selected = p;
        selectedMetric = metric;
        for j = 1:numel(model.grids{d})
            candidate = p;
            candidate(d) = model.grids{d}(j);
            if isequal(candidate,p)
                continue
            end
            [geom,valid] = fimGeometry('build',model,candidate,cfg);
            if ~valid
                continue
            end
            [H,info] = fimChannel('geometry',realization,pattern,geom,cfg);
            if ~info.cp_feasible
                continue
            end
            trial = uplinkOFDMRate(H,snr,cfg);
            score = pc_test_objective(trial,cfg.opt_mode);
            if score>best+cfg.objective_tol
                best = score;
                selected = candidate;
                selectedMetric = trial;
            end
        end
        if best>value+cfg.objective_tol
            value = best;
            p = selected;
            metric = selectedMetric;
            changed = true;
        end
    end
    if ~changed
        break
    end
end
end

function value = pc_test_objective(metric,mode)
if strcmp(mode,'sum_rate')
    value = metric.sumRate;
else
    value = metric.worstRate;
end
end

function pc_test_assertRejectedGeometry(model,p,cfg,label)
rejected = false;
try
    [~,valid] = fimGeometry('build',model,p,cfg);
    rejected = ~valid;
catch
    rejected = true;
end
assert(rejected,['Expected geometry rejection: ' label]);
end

function pc_test_assertRejectedConfig(cfg,label)
rejected = false;
try
    prepareFIMConfig(cfg);
catch
    rejected = true;
end
assert(rejected,['Expected config rejection: ' label]);
end

function pc_test_assertClose(actual,expected,tolerance,label)
assert(isequal(size(actual),size(expected)),['Shape mismatch: ' label]);
assert(all(isfinite(actual(:))) && all(isfinite(expected(:))),['Nonfinite comparison: ' label]);
if isempty(actual)
    return
end
errorSize = max(abs(actual(:)-expected(:)));
scale = max(1,max(abs(expected(:))));
assert(errorSize<=tolerance*scale,['Failed: ' label]);
end

function H=planarReference(realization,pattern,cfg)
% Direct formula independent of the consolidated geometry-channel dispatcher.
[X,Z]=ndgrid(((0:cfg.Nh-1)-(cfg.Nh-1)/2)*cfg.dh_lambda,((0:cfg.Nv-1)-(cfg.Nv-1)/2)*cfg.dv_lambda);
pos=[X(:),zeros(cfg.N,1),Z(:)]*(cfg.c0/cfg.fc_Hz);
H=complex(zeros(cfg.N,cfg.K,cfg.ofdm.num_tones));
for k=1:cfg.K
    p=realization.paths{k}; native=pattern.native_to_body'*p.U;
    phi=mod(atan2d(native(2,:),native(1,:)),360); theta=acosd(max(-1,min(1,native(3,:))));
    gain=interp2(pattern.phi,pattern.theta,pattern.G,phi,theta,'linear');
    if pattern.front_only, gain(p.U(2,:)<-1e-12)=0; end
    for f=1:cfg.ofdm.num_tones
        nu=cfg.ofdm.freq_offsets_Hz(f); freq=cfg.fc_Hz;
        if cfg.ofdm.frequency_dependent_steering, freq=freq+nu; end
        H(:,k,f)=exp(1i*2*pi*freq*(pos*p.U/cfg.c0))*(sqrt(gain(:)).*p.alpha(:).*exp(-1i*2*pi*nu*p.tau(:)));
    end
end
end

function writeTestPattern(path)
fid=fopen(path,'w'); assert(fid~=-1,'Cannot write temporary test pattern.');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'GAIN 0 dBi ; SYNTHETIC TEST ONLY\n');
fprintf(fid,'PATTERN Phi Theta Gain Phase ; TEST ONLY\n');
phi=(0:359).'; az=mod(phi-90+180,360)-180;
for theta=0:180
    el=90-theta; gain=-min(12*(az/120).^2+12*(el/12)^2,30);
    fprintf(fid,'%d %d %.9f 0\n',[phi,theta*ones(360,1),gain].');
end
end

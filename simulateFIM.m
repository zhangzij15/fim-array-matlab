function [result,figs,files] = simulateFIM(cfg,experiment)
%SIMULATEFIM 单个 K 的固定 / 子阵 / PCA 严格配对仿真。
% cdf 使用 cfg.snr_dB；sumrate / worstuser 扫描 cfg.snr_dB_list。
if nargin<1 || isempty(cfg), cfg = configFIMSimulation(); end
if nargin<2, experiment = 'cdf'; end
cfg.fim.mode = 'deformation'; cfg.deform.family = 'subarray';
cfg.comparison.include_pca = true;
switch experiment
    case 'cdf'
        snrs = cfg.snr_dB;
    case 'sumrate'
        cfg.opt_mode = 'sum_rate'; snrs = cfg.snr_dB_list;
    case 'worstuser'
        cfg.opt_mode = 'max_min_rate'; snrs = cfg.snr_dB_list;
    otherwise
        error('FIM:Experiment','Use cdf, sumrate or worstuser.');
end
cfg = prepareFIMConfig(cfg);
M = cfg.MC; K = cfg.K; F = cfg.ofdm.num_tones; S = numel(snrs);
model = fimGeometry('prepare',cfg);
pcaCfg = cfg; pcaCfg.deform.family = 'pca'; pcaCfg.comparison.include_pca = false;
pcaModel = fimGeometry('prepare',pcaCfg);
[pattern,profile,rawCDL] = fimChannel('prepare',cfg);
result.schema_version = 7;
result.cfg = cfg; result.experiment = experiment; result.snr_dB = snrs(:).';
result.created = datestr(now,30);
result.model = ['Paired fixed/PCA/pattern/subarray; ' profile.model '; perfect CSI; Shannon SE'];
result.pattern = pattern; result.profile = profile; result.rawCDL = rawCDL;
result.deformation_model = model; result.pca_model = pcaModel;
result.bs_deployment = fimScene('prepare',cfg); result.scene_cell = cell(M,1);
result.realizations = {};
if cfg.output.save_realizations, result.realizations = cell(M,1); end
result.comparison.labels = {'Fixed planar array','PCA','Pattern reconfigurable','Reconfigurable subarrays'};
result.comparison.common_realization = true;
result.comparison.note = 'Same UE, paths, fading, RBs, power and receiver; independent local searches.';
result.resource_allocation.user_active_mask = cfg.ofdm.user_active_mask;
result.resource_allocation.user_sample_mask = cfg.ofdm.user_sample_mask;
result.resource_allocation.user_num_active = cfg.ofdm.user_num_active;
for label = {'fixed','pca','orientation','fim'}
    result.(label{1}) = allocateStatistics(M,K,F,S,cfg);
end
[Y,T]=ndgrid(unique([0 cfg.orientation.yaw_grid_deg(:).']),unique([0 cfg.orientation.tilt_grid_deg(:).']));
states.yaw=Y(:); states.tilt=T(:); states.count=numel(Y);
states.zero=find(states.yaw==0 & states.tilt==0);
result.orientation_states=states;
result.orientation_selection=zeros(M,cfg.N,S,'uint32');
result.orientation_evaluations=zeros(M,S);
result.orientation_exact_evaluations=zeros(M,S);
result.deformation_parameters = zeros(M,model.parameter_count,S);
result.pca_parameters = zeros(M,pcaModel.parameter_count,S);
[flat,valid,reason] = fimGeometry('build',model,model.zero,cfg);
assert(valid,['Infeasible fixed reference: ' reason]);
[pcaFlat,valid,reason] = fimGeometry('build',pcaModel,pcaModel.zero,pcaCfg);
assert(valid,['Infeasible PCA reference: ' reason]);
assert(isequal(flat.positions_m,pcaFlat.positions_m) && isequal(flat.R,pcaFlat.R),'The planar references must match exactly.');
result.geometry.fixed = flat; result.geometry.pca_fixed = pcaFlat;
result.geometry.fim = {}; result.geometry.pca = {};
if cfg.output.save_geometry
    result.geometry.fim = cell(M,S); result.geometry.pca = cell(M,S);
end
counterFields = {'evaluations','exact_evaluations','gram_scores','constraint_rejections','cp_rejections',...
    'pca_evaluations','pca_exact_evaluations','pca_gram_scores','pca_constraint_rejections','pca_cp_rejections'};
for j = 1:numel(counterFields), result.(counterFields{j}) = zeros(M,S); end
result.requiredCP_s = zeros(M,1);
result.fimRequiredCP_s = zeros(M,S); result.pcaRequiredCP_s = zeros(M,S);
previous = rng; cleanup = onCleanup(@() rng(previous));
rng(cfg.seed,'twister'); started = tic;
fprintf('Four-array comparison: N=%d, K=%d, MC=%d, tones=%d, SNR points=%d.\n',cfg.N,K,M,F,S);
for mc = 1:M
    realization = fimChannel('realization',cfg,profile);
    result.scene_cell{mc} = realization.userScene;
    if cfg.output.save_realizations, result.realizations{mc} = realization; end
    work = fimDeformation('prepare',realization,pattern,model,cfg);
    pcaWork = fimDeformation('prepare',realization,pattern,pcaModel,pcaCfg);
    assert(isequal(work.Hfixed,pcaWork.Hfixed),'Baseline channel mismatch.');
    result.requiredCP_s(mc) = work.fixed_channel_info.delay_span_s;
    patternCandidates=fimChannel('candidates',realization,pattern,states,cfg);
    patternFixed=reshape(patternCandidates(:,:,states.zero,:),cfg.N,cfg.K,F);
    assert(norm(patternFixed(:)-work.Hfixed(:))<=1e-10*max(1,norm(work.Hfixed(:))),'Pattern and position baselines disagree.');
    for s = 1:S
        fixed = uplinkOFDMRate(work.Hfixed,snrs(s),cfg);
        [~,geom,sub,trace,work] = fimDeformation('optimize',work,snrs(s),cfg);
        [~,pgeom,pca,ptrace,pcaWork] = fimDeformation('optimize',pcaWork,snrs(s),pcaCfg);
        [~,selection,orientation,orientationTrace]=fimDeformation('orientation',patternCandidates,states,snrs(s),cfg);
        result.orientation_selection(mc,:,s)=uint32(selection.');
        result.orientation_evaluations(mc,s)=orientationTrace.evaluations;
        result.orientation_exact_evaluations(mc,s)=orientationTrace.exact_evaluations;
        result.orientation=recordMetric(result.orientation,orientation,mc,s,cfg);
        result.fixed = recordMetric(result.fixed,fixed,mc,s,cfg);
        result.fim = recordMetric(result.fim,sub,mc,s,cfg);
        result.pca = recordMetric(result.pca,pca,mc,s,cfg);
        result.deformation_parameters(mc,:,s) = geom.parameters;
        result.pca_parameters(mc,:,s) = pgeom.parameters;
        if cfg.output.save_geometry
            result.geometry.fim{mc,s} = geom; result.geometry.pca{mc,s} = pgeom;
        end
        result.fimRequiredCP_s(mc,s) = work.last_channel_info.delay_span_s;
        result.pcaRequiredCP_s(mc,s) = pcaWork.last_channel_info.delay_span_s;
        for name = {'evaluations','exact_evaluations','gram_scores','constraint_rejections','cp_rejections'}
            f = name{1}; result.(f)(mc,s) = trace.(f);
            result.(['pca_' f])(mc,s) = ptrace.(f);
        end
    end
    if mod(mc,cfg.progress_every)==0 || mc==M
        fprintf('  MC %d/%d complete, elapsed %.1f s.\n',mc,M,toc(started));
    end
end
result.elapsed_s = toc(started);
result.gain.sumRate = result.fim.sumRate-result.fixed.sumRate;
result.gain.worstRate = result.fim.worstRate-result.fixed.worstRate;
result.gain.objective = result.fim.objective-result.fixed.objective;
result.gain.orientation_sumRate=result.orientation.sumRate-result.fixed.sumRate;
result.gain.orientation_objective=result.orientation.objective-result.fixed.objective;
assert(all(result.gain.orientation_objective(:)>=-cfg.check.objective_tol),'Pattern optimization lost fixed baseline.');
result.gain.pca_sumRate = result.pca.sumRate-result.fixed.sumRate;
result.gain.pca_worstRate = result.pca.worstRate-result.fixed.worstRate;
result.gain.pca_objective = result.pca.objective-result.fixed.objective;
result.gain.subarray_minus_pca_sumRate = result.fim.sumRate-result.pca.sumRate;
result.gain.subarray_minus_pca_objective = result.fim.objective-result.pca.objective;
assert(all(result.gain.objective(:)>=-cfg.check.objective_tol),'Subarray lost fixed baseline.');
assert(all(result.gain.pca_objective(:)>=-cfg.check.objective_tol),'PCA lost fixed baseline.');
figs = {}; files = {};
if cfg.save_results
    if exist(cfg.output_dir,'dir')~=7, mkdir(cfg.output_dir); end
    stem = fullfile(cfg.output_dir,['comparison_' experiment '_' datestr(now,'yyyymmdd_HHMMSS')]);
    base = stem; n = 1;
    while exist([stem '.mat'],'file')==2, stem = [base '_' num2str(n)]; n=n+1; end
    if strcmp(cfg.output.mat_format,'default'), save([stem '.mat'],'result');
    else, save([stem '.mat'],'result',cfg.output.mat_format); end
    files = {[stem '.mat']};
    reportFIMParameters(result,stem,false,false);
    files = [files {[stem '_parameters.txt'],[stem '_parameters.csv']}];
end
if cfg.plot.enabled
    [figs,graphics] = plotFIMResults(result,cfg.output_dir);
    files = [files graphics.files];
end
end

function stats = allocateStatistics(M,K,F,S,cfg)
for name = {'sinr','snr','rate','rateNoCP','throughput_bps','rateScheduledBand','rateScheduledBandNoCP'}
    stats.(name{1}) = zeros(M,K,S);
end
for name = {'sumRate','worstRate','worstSINR','sumThroughput_bps','worstThroughput_bps','rank','rankDeficient','leakage','objective'}
    stats.(name{1}) = zeros(M,S);
end
stats.meanRatePerToneNoCP = zeros(K,F,S); stats.meanSINRPerTone = zeros(K,F,S);
if cfg.output.save_per_tone, stats.sinrPerTone = zeros(M,K,F,S); end
end

function stats = recordMetric(stats,metric,mc,s,cfg)
for name = {'sinr','snr','rate','rateNoCP','throughput_bps','rateScheduledBand','rateScheduledBandNoCP'}
    f = name{1}; stats.(f)(mc,:,s) = metric.(f).';
end
for name = {'sumRate','worstRate','worstSINR','sumThroughput_bps','worstThroughput_bps'}
    f = name{1}; stats.(f)(mc,s) = metric.(f);
end
used = metric.receiver.num_scheduled>0;
stats.rank(mc,s) = min(metric.receiver.rank(used));
stats.rankDeficient(mc,s) = any(metric.receiver.rank_deficient);
stats.leakage(mc,s) = max(metric.receiver.leakage_relative);
if strcmp(cfg.opt_mode,'sum_rate'), stats.objective(mc,s) = metric.sumRate;
else, stats.objective(mc,s) = metric.worstRate; end
stats.meanRatePerToneNoCP(:,:,s) = stats.meanRatePerToneNoCP(:,:,s)+metric.ratePerToneNoCP/cfg.MC;
stats.meanSINRPerTone(:,:,s) = stats.meanSINRPerTone(:,:,s)+metric.sinrPerTone/cfg.MC;
if cfg.output.save_per_tone
    stats.sinrPerTone(mc,:,:,s) = reshape(metric.sinrPerTone,1,cfg.K,cfg.ofdm.num_tones);
end
end

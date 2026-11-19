function [study,figs,files] = runFIMStudy(cfg,options)
%RUNFIMSTUDY 同一CDL模板、相同场景分布的K扫描；分批保存并可续算。
% 各K内四阵列严格配对。不同K使用同分布独立样本，不声称逐UE跨K配对。
% 不筛选样本，不平滑/缩放CDF，不更改优化器，不保证预设增益排序。
if nargin<1 || isempty(cfg), cfg = configFIMSimulation(); end
if nargin<2 || isempty(options), options = cfg.study; end
presentation = cfg;
cfg.study.resume_directory='';
options = defaults(options);
assert(strcmp(cfg.users.mode,'random'),'K-scan requires users.mode=random; use simulateFIM for specified users.');
assert(isempty(cfg.ofdm.rb_allocation),'Use shared RBs for this K comparison; leave ofdm.rb_allocation empty.');
assert(~cfg.fast_mode,'Disable fast_mode so the requested MC count is not silently reduced.');
cfg.K = options.K_values(1);
cfg.fim.mode = 'deformation'; cfg.deform.family = 'subarray';
cfg.comparison.include_pca = true;
cfg.MC = 1;
cfg.save_results = false; cfg.save_figures = false; cfg.plot.enabled = false;
cfg.output.save_per_tone = false;
cfg.output.save_realizations = false;
cfg.output.save_geometry = false;
cfg = prepareFIMConfig(cfg);
assert(all(options.K_values<=cfg.N),'Every K must be <= the number of receive elements.');
cfgTemplate = cfg;
figs = {}; files = {};
if isempty(options.resume_directory)
    runDir = newDirectory(cfg.output_dir);
    fprintf('New study directory: %s\n',runDir);
    fprintf('Generating/loading one CDL template for all K cases and batches.\n');
    [pattern,~,rawCDL] = fimChannel('prepare',cfg);
    patternChecked = true;
    plan.version = 3;
    plan.created = datestr(now,30);
    plan.cfg_template = cfgTemplate;
    plan.K_values = options.K_values;
    plan.MC_per_case = options.MC_per_case;
    plan.batch_MC = options.batch_MC;
    plan.pattern = pattern;
    plan.rawCDL = rawCDL;
    plan.note = 'All samples retained. Same fixed/PCA/pattern/subarray paths within K; independent same-distribution cohorts across K. One CDL template; exact configured tones.';
    saveNamed(fullfile(runDir,'plan.mat'),'plan',plan,cfg.output.mat_format);
else
    patternChecked = false;
    runDir = absolutePath(options.resume_directory);
    plan = readNamed(fullfile(runDir,'plan.mat'),'plan');
    if plan.version==3
        plan.cfg_template = restoreLegacyCDLConfig(plan.cfg_template);
        % Older plans keep their saved generator name and resume from the CDL snapshot.
        if strcmp(cfgTemplate.cdl.generator,'generateCDLBase') && ...
                strcmp(plan.cfg_template.cdl.generator,'generateCDL5G')
            cfgTemplate.cdl.generator = plan.cfg_template.cdl.generator;
        end
    end
    assert(plan.version==3 && isequaln(plan.cfg_template,cfgTemplate) && ...
        isequal(plan.K_values,options.K_values) && isequal(plan.MC_per_case,options.MC_per_case) && ...
        plan.batch_MC==options.batch_MC,...
        'Resume configuration differs from saved plan. Keep all simulation settings/MC/batch size unchanged, or start a new directory.');
    fprintf('Resuming study directory: %s\n',runDir);
end
% 该MAT是本任务真实CDL模板快照，不是合成数据或旧仿真后处理输入。
profilePath = fullfile(runDir,'cdl_snapshot.mat');
rawCDL = plan.rawCDL;
saveNamed(profilePath,'rawCDL',rawCDL,cfg.output.mat_format);
study.schema_version = 2;
study.created = plan.created;
study.directory = runDir;
study.K_values = plan.K_values;
study.MC_per_case = plan.MC_per_case;
study.plan = plan;
study.cases = cell(1,numel(plan.K_values));
study.batch_seeds = cell(1,numel(plan.K_values));
study.batch_files = cell(1,numel(plan.K_values));
study.note = plan.note;
fprintf('Requested MC per case: %s; K values: %s.\n',mat2str(plan.MC_per_case),mat2str(plan.K_values));
fprintf('Resume by setting options.resume_directory to the directory above. No sample or seed selection is performed.\n');
for c = 1:numel(plan.K_values)
    K = plan.K_values(c); target = plan.MC_per_case(c);
    batches = ceil(target/plan.batch_MC);
    combined = [];
    study.batch_seeds{c} = zeros(1,batches);
    study.batch_files{c} = cell(1,batches);
    for b = 1:batches
        n = min(plan.batch_MC,target-(b-1)*plan.batch_MC);
        % 唯一确定的case/batch seed；与是否续算无关。不根据结果挑选seed。
        seed = mod(double(cfgTemplate.seed)+1000003*(c-1)+104729*(b-1),2^32);
        batchPath = fullfile(runDir,sprintf('K%d_batch_%04d.mat',K,b));
        expected = struct('plan_created',plan.created,'case_index',c,'batch_index',b,'K',K,'MC',n,'seed',seed);
        batchCfg = cfgTemplate;
        batchCfg.K = K; batchCfg.MC = n; batchCfg.seed = seed;
        batchCfg.cdl.profile_mat = profilePath;
        batchCfg.output_dir = runDir;
        batchCfg = prepareFIMConfig(batchCfg);
        if exist(batchPath,'file')==2
            saved = readNamed(batchPath,'checkpoint');
            assert(isequal(saved.meta,expected),'Checkpoint metadata differs from this study plan.');
            batchResult = saved.result;
            batchResult.cfg = restoreLegacyCDLConfig(batchResult.cfg);
            fprintf('Reusing completed checkpoint: K=%d, batch %d/%d, MC=%d.\n',K,b,batches,n);
        else
            if ~patternChecked
                [currentPattern,~,currentCDL] = fimChannel('prepare',batchCfg);
                assert(isequaln(currentCDL,plan.rawCDL) && samePattern(currentPattern,plan.pattern),...
                    'Input pattern/CDL differs from this study plan; pending batches have not been run.');
                patternChecked = true;
            end
            fprintf('\nRunning K=%d, batch %d/%d, MC=%d, seed=%u.\n',K,b,batches,n,seed);
            started = tic;
            [batchResult,~,~] = simulateFIM(batchCfg,'cdf');
            validateBatch(batchResult,batchCfg,plan);
            checkpoint.meta = expected;
            checkpoint.result = batchResult;
            saveNamed(batchPath,'checkpoint',checkpoint,cfg.output.mat_format);
            fprintf('Saved %s (batch %.1f min).\n',batchPath,toc(started)/60);
        end
        validateBatch(batchResult,batchCfg,plan);
        combined = appendBatch(combined,batchResult);
        study.batch_seeds{c}(b) = seed;
        study.batch_files{c}{b} = batchPath;
        fprintf('K=%d: %d/%d MC complete, %d pooled UE samples.\n',K,combined.cfg.MC,target,combined.cfg.MC*K);
    end
    combined.cfg.plot = presentation.plot;
    combined.cfg.save_figures = presentation.save_figures;
    combined.batch_seeds = study.batch_seeds{c};
    combined.batch_MC = plan.batch_MC;
    combined.seed_note = 'cfg.seed is the first batch seed only; reproduce with batch_seeds and batch_MC, not one uninterrupted RNG stream.';
    combined.comparison.cross_K_note = 'Same configured user/path distribution and one CDL template; not paired UE realizations across K.';
    study.cases{c} = combined;
    printSummary(combined);
end
summaryFile = fullfile(runDir,'summary.csv');
writeSummary(study,summaryFile);
studyFile = fullfile(runDir,'study.mat');
saveNamed(studyFile,'study',study,cfg.output.mat_format);
files = {studyFile,summaryFile,fullfile(runDir,'plan.mat'),profilePath};
fprintf('\nCompleted study saved: %s\n',studyFile);
fprintf('Summary saved: %s\n',summaryFile);
fprintf('CDFs are empirical samples, with no smoothing. The requested SU/MU gain pattern is a hypothesis, not enforced.\n');
if options.make_plots
    [figs,graphics] = plotFIMResults(study,runDir);
    study.graphics = graphics;
    files = [files graphics.files];
    saveNamed(studyFile,'study',study,cfg.output.mat_format);
end
end

function o = defaults(o)
if ~isfield(o,'K_values'), o.K_values = [1 10]; end
if ~isfield(o,'MC_per_case'), o.MC_per_case = [200 200]; end
if ~isfield(o,'batch_MC'), o.batch_MC = 10; end
if ~isfield(o,'resume_directory'), o.resume_directory = ''; end
if ~isfield(o,'make_plots'), o.make_plots = true; end
assert(isnumeric(o.K_values) && isreal(o.K_values) && isvector(o.K_values) && ...
    all(isfinite(o.K_values(:))) && all(o.K_values(:)>=1 & o.K_values(:)==fix(o.K_values(:))) && ...
    numel(unique(o.K_values))==numel(o.K_values),'K_values must be distinct positive integers.');
o.K_values = o.K_values(:).';
assert(isnumeric(o.MC_per_case) && isreal(o.MC_per_case) && ~isempty(o.MC_per_case),'Invalid MC_per_case.');
if isscalar(o.MC_per_case), o.MC_per_case = repmat(o.MC_per_case,size(o.K_values)); end
o.MC_per_case = o.MC_per_case(:).';
assert(numel(o.MC_per_case)==numel(o.K_values) && all(isfinite(o.MC_per_case)) && ...
    all(o.MC_per_case>=1 & o.MC_per_case==fix(o.MC_per_case)),'MC_per_case must contain one positive integer per K.');
assert(isnumeric(o.batch_MC) && isscalar(o.batch_MC) && isfinite(o.batch_MC) && o.batch_MC>=1 && o.batch_MC==fix(o.batch_MC),'Invalid batch_MC.');
assert(ischar(o.resume_directory),'resume_directory must be a character-vector path.');
assert(isscalar(o.make_plots) && (islogical(o.make_plots) || isnumeric(o.make_plots)) && ...
    any(o.make_plots==[0 1]),'make_plots must be true/false.');
end

function yes = samePattern(a,b)
yes = isequaln(a.G,b.G) && isequaln(a.native_to_body,b.native_to_body) && ...
    isequaln(a.phi,b.phi) && isequaln(a.theta,b.theta) && a.front_only==b.front_only;
end

function cfg = restoreLegacyCDLConfig(cfg)
% Old v3 plans and checkpoints predate the ray settings and used cluster CDL-C.
if ~isfield(cfg.cdl,'model')
    defaults = configFIMSimulation();
    names = setdiff(fieldnames(defaults.cdl),fieldnames(cfg.cdl));
    for i = 1:numel(names)
        field = names{i};
        cfg.cdl.(field) = defaults.cdl.(field);
    end
    cfg.cdl.model = 'cluster';
    cfg.cdl.profile = 'CDL-C';
end
end

function validateBatch(r,cfg,plan)
assert(isequaln(r.cfg,cfg),'Checkpoint result configuration differs from the expected batch configuration.');
assert(isequaln(r.rawCDL,plan.rawCDL) && samePattern(r.pattern,plan.pattern),...
    'Checkpoint/current input pattern or CDL differs from the study plan.');
assert(strcmp(r.experiment,'cdf') && isscalar(r.snr_dB) && r.snr_dB==cfg.snr_dB,...
    'Checkpoint experiment or SNR differs from the study plan.');
for scheme = {'fixed','pca','orientation','fim'}
    d = r.(scheme{1});
    assert(isequal(size(d.rate),[cfg.MC cfg.K]) && all(isfinite(d.rate(:))) && all(d.rate(:)>=0),...
        'Invalid checkpoint UE rate samples.');
    assert(isequal(size(d.sumRate),[cfg.MC 1]) && all(isfinite(d.sumRate(:))),...
        'Invalid checkpoint sum-rate samples.');
end
end

function target = appendBatch(target,source)
if isempty(target), target = source; return; end
oldM = target.cfg.MC; newM = source.cfg.MC;
assert(target.cfg.K==source.cfg.K && isequal(target.snr_dB,source.snr_dB),'Cannot mix K or SNR in an empirical CDF.');
fields = {'sinr','snr','rate','rateNoCP','throughput_bps','rateScheduledBand','rateScheduledBandNoCP',...
    'sumRate','worstRate','worstSINR','sumThroughput_bps','worstThroughput_bps','rank','rankDeficient','leakage','objective'};
for label = {'fixed','pca','orientation','fim'}
    name = label{1};
    for i = 1:numel(fields)
        f = fields{i}; target.(name).(f) = cat(1,target.(name).(f),source.(name).(f));
    end
    for meanName = {'meanRatePerToneNoCP','meanSINRPerTone'}
        f = meanName{1};
        target.(name).(f) = (oldM*target.(name).(f)+newM*source.(name).(f))/(oldM+newM);
    end
end
fields = {'orientation_selection','orientation_evaluations','orientation_exact_evaluations','deformation_parameters','pca_parameters','evaluations','exact_evaluations','gram_scores',...
    'requiredCP_s','fimRequiredCP_s','constraint_rejections','cp_rejections','pca_evaluations',...
    'pca_exact_evaluations','pca_gram_scores','pca_constraint_rejections','pca_cp_rejections','pcaRequiredCP_s','scene_cell'};
for i = 1:numel(fields)
    f = fields{i}; target.(f) = cat(1,target.(f),source.(f));
end
gainFields = fieldnames(source.gain);
for i = 1:numel(gainFields)
    f = gainFields{i}; target.gain.(f) = cat(1,target.gain.(f),source.gain.(f));
end
target.cfg.MC = oldM+newM;
target.elapsed_s = target.elapsed_s+source.elapsed_s;
end

function printSummary(r)
fprintf('\nK=%d; MC=%d; pooled UE samples=%d; empirical CDF nominal step=%.6g.\n',r.cfg.K,r.cfg.MC,r.cfg.K*r.cfg.MC,1/(r.cfg.K*r.cfg.MC));
fprintf('Scheme       Mean UE SE     Mean sum SE    Gain over fixed     P5 UE rate(Mbps)\n');
for label = {'fixed','pca','orientation','fim'}
    name = label{1}; data = r.(name);
    fprintf('%-10s  %11.5f    %11.5f    %13.3f %%    %11.5f\n',name,mean(data.rate(:)),mean(data.sumRate(:)),...
        gainPercent(mean(data.sumRate(:)),mean(r.fixed.sumRate(:))),percentile5(data.throughput_bps)/1e6);
end
end

function writeSummary(study,path)
fid = fopen(path,'w'); assert(fid~=-1,'Cannot write summary CSV.');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'K,MC,pooled_UE_samples,scheme,mean_UE_SE_bit_s_Hz,mean_sum_SE_bit_s_Hz,mean_UE_rate_Mbps,P5_UE_SE_bit_s_Hz,P5_UE_rate_Mbps,relative_mean_sum_SE_gain_percent\n');
for c = 1:numel(study.cases)
    r = study.cases{c};
    for label = {'fixed','pca','orientation','fim'}
        name = label{1}; d = r.(name);
        fprintf(fid,'%d,%d,%d,%s,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n',...
            r.cfg.K,r.cfg.MC,r.cfg.K*r.cfg.MC,name,mean(d.rate(:)),mean(d.sumRate(:)),...
            mean(d.throughput_bps(:))/1e6,percentile5(d.rate),percentile5(d.throughput_bps)/1e6,...
            gainPercent(mean(d.sumRate(:)),mean(r.fixed.sumRate(:))));
    end
end
end

function v = percentile5(x)
x = sort(x(:)); v = x(max(1,ceil(0.05*numel(x))));
end

function g = gainPercent(value,baseline)
if baseline==0, g = NaN; else, g = 100*(value/baseline-1); end
end

function dirName = newDirectory(root)
root = absolutePath(root);
if exist(root,'dir')~=7, mkdir(root); end
base = fullfile(root,['fim_comparison_' datestr(now,'yyyymmdd_HHMMSS')]);
dirName = base; i = 1;
while exist(dirName,'dir')==7, dirName = [base '_' num2str(i)]; i = i+1; end
mkdir(dirName); assert(exist(dirName,'dir')==7,'Cannot create study directory.');
end

function p = absolutePath(p)
if isempty(p), p = pwd; end
if isempty(regexp(p,'^([A-Za-z]:[\\/]|[\\/])','once')), p = fullfile(pwd,p); end
end

function saveNamed(path,name,value,format)
% 在本地结果目录先写临时MAT再改名，避免中断把半写文件当成完成批次。
tmp = [path '.pending.mat'];
switch name
    case 'plan', plan = value;
    case 'rawCDL', rawCDL = value;
    case 'checkpoint', checkpoint = value;
    case 'study', study = value;
    otherwise, error('Unexpected saved variable.');
end
if strcmp(format,'default'), save(tmp,name); else, save(tmp,name,format); end
if exist(path,'file') == 2
    backup = [path '.previous'];
    if exist(backup,'file')==2, delete(backup); end
    movefile(path,backup);
    try
        movefile(tmp,path);
        delete(backup);
    catch err
        if exist(path,'file')~=2, movefile(backup,path); end
        rethrow(err);
    end
else
    movefile(tmp,path);
end
assert(exist(path,'file') == 2, ['Cannot finalize MAT file: ' path]);
end

function value = readNamed(path,name)
assert(exist(path,'file')==2,['Missing checkpoint: ' path]);
switch name
    case 'plan'
        load(path,'plan'); value = plan;
    case 'checkpoint'
        load(path,'checkpoint'); value = checkpoint;
    otherwise
        error('Unexpected checkpoint variable.');
end
end

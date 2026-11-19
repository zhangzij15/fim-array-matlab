function varargout = fimDeformation(action,varargin)
%FIMDEFORMATION 刚性子阵/PCA搜索；同一组参数覆盖全部 UE 和频点。
% work = fimDeformation('prepare',realization,pattern,model,cfg)
% [H,geom,metric,trace,work] = fimDeformation('optimize',work,snr_dB,cfg)
% work 中的缓存跨 SNR 复用几何、信道和 Gram；每个 SNR 仍独立从平面开始。
% 两种机构均从固定阵列开始，逐参数做离散坐标上升。
switch action
    case 'prepare'
        varargout{1} = prepareWork(varargin{:});
    case 'orientation'
        [varargout{1},varargout{2},varargout{3},varargout{4}] = orientation_optimizeFIM(varargin{:});
    case 'optimize'
        [varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}] = optimizeWork(varargin{:});
    otherwise
        error('Unknown fimDeformation action: %s',action);
end
end

function work = prepareWork(realization,pattern,model,cfg)
assert(any(strcmp(model.family,{'subarray','pca'})),'Unknown deformation family.');
assert(cfg.deform.cache_limit>=1 && cfg.deform.cache_limit==fix(cfg.deform.cache_limit),'cache_limit must be a positive integer.');
work.realization = realization;
work.pattern = pattern;
work.model = model;
work.cache = {};
work.cache_limit = cfg.deform.cache_limit;
work.cache_hits = 0;
work.cache_misses = 0;
work.cache_evictions = 0;
[flat,work] = getEntry(work,model.zero,cfg);
assert(flat.valid,['The fixed planar baseline is infeasible: ' flat.reason]);
work.Hfixed = flat.H;
work.fixed_geometry = flat.geometry;
work.fixed_channel_info = flat.channel_info;
work.last_channel_info = flat.channel_info;
end

function [H,geom,metric,trace,work] = optimizeWork(work,snr_dB,cfg)
if strcmp(work.model.family,'pca')
    [H,geom,metric,trace,work] = pcaOptimize(work,snr_dB,cfg);
    return
end
% 曲面改变阵元位置与坐标，必须重建整个 N x K x F 信道。
% 不采用 orientation 模式逐阵元独立候选、仅替换一行的假设。
initialHits = work.cache_hits;
initialMisses = work.cache_misses;
[current,work] = getEntry(work,work.model.zero,cfg);
metric = uplinkOFDMRate(current.H,snr_dB,cfg);
value = objective(metric,cfg.opt_mode);
fixedValue = value;
trace.fixed = value;
trace.values = value;
trace.evaluations = 1;
trace.exact_evaluations = 1;
trace.gram_scores = 0;
trace.accepted = 0;
trace.sweeps = 0;
trace.feasible = 1;
trace.invalid = 0;
trace.constraint_rejections = 0;
trace.cp_rejections = 0;
trace.search_method = '';
switch work.model.family
    case 'subarray'
        trace.search_method = 'coordinate ascent of rigid-subarray signed displacement, yaw and tilt';
        maxSweeps = cfg.subarray.sweeps;
        for sweep = 1:maxSweeps
            changed = false;
            for d = 1:work.model.parameter_count
                best = current;
                bestMetric = metric;
                bestValue = value;
                for i = 1:numel(work.model.grids{d})
                    parameters = current.parameters;
                    parameters(d) = work.model.grids{d}(i);
                    if isequal(parameters,current.parameters)
                        continue
                    end
                    [candidate,work] = getEntry(work,parameters,cfg);
                    [better,candidateMetric,actual,trace] = assessCandidate(candidate,bestValue,snr_dB,cfg,trace);
                    if better
                        best = candidate;
                        bestMetric = candidateMetric;
                        bestValue = actual;
                    end
                end
                if bestValue>value+cfg.objective_tol
                    current = best;
                    metric = bestMetric;
                    value = bestValue;
                    changed = true;
                    trace.accepted = trace.accepted+1;
                    trace.values(end+1,1) = value;
                end
            end
            trace.sweeps = sweep;
            if ~changed
                break
            end
        end
end
H = current.H;
geom = current.geometry;
metric = uplinkOFDMRate(H,snr_dB,cfg);
trace.exact_evaluations = trace.exact_evaluations+1;
finalValue = objective(metric,cfg.opt_mode);
assert(finalValue>=fixedValue-cfg.check.objective_tol*max(1,abs(fixedValue)),'Deformation objective fell below planar baseline.');
assert(all(diff(trace.values)>=-cfg.check.objective_tol),'Accepted deformation objective decreased.');
assert(current.channel_info.cp_feasible,'Selected geometry must satisfy CP.');
work.last_channel_info = current.channel_info;
trace.cache_hits = work.cache_hits-initialHits;
trace.cache_misses = work.cache_misses-initialMisses;
trace.cache_size = numel(work.cache);
trace.final_parameters = geom.parameters;
trace.final_value = finalValue;
end

function [better,metric,actual,trace] = assessCandidate(candidate,currentValue,snr_dB,cfg,trace)
better = false;
metric = [];
actual = -Inf;
trace.evaluations = trace.evaluations+1;
if ~candidate.valid
    trace.invalid = trace.invalid+1;
    if strcmp(candidate.rejection,'cp')
        trace.cp_rejections = trace.cp_rejections+1;
    else
        trace.constraint_rejections = trace.constraint_rejections+1;
    end
    return
end
trace.feasible = trace.feasible+1;
% 缓存 Gram 代表全部 K 个 UE；非全带共享分配需逐频点选择真实调度集合。
sharedAllocation = all(cfg.ofdm.user_sample_mask(:)) && all(cfg.ofdm.user_num_active(:)==cfg.ofdm.num_active);
% 下方闭式评分仅适用于白噪声 ZF；MMSE-IRC/外部相关干扰使用完整接收机。
gramMode = strcmp(cfg.receiver,'zf') && isempty(cfg.irc.external_covariance) && ...
    strcmp(cfg.opt.score_mode,'gram_verified') && cfg.check.svd_rel_tol==0 && candidate.gram_valid && sharedAllocation;
if gramMode
    rho = 10.^(snr_dB/10);
    if strcmp(cfg.power.snr_reference,'total_over_tone_noise')
        rho = rho/cfg.ofdm.num_active;
    end
    ratePerTone = log1p(rho./candidate.gram_inverse_diagonal)/log(2);
    rate = cfg.ofdm.cp_efficiency*(ratePerTone*cfg.ofdm.freq_weights(:));
    if strcmp(cfg.opt_mode,'sum_rate')
        score = sum(rate);
    else
        score = min(rate);
    end
    trace.gram_scores = trace.gram_scores+1;
    % 对评分差距接近 Gram 舍入误差的候选也做完整接收机复核。
    % Gram 筛选只用于良态满秩信道，不用于秩亏伪逆或显式 SVD 截断。
    margin = max(cfg.objective_tol,64*eps*max(1,abs(score))/candidate.gram_min_rcond);
    if isfinite(score) && score<currentValue-margin
        return
    end
end
metric = uplinkOFDMRate(candidate.H,snr_dB,cfg);
trace.exact_evaluations = trace.exact_evaluations+1;
actual = objective(metric,cfg.opt_mode);
better = actual>currentValue+cfg.objective_tol;
end

function [entry,work] = getEntry(work,parameters,cfg)
parameters = canonicalParameters(parameters,work.model);
for i = 1:numel(work.cache)
    if isequal(work.cache{i}.parameters,parameters)
        entry = work.cache{i};
        work.cache_hits = work.cache_hits+1;
        return
    end
end
work.cache_misses = work.cache_misses+1;
entry.parameters = parameters;
entry.H = [];
entry.gram = [];
entry.gram_inverse_diagonal = [];
entry.gram_valid = false;
entry.gram_min_rcond = 0;
entry.channel_info = struct();
[entry.geometry,valid,reason] = fimGeometry('build',work.model,parameters,cfg);
entry.valid = valid;
entry.reason = reason;
entry.rejection = '';
if ~valid
    entry.rejection = 'geometry';
else
    [entry.H,entry.channel_info] = fimChannel('geometry',work.realization,work.pattern,entry.geometry,cfg);
    if ~entry.channel_info.cp_feasible
        entry.valid = false;
        entry.reason = 'CP shorter than geometry-dependent arrival delay span';
        entry.rejection = 'cp';
    else
        if strcmp(cfg.receiver,'zf') && isempty(cfg.irc.external_covariance)
            [entry.gram,entry.gram_inverse_diagonal,entry.gram_valid,entry.gram_min_rcond] = gramCache(entry.H,cfg);
        end
    end
end
% 平面始终固定在首个缓存槽位，其他槽位按 FIFO 淘汰；不使用 containers.Map。
% cache_limit=1 时只缓存平面，仍允许直接计算其他形状。
if numel(work.cache)<work.cache_limit
    work.cache{end+1} = entry;
elseif work.cache_limit>1
    work.cache(2) = [];
    work.cache{end+1} = entry;
    work.cache_evictions = work.cache_evictions+1;
end
end

function parameters = canonicalParameters(parameters,model)
parameters = parameters(:).';
assert(numel(parameters)==model.parameter_count && all(isfinite(parameters)),'Invalid shape parameter vector.');
% 总弯曲角为零时，圆柱轴方向不影响物理平面，合并这些等价缓存项。
end

function [G,diagonal,valid,minRcond] = gramCache(H,cfg)
[~,K,F] = size(H);
G = complex(zeros(K,K,F));
diagonal = zeros(K,F);
valid = true;
minRcond = Inf;
I = eye(K);
for f = 1:F
    A = H(:,:,f)'*H(:,:,f);
    A = (A+A')/2;
    G(:,:,f) = A;
    rc = rcond(A);
    minRcond = min(minRcond,rc);
    [R,failed] = chol(A);
    if failed~=0 || ~isfinite(rc) || rc<cfg.opt.gram_rcond_min
        valid = false;
        continue
    end
    invR = R\I;
    % A=R'*R，逆矩阵的对角是 invR 各行的平方范数。
    diagonal(:,f) = sum(abs(invR).^2,2);
end
valid = valid && all(isfinite(diagonal(:))) && all(diagonal(:)>0);
end

function value = objective(metric,mode)
switch mode
    case 'sum_rate'
        value = metric.sumRate;
    case 'max_min_rate'
        value = metric.worstRate;
    otherwise
        error('Unknown objective: %s',mode);
end
end

function [H,geom,metric,trace,work] = pcaOptimize(work,snr_dB,cfg)
%FIMPCAOPTIMIZE 每个阵元在原阵面 x-z 内移动的离散坐标上升。
% 每次仅一个阵元的位置改变，因此仅重算这一行物理信道；其余行完全复用。
% 候选评分始终调用完整 uplinkOFDMRate，保留 MMSE-IRC/ZF、外部协方差、
% RB 分配及全频率权重。没有用 ZF Gram 近似代替 MMSE，也不逐频点选位置。
% 位置/间距约束用完整几何检查，CP 用所有阵元/UE/路径的共同到达时延跨度。
% 最终用完整 N 阵元信道重算验证行更新结果。
assert(strcmp(work.model.family,'pca'),'pcaOptimize requires a PCA model.');
model = work.model;
parameters = model.zero;
geom = work.fixed_geometry;
H = work.Hfixed;
metric = uplinkOFDMRate(H,snr_dB,cfg);
value = pcaObjective(metric,cfg.opt_mode);
fixedValue = value;
[delayMin,delayMax] = elementDelayBounds(geom.positions_m,work.realization,cfg);
trace.fixed = value;
trace.values = value;
trace.evaluations = 1;
trace.exact_evaluations = 1;
trace.gram_scores = 0;
trace.accepted = 0;
trace.sweeps = 0;
trace.feasible = 1;
trace.invalid = 0;
trace.constraint_rejections = 0;
trace.cp_rejections = 0;
trace.channel_row_evaluations = 0;
trace.full_channel_evaluations = 0;
trace.search_method = 'PCA coordinate ascent; exact single-row channel updates and full receiver scores';
% 此配置仅传给几何信道构建器，不进行 prepareFIMConfig，不能改变接收机维度。
% 接收机评分仍始终使用完整原 cfg 和 N x K x F 信道。
rowCfg = cfg;
rowCfg.N = 1;
searchClock = tic;
for sweep = 1:cfg.pca.sweeps
    changed = false;
    for d = 1:model.parameter_count
        element = ceil(d/2);
        bestValue = value;
        bestH = H;
        bestMetric = metric;
        bestGeometry = geom;
        bestParameters = parameters;
        bestMin = delayMin(element);
        bestMax = delayMax(element);
        for j = 1:numel(model.grids{d})
            trialParameters = parameters;
            trialParameters(d) = model.grids{d}(j);
            if isequal(trialParameters,parameters)
                continue
            end
            trace.evaluations = trace.evaluations+1;
            [trialGeometry,valid,~] = fimGeometry('build',model,trialParameters,cfg);
            if ~valid
                trace.invalid = trace.invalid+1;
                trace.constraint_rejections = trace.constraint_rejections+1;
                continue
            end
            rowGeometry.positions_m = trialGeometry.positions_m(element,:);
            rowGeometry.R = trialGeometry.R(:,:,element);
            [rowH,rowInfo] = fimChannel('geometry',work.realization,work.pattern,rowGeometry,rowCfg);
            trace.channel_row_evaluations = trace.channel_row_evaluations+1;
            trialMin = delayMin;
            trialMax = delayMax;
            trialMin(element) = rowInfo.delay_min_s;
            trialMax(element) = rowInfo.delay_max_s;
            span = max(trialMax)-min(trialMin);
            cpFeasible = ~cfg.ofdm.enabled || span<=cfg.ofdm.cp_duration_s+1e-15;
            if ~rowInfo.cp_feasible || ~cpFeasible
                trace.invalid = trace.invalid+1;
                trace.cp_rejections = trace.cp_rejections+1;
                continue
            end
            trace.feasible = trace.feasible+1;
            trialH = H;
            trialH(element,:,:) = rowH;
            trialMetric = uplinkOFDMRate(trialH,snr_dB,cfg);
            trace.exact_evaluations = trace.exact_evaluations+1;
            trialValue = pcaObjective(trialMetric,cfg.opt_mode);
            if trialValue>bestValue+cfg.objective_tol
                bestValue = trialValue;
                bestH = trialH;
                bestMetric = trialMetric;
                bestGeometry = trialGeometry;
                bestParameters = trialParameters;
                bestMin = rowInfo.delay_min_s;
                bestMax = rowInfo.delay_max_s;
            end
        end
        if bestValue>value+cfg.objective_tol
            value = bestValue;
            H = bestH;
            metric = bestMetric;
            geom = bestGeometry;
            parameters = bestParameters;
            delayMin(element) = bestMin;
            delayMax(element) = bestMax;
            changed = true;
            trace.accepted = trace.accepted+1;
            trace.values(end+1,1) = value;
        end
        if cfg.N>=64 && mod(d,2)==0 && (mod(element,32)==0 || element==cfg.N)
            fprintf('    PCA sweep %d/%d, element %d/%d, feasible scores=%d, elapsed %.1f s\n',...
                sweep,cfg.pca.sweeps,element,cfg.N,trace.exact_evaluations,toc(searchClock));
        end
    end
    trace.sweeps = sweep;
    if ~changed
        break
    end
end
% 统一返回完整信道构建器的结果和完整几何元数据，并检查不存在行更新偏差。
[fullH,info] = fimChannel('geometry',work.realization,work.pattern,geom,cfg);
trace.full_channel_evaluations = trace.full_channel_evaluations+1;
assert(info.cp_feasible,'Selected PCA geometry must satisfy the common CP constraint.');
relativeChannelError = norm(fullH(:)-H(:))/max(norm(fullH(:)),realmin);
assert(relativeChannelError<1e-10,'PCA row-update channel disagrees with full geometry rebuild.');
H = fullH;
metric = uplinkOFDMRate(H,snr_dB,cfg);
trace.exact_evaluations = trace.exact_evaluations+1;
finalValue = pcaObjective(metric,cfg.opt_mode);
assert(finalValue>=fixedValue-cfg.check.objective_tol,'PCA pcaObjective fell below the fixed baseline.');
assert(abs(finalValue-value)<=cfg.check.objective_tol*max(1,abs(finalValue)),...
    'PCA row-update and full-rebuild objectives differ.');
assert(all(diff(trace.values)>=-cfg.check.objective_tol),'Accepted PCA pcaObjective decreased.');
trace.final_parameters = parameters;
trace.final_value = finalValue;
trace.row_update_relative_error = relativeChannelError;
trace.cache_hits = 0;
trace.cache_misses = 0;
trace.cache_size = numel(work.cache);
work.last_channel_info = info;
end

function [earliest,latest] = elementDelayBounds(positions,realization,cfg)
% 与 fimChannel geometry 相同的有效时延；绝对阵面平移不被误当成延迟扩展。
N = size(positions,1);
earliest = Inf(N,1);
latest = -Inf(N,1);
for k = 1:cfg.K
    path = realization.paths{k};
    active = abs(path.alpha(:))>0;
    if any(active)
        referenceDelay = path.tau(:);
        delay = bsxfun(@minus,referenceDelay(active).',positions*path.U(:,active)/cfg.c0);
        earliest = min(earliest,min(delay,[],2));
        latest = max(latest,max(delay,[],2));
    end
end
if all(isinf(earliest))
    earliest = zeros(N,1);
    latest = zeros(N,1);
end
end

function value = pcaObjective(metric,mode)
if strcmp(mode,'sum_rate')
    value = metric.sumRate;
elseif strcmp(mode,'max_min_rate')
    value = metric.worstRate;
else
    error('Unknown PCA pcaObjective: %s',mode);
end
end

function [H,selection,metric,trace] = orientation_optimizeFIM(Hcand,states,snr_dB,cfg)
%OPTIMIZEFIM 全频带共用阵元角度的离散坐标上升。
% Hcand:N x K x C x F。每次改变一个阵元在所有用户、频点上的整行。
% 局部搜索不保证全局最优，但所有接受的更新均经完整接收机验证。
[N,K,C,F] = size(Hcand);
assert(N==cfg.N && K==cfg.K && F==cfg.ofdm.num_tones && C==states.count,'Invalid candidate shape.');
selection = repmat(states.zero,N,1);
H = orientation_selectChannel(Hcand,selection);
metric = uplinkOFDMRate(H,snr_dB,cfg);
fixedValue = orientation_objective(metric,cfg.opt_mode);
value = fixedValue;
trace.fixed = fixedValue;
trace.values = fixedValue;
trace.evaluations = 1;
trace.exact_evaluations = 1;
trace.gram_scores = 0;
trace.accepted = 0;
trace.sweeps = 0;
% 启发式只提供初值；必须真实带内速率优于固定阵列才采用。
power = zeros(N,K,C);
for f = 1:F
    % 启发式只累计该 UE 获分配资源上的信道能量。
    tonePower = bsxfun(@times,abs(Hcand(:,:,:,f)).^2,reshape(cfg.ofdm.user_sample_mask(:,f),1,K,1));
    power = power+cfg.ofdm.freq_weights(f)*tonePower;
end
if strcmp(cfg.opt_mode,'sum_rate')
    heuristic = reshape(sum(power,2),N,C);
else
    heuristic = reshape(min(power,[],2),N,C);
end
[~,init] = max(heuristic,[],2);
if any(init~=selection)
    initH = orientation_selectChannel(Hcand,init);
    initMetric = uplinkOFDMRate(initH,snr_dB,cfg);
    trace.exact_evaluations = trace.exact_evaluations+1;
    trace.evaluations = trace.evaluations+1;
    initValue = orientation_objective(initMetric,cfg.opt_mode);
    if initValue>value+cfg.objective_tol
        H = initH;
        selection = init;
        metric = initMetric;
        value = initValue;
        trace.values(end+1,1) = value;
    end
end
rho = 10.^(snr_dB/10);
if strcmp(cfg.power.snr_reference,'total_over_tone_noise')
    rho = rho/cfg.ofdm.num_active;
end
% 若人为设置 SVD 截断阈值，满秩 Gram 公式可能不再等于实际接收机。
% 此时自动用完整评分，保持两种模式的物理口径一致。
% 不同 RB 的用户集合或功率分配不同，旧的全 K 用户 Gram 公式不再适用。
% 这类分配自动使用逐频点调度感知的完整接收机评分。
sharedAllocation = all(cfg.ofdm.user_sample_mask(:)) && all(cfg.ofdm.user_num_active(:)==cfg.ofdm.num_active);
% 原 Gram 闭式评分只用于白噪声 ZF；其他接收机完整评估每个候选。
gramMode = strcmp(cfg.receiver,'zf') && isempty(cfg.irc.external_covariance) && ...
    strcmp(cfg.opt.score_mode,'gram_verified') && cfg.check.svd_rel_tol==0 && sharedAllocation;
G = [];
if gramMode
    G = orientation_channelGram(H);
end
if C>1
    for sweep = 1:cfg.orientation.sweeps
        changed = false;
        for n = 1:N
            scores = -Inf(C,1);
            scores(selection(n)) = value;
            % 缓存当前整行移除后的 K x K Gram，避免候选内 N x K QR。
            withoutRow = [];
            if gramMode
                withoutRow = G;
                for f = 1:F
                    oldRow = H(n,:,f);
                    withoutRow(:,:,f) = G(:,:,f)-oldRow'*oldRow;
                end
            end
            for c = 1:C
                if c==selection(n)
                    continue
                end
                row = reshape(Hcand(n,:,c,:),K,F);
                valid = false;
                if gramMode
                    [scores(c),valid] = orientation_gramScore(withoutRow,row,rho,cfg);
                    trace.gram_scores = trace.gram_scores+1;
                end
                if ~valid
                    trial = H;
                    trial(n,:,:) = reshape(row,1,K,F);
                    candidate = uplinkOFDMRate(trial,snr_dB,cfg);
                    scores(c) = orientation_objective(candidate,cfg.opt_mode);
                    trace.exact_evaluations = trace.exact_evaluations+1;
                end
                trace.evaluations = trace.evaluations+1;
            end
            % 即便评分使用 Gram，也按完整 QR/SVD 接收后的目标确认。
            % 若最高评分因数值误差没有改善，继续检查下一候选。
            [~,order] = sort(scores,'descend');
            for i = 1:C
                c = order(i);
                if scores(c)<=value+cfg.objective_tol
                    break
                end
                if c==selection(n)
                    continue
                end
                trial = H;
                trial(n,:,:) = reshape(Hcand(n,:,c,:),1,K,F);
                candidate = uplinkOFDMRate(trial,snr_dB,cfg);
                trace.exact_evaluations = trace.exact_evaluations+1;
                actual = orientation_objective(candidate,cfg.opt_mode);
                if actual>value+cfg.objective_tol
                    H = trial;
                    metric = candidate;
                    selection(n) = c;
                    value = actual;
                    if gramMode
                        G = orientation_channelGram(H);
                    end
                    changed = true;
                    trace.accepted = trace.accepted+1;
                    trace.values(end+1,1) = actual;
                    break
                end
            end
        end
        trace.sweeps = sweep;
        if ~changed
            break
        end
    end
end
metric = uplinkOFDMRate(H,snr_dB,cfg);
trace.exact_evaluations = trace.exact_evaluations+1;
assert(orientation_objective(metric,cfg.opt_mode)>=fixedValue-cfg.check.objective_tol*max(1,abs(fixedValue)),'FIM fell below fixed baseline.');
assert(all(diff(trace.values)>=-cfg.check.objective_tol),'Accepted orientation_objective decreased.');
end

function H = orientation_selectChannel(Hcand,selection)
[N,K,~,F] = size(Hcand);
H = complex(zeros(N,K,F));
for n = 1:N
    H(n,:,:) = reshape(Hcand(n,:,selection(n),:),1,K,F);
end
end

function G = orientation_channelGram(H)
[~,K,F] = size(H);
G = complex(zeros(K,K,F));
for f = 1:F
    G(:,:,f) = H(:,:,f)'*H(:,:,f);
end
end

function [value,valid] = orientation_gramScore(base,row,rho,cfg)
% 满秩 ZF: gamma_k=rho/[(H'*H)^(-1)]_kk。
% 近奇异候选交回完整 QR/SVD，避免将伪逆误当成无干扰 ZF。
[K,F] = size(row);
rates = zeros(K,F);
value = -Inf;
valid = false;
I = eye(K);
for f = 1:F
    r = row(:,f).';
    G = base(:,:,f)+r'*r;
    G = (G+G')/2;
    [R,failed] = chol(G);
    if failed~=0 || rcond(G)<cfg.opt.gram_rcond_min
        return
    end
    invR = R\I;
    % G=R'*R -> G^-1=inv(R)*inv(R)'，对角是 invR 的行平方和。
    diagonal = sum(abs(invR).^2,2);
    rates(:,f) = log1p(rho./diagonal)/log(2);
end
bandRate = cfg.ofdm.cp_efficiency*(rates*cfg.ofdm.freq_weights(:));
if strcmp(cfg.opt_mode,'sum_rate')
    value = sum(bandRate);
else
    value = min(bandRate);
end
valid = isfinite(value);
end

function value = orientation_objective(metric,mode)
if strcmp(mode,'sum_rate')
    value = metric.sumRate;
else
    value = metric.worstRate;
end
end

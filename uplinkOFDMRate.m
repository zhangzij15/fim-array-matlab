function metric = uplinkOFDMRate(H,snr_dB,cfg)
%UPLINKOFDMRATE 逐子载波接收后 SINR，以及频带内用户速率的唯一计算入口。
% cfg 需先经过 prepareFIMConfig；调用主引擎时会自动完成。
% H 的维度为 N x K x F；F 是实际评估频点数，不一定等于活动频点数。
% 每个频点只为获分配资源的 UE 构造 ZF/MMSE-IRC，不纳入未调度列。
% 先计算 log2(1+SINR)，再按全小区活动频带的频率权重积分。
[N,K,F] = size(H);
assert(N==cfg.N && K==cfg.K && F==cfg.ofdm.num_tones,'H must be N x K x F.');
assert(all(isfinite(H(:))),'Channel contains nonfinite values.');
scheduled = cfg.ofdm.user_sample_mask;
userCount = cfg.ofdm.user_num_active(:);
assert(isequal(size(scheduled),[K F]) && all(scheduled(:)==0 | scheduled(:)==1),'Invalid per-user frequency allocation mask.');
assert(numel(userCount)==K && all(isfinite(userCount)) && all(userCount>0) && all(userCount<=cfg.ofdm.num_active),'Each UE requires a positive allocated tone count.');
rho = 10.^(snr_dB/10)*ones(K,1);
if strcmp(cfg.power.snr_reference,'total_over_tone_noise')
    % 每个 UE 的总功率只分摊到它实际获分配的活动频点，不能除采样频点数 F。
    rho = rho./userCount;
end
P = cfg.noiseVar*rho;
assert(all(isfinite(P)) && all(P>0),'SNR outside supported numeric range.');
isMMSE = strcmp(cfg.receiver,'mmse-irc');
assert(isMMSE || strcmp(cfg.receiver,'zf'),'Unsupported receiver.');
% 外部协方差不含热噪声，也不含下方 H 中已显式计算的同频 UE 信号。
% 它在所有频点共用；白化因子只构造一次，避免每个频点做 N x N 分解。
externalCovariance = [];
if isfield(cfg,'irc') && isfield(cfg.irc,'external_covariance')
    externalCovariance = cfg.irc.external_covariance;
end
hasExternal = ~isempty(externalCovariance) && any(externalCovariance(:)~=0);
whiteningL = [];
if isMMSE && hasExternal
    C = eye(N)+externalCovariance/cfg.noiseVar;
    C = (C+C')/2;
    [whiteningR,failed] = chol(C);
    assert(failed==0,'Thermal noise plus external interference must be positive definite.');
    whiteningL = whiteningR';
end
metric.sinrPerTone = zeros(K,F);
metric.desired = zeros(K,F);
metric.interference = zeros(K,F);
metric.noise = zeros(K,F);
metric.externalInterference = zeros(K,F);
metric.effectiveChannel = complex(zeros(K,K,F));
receiver.rank = zeros(1,F);
receiver.rank_deficient = false(1,F);
receiver.leakage_relative = zeros(1,F);
receiver.rcondR = zeros(1,F);
receiver.method = cell(1,F);
receiver.num_scheduled = sum(scheduled,1);
for f = 1:F
    users = find(scheduled(:,f));
    U = numel(users);
    if U==0
        receiver.method{f} = 'inactive';
        continue
    end
    Hf = H(:,users,f);
    if isMMSE
        [W,info] = stableMMSEIRC(Hf,rho(users),whiteningL,cfg);
    else
        [W,info] = stableZF(Hf,cfg);
    end
    E = W'*Hf;
    powerMatrix = bsxfun(@times,abs(E).^2,P(users).');
    desired = real(diag(powerMatrix));
    powerMatrix(1:U+1:end) = 0;
    interference = sum(powerMatrix,2);
    noise = cfg.noiseVar*sum(abs(W).^2,1).';
    denominator = interference+noise;
    if hasExternal
        external = max(real(sum(conj(W).*(externalCovariance*W),1)).',0);
        denominator = denominator+external;
        metric.externalInterference(users,f) = external;
    end
    active = denominator>0;
    % MLang 对空数组广播支持有限；全零信道时避免对空掩码执行赋值运算。
    for u = 1:U
        if active(u)
            metric.sinrPerTone(users(u),f) = desired(u)/denominator(u);
        end
    end
    metric.desired(users,f) = desired;
    metric.interference(users,f) = interference;
    metric.noise(users,f) = noise;
    metric.effectiveChannel(users,users,f) = E;
    receiver.rank(f) = info.rank;
    receiver.rank_deficient(f) = info.rank<U;
    receiver.leakage_relative(f) = info.leakage_relative;
    receiver.rcondR(f) = info.rcondR;
    receiver.method{f} = info.method;
end
weights = cfg.ofdm.freq_weights(:);
metric.snrPerTone = zeros(K,F);
hasNoise = metric.noise>0;
metric.snrPerTone(hasNoise) = metric.desired(hasNoise)./metric.noise(hasNoise);
metric.snr = expm1(log(2)*(log1p(metric.snrPerTone)/log(2)*weights));
metric.ratePerToneNoCP = log1p(metric.sinrPerTone)/log(2);
metric.rateNoCP = metric.ratePerToneNoCP*weights;
metric.rate = cfg.ofdm.cp_efficiency*metric.rateNoCP;
% rate 始终以全小区活动频带归一化；未调度频点贡献零，便于比较不同分配。
allocationFraction = userCount/cfg.ofdm.num_active;
metric.allocatedBandwidth_Hz = cfg.ofdm.active_bandwidth_Hz*allocationFraction;
metric.rateScheduledBandNoCP = metric.rateNoCP./allocationFraction;
metric.rateScheduledBand = metric.rate./allocationFraction;
% 全小区活动频带上的速率等效 SINR，不能称为某个子载波的物理 SINR。
% CP 降低净速率，但不改变此定义下的速率等效 SINR。
metric.sinr = expm1(log(2)*metric.rateNoCP);
metric.sumRate = sum(metric.rate);
metric.worstRate = min(metric.rate);
metric.worstSINR = min(metric.sinr);
% 带宽按活动频点定义，不能再额外乘活动频点/FFT 比例，否则重复扣减。
metric.throughput_bps = cfg.ofdm.active_bandwidth_Hz*metric.rate;
metric.sumThroughput_bps = sum(metric.throughput_bps);
metric.worstThroughput_bps = min(metric.throughput_bps);
metric.txPowerPerTone = P;
metric.txPower = P;
% 保留上面 K x 1 的旧字段：每个 UE 在获分配频点上的功率。
% 新字段显式给出 K x F 实际功率矩阵；未调度频点功率严格为零。
metric.txPowerPerUserPerTone = bsxfun(@times,double(scheduled),P);
metric.txPowerTotal = P.*userCount;
metric.scheduledMask = logical(scheduled);
metric.userNumActive = userCount;
metric.receiver = receiver;
assert(all(isfinite(metric.rate)) && all(metric.rate>=0),'Invalid user rates.');
end

function [W,info] = stableMMSEIRC(H,rho,whiteningL,cfg)
% A=L\H*diag(sqrt(P/noiseVar)); L*L'=I+Rexternal/noiseVar。
% W 与 (H*diag(P)*H'+noiseVar*I+Rexternal)\H*diag(sqrt(P))
% 逐列成比例，因而给出相同 SINR。利用增广 QR 求解小型用户维问题。
[N,K] = size(H);
if isempty(whiteningL)
    Hw = H;
else
    Hw = whiteningL\H;
end
A = bsxfun(@times,Hw,sqrt(rho(:)).');
assert(all(isfinite(A(:))),'MMSE-IRC scaled channel exceeds numeric range.');
% 增广单位阵提供真实热噪声正则项；没有人为 SVD 截断或零迫近似。
[~,R] = qr([A;eye(K)],0);
% 不从 Q 的顶端切片恢复 A/R：极弱信道下 Householder 消减会丢失分量。
% 良态时直接从 A 回代；很强且相关的信道用 SVD，避免两次回代放大误差。
regularizedRcond = rcond(R);
method = 'mmse-irc-qr';
if isfinite(regularizedRcond) && regularizedRcond>=1e-4
    W = (A/R)/R';
else
    [Ua,Sa,Va] = svd(A,'econ');
    singular = diag(Sa);
    factors = zeros(size(singular));
    for k = 1:numel(singular)
        if singular(k)<=1
            factors(k) = singular(k)/(1+singular(k)^2);
        else
            factors(k) = 1/(singular(k)+1/singular(k));
        end
    end
    r = numel(singular);
    W = bsxfun(@times,Ua(:,1:r),factors.')*Va(:,1:r)';
    method = 'mmse-irc-svd';
end
if ~isempty(whiteningL)
    W = whiteningL'\W;
end
% 列缩放不改变 SINR；先按最大模缩放，避免极端功率下范数溢出/下溢。
for k = 1:K
    scale = max(abs(W(:,k)));
    if scale>0
        W(:,k) = W(:,k)/scale;
        W(:,k) = W(:,k)/sqrt(sum(abs(W(:,k)).^2));
    end
end
assert(all(isfinite(W(:))),'MMSE-IRC weights exceed numeric range.');
% 诊断原始信道秩，不能用带单位正则项的增广矩阵秩掩盖信道秩亏。
[~,channelR] = qr(H,0);
rc = 0;
if N>=K
    rc = rcond(channelR);
end
if isfinite(rc) && rc>=cfg.check.zf_rcond_switch
    info.rank = K;
else
    s = svd(H,'econ');
    if isempty(s) || max(s)==0
        info.rank = 0;
    else
        info.rank = sum(s>max(N,K)*eps(max(s)));
    end
end
E = W'*H;
off = E;
off(1:K+1:end) = 0;
info.rcondR = rc;
info.leakage_relative = norm(off,'fro')/max(norm(E,'fro'),realmin);
info.method = method;
end

function [W,info] = stableZF(H,cfg)
% W:N x K，判决 z=W'*y。用 QR 求解，避免显式逆矩阵。
% 近奇异时使用 SVD 伪逆，并保留其真实残余干扰。
[N,K] = size(H);
[Q,R] = qr(H,0);
rc = rcond(R);
if isfinite(rc) && rc>=cfg.check.zf_rcond_switch
    W = Q/R';
    info.rank = K;
    info.method = 'qr';
else
    [U,S,V] = svd(H,'econ');
    s = diag(S);
    if isempty(s) || max(s)==0
        keep = false(size(s));
    elseif cfg.check.svd_rel_tol>0
        keep = s>cfg.check.svd_rel_tol*max(s);
    else
        keep = s>max(N,K)*eps(max(s));
    end
    W = complex(zeros(N,K));
    info.rank = sum(keep);
    info.method = 'svd';
    if info.rank>0
        W = bsxfun(@rdivide,U(:,keep),s(keep).')*V(:,keep)';
    end
end
% 逐列单位范数缩放仅用于数值稳定；这不是下行发射功率归一化。
scale = sqrt(sum(abs(W).^2,1));
for k = 1:K
    if scale(k)>0
        W(:,k) = W(:,k)/scale(k);
    end
end
E = W'*H;
off = E;
off(1:K+1:end) = 0;
info.rcondR = rc;
info.leakage_relative = norm(off,'fro')/max(norm(E,'fro'),realmin);
end

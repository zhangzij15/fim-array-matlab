function [snapshots,files] = run_fim_angular_sinr(patternTxt,options)
%RUN_FIM_ANGULAR_SINR Simulate and export directional SINR for four arrays.
%   run_fim_angular_sinr()
%   run_fim_angular_sinr("D:\patterns\measured.txt")
%   run_fim_angular_sinr("D:\patterns\measured.txt",KValues=1,NumRB=24,Sweeps=1)
%
% The default reproduces the 10-degree 288-RB snapshots, one realization per K.
% An empty patternTxt uses the bundled pattern. A TXT must have two header
% lines followed by 360*181 Phi Theta gain_dB phase rows (Phi varies fastest).
% The project uses the gain column, normalizes power over the sphere, and
% aligns the pattern peak to local +y unless PatternAlignment="manual".
% Set PatternNativeToBody for a calibrated 3-by-3 rotation. FrontOnly=false
% retains the supplied back lobes. Results go into a timestamped local folder.
arguments
    patternTxt (1,1) string = ""
    options.KValues (1,:) double {mustBeMember(options.KValues,[1 10])} = [1 10]
    options.NumRB (1,1) double {mustBeInteger,mustBePositive} = 288
    options.Sweeps (1,1) double {mustBeInteger,mustBeNonnegative} = 3
    options.OutputDir (1,1) string = ""
    options.ShowFigures (1,1) logical = false
    options.PatternAlignment (1,1) string {mustBeMember(options.PatternAlignment,["peak","manual"])} = "peak"
    options.PatternNativeToBody (3,3) double {mustBeReal,mustBeFinite} = eye(3)
    options.FrontOnly (1,1) logical = true
end

assert(~isempty(options.KValues) && numel(unique(options.KValues))==numel(options.KValues), ...
    'FIM:AngularSINRKValues','KValues must contain 1, 10, or both without repeats.');
cfg = angularStudyConfig(options);
if strlength(patternTxt)>0
    assert(isfile(patternTxt),'FIM:AngularSINRPatternNotFound', ...
        'Antenna pattern TXT not found: %s',patternTxt);
    cfg.pattern.file = char(patternTxt);
end
outputRoot = fullfile(cfg.output_dir,'angular_sinr');
if strlength(options.OutputDir)>0
    outputRoot = char(options.OutputDir);
end
stamp = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
outputDir = fullfile(outputRoot,['scan_' stamp]);
[created,message] = mkdir(outputDir);
assert(created,'FIM:AngularSINROutput','Cannot create output folder: %s',message);

azimuth_deg = -134:2:134;
elevation_deg = -40:2:34;
[az,el] = meshgrid(azimuth_deg,elevation_deg);
directions = [cosd(el(:)).*sind(az(:)),cosd(el(:)).*cosd(az(:)),sind(el(:))].';
schemes = {'fixed','pca','orientation','fim'};
collected = cell(1,numel(options.KValues));
summary = table();
files = {};
fprintf('Pattern TXT: %s\n',cfg.pattern.file);
for i = 1:numel(options.KValues)
    cfg.K = options.KValues(i);
    if cfg.K==1
        cfg.seed = 20365657;
    else
        cfg.seed = 21575118;
    end
    % The simulator samples once and optimizes all four arrays on that realization.
    result = simulateFIM(cfg,'cdf');
    snapshot = angularScan(result,directions,azimuth_deg,elevation_deg,schemes);
    collected{i} = snapshot;
    rows = table(repmat(snapshot.K,4,1),repmat(snapshot.Seed,4,1), ...
        repmat(snapshot.TargetUE,4,1),string(schemes(:)),snapshot.MMSE_IRC_dB(:), ...
        snapshot.PeakScan_dB(:),snapshot.PeakAzimuth_deg(:),snapshot.PeakElevation_deg(:), ...
        snapshot.SumSE_bit_s_Hz(:), ...
        'VariableNames',{'K','Seed','TargetUE','Scheme','MMSE_IRC_dB','PeakScan_dB', ...
        'PeakAzimuth_deg','PeakElevation_deg','SumSE_bit_s_Hz'});
    summary = [summary; rows]; %#ok<AGROW>
    figures = renderAngularSINR(snapshot,outputDir,options.ShowFigures);
    files = [files figures]; %#ok<AGROW>
end
snapshots = [collected{:}];
matFile = fullfile(outputDir,'angular_sinr_data.mat');
csvFile = fullfile(outputDir,'angular_sinr_summary.csv');
save(matFile,'snapshots','summary');
writetable(summary,csvFile);
files = [files {matFile,csvFile}];
fprintf('Angular SINR figures and numerical results saved in %s\n',outputDir);
end

function cfg = angularStudyConfig(options)
cfg = configFIMSimulation();
cfg.bs.downtilt_deg = 10;
cfg.dv_lambda = 1.5;
cfg.users.radius_range_m = [20 200];
cfg.ofdm.num_rb = options.NumRB;
cfg.ofdm.nfft = max(512,2^nextpow2(cfg.ofdm.subcarriers_per_rb*cfg.ofdm.num_rb));
cfg.ofdm.frequency_mode = 'grouped';
cfg.ofdm.num_frequency_samples = cfg.ofdm.num_rb;
cfg.snr_dB_list = cfg.snr_dB;
cfg.subarray.displacement_grid_lambda = 0;
cfg.subarray.yaw_grid_deg = [-30 0 30];
cfg.subarray.tilt_grid_deg = -40:10:40;
cfg.subarray.min_spacing_lambda = 0;
cfg.subarray.sweeps = options.Sweeps;
cfg.orientation.yaw_grid_deg = cfg.subarray.yaw_grid_deg;
cfg.orientation.tilt_grid_deg = cfg.subarray.tilt_grid_deg;
cfg.orientation.sweeps = options.Sweeps;
cfg.pca.z_grid_lambda = [-1.5 0 1.5];
cfg.pca.min_spacing_lambda = 0;
cfg.pca.sweeps = options.Sweeps;
cfg.pattern.alignment = char(options.PatternAlignment);
cfg.pattern.native_to_body = options.PatternNativeToBody;
cfg.front_only = options.FrontOnly;
cfg.MC = 1;
cfg.plot.enabled = false;
cfg.save_results = false;
cfg.save_figures = false;
cfg.output.save_realizations = true;
cfg.output.save_geometry = true;
end

function snapshot = angularScan(result,Uglobal,azimuth_deg,elevation_deg,schemes)
cfg = result.cfg;
realization = result.realizations{1};
[~,order] = sort(result.fixed.sinr(1,:,1),'ascend');
target = order(ceil(cfg.K/2));  % Lower median Fixed SINR in the same snapshot.
interferers = setdiff(1:cfg.K,target);
U = realization.bs_frame.R_global_from_array'*Uglobal;
shape = [numel(elevation_deg) numel(azimuth_deg)];
snapshot = struct('K',cfg.K,'Seed',cfg.seed,'TargetUE',target, ...
    'BSDowntilt_deg',cfg.bs.downtilt_deg,'PatternFile',cfg.pattern.file, ...
    'PatternDiagnostics',result.pattern.diagnostics,'Azimuth_deg',azimuth_deg, ...
    'Elevation_deg',elevation_deg,'SINR_dB',NaN([shape 4]),'Config',cfg);

geometry = {result.geometry.fixed,result.geometry.pca{1}, ...
    fimGeometry('orientation',result.geometry.fixed,result.orientation_states, ...
    result.orientation_selection(1,:,1),cfg),result.geometry.fim{1}};
for j = 1:4
    [H,info] = fimChannel('geometry',realization,result.pattern,geometry{j},cfg);
    assert(info.cp_feasible,'FIM:AngularSINRInfeasible','Selected %s geometry violates the cyclic prefix.',schemes{j});
    metric = uplinkOFDMRate(H,result.snr_dB(1),cfg);
    [amplitude,phase,active] = scanSteering(geometry{j},result.pattern,U,cfg);
    assert(any(active),'FIM:AngularSINRNoResponse','%s has no antenna response on the scan grid.',schemes{j});
    scan = widebandAngularSINR(amplitude,phase,H,cfg,metric,target,interferers);
    map = NaN(shape);
    map(active) = 10*log10(max(scan,realmin));
    snapshot.SINR_dB(:,:,j) = map;
    [peak,peakIndex] = max(scan);
    angleIndex = find(active);
    [iEl,iAz] = ind2sub(shape,angleIndex(peakIndex));
    snapshot.MMSE_IRC_dB(j) = 10*log10(metric.sinr(target));
    snapshot.PeakScan_dB(j) = 10*log10(max(peak,realmin));
    snapshot.PeakAzimuth_deg(j) = azimuth_deg(iAz);
    snapshot.PeakElevation_deg(j) = elevation_deg(iEl);
    snapshot.SumSE_bit_s_Hz(j) = metric.sumRate;
end
end

function [amplitude,phase,active] = scanSteering(geometry,pattern,U,cfg)
% Same pattern coordinates and front mask as the physical channel builder.
amplitude = zeros(cfg.N,size(U,2));
for n = 1:cfg.N
    local = geometry.R(:,:,n)'*U;
    native = pattern.native_to_body'*local;
    phi = mod(atan2d(native(2,:),native(1,:)),360);
    theta = acosd(max(-1,min(1,native(3,:))));
    gain = max(0,interp2(pattern.phi,pattern.theta,pattern.G,phi,theta,'linear'));
    if pattern.front_only, gain(local(2,:)<-1e-12) = 0; end
    amplitude(n,:) = sqrt(gain);
end
power = sum(amplitude.^2,1);
active = power>0;
amplitude = amplitude(:,active)./sqrt(power(active));
phase = (2*pi/cfg.c0)*(geometry.positions_m*U(:,active));
end

function bandSINR = widebandAngularSINR(amplitude,phase,H,cfg,metric,target,interferers)
% For unit-norm steering a, use w=C^(-1)*a in normalized noise units.
% w'*C*w = a'*C^(-1)*a. Woodbury uses only the interfering-user dimension.
offsets = cfg.ofdm.freq_offsets_Hz(:).';
if ~cfg.ofdm.frequency_dependent_steering, offsets(:) = 0; end
intervals = diff(offsets);
uniform = ~isempty(intervals) && all(abs(intervals-intervals(1))<1e-6);
if uniform, advance = exp(1i*intervals(1)*phase); end
rho = metric.txPowerPerTone(:)/cfg.noiseVar;
averageLog = zeros(1,size(amplitude,2));
for f = 1:numel(offsets)
    if ~uniform || mod(f-1,32)==0
        A = amplitude.*exp(1i*(cfg.fc_Hz+offsets(f))*phase);
    else
        A = A.*advance;
    end
    h = H(:,target,f);
    denominator = ones(1,size(A,2));
    if ~isempty(interferers)
        S = H(:,interferers,f).*sqrt(rho(interferers)).';
        covariance = eye(numel(interferers))+S'*S;
        h = h-S*(covariance\(S'*h));
        projected = S'*A;
        denominator = 1-real(sum(conj(projected).*(covariance\projected),1));
    end
    assert(all(denominator>0),'FIM:AngularSINRDenominator','IRC scan covariance produced an invalid denominator.');
    sinr = rho(target)*abs(h'*A).^2./denominator;
    assert(max(sinr)<=metric.sinrPerTone(target,f)*(1+1e-8)+1e-10, ...
        'FIM:AngularSINRUpperBound','The directional scan exceeded unrestricted MMSE-IRC.');
    averageLog = averageLog+cfg.ofdm.freq_weights(f)*log1p(sinr);
end
bandSINR = expm1(averageLog);
end

function files = renderAngularSINR(snapshot,outputDir,visible)
%RENDERANGULARSINR Export the four scan maps, a paired difference and a cut.
if islogical(visible)
    states = {'off','on'};
    visible = states{1+visible};
end
names = {'Fixed array','PCA','Pattern reconfigurable','FIM'};
keys = {'fixed','pca','orientation','fim'};
az = snapshot.Azimuth_deg;
el = snapshot.Elevation_deg;
maps = snapshot.SINR_dB;
finiteValues = sort(maps(isfinite(maps)));
assert(~isempty(finiteValues),'FIM:EmptyAngularScan','No finite angular scan SINR.');
upper = 5*ceil(finiteValues(end)/5);
lower = 5*floor(finiteValues(max(1,ceil(.01*numel(finiteValues))))/5);
mapLimits = [min(lower,upper-5) upper];
files = cell(1,6);
context = sprintf('K=%d, UE %d, BS tilt %.1f%c',snapshot.K,snapshot.TargetUE, ...
    snapshot.BSDowntilt_deg,char(176));
prefix = sprintf('angular_sinr_K%d_',snapshot.K);

for j = 1:4
    fig = figure('Visible',visible,'Color','white','Position',[80 100 880 610]);
    cleanup = onCleanup(@() closeAngularSINRHidden(fig));
    ax = axes(fig);
    angularSINRMap(ax,az,el,maps(:,:,j),mapLimits,parula(256),'Scan SINR (dB)');
    title(ax,sprintf('%s directional IRC scan',names{j}));
    description = sprintf('%s | peak scan %.2f dB | unrestricted MMSE-IRC %.2f dB', ...
        context,snapshot.PeakScan_dB(j),snapshot.MMSE_IRC_dB(j));
    if any(maps(:,:,j)<mapLimits(1),'all')
        description = sprintf('%s | color floor %g dB',description,mapLimits(1));
    end
    subtitle(ax,description);
    files{j} = fullfile(outputDir,[prefix keys{j} '.png']);
    print(fig,files{j},'-dpng','-r160');
    clear cleanup
end

fixed = maps(:,:,1);
fim = maps(:,:,4);
difference = fim-fixed;
threshold = -6;
common = isfinite(fim) & isfinite(fixed) & max(fim,fixed)>=threshold;
difference(~common) = NaN;
fig = figure('Visible',visible,'Color','white','Position',[80 100 880 610]);
cleanup = onCleanup(@() closeAngularSINRHidden(fig));
ax = axes(fig);
differenceColors = interp1([0 .5 1],[.13 .33 .68;.98 .98 .98;.79 .22 .14],linspace(0,1,257));
angularSINRMap(ax,az,el,difference,[-6 6],differenceColors,'FIM - Fixed (dB)');
title(ax,'Paired directional IRC scan gain');
subtitle(ax,sprintf('%s | both valid, at least one %g dB or higher | colors clipped at %c6 dB', ...
    context,threshold,char(177)));
files{5} = fullfile(outputDir,[prefix 'fim_minus_fixed.png']);
print(fig,files{5},'-dpng','-r160');
clear cleanup

fim(~isfinite(fim)) = -Inf;
[~,peakIndex] = max(fim(:));
[peakEl,peakAz] = ind2sub(size(fim),peakIndex);
curves = reshape(maps(:,peakAz,:),numel(el),4);
curves(~isfinite(curves)) = NaN;
fig = figure('Visible',visible,'Color','white','Position',[80 100 880 610]);
cleanup = onCleanup(@() closeAngularSINRHidden(fig));
ax = axes(fig,'FontSize',11);
hold(ax,'on');
curveColors = [.20 .23 .29;.23 .52 .69;.13 .54 .36;.76 .22 .16];
styles = {'--',':','-.','-'};
for j = 1:4
    plot(ax,el,curves(:,j),'Color',curveColors(j,:),'LineStyle',styles{j},'LineWidth',1.8);
end
xlabel(ax,['Global elevation (' char(176) ')']);
ylabel(ax,'Scan SINR (dB)');
title(ax,sprintf('Elevation cut at global azimuth %g%c',az(peakAz),char(176)));
subtitle(ax,sprintf('%s | azimuth selected at the FIM full-scan peak',context));
if isfinite(curves(peakEl,4))
    plot(ax,el(peakEl),curves(peakEl,4),'o','MarkerSize',6, ...
        'Color',curveColors(4,:),'HandleVisibility','off');
end
legend(ax,names,'Location','best');
xlim(ax,[min(el)-0.5*isscalar(el) max(el)+0.5*isscalar(el)]);
grid(ax,'on');
box(ax,'on');
files{6} = fullfile(outputDir,[prefix 'elevation_cut.png']);
print(fig,files{6},'-dpng','-r160');
end

function angularSINRMap(ax,az,el,values,limits,colors,colorLabel)
azEdges = angularSINREdges(az);
elEdges = angularSINREdges(el);
[x,y] = meshgrid(azEdges,elEdges);
values = values([1:end end],[1:end end]);
valid = isfinite(values);
values(~valid) = limits(1);
surface(ax,x,y,zeros(size(x)),values,'EdgeColor','none','FaceColor','flat', ...
    'FaceAlpha','flat','AlphaData',double(valid),'AlphaDataMapping','none');
view(ax,2);
set(ax,'XDir','normal','YDir','normal','Color',[.90 .90 .90], ...
    'Layer','top','FontSize',11,'Box','on','TickDir','out');
xlim(ax,[azEdges(1) azEdges(end)]);
ylim(ax,[elEdges(1) elEdges(end)]);
colormap(ax,colors);
clim(ax,limits);
key = colorbar(ax);
key.Label.String = colorLabel;
xlabel(ax,['Global azimuth (' char(176) ')']);
ylabel(ax,['Global elevation (' char(176) ')']);
end

function edges = angularSINREdges(values)
values = values(:).';
if isscalar(values)
    edges = values+[-.5 .5];
    return
end
halfStep = diff(values)/2;
edges = [values(1)-halfStep(1),values(1:end-1)+halfStep,values(end)+halfStep(end)];
end

function closeAngularSINRHidden(fig)
if isgraphics(fig) && strcmp(fig.Visible,'off')
    close(fig);
end
end

# FIM 四种阵列形态上行对比：精简工程

本项目包含四种阵列的 MATLAB 仿真源码和一份合成方向图。`optional/` 中提供 PCA 无间距轻量入口。已在基础 MATLAB R2024b 中验证，无需 5G Toolbox。真实方向图和外部 CDL MAT 可按需指定。

## 1. 快速开始

1. 将 MATLAB 当前文件夹切换到包含 `configFIMSimulation.m` 的仓库根目录。
2. 运行以下自检，确认当前环境可以使用内置方向图、CDL 模板与四方案仿真。验证设置见 [VALIDATION_zh.md](VALIDATION_zh.md)。

```matlab
run_fim_validation
run_fim_cdl_validation
```

3. 在 `run_24rb_four_array_comparison.m` 顶部设置 `UE_CASE=1`、`MC_PER_CASE=1`、`SAMPLES_PER_RB=1`，然后运行主入口。

```matlab
run_24rb_four_array_comparison
```

`run_fim_comparison` 是简短别名。`PATTERN_TXT` 和 `CDL_PROFILE_MAT` 留空时，使用随仓库提供的 6.7 GHz 合成方向图和内置 CDL 参数表。用 `CDL_PROFILE` 选择 CDL-A 至 CDL-E。`CDL_MODEL='cluster'` 复现簇级研究，`CDL_MODEL='rays'` 启用标准子径。

## 2. 本包主入口使用的配置

主入口有意与 `configFIMSimulation.m` 的通用默认值分开：它会覆盖以下参数。

| 项目 | 主入口设置 |
|---|---|
| 阵列规模 | 16×8=128 阵元 |
| 阵元间距 | 水平0.5λ，竖直1.5λ |
| 载频 | 6.7 GHz |
| UE | 面积均匀撒点，水平距离20–200m，扇区角±60°，高度1.5m |
| OFDM | 24 RB×12子载波=288活动子载波；30kHz SCS；8.64MHz |
| FFT | **512点**；不能将 `nfft` 设为51 |
| 频率计算 | 默认每RB一个代表频点，共24个；可改为48或完整288个 |
| CP | CP/Tu=0.12，并计入速率 |
| CDL | 默认 CDL-C、100 ns 名义 RMS 时延扩展、簇级模式 |
| 发射功率 | 每UE总功率23dBm |
| 噪声 | 290K、7dB噪声系数 |
| 接收机/速率 | MMSE-IRC、Shannon频谱效率 |
| MC | K=1和K=10，每组默认20 MC，每批1 MC |
| 子阵 | 水平4块×竖直2块，共8个子阵；每块4×4阵元 |
| 子阵搜索 | 法向位移固定为0；yaw={-30°,0°,30°}，tilt={-20°,0°,20°} |
| PCA | 每阵元独立x/z移动；最小中心间距设为0 |
| 方向图重构 | 阵元位置固定；每阵元独立选择yaw/tilt状态 |

**阵列维度约定：** `cfg.Nh` 是水平方向阵元数，`cfg.Nv` 是竖直方向阵元数，总阵元数 `N=Nh*Nv`。本包主入口设为16×8=128。子阵划分参数 `cfg.subarray.m=4`、`cfg.subarray.n=2` 表示水平4块、竖直2块，不是每块阵元数；所以每块包含(16/4)×(8/2)=4×4个阵元。

`cfg.ofdm.nfft=512` 是FFT长度，不是活动子载波数。288个活动子载波映射在512点FFT的有效有符号频率索引范围内。`SAMPLES_PER_RB=1` 仅减少信道评估频点，实际RB数量、288个活动子载波、8.64MHz物理带宽和每UE总功率均不改变。

## 3. 正式输入

### 方向图 TXT

`patterns/Synthetic_6p7GHz_H120_V24.txt` 是解析合成的功率方向图。主瓣指向 `Phi=90°`、`Theta=90°`，水平和俯仰切面的半功率全宽分别为 120° 和 24°。使用真实方向图时，在 `PATTERN_TXT` 中填写文件路径。从通用配置直接运行时，设置 `cfg.pattern.file`。

格式要求：两行表头，随后360×181行、四列Phi / Theta / 功率增益(dB) / 相位；Phi先变化。计算实际使用功率方向图。方向图归一化、坐标对齐和格式检查由 `fimChannel.m` 完成。

### CDL

将 `CDL_PROFILE_MAT` 留空时，`generateCDLBase.m` 按 CDL-A 至 CDL-E 的标准参数表生成 `rawCDL`。CDL-A/B/C 描述非视距场景，CDL-D/E 包含视距分量。也可通过 `cfg.cdl.generator` 选择符合 `rawCDL` 接口的自定义生成器。[3GPP TR 38.901 V18.0.0 §7.7.1](https://www.etsi.org/deliver/etsi_tr/138900_138999/138901/18.00.00_60/tr_138901v180000p.pdf)

| `CDL_MODEL` | 信道实现 |
|---|---|
| `'cluster'` | 每簇一条标量等效路径，用于复现此前研究。 |
| `'rays'` | 每个散射簇按 20 条子径展开，随机耦合收发角度并生成极化初相。CDL-D/E 另有相干视距分量。四种阵列共用同一 UE 的子径及随机初相，各自的阵元坐标、朝向和方向图参与信道计算。 |

运行 CDL-D 子径示例时，可在 `run_24rb_four_array_comparison.m` 顶部设置以下选项，再运行入口。其他模板只需替换 `CDL_PROFILE`。

```matlab
UE_CASE = 1;
MC_PER_CASE = 1;
SAMPLES_PER_RB = 1;
CDL_PROFILE = 'CDL-D';
CDL_MODEL = 'rays';
CDL_PROFILE_MAT = '';
```

直接调用 `simulateFIM(cfg,'cdf')` 或 `runFIMStudy(cfg)` 时，在 `configFIMSimulation()` 返回的结构中设置 `cfg.cdl.profile` 和 `cfg.cdl.model`。子径模式可用以下选项。

| 配置 | 含义 |
|---|---|
| `cfg.cdl.extra_args` | 默认 `{'DelaySpread',100}`，时延扩展单位为 ns。 |
| `cfg.cdl.align_to_ue` | 默认按 UE 几何方向旋转散射角谱。设为 `false` 时按模板全局角度生成，仍可指定角度缩放。 |
| `cfg.cdl.rx_slant_deg`、`cfg.cdl.tx_slant_deg` | 接收阵元、UE 端口的极化倾角，单位为度，默认均为 0。 |
| `cfg.cdl.velocity_mps`、`cfg.cdl.time_s` | 所有 UE 共用的全局速度向量 `[vx vy vz]` 和快照时刻，单位分别为 m/s 和秒。默认零速度和零时刻。 |
| `cfg.cdl.k_factor_dB` | CDL-D/E 可指定视距分量相对所有散射分量的总 K 因子。留空使用表中功率。 |
| `cfg.cdl.angle_scaling_deg` | 可选 2×4 矩阵，列依次为 AoD、AoA、ZoD、ZoA，第一行是目标角度均值，第二行是目标 RMS 角扩展，单位为度。一列两个 `NaN` 表示保留该维参数。缩放在 UE 方向对齐前完成。 |
| `cfg.cdl.ray_coupling` | 默认 `'random'` 进行随机耦合。`'aligned'` 供固定子径配对核对。 |

若已有包含 `rawCDL` 或 `result.rawCDL` 的 MAT 文件，填写 `CDL_PROFILE_MAT`。输入应包含线性非负功率 `pow` 和秒单位时延 `taus`。子径模式要求 `name` 标识 CDL-A 至 CDL-E，按该标准表的行序提供 `AoA`、`AoD`、`EoA`、`EoD`（弧度），并保留 CDL-D/E 最前面的视距及同簇散射行。默认 `EoA/EoD` 为仰角，旧文件采用天顶角时设置 `cfg.cdl.eoa_convention='zenith'`。

#### 信道时间序列

使用 `fimChannel('series',realization,pattern,geom,cfg,time_s)` 计算同一次子径实现的多个频域快照。例如，按 1 ms 间隔采样 100 ms 的 CDL-D 信道。

```matlab
cfg = configFIMSimulation();
cfg.cdl.model = 'rays';
cfg.cdl.profile = 'CDL-D';
cfg.cdl.velocity_mps = [1.5 0 0];
cfg = prepareFIMConfig(cfg);

[pattern,profile] = fimChannel('prepare',cfg);
realization = fimChannel('realization',cfg,profile);
model = fimGeometry('prepare',cfg);
[geom,valid] = fimGeometry('build',model,model.zero,cfg);
assert(valid);
[H,info] = fimChannel('series',realization,pattern,geom,cfg,0:1e-3:0.1);
```

`H` 的四个维度依次为接收阵元、UE、实际评估频点和时间，形状为 `N × K × F × T`。`info.time_s` 为时间行向量。时间点可以非均匀，单位为秒，从该实现的随机初相时刻起算，取非负值。输出快照按输入时间点的次序排列。

整段序列共享阵列几何、UE 位置、角度、时延和路径损耗，子径按固定多普勒频移推进相位。默认零速度产生恒定信道。分析多普勒谱时，时间采样率应高于最大多普勒绝对值的两倍。长序列可用同一个 `realization` 分段计算，分段时继续使用相对于同一起点的时间。也可用 `fimChannel('time',realization,t_s)` 获取单个时刻的实现，再用 `fimChannel('geometry',snapshot,pattern,geom,cfg)` 计算其频域信道。

## 4. 首次运行和续算

首次建议在入口顶部修改：

```matlab
UE_CASE = 1;
MC_PER_CASE = 1;
SAMPLES_PER_RB = 1;
```

确认方向图和CDL输入可用、数值结果正常保存后，再增加MC数量。`SAMPLES_PER_RB` 支持1、2、3、4、6、12；12代表每RB 12个子载波，即完整288频点。

续算时，将 `RESUME_DIRECTORY` 设置为原任务目录，并保持UE数量、MC数量、批次大小、频点数、阵列、网格、方向图和CDL模板等配置一致。改变这些设置请启动新任务，避免混用检查点。

早期使用 `generateCDL5G` 的研究按原任务保存的 `rawCDL` 快照续算。

结果默认写入：

```text
results/four_array_24rb_lowcost/fim_comparison_时间戳/
```

主要输出包括 `study.mat`、`summary.csv`、`plan.mat`、`cdl_snapshot.mat`、逐批checkpoint、THP/SNR/SINR及增益CDF、参数表、四种阵列几何图和BS/UE部署图。

如果数值仿真已完成但图形导出失败，`study.mat` 和批次检查点仍会保留。可单独重绘：

```matlab
[figs,graphics] = plotFIMResults('完整路径/study.mat');
```

SYS/MLang等兼容环境的图形导出行为可能与标准MATLAB不同；首次运行建议先检查图形输出。

## 5. 文件结构

核心MATLAB文件：

| 文件 | 功能 |
|---|---|
| `run_24rb_four_array_comparison.m` | 24 RB四方案主入口 |
| `run_fim_comparison.m` | 主入口别名 |
| `configFIMSimulation.m` | 通用配置默认值 |
| `prepareFIMConfig.m` | 参数、资源、CP、功率和几何配置检查 |
| `runFIMStudy.m` | K扫描、MC分批、断点续算和统计汇总 |
| `simulateFIM.m` | 单个K下四方案配对仿真 |
| `fimScene.m` | BS/UE部署、坐标和路径损耗 |
| `fimChannel.m` | 方向图、CDL适配、随机实现和宽带物理信道 |
| `fimGeometry.m` | 固定、PCA、子阵及方向图重构几何 |
| `fimDeformation.m` | 子阵、PCA和方向图重构优化器 |
| `uplinkOFDMRate.m` | ZF/MMSE-IRC、SINR、频谱效率和吞吐率 |
| `plotFIMResults.m` | CDF、几何、部署图和结果重绘 |
| `reportFIMParameters.m` | 参数表TXT/CSV/PNG |
| `getAntennaPattern.m` | 方向图读取兼容接口 |
| `generateCDLBase.m`、`cdlProfileTable.m` | 基础 MATLAB CDL 簇模板生成器及标准参数表 |
| `expandCDLRays.m` | 子径角度展开、随机耦合及可选的角度、K 因子缩放 |
| `cdlRayCoefficients.m`、`cdlPolarizationField.m` | 子径复系数、多普勒和极化场投影 |
| `generateCDL5G.m` | 旧生成器名称的兼容入口 |
| `run_fim_validation.m` | 数值自检和集成验证 |
| `run_fim_cdl_validation.m` | 子径模式与五种 CDL 模板自检 |

`optional/` 中保留 `run_24rb_pca_no_spacing.m` 和其说明，作为PCA无间距轻量运行入口。需要时先执行 `addpath(fullfile(pwd,'optional'))`。核心主流程不依赖该可选脚本。

## 6. 模型边界

子径模式覆盖陆地 CDL-A 至 CDL-E 的频域链路快照，载频范围为 0.5–100 GHz，最大频率跨度为 2 GHz。每个接收阵元采用单馈端口，由 TXT 功率方向图和 3GPP Model-2 极化倾角计算电场。每个 UE 采用单个各向同性端口。模型加入子径多普勒相位，省略时域波形和由多普勒引起的载波间干扰。

工程以对数距离模型计算路径损耗，以完美 CSI 计算 ZF 或 MMSE-IRC，以 Shannon 公式估算频谱效率。模型省略 SRS/DMRS 和信道估计、MCS/BLER、实测复方向图相位、阵元互耦、实体碰撞、结构应变和执行器动态。

因此，轻量频点和粗搜索结果适合检查流程、观察趋势，不应直接作为完整精度的最终性能结论。

## 7. 本地文件与版本管理

仓库跟踪仿真源码、自检入口、说明和内置合成方向图。仿真输出保存在本地 `results/`。将额外的方向图 TXT、外部 CDL MAT 及本地实验数据放在 `local/`，并通过配置指向所需文件。这两个目录和临时导出文件由 `.gitignore` 排除。

旧检查点保存了绝对路径。迁移至其他目录或电脑后，续算前需核对原配置和输入路径。已有 `study.mat` 可以用 `plotFIMResults` 在新位置重绘。维护约定见 [AGENTS.md](AGENTS.md)。

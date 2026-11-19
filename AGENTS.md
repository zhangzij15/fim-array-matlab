# FIM 阵列仿真工程交接

## 使用环境和启动

- 在 MATLAB 中进入包含 `configFIMSimulation.m` 的仓库根目录。基础 MATLAB R2024b 已验证，合成方向图和标准 CDL 参数表随仓库提供。
- 运行 `run_fim_validation` 和 `run_fim_cdl_validation`。前者应通过 39 组自检，后者应通过 CDL-A 至 CDL-E、K=1/K=10 的四方案测试和时间序列检查。
- 首次试跑时，在 `run_24rb_four_array_comparison.m` 顶部设置 `UE_CASE=1`、`MC_PER_CASE=1`、`SAMPLES_PER_RB=1`，再运行 `run_fim_comparison`。用 `CDL_PROFILE` 选择剖面，用 `CDL_MODEL='rays'` 启用标准子径。`CDL_PROFILE_MAT=''` 使用内置标准模板。结果默认写入 `results/`。
- 完整用法见 [README.md](README.md)，验证记录见 [VALIDATION_zh.md](VALIDATION_zh.md)。`optional/` 提供 PCA 无间距约束的入口与说明。

## 代码入口和数据流

| 文件 | 职责 |
|---|---|
| `configFIMSimulation.m`、`prepareFIMConfig.m` | 公共默认值、参数检查、频率和资源网格 |
| `runFIMStudy.m`、`simulateFIM.m` | K/MC 批量、权重和续算，以及单批次四方案配对仿真 |
| `fimScene.m`、`fimGeometry.m` | BS/UE 部署、路径损耗、阵元及子阵几何 |
| `fimChannel.m` | 读取方向图和 CDL，生成共用的传播实现并构建宽带信道 |
| `generateCDLBase.m`、`cdlProfileTable.m`、`expandCDLRays.m` | 基础 MATLAB 标准 CDL-A 至 CDL-E 参数和子径 |
| `fimDeformation.m`、`uplinkOFDMRate.m` | 有限网格阵列优化、ZF/MMSE-IRC 接收和频谱效率 |
| `plotFIMResults.m`、`reportFIMParameters.m` | 重绘经验 CDF、部署与几何图，导出参数和统计 |

`cfg.cdl.model='cluster'` 用于复现簇级近似结果，`'rays'` 使用标准子径，支持交叉极化、D/E 视距分量和多普勒。可用 `cfg.cdl.profile_mat` 指定外部 `rawCDL` MAT。`generateCDL5G.m` 是旧名称的兼容入口，实际调用基础 MATLAB 实现。内置 6.7 GHz 合成方向图的水平和俯仰半功率宽度分别为 120° 和 24°。可用 `cfg.pattern.file` 或脚本中的 `PATTERN_TXT` 指定其他方向图。

`fimChannel('series',realization,pattern,geom,cfg,time_s)` 可复用同一次随机实现，生成 `N × K × F × T` 的时间序列。在创建 `realization` 前设置 `cfg.cdl.velocity_mps`，各时刻共享几何、路径和初相。示例见 [README.md](README.md#信道时间序列)。

## 对比试验的约定

- Fixed 的阵元位置和朝向固定。PCA 逐阵元在原阵面内移动，各阵元朝向固定。方向图重构固定阵元位置，只调整各阵元方向图的朝向。FIM 对刚性子阵做位姿调整，子阵内阵元坐标与朝向一起变换。
- 同一次 MC 的四种方案共享 UE、路径、时延、衰落、RB、功率及接收机配置，各自独立搜索。跨 K 的样本来自同一设定分布。复现同一 K 的逐场景配对时，使用同一个研究计划与批次种子。
- 主入口默认 24 RB，即 288 个活动子载波。若研究设置为 288 RB，则共有 3456 个活动子载波。在 `cfg = configFIMSimulation()` 后设 `cfg.ofdm.num_rb=288`、`cfg.ofdm.nfft=4096`、`cfg.ofdm.scs_Hz=30e3`，再调用 `runFIMStudy(cfg)`。分组模式每 RB 至少一个代表频点，可设 `cfg.ofdm.frequency_mode='grouped'` 和 `cfg.ofdm.num_frequency_samples=288`。
- 开始新研究时记录种子、CDL 模式、机械下倾、搜索网格、扫描轮数、MC 和实际评估频点。迁移历史结果时另行复制完整任务目录，续算前核对检查点保存的绝对路径及原配置。可用 `plotFIMResults` 读取 `study.mat` 重绘。

接收机在频域链路模型下采用完美 CSI，速率按 Shannon 公式计算。工程尚未建模 SRS、DMRS、时域波形和载波间干扰。性能结论以实际配对结果和足量 MC 为依据。

## 本地文件与 Git

- 仅跟踪源码、自检入口、说明和内置合成方向图。将生成结果留在 `results/`，将外部方向图、CDL MAT 等本地输入放在 `local/`，本地 ZIP 和历史 Julia 包装器也留在版本管理之外。
- 本地提交与远端操作分别授权。获得用户授权后再创建或连接远端、推送或设置自动推送。

# PCA无间距约束：可选轻量入口

本目录中的 `run_24rb_pca_no_spacing.m` 是主工程的可选入口，不是核心四方案主流程的依赖。先在工程根目录运行：

```matlab
addpath(fullfile(pwd,'optional'));
[study,figs,files] = run_24rb_pca_no_spacing();
```

支持 `smoke`、`light`、`balanced`、`full` 四种档位。默认 `light`：每个K最多5 MC、每RB一个代表频点、候选网格粗化为两端和零点；`smoke`最多1 MC；`balanced`最多20 MC、每RB两个代表频点；`full`沿用 `configFIMSimulation.m` 当前配置（本包默认20 MC/组、全288频点和配置中的完整搜索网格）。

该入口将PCA最小中心间距设为0，但仍保留PCA移动边界、固定阵元朝向和CP可行性检查。它会对四种阵列方案一起重算，不应把不同档位的结果直接拼接比较。

```matlab
run_24rb_pca_no_spacing([],[],'smoke');
run_24rb_pca_no_spacing([],[],'balanced');
```

续算时必须显式指定原任务目录，并保持原档位及其余配置不变。入口使用根目录 `configFIMSimulation.m` 配置的方向图和 CDL 生成器，默认使用合成方向图与 `generateCDLBase.m`。

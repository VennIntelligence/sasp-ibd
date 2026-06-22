# 通宵自主运行日志（Claude 主控，2026-06-23 夜 → 06-24 晚）

> 用户睡觉期间，Claude 自主：codex 实验结束→评估→据结果布置下一个验证/探索实验→长等待→循环。
> 所有产物进 git，每轮在此追加一条。**明早/明晚看本文件 + git 历史即可掌握全部。**
> 等待策略：后台等待器仅在每个实验**完成时**唤醒 Claude 一次（中途 10 分钟一次廉价 ssh 探测，不耗 token）；不做中途进度汇报。
> 主线（PIVOT）：难治性 IBD 炎症模块的**因果解剖**；当前因果锚点 = CCL8(血液,FinnGen复制) + CXCR2(血液,方向性)。

## 实验队列与决策树

- **#1 肠道 eQTL × 模块**（已完成）：GTEx colon 健康肠 → **0 因果命中**，28/38 基因无肠工具（免疫诱导基因不表达）。结论：组织 eQTL 是错工具，需细胞类型匹配 eQTL。
- **#2 多语境因果地图**（运行中）：全量 allpairs 升级肠 coloc + **免疫细胞 eQTL（单核/中性粒/刺激态）** × 模块。决定因果地图能否从 2 个扩大。
- **#3+ 据 #2 结果二选一**：
  - 若**因果地图扩大**（OSM/TREM1/IL13RA2 等新命中）→ 对新命中做：FinnGen R12 复制、reverse-MR（排反向因果）、MVMR 调炎症/CRP（排泛炎症）、Steiger filtering、MR 敏感性（Egger/weighted-median/LOO）、pQTL coloc（若有蛋白层验证更硬）。
  - 若**仍只有 CCL8/CXCR2** → 集中加固这两个：CXCR2 多 GWAS/FinnGen 复制、两者 reverse-MR + MVMR-炎症 + Steiger + 敏感性 + pQTL(deCODE/UKB-PPP) coloc；并探索把基因集扩到更广 IBD 炎症基因看方法能否多挖出因果点。
- **#N 通用补强（任意时机有空槽就做）**：单细胞反卷积 bulk 应答队列（哪类细胞丰度预测无应答，把单细胞接回应答）；UST/VDZ 队列拓宽描述性预测。

## codex brief 标准运行准则（今后每个 brief 都带）
1. **吃满主机 16 核/32 线程**：CPU 密集活用 `joblib.Parallel(n_jobs=30)`/`multiprocessing` 并行，别串行干等。
2. **耐心长等待、绝不中断**：下载/流式/API 限速可能耗时几小时是正常的；用大超时 + sleep 指数退避重试；让长命令跑到底，不因"慢"提前放弃或降级。
3. **GPU 安全（硬约束）**：默认纯 CPU。若必须上 GPU，**只准 `CUDA_VISIBLE_DEVICES=1`**；**绝不碰 0 号卡（已烧坏，碰了会主机断电/掉驱动）**。运行中 `nvidia-smi` 失败或 CUDA 驱动错 = 掉驱动 → **立即终止任务、写 STOP、不重试**，等用户次日处理。

## 运行记录

- **2026-06-23 ~00:55 JST** — #2 派出（tmux `cm2`，codex，纯 CPU）。
- **2026-06-23 ~01:1x JST** — 据用户要求给 #2 brief 加"运行准则"（并行+耐心）并**重启** `cm2`。
- **2026-06-23 ~01:3x JST** — #2 卡在 eQTL Catalogue API（按 gene_id 查回 400）。实测出可用法并写进 brief 重启：
  **associations 端点必须按区域查** `…/datasets/{QTD}/associations?pos={chr}:{start}-{end}&size=1000` 再过滤 molecular_trait_id；
  免疫数据集 id（已验证）：单核 QTD000021(BLUEPRINT)/QTD000504(DICE)、**中性粒 QTD000026(BLUEPRINT)**、刺激态 QTD000414(Quach LPS)。
- **2026-06-23 ~02:35 JST — #2 完成。** 多语境因果地图（`results/tables/module_causal_map_multicontext.tsv`，`src/26`）。
  **正面**：CXCR2 在**中性粒**独立确认因果保护（OR 0.851, FDR 1.05e-07, **coloc PP4 0.935**），加血液=双语境双确认。
  **阴性/路人**：CCL8 单核无可用工具(仍 blood-only)；OSM/OSMR/IL13RA2 无信号；**TREM1=路人**(中性粒 eQTL p=2.4e-27 但 MR OR 1.002 p=0.92)。
  因果地图未扩大（仍 CCL8+CXCR2）→ 走"加固核心 + 坐实路人"分支。技术债：肠 allpairs phenotype-id 映射没对上，肠 coloc 仍近似(低优先后补)。
- **2026-06-23 ~02:40 JST — #3 派出**（会话 `cm3`，纯 CPU）：CCL8/CXCR2 加固。
- **2026-06-23 ~03:20 JST — #3 完成**（`src/27`，`results/tables/{reverse_mr,steiger,mvmr_crp,finngen_cxcr2,mr_sensitivity,pqtl_ccl8}.tsv`）。
  **CXCR2 强**：FinnGen 复制 IBD OR 0.796(p1.3e-5)/UC 0.736(p1.6e-6)✓；多工具敏感性 IVW/Egger/wmedian/LOO 全一致~0.84 无多效性✓；MVMR-CRP 方向稳(IBD OR0.749,p1.1e-10, 但 CRP F弱)✓；Steiger✓。**疑点**：单SNP反向MR IBD→CXCR2 p=3e-200（疑共定位伪反向）。
  **CCL8 部分**：血液+FinnGen+Steiger✓；MVMR不可估(1工具)、pQTL被deCODE墙挡、反向MR同疑点。
  **路人坐实**：OSM/OSMR/TREM1/IL13RA2 全 bystander。
- **2026-06-23 ~03:25 JST — #4 派出**（会话 `cm4`，纯 CPU）：补缺口——正规双向反向MR(排cis区澄清p=3e-200)+CCL8多工具MVMR/敏感性+换可及pQTL源(SCALLOP)。等待中。

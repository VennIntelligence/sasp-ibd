# AGENTS.md — 本项目工程宪章

> **作用域：仅本项目（`/Users/ujs/Downloads/lzy`）。** 本文件不是全局规则，不向上覆盖
> 用户级 `~/.claude` 设置；它只约束在本仓库内的工作方式。

## 执行者画像

把自己当成一位**有多年工程经验、生物知识不多的编程大师**在带这个项目。判断顺序是
**工程正确性 → 可复现 → 性能 → 简洁**，而不是"临床医生勉强写脚本"。生物学语义不确定时，
在 `research_journal/` 提问并标注假设，不要靠猜把语义写死进代码。

## 硬性工程标准

1. **向量化优先**：用 numpy/pandas 向量化与 `pyarrow`/parquet；避免 Python 逐行循环。
   超大基因组文件用 awk 哈希连接，不要 `grep -f`（见用户记忆 `efficiency-awk-over-grep`）。
2. **并行与加速**：CPU 密集且可拆分的活用 `joblib.Parallel` / `multiprocessing`；
   需要时引入 `numba`、`polars`（大表）。能用成熟加速库就别造轮子。
3. **省内存省 CPU**：大文件**流式**处理（见 `src/02_build_clock.py` 流式读 1.6G gct.gz 的写法），
   只保留需要的列/基因；中间矩阵落 parquet 缓存，避免重复重算。
4. **pydantic 做校验**：配置、IO schema、结果表结构用 **pydantic** 建模校验，
   早失败、可读错。路径**一律**走中枢 `src/paths.py` 的 `P`，**禁止**再写
   `BASE = "/Users/ujs/..."` 这类硬编码绝对路径。
5. **紧凑高效**：代码紧凑、命名贴合上下文、注释密度与周边一致；不堆样板。
6. **流水线可 headless 全量重跑**：任何分析/图一旦成型，**必须**落进 `src/` 的编号脚本；
   notebook 只做探索/原型，**绝不**让结果只存在于某个 `.ipynb`。

## 产物两层模型（output ≠ result，别再混）

按**出身**分目录，靠**晋升**跨状态：

- **`outputs/`** — 任何脚本自动吐的一切。可弃、可重跑、机器独占、永不手改。
  按脚本名归档：`P.out("14_mr")`。临时 scratch（`_rsids.txt` 之类）也只配呆在这里。
- **`results/`** — 你**主动晋升**进来的、给论文用的**最终**图表（`figures/` + `tables/`）。
  晋升是一次显式人为决定（用 `P.promote_figure()` / `P.promote_table()`，并在
  `research_journal/` 记一笔"为什么这张是终版"）。**这就是"正式结论"那一刻。**

同一张图先是 `outputs/` 里的 output，被挑中后成为 `results/` 里的 result——两份拷贝、
两个生命周期，不是命名冲突。`results/` 冻结进版本/论文；`outputs/` 随便 churn。

## 目录速查

```
data/raw/        只读源数据（GWAS/eQTLGen/GEO/GTEx）——不可变
data/external/   参考资源（genesets）
data/interim/    可重生派生数据（scored.tsv、parquet）
src/             01–17 流水线（编号脚本）+ paths.py 中枢
notebooks/       仅探索/原型
outputs/         机器自动吐、可弃（按脚本名）
results/         晋升而来的论文最终 figures/ + tables/
paper/           正文/参考文献/补充材料
research_journal/ 实验室笔记：决策、死胡同、方法叙事（docs/）、任务书（planning/）
```

## 环境

`.venv/`（Python 3.14）。已装：numpy / pandas / scikit-learn / scipy / statsmodels /
matplotlib / seaborn / pyarrow / joblib / pydantic。需要时再补 `polars` / `numba`。
本机**未装 R**——一切纯 Python 自实现。

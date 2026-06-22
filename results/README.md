# results/ — 论文最终图表（晋升而来）

只放你**主动挑中**为终版的 `figures/` 与 `tables/`。这里的每个文件都是一次显式人为
决定的产物（"这张/这表是正式结论"），从 `../outputs/` 拷贝晋升而来，应当冻结、进论文。

- `figures/` 主图 + 因果整合图 + 时钟精度图 + MR 森林图 + 总览图
- `tables/` instruments / mr_{IBD,CD,UC} / coloc_IBD / triangulation

晋升用 `P.promote_figure()` / `P.promote_table()`，并在 `../research_journal/` 记一笔
"为什么这张是终版"。机器随手吐的东西不要直接丢进来——那是 `../outputs/` 的事。

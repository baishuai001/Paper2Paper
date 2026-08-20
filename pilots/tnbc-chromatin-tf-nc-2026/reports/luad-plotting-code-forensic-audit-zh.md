# LUAD 图形实现取证审计与重建决定

日期：2026-08-20  
状态：旧版成图撤回；正在按 TNBC Code Ocean v1（commit `edf5314`）重建

## 结论

用户指出的“字号、配色廉价，Figure 3B/3D 像涂鸦”成立。问题不是 LUAD
数据天然只能画成这样，也不是 TNBC 官方代码无法复用，而是旧实现走错了路线：

1. 统计分析与可视化被分开重写，只复用了部分阈值和输入概念，没有逐段移植
   TNBC 官方绘图构造器；
2. 多个面板改变了原文的统计对象或信息粒度；
3. 所有原子图又被塞入固定 `13 × 10 in` 横向报告模板，附加页眉、副标题和页脚，
   造成关键面板被二次缩小；
4. 一个通用 `theme_pub(base_size = 8.2)` 被用于性质不同的热图、森林图、网络图和
   药敏图，抹掉了 TNBC 官方针对每类图形设置的视觉层级；
5. 动态 HCL 色板、过饱和模块色、圆角卡片和 UI 风格说明框替代了原文固定、低饱和、
   全篇一致的颜色语义。

因此，旧版 Figure 1–5 和 Supplementary Figure 1–9 不能作为投稿版成图，必须先恢复
科学对象，再恢复官方视觉语法，最后重新拼版。

## 两处需要立即撤回的科学性错误

### Figure 3B

TNBC 官方代码只把满足“靶基因对共享至少 3 个 HC-TF”的边及其端点交给
`graph_from_data_frame()`。LUAD 在同一规则下只有 3 条边和 6 个端点：

| 靶基因对 | 共同 HC-TF |
|---|---|
| NDNF–ADGRF5 | CSRNP1、ETV1、NR3C2 |
| ADGRD1–SELENBP1 | CRY2、MXD4、ZBTB18 |
| YPEL3–CACFD1 | MAGED2、MXD4、ZNF444 |

旧代码却先把所有“被至少 3 个 TF 调控”的 30 个靶基因作为固定 vertices，再加 3 条边，
因而额外显示了 24 个原文代码不会显示的无边节点。Fruchterman–Reingold 布局把这些
节点随机撒满画布，再强制标注全部名称，直接制造了“散点涂鸦”。

修订决定：审计版严格显示 6 节点/3 边；解释版把 3 个 dyad 做成紧凑小倍图并列出
共同 TF。24 个无边候选只能进入补充表，不能继续出现在主图网络中。结果应解释为：
在锚点阈值下未形成高阶靶基因群落，而不是把随机散开的孤点当成群落结构。

### Figure 3D

LUAD 的完整 TF 协作网络共有 31 个节点、134 条正边、1 个连通分量、0 个孤立节点，
密度为 0.288。旧代码使用边权 80% 分位数作为未预注册的显示阈值，仅留下 33 条边，
删除 101/134（75.4%）有效边，制造出 11 个组分和 10 个视觉孤立点。图中没有披露该
裁边规则；节点大小还错误映射为 regulon 靶基因数，而 TNBC 官方代码映射网络 degree。

修订决定：直接使用全部 134 条正边；节点大小=degree，节点统一暗红 `#981111`，
灰色半透明边，边宽=共享激活靶基因数，标签 repel，同时给出节点 degree 与边权图例。
这会与 Figure 3E 的 hub 定义保持一致。真实 degree 最高的 TF 为 CRY2（19）、
MAGED2（18）、ZNF19（16）、NR3C2 和 ZBTB18（各 13）。

## 其他主图中发现的对象偏移

- Figure 1B–E：官方使用逐 TF 行标准化、ComplexHeatmap 和固定临床/分子注释色；
  旧版存在动态色板、跨面板颜色语义不一致和部分函数未在绘图内执行行标准化的问题。
- Supplementary Figure 2B：TNBC 为两组 Euler/Venn，旧版换成了分面柱状图。
- Supplementary Figure 2F：TNBC 为全部发现 TF 在四个系统中的逐样本 NES 点分布，
  旧版换成了 cohort-effect 热图。
- Figure 2B/2D：TNBC 使用 ComplexUpset；旧版是自制 UpSet，motif 数据库层信息丢失。
- Figure 2C：TNBC 展示全部发现 TF 的 sample-level 启动子、NES 和 motif 证据；旧版压缩为
  31 个 HC-TF × 3 个系统摘要，不能作为同构主图。
- Figure 2E：TNBC 展示 mean LOR、跨样本支持比例及 motif case examples；旧版换成了
  primary/sensitivity prevalence 气泡图。
- Figure 5E：TNBC 是 OS × RFS 的 3×3 预后分类计数图；旧版画成 replicated drug-effect
  heatmap，生物学对象不相同。
- Figure 5H：TNBC 是已验证 drug–TF 对的 PDX waterfall。LUAD 目前没有通过冻结规则的
  PDX 验证对，因此不能用其他图形冒充；只能明确标记“不支持生成等价面板”。
- Supplementary Figure 5–7：旧版分别把 skewness、相关结构、成对置换、Euler/散点和
  pathway–drug 矩阵换成了通用热图/UpSet/巨型气泡图，必须恢复原图语法。

## 重建规范

1. 以 TNBC 官方代码块为骨架，LUAD wrapper 仅适配列名、队列名和 endpoint；不得重新发明
   图形类型或显示阈值。
2. 使用固定锚点色板。LUAD/LUSC 对应官方 TNBC/non-TNBC 的 `#435773/#9bc1e5`；
   患者/PDX/细胞系为 `#7b007c/#6c8438/#147574`；风险/保护为
   `#A60311/#04588C`。
3. 服务器没有 Helvetica 授权字体；全篇统一使用度量兼容的 Liberation Sans，避免 Cairo
   自动混入 Nimbus Roman。若以后提供合法 Helvetica/Arial，可通过环境变量统一替换。
4. 原子图和主图均输出 Cairo 矢量 PDF；KM 图不得再使用 PNG 回退。
5. 取消固定 13×10 英寸横向报告模板、页眉、副标题和页脚。主图采用接近期刊版心的纵向
   尺寸和 TNBC 的非等宽面板层级，关键网络获得足够空间。
6. 最终 panel letter 为 14–16 pt；最终可见标签不得低于约 5.5–6 pt。
7. 没有数据支持的阳性面板必须明确显示负结果或不支持，不得为了视觉对齐降低阈值或
   伪造等价证据。
8. 每个 PDF 必须经过逐页渲染、字体清单、裁切、重叠、空白与最小字号检查；对照图册的
   TNBC crop 必须重新校准。

## 代码路线

- 共享锚点样式：`code/publication_rebuild/tnbc_anchor_theme.R`
- Figure 1–2 / Supplementary 1–4：独立官方风格 wrapper
- Figure 3：独立官方网络 wrapper，并修订协议与 Figure Legend
- Figure 4–5 / Supplementary 5–7：独立官方风格 wrapper
- 最终拼版：仅使用矢量面板，取消报告式装饰；逐图保留 TNBC 的视觉层级

旧版仅作为审计历史保留，不再被称为投稿版或严格复刻版。

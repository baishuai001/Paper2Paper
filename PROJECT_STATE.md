# Paper2Paper current state

Updated: 2026-08-08

## Product boundary

Paper2Paper is an alpha workflow for a beginner and AI to audit a real anchor
paper, compare bounded imitation routes, verify data and code with real inputs,
and create a separate formal manuscript project only after user approval. Its
success criterion is a scientifically defensible paper package, not continual
platform expansion.

The current schema is 2.0.0 and the package version is 0.5.0-alpha.1. Core,
reusable modules, real-paper Pilots and formal manuscript projects are separate
layers. Local issue resolution, module maturity and cross-paper rule promotion
are separate records.

## Canonical server workspace

The server-local project working directory follows:

`/media/desk16/<server-account>/projects/Paper2Paper`

The concrete account path is kept outside this public repository. Credentials,
hostnames and tokens must never be stored here.

The core workflow currently passes its repository suite on the canonical
server's Python 3.10 runtime. The gastric scientific module was reference-run
under its recorded Python 3.12 environment; its exact lock is not currently
resolvable from the server package index because `lifelines==0.30.3` is absent.
The failed isolated install was not replaced with a nearby version, so the
server must not be described as module-verified until the exact artifacts are
available or a new environment release is rerun and reviewed.

## Real-paper Pilots

- `P2P-GASTRIC-NRRS`: a bounded training/reference reconstruction. GN-R01 has
  a passed 300-patient real-data run and a separate focused-test run. It is not
  a manuscript candidate and does not test other modalities, other papers or
  the whole workflow. It is awaiting a user retain/close decision; none is inferred.

## Reusable evidence

- Resource integrity and patient-ID linkage each have one real-paper reference;
  they are not cross-paper confirmed.
- The fixed bulk-signature implementation is reference-verified on GN-R01 and
  has focused unit/adversarial coverage; it is not transfer-verified.
- Frozen GN-R01 values are exact-input change detectors only.
- The dedicated negative identity-linkage case is registered but not run; that
  promotion remains observed rather than being overstated as provisional.
- The gastric patient-level derived table remains local-only pending a
  redistribution review. A tracked receipt records its byte size, SHA256,
  dimensions and columns; repository-visible claims point to tracked aggregate
  source tables instead of depending on that ignored file.

## Formal manuscript projects

None exists. A formal project must be created with `paper2paper promote` after
a manuscript-candidate route has complete evidence review and explicit user
approval. Pilot results are not copied as formal manuscript results.

In this alpha, `run_kind=real_data` remains a declared provenance label rather
than a cryptographic proof of dataset identity. Reviewers must compare each run
manifest with the registered candidate, resource and cohort-usage records.
Likewise, registry case-to-run and assertion records are auditable declarations,
not cryptographic proof that every asserted property was executed by the bound
run. Automated checks enforce identities, statuses, real-data command and
manifest equality, independent anchor/dataset keys and declared unit-test
coverage; a reviewer must still inspect native logs and artifacts before
accepting transfer evidence.

## Immediate boundary

The schema/registry migration is complete. Keep the automated checks green and
use the canonical server directory for continued work. No active target-paper
Pilot or formal manuscript project currently exists. Do not add generic
features until a newly approved real-paper Pilot exposes a concrete blocker.

# Paper2Paper 当前状态

更新日期：2026-08-08

## 产品边界

Paper2Paper 是一个处于 Alpha 阶段的工作流，面向初学者并由 AI 协助，用于审计真实的锚点论文、比较范围明确的模仿路线、使用真实输入验证数据与代码，并且只在用户批准后创建独立的正式论文项目。

它的成功标准是形成一套在科学上经得起论证的论文材料，而不是持续扩建平台。

当前 schema 版本为 `2.0.0`，软件包版本为 `0.5.0-alpha.1`。核心工作流、可复用模块、真实论文 Pilot 和正式论文项目属于相互分离的层次。局部问题是否解决、模块成熟度，以及跨论文规则是否晋升，分别使用独立记录管理。

## 规范服务器工作目录

服务器上的项目工作目录遵循以下形式：

`/media/desk16/<服务器账号>/projects/Paper2Paper`

具体账号路径保存在该公开仓库之外。登录凭据、主机名和令牌绝不能存入这里。

目前，核心工作流已经在规范服务器的 Python 3.10 运行环境中通过仓库测试套件。

胃癌科学模块是在其记录的 Python 3.12 环境中完成参考案例运行的。由于服务器当前的软件包源中缺少 `lifelines==0.30.3`，该模块的精确锁定环境目前无法在服务器上解析和安装。

失败的隔离安装没有通过替换成相近版本来绕过。因此，在获得精确依赖安装文件，或者建立新的环境 release 并重新运行、审查之前，不能把服务器描述为“胃癌模块已验证”。

## 真实论文 Pilot

- `P2P-GASTRIC-NRRS`：一个范围受限的训练/参考重建项目。GN-R01 已完成一次包含300名患者的真实数据运行，并且另有一次独立记录的聚焦测试运行。它不是正式论文候选路线，也没有检验其他模态、其他论文或整个 Paper2Paper 工作流。目前仍在等待用户决定将其保留为训练案例还是关闭；系统没有擅自推断用户决定。

## 可复用证据

- 资源完整性检查和患者 ID 连接各自拥有一个真实论文参考案例，但尚未得到跨论文确认。

- 固定 bulk signature 实现已经在 GN-R01 中达到“参考案例验证”级别，并具有聚焦的单元测试和反向测试覆盖；但尚未达到“独立迁移验证”级别。

- GN-R01 的冻结数值只能用于检测相同输入条件下是否发生变化。

- 专门针对患者身份连接错误的反向案例已经登记，但尚未实际运行。因此，相应规则的晋升成熟度仍保持为 `observed`，不能夸大为 `provisional`。

- 胃癌患者级衍生表在完成再分发审查前，只保存在本地。一个已纳入 Git 的凭证文件记录了该衍生表的字节数、SHA256、行列维度和字段名称。仓库中可以查看的主张指向已纳入 Git 的汇总 source table，而不是依赖这个被 Git 忽略的患者级文件。

## 正式论文项目

目前不存在正式论文项目。

只有当一条 `manuscript_candidate` 路线完成全部证据审查，并且获得用户明确批准之后，才能通过 `paper2paper promote` 创建正式项目。Pilot 的结果不会被直接复制并冒充正式论文结果。

在当前 Alpha 版本中，`run_kind=real_data` 仍然只是一个申报的数据来源标签，并不是对数据集身份的密码学证明。审查者必须将每次运行的 manifest 与已登记的数据候选、资源记录和队列用途记录逐项比较。

同样，registry 中案例与运行之间的关联记录，以及 assertion 记录，都是可以审计的声明；它们并不能从密码学层面证明，被登记的运行确实执行了每一项所声明的性质检查。

自动化检查目前能够强制校验：

- 各类 ID 和引用关系；
- 状态是否合法；
- 真实数据案例的命令和 manifest 是否与对应运行一致；
- 锚点论文和数据集标识是否满足独立性要求；
- 声明的单元测试是否具有相应覆盖记录。

但是，在接受独立迁移证据之前，审查者仍然必须检查原始运行日志和实际产物。

## 当前工作边界

Schema 和 registry 迁移已经完成。后续工作应保持所有自动化检查持续通过，并使用规范服务器目录继续推进。

目前没有活动中的目标论文 Pilot，也没有正式论文项目。在用户批准新的真实论文 Pilot 并由其暴露出具体阻断问题之前，不应继续增加通用功能。

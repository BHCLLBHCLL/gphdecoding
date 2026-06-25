# GPH 二进制格式描述（逆向工程）

## 1. 概述

GPH 是一种几何/多面体网格文件格式，可能来自 SCTpre 或 CGNS 相关 CFD 工具。

- **魔数标识**: `CRDL-FLD`（8 字节 ASCII）
- **字节序**: **大端序**（Big-Endian），所有多字节整数和浮点数
- **对齐**: 4 字节对齐

## 2. 文件布局

| 偏移 (hex) | 大小 | 节名称 | 描述 |
|------------|------|--------|------|
| 0x0000-0x001C | 28 B | 文件头 | CRDL-FLD 标识 + 维度 |
| 0x001C-0x0078 | 92 B | FileRevision | I4 标量（如 2025） |
| 0x0078-0x00D8 | 96 B | Application | C1[8]，如 "SCTpre" |
| 0x00D8-0x0134 | 92 B | ApplicationVersion | I4 |
| 0x0134-0x0190 | 92 B | ReleaseDate | 日期字符串 |
| 0x0190-0x01EC | 92 B | GridType | 字符串 |
| 0x01EC-0x0248 | 92 B | Dimension | I4（空间维度） |
| 0x0248-0x02A4 | 92 B | Bias | I4 |
| 0x02A4-0x0300 | 92 B | Date | 日期字符串 |
| 0x0300-0x03A8 | 168 B | Comments | 注释字符串 |
| 0x03A8-0x04E0 | 312 B | Cycle | I4 + 单位信息 |
| 0x04E0-0x0560 | 128 B | Unused | 保留 |
| 0x0560-0x05D8 | 120 B | Encoding | 如 "UTF-8" |
| 0x05D8-0x0600 | 40 B | HeaderDataEnd | 头部结束标记 |
| 0x0600-0x0628 | 40 B | OverlapStart_0 | 重叠区开始 |
| 0x0628-0x08DC | 692 B | LS_CvolIdOfElements | I4[n_cells] 每单元 cvol_id（见 §5.1）|
| 0x08DC-0x26B0 | 变长 | LS_Links | 面拓扑（owner/neighbor/npe/conn，§5.2）|
| 0x26B0-0x2C60 | 变长 | LS_Nodes | R4/R8[n,3] 顶点坐标（三轴块，§5.3）|
| … | 变长 | LS_SurfaceRegions | 命名边界面区域（§5.4）|
| … | 变长 | LS_SolverUnusedRegions | 求解器内部区域名（可选）|
| … | 变长 | LS_VolumeRegions | 体区域名 → CGNS Zone（§5.5）|
| … | 变长 | LS_Parts | Part 名 + cvol_id 描述符（§5.6）|
| … | 变长 | LS_Assemblies | XML 装配树 → Zone 命名（§5.7）|
| 0x42B0-0x45A4 | 756 B | Element_InformationFlag | 单元标志 |
| 0x45A4-EOF | - | OverlapEnd | 文件尾 |

> **注意**：上表偏移仅适用于 `box.gph` 等早期样例；`tr03.gph`、`laptop_*.gph` 等文件体积更大，须用 `gph_parser.py` 或 `gph_model.find_section()` **动态定位**各节，不可硬编码偏移。

## 3. 记录格式

每个命名字段遵循以下模式：

```
[长度 I4] = 0x20 (32)
[名称 C1[32]] = 32 字节 ASCII，空格填充
[长度 I4] = 0x20 (32)   ; 可选，节标记
[描述符 变长]
  - 0x0C 0x04 [dim0] [dim1]  : 维度，类型 04 = I4
  - 0x0C 0x08 [dim0] [dim1]  : 类型 08 = R4 或 I8
[值 / 数组数据]
```

### 数据类型代码（推断）

| 代码 | 类型 | 说明 |
|------|------|------|
| 0x04 | I4 | 32 位有符号整数 |
| 0x08 | R4 / I8 | 32 位浮点或 64 位整数（视上下文） |

## 4. 元数据字段

| 字段 | 类型 | 示例 |
|------|------|------|
| FileRevision | I4 | 2025 |
| Application | C1[8] | "SCTpre" |
| ApplicationVersion | I4 | 1 |
| ReleaseDate | string | 日期 |
| GridType | string | 网格类型 |
| Dimension | I4 | 3 |
| Bias | I4 | - |
| Date | string | - |
| Comments | string | - |
| Cycle | I4 | 仿真步 |
| Encoding | C1 | "UTF-8" |

## 5. 数据数组

| 数组名 | 类型 | 说明 |
|--------|------|------|
| LS_CvolIdOfElements | I4[n_cells] | 每单元的 **cvol_id**（与 LS_Parts 描述符对应，§5.1）|
| LS_Links | I4[] | 面拓扑（owner / neighbor / npe / [face_type] / conn，§5.2）|
| LS_Nodes | R4/R8[n,3] | 顶点坐标（三轴 float32 或 float64 块，§5.3）|
| LS_SurfaceRegions | 变长 | 命名边界面区域 → CGNS ZoneBC（§5.4）|
| LS_SolverUnusedRegions | 变长 | 求解器内部区域名（可选）|
| LS_VolumeRegions | 变长 | 体区域名列表 → CGNS Zone（§5.5）|
| LS_Parts | 变长 | Part 定义（名 + cvol_id，§5.6）|
| LS_Assemblies | XML | 装配层次 → Zone 命名规则（§5.7）|
| Element_InformationFlag | 变长 | 单元信息标志 |

### 5.1 LS_CvolIdOfElements

每个体单元一条 I4，值为该单元所属 Part 的 **cvol_id**（scFLOW 内部分区标识）。

- **不是** `LS_Parts` 列表中的 1-based 序号。
- `tr03.gph`：63 882 单元，cvol_id ∈ `{1, 2}` 对应 Case / Rotate。
- `laptop_simplified_voxel_less.gph`：cvol_id ∈ `{1, 9, 11}`，与 Part 列表顺序无关。
- `laptop_simplified_more_regions.gph`：173 094 单元，**70** 个不同 cvol_id；多数几何子块 id 由复合 Part `air_domain` 的成员列表统辖（§5.6）。

**大文件前缀块**：`laptop_simplified_voxel_v4.gph` 等超大网格的 `LS_CvolIdOfElements` 节在真正的 `I4[n_cells]` 数组之前可能有一个 **4 字节的 metadata 块**（值为 cell 总数，如 `27553410`）。解析时必须取节内**最大的** I4 数据块，而非第一个块——误读会导致所有 Zone 的 cell 掩码均为全网格（进而 ZoneBC 的 PointList 也全部错误）。

`gph_model.parse_ls_parts()` 扫描每个 Part 名块后的描述符链，以 `LS_CvolIdOfElements` 的唯一值集合为权威来源建立 `part_name → cvol_id` 映射，再与 per-cell 数组匹配生成各 Zone 的 cell 子集（见 §5.6）。

### 5.2 LS_Links（含多面体网格）

按数据块顺序（经通用块扫描器定位）：

| 块 | 大小 | 含义 |
|----|------|------|
| owner | n_faces × I4 | 面的 owner 单元（0-indexed）|
| neighbor | n_faces × I4 | 对侧单元；`0xFFFFFFFF` = 边界面 |
| npe | n_faces × I4 | 每面节点数（纯三角时全为 3；`tr03` 为 3..11）|
| face_type | 1 × I4 | **仅部分旧文件**：单元类型标志 |
| conn | sum(npe) × I4 | 0-based 顶点索引，**CSR 布局** |

CSR 索引：`face_offsets[0]=0`，`face_offsets[i+1]=face_offsets[i]+npe[i]`；面 `i` 的节点为 `conn[face_offsets[i]:face_offsets[i+1]]`。

#### 超大 conn 分块（>~1 GiB）

当 `sum(npe) × 4` 超过单个数据块约 **1 GiB**（`1073741824` 字节）时，scFLOW 将 `conn` 拆成**多段**：

```
[12, bc1][conn 第 1 段 payload][bc1]   ← 标准块（bc1 常为 1073741824）
[I4=1073741824][conn 第 2 段 raw...]  ← 裸 byte_count + 数据，无 [I4=12] 头
... 可重复多个完整 1 GiB 裸块 ...
[I4=1073741824][I4=bcN][conn 末段 raw...]  ← 末段：重复 1 GiB 标记 + 实际 payload 字节数 + 数据
```

两段 conn（如 `laptop_simplified_voxel_v4.gph`）的续接块为单头裸块 ``[I4=bcN][payload]``（``bcN ==`` 剩余字节数）。**三段及以上**时，末段使用双头格式（见上行）；若误将 ``bcN`` 当作首个顶点索引读入，会在第二段 1 GiB 边界出现单点索引等于末段字节数（如 ``tests/box.gph`` 的 ``173846272``），导致面翘曲。

**conn 块选择**：若无块字节数恰好等于 `sum(npe)×4`，取除 owner/neighbor/npe 三数组外 **byte_count 最大且 ≥ 12** 的 I4 块（勿用 `3×n_faces×4` 作下界——多面体网格首段 conn 常被 cap 在 1 GiB）。

实测 `laptop_simplified_voxel_v4.gph`（~3.7 GiB，**2 段** conn）：

| 项目 | 数值 |
|------|------|
| 面数 | 88,833,031 |
| 单元数 | 27,553,410 |
| `sum(npe)` | 360,934,738 |
| conn 第 1 段 | 268,435,456 条目（1 GiB） |
| conn 第 2 段 | 92,499,282 条目（~370 MiB） |

实测 `laptop_simplified_denser_v2_gph.gph`（~5.9 GiB，**3 段** conn）：

| 项目 | 数值 |
|------|------|
| 面数 | 114,039,102 |
| 单元数 | 20,687,038 |
| 顶点数 | 83,664,081 |
| `sum(npe)` | 549,000,094 |
| conn 第 1 段 | 268,435,456 条目（1 GiB，标准块） |
| conn 第 2 段 | 268,435,456 条目（1 GiB，裸续接） |
| conn 第 3 段 | 12,129,182 条目（~48 MiB，末段） |

实测 `laptop_simplified_voxel_v6.gph`（~4.9 GiB，**2 段** conn）：

| 项目 | 数值 |
|------|------|
| 面数 | 126,318,473 |
| 单元数 | 38,895,916 |
| 顶点数 | 48,526,564 |
| `sum(npe)` | 513,041,554 |
| conn 第 1 段 | 268,435,456 条目（1 GiB，标准块） |
| conn 第 2 段 | 244,606,098 条目（~954 MiB，裸续接） |

`gph2cgns.py` 从 `gph_model._read_conn_continuations` 导入续接逻辑（返回 `(got, pos, n_continuations)` 三元组）。`gph2cgns.py` / `gph_model.py` 在主块之后循环读取裸 `byte_count` 续接块（亦支持标准 `[12,bc]` 块）并拼接后再做 CSR 索引。`gph_parser.py` / `gphviewer.py` 通过 `parse_ls_links_summary` 报告 `conn_got`、`conn_chunks`、`conn_complete`。**旧版解析器在 conn 块选错或续接不完整时直接失败**（报 `LS_Links parse failed`）。

纯三角 legacy 文件可能使用列主序 conn；`gph2cgns` 根据块大小与 `sum(npe)` 自动判别。

### 5.3 LS_Nodes

节点坐标按 X / Y / Z 三个等长轴块存储（经块扫描取三个最大且等长的块）。描述符 `type` 字段区分编码：**4 = float32（R4）**，**8 = float64（R8）**。FPH 求解结果文件（如 `tests/tr03_9.fph`）使用 float32；ANSA / 多数 GPH 导出使用 float64。

#### float64 块布局（标准）

```
[16B 描述符] 00 00 00 0C / 00 00 00 08 / n_verts / 00 00 00 01
[8B  块头]   00 00 00 0C / byte_count
[n_verts × 8B]  标准大端 IEEE-754 float64
[4B  尾部]   byte_count（哨兵，等于块字节数）
```

- 轴顺序：**X, Y, Z**（与 CGNS 一致，无需重排）
- 每块末尾有 4 字节的 `byte_count` 哨兵

#### float32 块布局（FPH）

```
[16B 描述符] 00 00 00 0C / 00 00 00 04 / n_verts / 00 00 00 01
[8B  块头]   00 00 00 0C / byte_count   (byte_count = n_verts × 4)
[n_verts × 4B]  标准大端 IEEE-754 float32
[4B  尾部]   byte_count（哨兵）
```

解析时升宽为 float64 写入 CGNS。若误将 float32 payload 按 float64 读取，坐标幅值会落在 ~1e-13 量级（denormal），旧版 ``1e-30`` 阈值会误判为“合理”并导致顶点数减半（``n = bc // 8``）。

示例（来自 `box_ansa.gph`）：坐标值 `float32(0.01)` 加宽为 float64 = `0x3F847AE140000000`，在文件中按标准大端写为 `3F 84 7A E1 40 00 00 00`。

示例（来自旧 `box.gph`）：坐标值 0.01（精确）= `0x3F847AE147AE147B`，在文件中按标准大端写为 `3F 84 7A E1 47 AE 14 7B`。

#### 历史"字反转 float64"误读

原始 `gph2cgns.py` 跳过了 12 字节块头（而正确应该跳过 8 字节），并使用字反转读法。对**特定**字节序列，这两个错误恰好相互抵消，使少数顶点产生看似合理的数值——但大多数顶点会得到 1e37 等夸张幅值。当前实现已修正为正确的"跳 8 字节 + 标准大端"读法。

#### 防御性方言自检（`parse_ls_nodes_xyz`）

`gph_model.parse_ls_nodes_xyz()` 为 `gph2cgns.py` / `fph2cgns.py` / `gph_parser.py` / `gphviewer.py` 的统一入口。对候选编码（float32 BE、float64 BE、word-reversed float64）采样打分（``_score_coord_axes``），并参考描述符 `type` 与 `n_verts`：

- 坐标幅值合理区间：``[_COORD_MIN_ABSMAX, _COORD_MAX_ABSMAX]`` = ``[1e-4, 1e6]`` 米量级
- 顶点数取自描述符 ``dim0``（``dim0 > 1``），而非 ``byte_count // 8``
- word-reversed 选中时对磁盘轴序 X,Z,Y 做 Y/Z 交换

### 5.4 LS_SurfaceRegions

每个表面区域由 **三个连续数据块** 组成：

1. **name** — ASCII，NUL/空格填充  
2. **face_ids** — I4[m]，全局面索引（0-based，对应 `LS_Links` 面数组下标）  
3. **weights** — I4[m]，与 face_ids 等长（目前观测值均为 1）

`gph2cgns.py` 为每个 Zone 投影这些区域，在 `ZoneBC` 下为每个区域名写一个 `BC_t`（`PointList` 为 zone 内 1-based 面号；该 zone 无此区域面时省略 ` data` 数据集，与 `tr03_orig.cgns` 一致）。

### 5.5 LS_VolumeRegions

ASCII 字符串块列表，顺序即 CGNS 中体区域 Zone 的生成顺序（如 `FluidRegion`、`Rotate_MovingVolumeRegion` 等）。`gph2cgns` 按名称与 `LS_CvolIdOfElements` / Part 映射匹配 cell 子集。

### 5.6 LS_Parts

每个 Part 记录布局（节选）：

```
[12,4,1,1] … [12,4,255,4] [12,1,255,1]
[<part name>, 255 B ASCII]
[trailer=255]
… post-name 描述符 / 数据块 …
```

#### 5.6.1 简单 Part（单一 cvol_id）

post-name 区通常为 **`[1, cvol_id]`** 形式的 `[12,4,X,4]` 链——前导 `1` 为 marker，末尾值为不透明 Part id：

| 文件 | Part | 描述符链 | cvol_id |
|------|------|----------|---------|
| `tr03.gph` | Case / Rotate | `[1,1]` / `[1,2]` | 1 / 2 |
| `laptop_simplified_voxel_less.gph` | air_domain / rotation1 / rotation2 | `[1,1]` / `[1,9]` / `[1,11]` | 1 / 9 / 11 |
| `laptop_simplified_more_regions.gph` | outlet11 / rotation1 | `[1,2]` / `[1,7]` | 2 / 7 |

#### 5.6.2 复合 Part（cvol_id 成员列表）

多区域 laptop 模型中，背景流体 Part（如 `air_domain`）往往**不**对应单一 cvol_id，而携带显式成员列表：

```
[12, 4, N, 4]     ← N = 后续 I4 列表长度（不是 cvol_id！）
I4[N]             ← 该 Part 拥有的全部 cvol_id
```

实测 `laptop_simplified_more_regions.gph` 中 `air_domain`：

- `[12,4,66,4]` + `I4[66]` = `{1,8,9,…,72}` \ `{2,3,4,7}`（66 个 id）
- 151 375 cells（= 全网格 173 094 减去 outlet/rotation 四个 Part 的 21 719 cells）
- 误把 `N=66` 当作 cvol_id 只会得到 30 cells（cvol_id 66 在 mesh 中仅 30 单元）

#### 5.6.3 解析与 Zone 掩码

**唯一实现**：`gph_model.parse_ls_parts(data, cvol_id=…)` → `[(name, PartCvolSpec), …]`，其中 `PartCvolSpec = int | frozenset[int]`。`gph2cgns.py` 直接导入；`gph_parser.py` / `gphviewer.py` 共用。

1. 先解析 `LS_CvolIdOfElements`，唯一值集合 $S$ 为权威来源。
2. 对每个 Part post-name 区：若存在 `[12,4,N,4]` 且紧跟 `I4[N]`（值均 ∈ $S$），返回 `frozenset` 成员列表。
3. 否则按简单 Part：取 post-name 链中**最后一个属于 $S$ 的值**。
4. Zone cell 掩码：`part_cvol_cell_mask(cvol_id, spec)` — 单 id 用 `==`，集合用 `np.isin`。

> **历史**：§11.9–§11.11 演进见 `DEV_SUMMARY.md`；§11.12 增加复合 Part 成员列表（#18）。

Zone 命名规则（与 FLDUTIL 对齐，见 `DEV_SUMMARY.md` §10）：

- 路径深度 ≥ 2（含多个 `.`）→ 直接用完整路径（如 `laptop_3d_geom.____.air_domain`）  
- 否则加 `FPHPARTS.` 前缀（如 `FPHPARTS.tr03.Case`）  
- 根级 Part + 空 assembly 前缀 → `fan2.fan1.rotation1` 等

### 5.7 LS_Assemblies

UTF-8 XML，描述 `<assembly>` / `<part>` 层次。`gph2cgns` 解析 `part_paths` 与 `root_empty_prefix`（首个顶层 assembly 下、与根级 part 数量相同的**空** assembly 名拼接），用于 Part Zone 命名。

## 6. CGNS 多 Zone 输出概要

`gph2cgns.py` 根据 `LS_VolumeRegions` + `LS_Parts` + `LS_Assemblies` 生成多个 `Zone_t`，每个 Zone 为全局网格的子集（顶点/面重编号），并附带：

- `GridElements_Faces`（NGON_n=22）  
- 与 Zone 同名的 `NFACE_n`（带符号面索引）  
- `ZoneBC`：来自 `LS_SurfaceRegions` 的命名 BC 族  
- `FlowSolution` 占位  

无分区元数据时回退为 `FluidRegion` + `FPHPARTS.box_vol` 两 Zone（`box_ansa.gph` 行为）。

### 6.1 超大 GPH 文件（>512 MiB）

`gph2cgns.py`、`gph_parser.py`、`gphviewer.py`（`GphDocument.load`）对超过 **512 MiB** 的文件使用 **内存映射（mmap）**，避免整文件读入 RAM。坐标与 conn 数组通过 `numpy.frombuffer` 批量读取。

`laptop_simplified_voxel_v4.gph` 转换参考：解析 ~5 分钟，完整 CGNS 写出 ~16 分钟（输出约 14 GiB）。

`laptop_simplified_denser_v2_gph.gph`（~5.9 GiB，conn 三段 ~2.05 GiB）：解析 ~9 分钟。

`laptop_simplified_voxel_v6.gph`（~4.9 GiB，conn 两段 ~1.91 GiB）：解析 ~8 分钟。

## 7. 使用 Python 解析与验证

解析 GPH 结构：

```bash
python gph_parser.py [file.gph]    # 默认 tests/box_ansa.gph
```

对比体区域 Zone cell 数（需 `tests/{name}_orig.cgns` 参考文件）：

```bash
python tests/test_volume_zone_cells.py          # 简洁对比
python tests/test_volume_zone_cells.py -v       # 含 LS_Parts 链 / cvol_id / 复合 Part 成员列表
python tests/test_volume_zone_cells.py -v tests/laptop_simplified_more_regions.gph
```

将输出节布局、数据采样和完整格式说明；`-v` 模式便于核对 `LS_Parts` 与 `LS_CvolIdOfElements` 映射是否正确。

## 8. 参考

- 与 CGNS 的 ADF（Advanced Data Format）相似
- 32 字符标签符合 ADF 节点标签约定
- CRDL-FLD 可能表示 "Card/Record Field" 或厂商自定义格式

---

## 9. CGNS 输出 HDF5 格式约束（重要）

> 本节描述的不是 GPH 输入格式，而是 `gph2cgns.py` 写出的 **CGNS 输出文件**在底层 HDF5 层面必须满足的兼容性约束。详细推导见 `DEV_SUMMARY.md` 第 9 章。

### 9.1 HDF5 superblock 版本对下游 CFD 工具的影响

HDF5 文件起始 8 字节是固定魔数 `\x89HDF\r\n\x1a\n`，紧接其后的第 9 字节（offset 8）是 **superblock 版本号**。常见的三种版本对应不同的 HDF5 文件格式特性：

| superblock 版本 | 触发条件 | object header 格式 | 组存储 | 典型读取器兼容性 |
|----------------|----------|---------------------|--------|------------------|
| **v0** | `libver_low='earliest'`（h5py 默认） | v1（无 4 字节签名） | v1.6 sym-table（SNOD + TREE + HEAP） | **所有** CGNS 工具：ANSA、Tecplot、ParaView、CGNS 官方库 |
| **v1** | 非默认 `H5Pset_sym_k` 或 `H5Pset_istore_k` | v1 | sym-table，但 B-tree 参数可调 | 同 v0，绝大多数工具 |
| **v2** | `libver_low='v108'` 或更高 | v2（带 `OHDR` 4 字节签名） | compact link 或 dense（fractal heap） | ⚠️ **部分老旧 CGNS 模块不支持**，例如 **ANSA** 会报 `No bases found!` |

### 9.2 已知不兼容案例：ANSA "No bases found!"

**症状**：用 ANSA 导入 `gph2cgns.py` 生成的 CGNS 文件时报错：

```
No bases found!
```

**根因**：之前的提交曾设置 `h5py.File(libver=("v108","v108"))` 以启用 HDF5 1.8 紧凑链接存储（compact link storage），把输出文件从 73 KB 缩到 31 KB。该设置生成 **v2 superblock**，而 ANSA 的 CGNS 读取器只接受 v0 / v1 superblock，遇到 v2 时直接判定文件无 Base 节点。

**修复**：固定使用 `libver=('earliest', 'v108')`（h5py 默认），即 **v0 superblock + v1.6 sym-table 布局**。代价是文件比官网 vendor 参考大约 2 倍（73 KB vs 35 KB），但保证所有目标 CFD 工具都能读取。

### 9.3 vendor (`box_ansa_orig.cgns`) 的特殊紧凑布局

官网 `FLDUTIL` 工具产生的参考文件 `box_ansa_orig.cgns` 兼具两者之长：

- superblock **v0**（ANSA 可读）
- 但**无** SNOD / TREE / HEAP（即不是 v1.6 sym-table）
- 也**无** OHDR 签名（即不是 v2 object header）

逐字节解码根组 object header 可见使用了 **v1 object header 内嵌 HDF5 1.8 Link Info / Group Info 消息**：

```
Object header @ 0x60 (root group):
  version: 1                                    ← v1 object header（无签名）
  msg 0: type=0x02 (Link Info), size=40        ← HDF5 1.8 紧凑链接
  msg 1: type=0x0A (Group Info), size=8        ← HDF5 1.8 组信息
  msg 6: type=0x0C (Attribute), size=40
  ...
```

**该布局当前无法通过 h5py 复刻**：

- `h5py.h5p.PropGCID` 未导出 `set_link_phase_change` Python 包装。
- 即使通过 `ctypes` 直接调用 C 层 `H5Pset_link_phase_change` 在 GCPL / FCPL 上设置成功，HDF5 在 `libver_low='earliest'` 时仍会回退到 v1.6 sym-table 存储，无法触发 1.8 紧凑链接消息。

因此 `gph2cgns.py` 选择最稳妥的 **v0 superblock + v1.6 sym-table**（73 KB），与 vendor 内容字段级 `PERFECT MATCH`，仅 HDF5 元数据布局不同。

### 9.4 HDF5 签名速查表（用于逐字节诊断 CGNS 输出）

| 签名（4 字节 ASCII） | 含义 | 出现意味着 |
|---------------------|------|------------|
| `\x89HDF\r\n\x1a\n` | HDF5 魔数（offset 0） | 文件是 HDF5 容器 |
| `OHDR` | v2 object header | superblock 必为 v2，ANSA 等老读取器可能拒绝 |
| `OCHK` | object header continuation chunk (v2) | 同上 |
| `SNOD` | v1.6 symbol-table node | v0/v1 superblock + v1.6 组存储 |
| `TREE` | v1 B-tree | 同上，每组一个 |
| `HEAP` | v1 local heap | 同上，每组一个 |
| `FRHP` | fractal heap header | dense link 或 attribute 存储，HDF5 1.8+ |
| `BTHD` | v2 B-tree header | dense link 索引，HDF5 1.8+ |
| `GCOL` | global heap collection | 变长数据集（如 string）共享存储 |

诊断命令示例：

```bash
python3 -c "
data = open('your.cgns','rb').read()
print(f'size={len(data)} sb_v={data[8]}')
for sig in (b'OHDR', b'SNOD', b'TREE', b'HEAP', b'FRHP', b'BTHD'):
    print(f'  {sig.decode()}: {data.count(sig)}')
"
```

### 9.5 输出文件 sanity check

`gph2cgns.py` 生成的 CGNS 应满足：

- `data[8] == 0` （superblock v0）
- `OHDR` count = 0
- `SNOD` count > 0（与组数量相当）
- 通过 `h5py.File(path, 'r')` 可正常打开并读取 `Base` / `Base/<Zone>` / 各 `Elements_t` 子节点
- ZoneType 数据集字节序列为 `Unstructured`

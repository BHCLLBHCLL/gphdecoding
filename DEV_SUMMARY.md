# gph2cgns 开发总结

## 目录

1. [项目背景](#1-项目背景)
2. [第一阶段：坐标解析 Bug 修复](#2-第一阶段坐标解析-bug-修复)
3. [第二阶段：面单元与边界条件完善](#3-第二阶段面单元与边界条件完善)
4. [GPH 格式逆向工程结论](#4-gph-格式逆向工程结论)
5. [最终输出格式说明](#5-最终输出格式说明)
6. [改动文件清单](#6-改动文件清单)
7. [经验总结](#7-经验总结)
8. [第三阶段：ANSA 方言适配 + FLDUTIL 输出对齐](#8-第三阶段ansa-方言适配--fldutil-输出对齐)
9. [第四阶段：HDF5 superblock 版本与 ANSA 兼容性](#9-第四阶段hdf5-superblock-版本与-ansa-兼容性)
10. [第五阶段：多 Zone 分区、多面体网格与命名 BC](#10-第五阶段多-zone-分区多面体网格与命名-bc)
11. [第六阶段：超大 GPH（conn 分块 + mmap）](#11-第六阶段超大-gphconn-分块--mmap)

---

## 1. 项目背景

`gph2cgns` 工具用于将 SCTpre 生成的专有二进制格式 `.gph` 文件
（魔数 `CRDL-FLD`，大端序）转换为工业标准 CFD 网格格式 CGNS/HDF5。

测试文件 `box.gph`：10mm × 10mm × 10mm 立方体网格，52 个顶点，
135 个四面体单元，307 个三角面。

---

## 2. 第一阶段：坐标解析 Bug 修复

### 2.1 问题现象

`CoordinateX` 应在 `[0, 0.01]` 范围内，实际解析结果出现大量 **89000+** 异常值。

| 项目 | 期望值 | 实际值（修复前） |
|------|--------|-----------------|
| 顶点数 | 52 | 108（错误） |
| CoordinateX 范围 | [0, 0.01] | [-6.8e8, 3.5e34] |
| CoordinateY 范围 | [0, 0.01] | [-1.5e29, 7.9e34] |
| CoordinateZ 范围 | [0, 0.01] | [-4.3e15, 4.2e36] |
| 含 89000+ 的分量数 | 0 | 22 |

### 2.2 调查过程

**步骤一：排查硬编码偏移量**

原代码使用硬编码参数读取 float32 坐标：

```python
# 原代码（错误）
nodes_data_start = 0x2750
n_vertices = (0x2C60 - 0x2750) // 12  # = 108（错误）
xyz[i, 0] = read_f32_be(data, base)   # float32 误读
```

调整偏移和步长后，89000+ 问题依然存在，说明根因不是偏移错误。

**步骤二：定位 89000+ 的来源**

扫描文件中所有 float32 > 80000 的位置，发现异常值 `89128.96`
对应 4 字节 `47 AE 14 7B`，集中且规律出现于 `LS_Nodes` 节。

**步骤三：关键发现——坐标实为 float64**

同区域另一高频值 `1.035` 对应 4 字节 `3F 84 7A E1`，拼接后：

```
字节序列: 3F 84 7A E1  47 AE 14 7B
↓ 标准大端 float64
0x3F847AE147AE147B  =  0.01（精确值）
```

验证：
```python
struct.unpack(">d", bytes([0x3F,0x84,0x7A,0xE1,0x47,0xAE,0x14,0x7B]))[0]
# → 0.01000000000000000021  ✓
```

| 读法 | 字节 | 错误结果 |
|------|------|---------|
| float32（错误） | `3F 84 7A E1` | 1.035 |
| float32（错误） | `47 AE 14 7B` | **89128.96** ← 异常值根源 |
| float64（正确） | 8字节合并 | **0.01** ✓ |

**步骤四：发现字反转存储**

坐标以 **字反转（word-reversed / middle-endian）** 格式存储：

```
标准大端 float64:   [3F 84 7A E1][47 AE 14 7B]  (高32位在前)
GPH 文件实际存储:  [47 AE 14 7B][3F 84 7A E1]  (低32位在前)
```

**步骤五：揭示 LS_Nodes 三轴块结构**

`LS_Nodes` 节将 X、Z、Y 三个坐标轴分别存储为独立数据块：

```
LS_Nodes 节 (0x26B0–0x2C60)
├── 40B  节标签头
├── 4×16B 元数据描述符
├── [16B 描述符: 12/8/52/1] + [12B 块头] + 52×8B X 坐标（字反转 float64）
├── [16B 描述符: 12/8/52/1] + [12B 块头] + 52×8B Z 坐标（字反转 float64）
└── [16B 描述符: 12/8/52/1] + [12B 块头] + 52×8B Y 坐标（字反转 float64）
```

文件内轴块顺序为 **X → Z → Y**，读取后需重排为 X, Y, Z。

### 2.3 根本原因

原代码存在三个层叠错误：

| # | 错误 | 原始 | 正确 |
|---|------|------|------|
| 1 | 数据类型 | float32（4字节） | float64（8字节，字反转） |
| 2 | 起始偏移 | 硬编码 `0x2750` | 动态扫描描述符 `[12/8/n/1]` |
| 3 | 顶点计数 | `(0x2C60-0x2750)//12 = 108` | 从描述符读取 `dim0 = 52` |

### 2.4 修复方案

新增字反转读取函数、动态节定位和三轴块解析：

```python
def read_f64_wr(data: bytes, pos: int) -> float:
    """读取字反转 float64：存储顺序 [低32位大端][高32位大端]。"""
    lower = int.from_bytes(data[pos:pos+4], "big")
    upper = int.from_bytes(data[pos+4:pos+8], "big")
    return struct.unpack(">d", ((upper << 32) | lower).to_bytes(8, "big"))[0]
```

### 2.5 修复效果

| 项目 | 修复前 | 修复后 |
|------|--------|--------|
| 顶点数 | 108 | **52** |
| X/Y/Z 范围 | [-∞, +∞]（含89000+） | **[0, 0.01]** |
| 89000+ 异常值 | 22 | **0** |
| 坐标精度 | float32 | **float64** |

---

## 3. 第二阶段：面单元与边界条件完善

### 3.1 需求

测试发现转换后的 CGNS 文件：
1. **只有体单元，缺少面单元数据**（无法在 CFD 软件中定义边界条件）
2. **命名不符合规范**（参考文件 `box_ngons.cgns` 的命名约定）

### 3.2 参考文件分析

`box_ngons.cgns` 的 CGNS 结构：

```
CGNSLibraryVersion (R4, 4.2)
Base (CGNSBase_t, [3,3])
  box_vol (Zone_t, [[n_verts],[n_cells],[0]])
    ZoneType          → "Unstructured"
    GridCoordinates   → CoordinateX/Y/Z (R8)
    NGonElements      → 面元素 (NGON_n, type=22) + ElementStartOffset
    NFaceElements     → 体元素 (NFACE_n, type=23) + ElementStartOffset
    ZoneBC
      box_surfs       → BCWall, GridLocation=FaceCenter, PointList
```

### 3.3 LS_Links 节逆向解析

通过系统分析发现 `LS_Links` 节含 **5 个数据块**，存储 307 个三角面的拓扑信息：

| 块 | 大小 | 含义 |
|----|------|------|
| owner | 307×I4 | 每个面的拥有单元 ID（0-indexed） |
| neighbor | 307×I4 | 对侧单元 ID，`0xFFFFFFFF` = 边界面 |
| npe | 307×I4 | 每面节点数（全为 3 = 三角面） |
| face_type | 1×I4 | 值=4（四面体单元） |
| conn | 921×I4 | 节点索引，**列主序**存储 |

**关键发现：节点索引列主序存储**

```
conn 数组布局（共 307面 × 3节点 = 921 个值）:
[node0_face0, node0_face1, ..., node0_face306,   ← 307个 node[0]
 node1_face0, node1_face1, ..., node1_face306,   ← 307个 node[1]
 node2_face0, node2_face1, ..., node2_face306]   ← 307个 node[2]
```

按行主序误读则每个面会出现重复节点，这是面元素解析错位的原因。

**网格拓扑验证：**

```
307 面中：  边界面 75 个，内部面 232 个
验证：2 × 232 + 75 = 539 ≈ 135 × 4 = 540 (四面体，每个4个面) ✓
```

### 3.4 哨兵值过滤

块的 byte_count 值（1228 = 307×4）会作为哨兵值泄漏进 owner/neighbor 数组末尾。过滤条件：

```python
# 任何 cell index > n_faces 均为非法哨兵值
owner[owner > n_faces] = -1
neigh[(neigh > n_faces) & (neigh != -1)] = -1
# 节点索引越界截断
face_nodes[face_nodes >= n_vertices] = n_vertices - 1
```

### 3.5 CGNS 输出重写

按 `box_ngons.cgns` 规范重写输出，命名对比：

| 项目 | 修改前 | 修改后（对齐参考） |
|------|--------|--------------------|
| Zone 名称 | `Zone1` | `box_vol` |
| Zone data shape | `(1,3)` | `(3,1)` |
| CGNSLibraryVersion 类型 | C1 string | **R4 float32** |
| 体单元 | `Hexahedra`（HEXA_8，解析错误） | `NFaceElements`（NFACE_n=23） |
| 面单元 | ❌ 无 | ✅ `NGonElements`（NGON_n=22） |
| 边界条件 | ❌ 无 | ✅ `ZoneBC/box_surfs`（BCWall，FaceCenter） |

### 3.6 最终验证结果

```
Reading: box.gph
  Vertices : 52
  Faces    : 307  (triangular, NGON_n)
  Cells    : 135  (tetrahedral, NFACE_n)
  BC faces : 75
```

| 检查项 | 结果 |
|--------|------|
| NGonElements 节点索引范围 | [1, 52] ✓ |
| NFaceElements 面索引范围 | [1, 307] ✓ |
| ZoneBC PointList 边界面数 | 75 ✓ |
| 坐标范围 X/Y/Z | [0, 0.01] ✓ |
| 89000+ 异常值 | 0 ✓ |

---

## 4. GPH 格式逆向工程结论

### 4.1 文件整体布局

| 偏移 | 大小 | 节名 | 内容 |
|------|------|------|------|
| 0x0000 | 28B | file_header | `CRDL-FLD` 魔数 + 维度 |
| 0x001C–0x05D8 | ~1.2KB | 元数据 | FileRevision/Application/Dimension 等 |
| 0x0628–0x08DC | 692B | LS_CvolIdOfElements | I4[135] 控制体ID |
| 0x08DC–0x26B0 | 7.6KB | LS_Links | 面连通性（5块数据） |
| 0x26B0–0x2C60 | 1.4KB | LS_Nodes | 顶点坐标（字反转 float64，三轴块） |
| 0x2C60–0x42B0 | 5.7KB | LS_SurfaceRegions | 面区域数据 |
| 0x42B0–0x45A4 | 756B | Element_InformationFlag | 单元标志 |

### 4.2 通用节记录格式

```
[I4=32][C1[32]=节名称][I4=32]   ← 40B 标签头
{描述符块} × N                   ← 每块 16B：[I4=12][I4=type][I4=dim0][I4=dim1]
{数据块} × M                     ← 每块：[12B 块头][payload]
```

### 4.3 字反转 float64（Word-Reversed Float64）

```
标准大端 float64:   [高32位][低32位]
GPH 文件存储:       [低32位][高32位]  ← PDP-endian / 混合端序
```

示例：坐标 `0.01` = `0x3F847AE147AE147B`，文件中存为 `47 AE 14 7B 3F 84 7A E1`。

```python
def read_f64_wr(data, pos):
    lower = int.from_bytes(data[pos:pos+4], 'big')
    upper = int.from_bytes(data[pos+4:pos+8], 'big')
    return struct.unpack(">d", ((upper<<32)|lower).to_bytes(8,'big'))[0]
```

### 4.4 LS_Links 五块数据结构

面连通性节点索引以**列主序（column-major）**存储：

```
[块1: owner 307×I4]   [块2: neighbor 307×I4]   [块3: npe 307×I4]
[块4: face_type 1×I4]
[块5: conn 921×I4]
  = [node0×307] + [node1×307] + [node2×307]
```

获取第 i 个面的 3 个节点（0-indexed）：
```python
node0 = conn[i]
node1 = conn[307 + i]
node2 = conn[307*2 + i]
```

---

## 5. 最终输出格式说明

生成的 `box.cgns` 结构（对齐 `box_ngons.cgns`）：

```
CGNSLibraryVersion  (R4, 4.2)
Base                (CGNSBase_t, I4 [3,3])
  box_vol           (Zone_t, I8 [[52],[135],[0]])
    ZoneType            → "Unstructured"
    GridCoordinates
      CoordinateX       (R8, 52个顶点, 范围 [0, 0.01])
      CoordinateY       (R8, 52个顶点, 范围 [0, 0.01])
      CoordinateZ       (R8, 52个顶点, 范围 [0, 0.01])
    NGonElements        (Elements_t, type=[22,0])
      ElementRange          [1, 307]    ← 307个三角面
      ElementStartOffset    shape(308,) ← 等步长3
      ElementConnectivity   shape(921,) ← 1-indexed顶点
    NFaceElements       (Elements_t, type=[23,0])
      ElementRange          [308, 442]  ← 135个四面体单元
      ElementStartOffset    shape(136,)
      ElementConnectivity   shape(538,) ← 1-indexed面索引
    ZoneBC
      box_surfs         (BC_t, "BCWall")
        PointList           shape(75,1) ← 75个边界面，1-indexed
        GridLocation        "FaceCenter"
```

---

## 6. 改动文件清单

| 文件 | 改动类型 | 主要内容 |
|------|----------|----------|
| `gph2cgns.py` | **核心重写** | 新增 `read_f64_wr()`、`_find_section()`、`_parse_ls_nodes()`、`_parse_ls_links()`；重写 `parse_gph_mesh()`；完全重写 CGNS 输出层（NGon/NFace/ZoneBC） |
| `gph_model.py` | 同步修复 | LS_Nodes 解析使用字反转 float64 |
| `gph_parser.py` | 文档+修复 | 更正节采样代码；补充字反转 float64 格式说明 |
| `GPH_FORMAT_SPEC.md` | 文档更新 | 补充 LS_Nodes 字反转 float64 规范；更新 LS_Links 五块结构说明 |
| `box.cgns` | 重新生成 | 使用修复后脚本输出的正确 CGNS 文件 |
| `box_ngons.cgns` | 新增（参考） | 从 main 分支拉取的参考文件，用于验证输出格式 |
| `DEV_SUMMARY.md` | 文档 | 本总结文档 |

**提交历史：**

```
15e1bf5  feat: 添加面单元(NGonElements)、体单元(NFaceElements)和边界条件(ZoneBC)
9db2c4c  regen: 使用修复后的 gph2cgns 重新生成 box.cgns
2085c5c  docs: 添加 gph2cgns 坐标解析 Bug 修复开发总结
e25a693  fix: gph2cgns逆向解析坐标值错误 - 坐标为字反转float64非float32
```

---

## 7. 经验总结

1. **"看起来合理"的值可能是最危险的陷阱**。float32 误读得到的
   `1.035` 在 [0,1] 范围内，伪装成了正常坐标，掩盖了真正的错误
   （实为 float64 坐标 `0.01` 的高32位误读）。

2. **规律性出现的"魔法常量"是关键线索**。`89128.96`（`0x47AE147B`）
   在文件中高频出现，最终证明是 float64 坐标 `0.01` 的低32位被
   float32 误读的产物，直接指向了格式解析的根本错误。

3. **数据布局顺序不可想当然**。LS_Links 节的节点索引以列主序
   （column-major）存储，而非通常预期的行主序（row-major），导致
   按面逐行读取时出现节点重复。

4. **哨兵值会污染相邻数据**。块的 `byte_count` 值（1228）因内存
   布局关系泄漏进 `owner/neighbor` 数组末尾，需要主动识别并过滤。

5. **混合端序真实存在于工业软件**。PDP-endian（字反转 float64）虽
   罕见，但确实出现在历史遗留格式中，逆向时需将其纳入候选假设。

6. **动态解析优于硬编码偏移**。通过扫描描述符定位数据块，使代码
   对不同大小的 GPH 文件具有健壮性，而不是依赖特定测试文件的固定偏移。

---

## 8. 第三阶段：ANSA 方言适配 + FLDUTIL 输出对齐

### 8.1 新输入：`box_ansa.gph`

测试文件由 ANSA / 新版 scFLOW 导出，是一个 8 顶点立方体、6 个四面体、
18 个三角面、12 个边界面。同目录提供官网 `FLDUTIL` 工具生成的参考文件
`box_ansa_orig.cgns`，作为"金标准"。

### 8.2 GPH 方言差异（旧 SCTpre vs 新 ANSA）

| 维度 | 旧 `box.gph` (SCTpre) | 新 `box_ansa.gph` (ANSA) |
|---|---|---|
| `LS_Nodes` 浮点编码 | **标准大端 float64**（与 ANSA 一致） | 标准大端 float64 |
| `LS_Nodes` 块头 | **8 字节 + 4 字节尾部哨兵** | 8 字节 + 4 字节尾部哨兵 |
| 轴顺序 | **X, Y, Z**（与 ANSA 一致） | X, Y, Z |
| `LS_Links` 块数 | 5（含 face_type） | **4（无 face_type）** |
| `conn` 存储 | 列主序 | **行主序** |

#### 重要修正：原"字反转 float64"理论是误判

第二阶段总结中曾断言旧 `box.gph` 使用字反转 float64 + 12 字节块头。重新验证发现这是**错误推论**：

- 字节序列在 `0x2730` 处实际为 `3f 84 7a e1 47 ae 14 7b ...`，按标准大端 float64 读取得 `0.01` ✓
- 用字反转读取同样字节得 `1.999e+37` ✗（明显非物理）

历史误判源于原代码 off-by-4 字节偏移（跳过 12 字节而非 8 字节），与字反转读法两次错误"互相抵消"，使**部分**顶点产生看似合理的数值——但大多数顶点本应是 1e37+ 量级。当前实现统一采用"跳 8 字节 + 标准大端"读法。

`gph2cgns.py` 现在通过通用的"数据块扫描器"（识别 `[I4=12][I4=bc][bc 字节][I4=bc 尾部]` 模式）定位每个坐标块，**完全不依赖固定的块头长度**，且仍同时尝试两种浮点编码作为防御性兜底。

### 8.3 FLDUTIL 输出约定（参考 `box_ansa_orig.cgns` 逆向）

| 约定 | 值 |
|------|----|
| 根节点属性 | `name="HDF5 MotherNode"`, `label="Root Node of HDF5 File"`, `type="MT"` |
| 根数据集 | `' hdf5version'`（33 字节零）、`' format'`（`IEEE_LITTLE_32`） |
| `CGNSLibraryVersion` | R4，值 = **3.21** |
| 描述符 | `Base/ReferenceState/ReferenceStateDescription` = `"Software Cradle FLDUTIL"` |
| Zone 数 | **2 个同内容 Zone**：`FluidRegion` 与 `FPHPARTS.box_vol` |
| Zone 数据类型 | **I4**（int32），不是 I8 |
| 面单元名 | **`GridElements_Faces`**（NGON_n=22），不是 `NGonElements` |
| 体单元名 | **与父 Zone 同名**（NFACE_n=23） |
| 元素/索引数组 dtype | 全部 I4 |
| NFace 面索引 | **有符号**：cell 为 neighbor 时取负 |
| NFace 单 cell 排序 | 按 `abs(face_id)` 升序 |
| BC family 名称 | `box_surfs`，值 = `"Null"`（不是 `"BCWall"`） |
| `BC_t` 数据 dtype | `int8`，长度 4 |
| `PointList` 形状 | `(n,)`（一维），不是 `(n, 1)` |
| `FlowSolution` | 占位节点，仅含 `GridLocation = "CellCenter"` |
| 组属性 | `flags`(int32×1)、`label`(\|S33)、`name`(\|S33)、`type`(\|S3) |

### 8.4 顶点重排序（first-use ordering）

FLDUTIL 不直接使用 GPH 中的顶点顺序，而是**按面连通性中首次出现顺序**
重新编号顶点。例如：

```
GPH conn (row-major, 0-based):
  f0: [1, 2, 3]   ← 新编号 1, 2, 3
  f1: [2, 0, 3]   ← 0 首次出现 → 新编号 4
  f2: [0, 1, 3]   ← 已编号
  f3: [1, 0, 2]
  f4: [5, 1, 2]   ← 5 首次出现 → 新编号 5
  f5: [1, 4, 2]   ← 4 首次出现 → 新编号 6
  ...
```

由此得到的顶点排列恰好与官网工具完全一致：

```
GPH index : 1  2  3  0  5  4  6  7
CGNS 1-id : 1  2  3  4  5  6  7  8
```

### 8.5 带符号 NFACE_n 连通性

CGNS NFACE_n 规定：

- 单元拥有（owner）该面 → 正面 ID
- 单元在该面的另一侧（neighbor）→ 负面 ID

`box_ansa_orig.cgns` 中典型 cell 3 = `[-8, 9, 10, 11]`：

- face 8（gph owner=1, neigh=2）→ cell 2(0-idx) 是 neighbor → `-8`
- face 9, 10, 11（gph owner=2）→ cell 2 是 owner → 正

每个 cell 的面列表按 `abs(face_id)` 升序排列。

### 8.6 验证结果

```
Reading: box_ansa.gph
  Vertices : 8
  Faces    : 18  (triangular, NGON_n)
  Cells    : 6   (NFACE_n)
  BC faces : 12
Writing: box_ansa.cgns
Done.
PERFECT MATCH against box_ansa_orig.cgns
```

逐字段（`group / dataset / attr / dtype / shape / data`）比较，全部一致。

旧 `box.gph` 回归测试：仍可正确解析 52 顶点、307 面、135 单元；并顺带
修正了一处边界面 off-by-one（74 而非旧版 75，满足
`2·233 + 74 = 540 = 4·135` 拓扑恒等式）。

### 8.7 改动文件清单

| 文件 | 改动 |
|------|------|
| `gph2cgns.py` | 核心重写：方言自检 + first-use 重排 + 带符号 NFACE_n + 完整 FLDUTIL CGNS writer |
| `gph_model.py` | `LS_Nodes` 解析改为双方言自动判别（保留 word-reversed 兼容） |
| `gph_parser.py` | 同步双方言支持；输出新增 `LS_Nodes_dialect` 字段 |
| `gphviewer.py` | 顶点数据类型由 `R4[n,3]` 改为 `R8[n,3]`；禁用历史上不正确的逐 cell 顶点编辑 |
| `GPH_FORMAT_SPEC.md` | 增加方言 A/B 对比与 `LS_Links` 详细块布局 |
| `DEV_SUMMARY.md` | 本节（第三阶段）|

---

## 9. 第四阶段：HDF5 superblock 版本与 ANSA 兼容性

### 9.1 问题现象

第三阶段完成后，曾尝试通过 `h5py.File(libver=("v108", "v108"))` 启用 HDF5 1.8 紧凑链接存储（compact link storage），把输出文件从 73 KB 缩到 31 KB（甚至比 35 KB 的官网参考更小），并通过逐字段比对验证 `PERFECT MATCH`。但用户用 **ANSA** 导入这个 31 KB 文件时，前置网格器报错：

> **No bases found!**

ANSA 拒绝识别文件中的 `CGNSBase_t` 节点，即使该节点结构上确实存在。

### 9.2 根因：HDF5 superblock 版本不兼容

HDF5 文件起始 8 字节是魔数 `\x89HDF\r\n\x1a\n`，紧接着第 9 字节（offset 8）就是 **superblock 版本号**。三个候选格式的字节级差异：

| superblock 版本 | 触发条件 | object header | 组存储 | OHDR 签名 | SNOD/TREE/HEAP |
|----------------|----------|---------------|-------|-----------|----------------|
| **v0** | `libver_low='earliest'` | v1（无签名） | sym-table（v1.6） | 0 | 45 / 45 / 45 |
| **v1** | 非默认 `sym_k` / `istore_k` | v1 | sym-table | 0 | 45+ |
| **v2** | `libver_low='v108'+` | v2 OHDR | compact link / dense | 84 | 0 / 0 / 0 |

逐字节确认三个文件：

```
box_ansa_orig.cgns (vendor)   : 35 448 B  superblock=0  OHDR=0   SNOD/TREE/HEAP=0/0/0  ← 特殊!
我们 libver=v108 (31 KB)      : 31 260 B  superblock=2  OHDR=84  SNOD/TREE/HEAP=0/0/0  ← ANSA 拒绝
我们 libver=default (73 KB)   : 73 360 B  superblock=0  OHDR=0   SNOD/TREE/HEAP=45/45/45 ← ANSA 可读
```

**ANSA 的 CGNS 读取器只接受 superblock v0 / v1**，遇到 v2 superblock 直接判定文件为空，抛出 "No bases found!"。

### 9.3 vendor 的特殊布局（v0 + Link Info 消息）

官网参考 `box_ansa_orig.cgns` 兼具两者之长：

- superblock **v0**（ANSA 可读）
- 但**没有** SNOD / TREE / HEAP（即不是 v1.6 sym-table 布局）
- 也**没有** OHDR 签名（即不是 v2 object header）

逐字节解码根组 object header（offset 0x60）发现：

```
=== Object header @ 0x60 ===
  version: 1                    ← v1 object header（无签名）
  num_msgs=15, ref_count=1
  msg  0: type=0x02 (Link Info), size=40    ← HDF5 1.8 特性
  msg  1: type=0x0A (Group Info), size=8    ← HDF5 1.8 特性
  msg  2: type=0x10 (Object Header Continuation), size=16
  msg  6: type=0x0C (Attribute), size=40
  ...
```

即 **"v0 superblock + v1 object header 内嵌 HDF5 1.8 Link Info / Group Info 消息"** 的混合布局——既保持 v0 superblock 的向后兼容，又通过 1.8 的紧凑链接消息避免 SNOD / TREE / HEAP 开销。

### 9.4 为什么 h5py 写不出 vendor 的紧凑 v0 布局

触发该布局需在 GCPL 上调用 `H5Pset_link_phase_change(max_compact, min_dense)`，把 `min_dense` 设得足够大以禁用 dense storage，配合 `libver=('earliest', 'v108')`。但实测：

1. **h5py 的 `PropGCID` 不导出** `set_link_phase_change` 方法（仅有 `set_attr_phase_change` 等近邻 API）。
2. **通过 `ctypes` 绕过 h5py 直接调用 C 层** `H5Pset_link_phase_change` 虽然返回成功（`herr_t = 0`），且 `H5Pget_link_phase_change` 能读出设置的值，但 `H5Gcreate2` 在 `libver_low='earliest'` 时仍回退到 v1.6 sym-table 布局，生成的文件依然带 SNOD / TREE / HEAP。
3. 类似地尝试 `H5Pset_sym_k(1, 1)` / `H5Pset_istore_k(1)` 缩小 B-tree 节点，文件可压到 44 KB，但 superblock 会被 HDF5 自动从 v0 升到 **v1**（因为非默认 `sym_k` 值需要 v1 superblock 来存储）——这也偏离 vendor 的 v0 格式。

结论：**在当前 h5py 3.x（HDF5 lib 2.0.0）版本下无法精确复刻 vendor 的 35 KB 紧凑布局**。要么 73 KB v0 / sym-table（最稳），要么 31 KB v2 / compact（ANSA 拒绝）。

### 9.5 取舍与最终选择

| 候选 | superblock | ANSA 可读 | 大小 | 内容 vs 参考 |
|------|-----------|-----------|------|--------------|
| `libver=(earliest, v108)` ← **采纳** | **v0** | ✅ | 73 KB | PERFECT MATCH |
| `libver=(v108, v108)` | v2 | ❌ "No bases found!" | 31 KB | PERFECT MATCH |
| ctypes `set_sym_k(1,1) + set_istore_k(1)` | v1 | ⚠️ 未验证 ANSA | 44 KB | PERFECT MATCH |
| vendor 参考（无法复刻） | v0 + 1.8 LinkInfo | ✅ | 35 KB | (基准) |

**最终选择**：`libver=('earliest', 'v108')`（即 h5py 默认），写出 v0 superblock + v1.6 sym-table 布局。代价是文件比 vendor 大约 2 倍（多出的 ~42 KB 全部是 HDF5 元数据 SNOD/TREE/HEAP，与 CGNS 内容无关），但换来**与 ANSA 等所有 CGNS 工具的完全兼容性**。

### 9.6 经验总结

1. **`PERFECT MATCH` 是字段级等价，不保证字节级等价**。两个 CGNS 文件在 HDF5 层可以有完全不同的元数据编码（superblock 版本、对象头版本、链接存储模式）但保持完全相同的 CGNS 树内容。验证 CGNS 输出时不能只比内容，还要确认 superblock 版本符合目标读取器期望。
2. **第三方 CGNS 读取器对 HDF5 文件格式版本的兼容范围差异巨大**。Tecplot、ParaView、官方 CGNS lib 通常都支持 v2 superblock；但部分 GUI 前置网格器（如 ANSA）的 CGNS 模块基于较旧的 HDF5 库，只接受 v0 / v1 superblock。**默认输出最老的 superblock 版本**是工业 CFD 互操作性的安全选择。
3. **h5py 的 Python API 是 HDF5 C API 的精选子集**。GCPL 的 `H5Pset_link_phase_change`、FCPL 的 `H5Pset_sym_k` / `H5Pset_istore_k` 等关键 API 没有 Python 包装。虽然可以通过 ctypes 直接调用，但 HDF5 内部还有"libver 控制是否使用 1.8 特性"的硬约束，并非所有设置都能在所有 libver 下生效。
4. **逐字节验证（superblock 版本、特征签名计数）**比仅依赖 h5py 高层 API 自检更可靠。本次 bug 之所以能定位，关键就在于对 `box_ansa_orig.cgns` 做了 superblock 字段解析和 OHDR / SNOD / TREE / HEAP 签名扫描，把"看起来一样的两个 CGNS"还原为字节级真实差异。
5. **修复决策应优先正确性 > 文件大小**。31 KB → 73 KB 是 130% 的体积膨胀，但只要内容不变、所有目标工具都能读取，体积本身就是非功能性指标，让位于功能正确性。

### 9.7 改动文件清单

| 文件 | 改动 |
|------|------|
| `gph2cgns.py` | 把 `h5py.File(..., libver=("v108","v108"))` 改回 `libver=("earliest","v108")`；更新模块级注释解释 superblock 版本选择 |
| `DEV_SUMMARY.md` | 本节（第四阶段） |
| `GPH_FORMAT_SPEC.md` | 新增"CGNS 输出 HDF5 格式约束"章节 |

---

## 10. 第五阶段：多 Zone 分区、多面体网格与命名 BC

### 10.1 新能力概览

在 ANSA / HDF5 superblock 对齐之后，`gph2cgns.py` 又扩展了三类能力，并同步到 `gph_model.py`、`gph_parser.py`、`gphviewer.py` 与格式文档：

| 能力 | GPH 来源 | CGNS 输出 |
|------|----------|-----------|
| **多 Zone 分区** | `LS_VolumeRegions` + `LS_Parts` + `LS_Assemblies` | 每个区域/Part 一个 `Zone_t` 子网格 |
| **多面体网格** | `LS_Links` 可变 `npe` + CSR `conn` | `NGON_n` 变长面 + `NFACE_n` |
| **命名边界条件** | `LS_SurfaceRegions` | 每 Zone 下多个 `BC_t` 族（空 PointList 亦保留）|

测试样例：`tr03.gph`（6 万+ 四面体/多面体混合）、`laptop_simplified_voxel_less.gph`（三 Part + 空 assembly 前缀）。

### 10.2 多 Zone 分区与 Zone 命名

**体区域**：`LS_VolumeRegions` 中的 ASCII 名按文件顺序各生成一个 Zone（如 `FluidRegion`、`Rotate_MovingVolumeRegion`），cell 归属由名称与 `LS_CvolIdOfElements` 启发式匹配。

**Part 区域**：`LS_Parts` 中每个 Part 再生成一个 Zone。命名规则（与 FLDUTIL 参考对齐）：

| 样例文件 | Part | 生成的 Zone 名 |
|----------|------|----------------|
| `box_ansa.gph` | `box_vol` | `FPHPARTS.box_vol` |
| `tr03.gph` | `Case`（嵌套 1 层 assembly） | `FPHPARTS.tr03.Case` |
| `tr03.gph` | `Rotate`（根级） | `FPHPARTS.Rotate` |
| `laptop_*.gph` | `air_domain`（深度 2） | `laptop_3d_geom.____.air_domain` |
| `laptop_*.gph` | `rotation1`（根级 + 空前缀） | `fan2.fan1.rotation1` |

`LS_Assemblies` XML 提供 `part_paths`；对根级 Part 若存在与根 part 数量相同的**空** `<assembly>` 兄弟，其名称拼接为 `root_empty_prefix`（laptop 的 `fan2.fan1`）。

### 10.3 cvol_id ↔ Part 映射修复（#13）

**问题**：`laptop_simplified_voxel_less.gph` 转换后 `rotation1` / `rotation2` Zone 为 **0 单元**，ANSA 无法使用。

**根因**：曾错误假定 `LS_CvolIdOfElements[i]` 等于 Part 在 `LS_Parts` 列表中的 1-based 序号。实际上 cvol_id 是写在每个 Part 名块**之后**的描述符 `[12, 4, cvol_id, 4]` 中的不透明 ID：

- `tr03`：Part 顺序与 cvol_id `{1,2}` 巧合一致  
- `laptop`：Part 列表顺序对应 cvol_id **`{1, 9, 11}`**，而非 `{1, 2, 3}`

**修复**：`_parse_ls_parts_with_cvol_ids()` 扫描每个 255 字节名块后的最后一个 `[12,4,X,4]` 描述符；Zone cell 掩码改为 `cvol_id == X`。

### 10.4 多面体（混合单元）网格

`LS_Links` 的 `npe` 数组不再假定全为 3。`tr03.gph` 中面节点数为 3–11 不等；连通性按 **CSR** 存储：

```
face_offsets[i+1] = face_offsets[i] + npe[i]
face i 的节点 = conn[face_offsets[i] : face_offsets[i+1]]
```

`gph2cgns` 写 CGNS 时使用 `ElementStartOffset` + 扁平 `ElementConnectivity`（NGON_n），NFACE_n 仍用带符号面 ID（owner 正 / neighbor 负）。

### 10.5 命名 BC（`LS_SurfaceRegions`，#12）

每个表面区域 = 三个连续块：`name` | `face_ids[]` | `weights[]`（0-based 全局面号）。

对每个 CGNS Zone：

- 将全局 face id 映射到 zone 局部 1-based 面号  
- 每个区域名一个 `BC_t`，`GridLocation=FaceCenter`，值为 `"Null"`（FLDUTIL 约定）  
- 该 zone 不含此区域的面时，仍保留 `BC_t` 组但**不写** `PointList/ data`（与 `tr03_orig.cgns` 一致）

`tr03` 每 zone 约 15 个 BC 族（`inlet`、`outlet`、`Rotate_Plane`、`*_Moving`/`*_Static` 等）。

### 10.6 配套工具与文档同步

| 文件 | 同步内容 |
|------|----------|
| `gph_model.py` | 通用块扫描器；`LS_Parts`/`LS_VolumeRegions`/`LS_Assemblies`/`LS_SurfaceRegions` 解析；`LS_Nodes`/`LS_Links` 与 gph2cgns 对齐 |
| `gph_parser.py` | 动态节布局；mesh/partition 摘要输出；更新 `format_description()` |
| `gphviewer.py` | 树中显示新节；状态栏 mesh 摘要；分区/拓扑只读 |
| `GPH_FORMAT_SPEC.md` | §5.1–5.7 分区节；§6 多 Zone CGNS；§9 HDF5 约束 |
| `DEV_SUMMARY.md` | 本节；超大文件见 §11 |

**相关提交**（main 分支近期）：

```
f93f997  fix(gph2cgns): correct cvol_id ↔ part mapping (#13)
e56e6d7  feat(gph2cgns): multi-zone partitioning + polyhedral mesh (#11)
dcc8bce  feat(gph2cgns): per-zone named BC families from LS_SurfaceRegions (#12)
```

---

## 11. 第六阶段：超大 GPH（conn 分块 + mmap）

### 11.1 问题现象

`laptop_simplified_voxel_v4.gph`（磁盘 **~3.7 GiB**）转换时报：

```
Vertices : 0
Cells    : 0  (LS_Links parse failed)
Error: could not extract mesh data from GPH.
```

### 11.2 根因

1. **conn 超 1 GiB 被拆分**：88M 面、`sum(npe)=360,934,738` 的 conn 需 ~1.44 GiB。文件里第一段 conn 块为 **1073741824 B**（268,435,456 个 I4），trailer 后紧跟 **`[I4=byte_count][raw 续接 payload]`**（无 `[I4=12]` 头），第二段 ~370 MiB。
2. **旧逻辑直接失败**：`_parse_ls_links` 在 `conn_block_entries < sum(npe)` 时 `return None`。
3. **性能/内存**：3370 万顶点用 Python 循环读 float64 极慢；整文件 `read()` 占满 3.7 GB RAM。

### 11.3 修复

| 改动 | 说明 |
|------|------|
| conn 续接 | 主块不足时读取裸 `byte_count` + raw 续块（亦支持标准 `[12,bc]` 块） |
| 坐标读取 | `np.frombuffer` 批量读三轴 float64 |
| 大文件 I/O | `>512 MiB` 使用 **mmap**（`parse_gph_mesh` / `open_gph_buffer` / `GphDocument.load`） |
| 面→单元映射 | 88M 面 `argsort` 分组，替代 Python 循环 |

### 11.4 验证（v4）

```
Vertices : 33723176
Faces    : 88833031  (polyhedral, 4..7 nodes per face)
Cells    : 27553410
Zones    : 4  (FluidRegion + out_air + rotation1 + rotation2)
```

解析 ~5 min，完整 CGNS ~16 min（输出 ~14 GiB）。

### 11.5 配套同步

| 文件 | 同步内容 |
|------|----------|
| `gph2cgns.py` | conn 续接、mmap、frombuffer、向量化 cell-face 映射 |
| `gph_model.py` | `open_gph_buffer`、conn 分块检测、`parse_ls_links_summary` / `parse_ls_nodes_vertices` 对齐 |
| `gph_parser.py` | mmap 解析、format 文本补充 conn 分块 |
| `gphviewer.py` | 大文件 mmap 打开、只读提示、状态栏 mesh 摘要 |
| `GPH_FORMAT_SPEC.md` | §5.2 conn 分块、§6.1 超大文件 |
| `DEV_SUMMARY.md` | 本节 |

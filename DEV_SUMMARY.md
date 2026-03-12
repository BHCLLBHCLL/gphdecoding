# gph2cgns 开发总结

## 目录

1. [项目背景](#1-项目背景)
2. [第一阶段：坐标解析 Bug 修复](#2-第一阶段坐标解析-bug-修复)
3. [第二阶段：面单元与边界条件完善](#3-第二阶段面单元与边界条件完善)
4. [GPH 格式逆向工程结论](#4-gph-格式逆向工程结论)
5. [最终输出格式说明](#5-最终输出格式说明)
6. [改动文件清单](#6-改动文件清单)
7. [经验总结](#7-经验总结)

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

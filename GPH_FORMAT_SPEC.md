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
| 0x0628-0x08DC | 692 B | LS_CvolIdOfElements | I4[135] 控制体 ID |
| 0x08DC-0x26B0 | 7636 B | LS_Links | I4[] 单元连接关系 |
| 0x26B0-0x2C60 | 1456 B | LS_Nodes | R8[n,3] 顶点坐标（字反转 float64，三轴块）|
| 0x2C60-0x42B0 | 5712 B | LS_SurfaceRegions | 面区域数据 |
| 0x42B0-0x45A4 | 756 B | Element_InformationFlag | 单元标志 |
| 0x45A4-EOF | - | OverlapEnd | 文件尾 |

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
| LS_CvolIdOfElements | I4[n] | 每个单元的控制体 ID，n 由描述符决定 |
| LS_Links | I4[] | 面拓扑（owner / neighbor / npe / [face_type] / conn） |
| LS_Nodes | R8[n,3] | 顶点坐标（双方言：字反转 float64 或标准大端 float64） |
| LS_SurfaceRegions | 变长 | 表面区域定义 |
| Element_InformationFlag | 变长 | 单元信息标志 |

### LS_Links 详细格式

按数据块顺序存储以下数组（每块前后是描述符 / 8B 或 12B 块头 / 4B 哨兵）：

| 块 | 元素数 × 字节 | 含义 |
|----|---------------|------|
| owner    | n_faces × I4 | 每个面所属（owner）单元 ID（0-indexed） |
| neighbor | n_faces × I4 | 对侧（neighbor）单元 ID；`0xFFFFFFFF` = 边界面 |
| npe      | n_faces × I4 | 每面节点数（三角面时全 = 3） |
| face_type | 1 × I4      | **仅旧 SCTpre 文件**：单元类型标志（如 4 = 四面体） |
| conn     | (n_faces × npe) × I4 | 面-节点连通性，0-indexed |

`conn` 的内部布局也有两种方言：

- **旧 SCTpre**：列主序，`conn[k * n_faces + i]` 表示面 `i` 的第 `k` 个节点
- **新 ANSA**：行主序，`conn[i * npe + k]` 表示面 `i` 的第 `k` 个节点

`gph2cgns.py` 根据是否存在 `face_type` 块以及 `conn` 总长度 / `n_faces` 自动判别。

### LS_Nodes 详细格式

节点坐标按 X / Y / Z 三个独立轴块存储。**所有已观察到的 GPH 文件均使用相同的标准大端 float64 编码**——历史上猜测的"字反转 float64"是早期解析代码 off-by-4 字节偏移导致的虚假现象。

#### 实际编码（统一格式）

```
[16B 描述符] 00 00 00 0C / 00 00 00 08 / n_verts / 00 00 00 01
[8B  块头]   00 00 00 0C / byte_count
[n_verts × 8B]  标准大端 IEEE-754 float64
[4B  尾部]   byte_count（哨兵，等于块字节数）
```

- 轴顺序：**X, Y, Z**（与 CGNS 一致，无需重排）
- 每块末尾有 4 字节的 `byte_count` 哨兵

示例（来自 `box_ansa.gph`）：坐标值 `float32(0.01)` 加宽为 float64 = `0x3F847AE140000000`，在文件中按标准大端写为 `3F 84 7A E1 40 00 00 00`。

示例（来自旧 `box.gph`）：坐标值 0.01（精确）= `0x3F847AE147AE147B`，在文件中按标准大端写为 `3F 84 7A E1 47 AE 14 7B`。

#### 历史"字反转 float64"误读

原始 `gph2cgns.py` 跳过了 12 字节块头（而正确应该跳过 8 字节），并使用字反转读法。对**特定**字节序列，这两个错误恰好相互抵消，使少数顶点产生看似合理的数值——但大多数顶点会得到 1e37 等夸张幅值。当前实现已修正为正确的"跳 8 字节 + 标准大端"读法。

#### 防御性方言自检

`gph2cgns.py` / `gph_model.py` / `gph_parser.py` 仍同时尝试标准大端与字反转两种读法，选取幅值物理合理（`1e-30 ≤ |v| ≤ 1e6` 或 `v == 0`）的结果，以防遇到其他未知的 GPH 变体。

## 6. 使用 Python 解析

运行：

```bash
python gph_parser.py [box.gph]
```

将输出节布局、数据采样和完整格式说明。

## 7. 参考

- 与 CGNS 的 ADF（Advanced Data Format）相似
- 32 字符标签符合 ADF 节点标签约定
- CRDL-FLD 可能表示 "Card/Record Field" 或厂商自定义格式

---

## 8. CGNS 输出 HDF5 格式约束（重要）

> 本节描述的不是 GPH 输入格式，而是 `gph2cgns.py` 写出的 **CGNS 输出文件**在底层 HDF5 层面必须满足的兼容性约束。详细推导见 `DEV_SUMMARY.md` 第 9 章。

### 8.1 HDF5 superblock 版本对下游 CFD 工具的影响

HDF5 文件起始 8 字节是固定魔数 `\x89HDF\r\n\x1a\n`，紧接其后的第 9 字节（offset 8）是 **superblock 版本号**。常见的三种版本对应不同的 HDF5 文件格式特性：

| superblock 版本 | 触发条件 | object header 格式 | 组存储 | 典型读取器兼容性 |
|----------------|----------|---------------------|--------|------------------|
| **v0** | `libver_low='earliest'`（h5py 默认） | v1（无 4 字节签名） | v1.6 sym-table（SNOD + TREE + HEAP） | **所有** CGNS 工具：ANSA、Tecplot、ParaView、CGNS 官方库 |
| **v1** | 非默认 `H5Pset_sym_k` 或 `H5Pset_istore_k` | v1 | sym-table，但 B-tree 参数可调 | 同 v0，绝大多数工具 |
| **v2** | `libver_low='v108'` 或更高 | v2（带 `OHDR` 4 字节签名） | compact link 或 dense（fractal heap） | ⚠️ **部分老旧 CGNS 模块不支持**，例如 **ANSA** 会报 `No bases found!` |

### 8.2 已知不兼容案例：ANSA "No bases found!"

**症状**：用 ANSA 导入 `gph2cgns.py` 生成的 CGNS 文件时报错：

```
No bases found!
```

**根因**：之前的提交曾设置 `h5py.File(libver=("v108","v108"))` 以启用 HDF5 1.8 紧凑链接存储（compact link storage），把输出文件从 73 KB 缩到 31 KB。该设置生成 **v2 superblock**，而 ANSA 的 CGNS 读取器只接受 v0 / v1 superblock，遇到 v2 时直接判定文件无 Base 节点。

**修复**：固定使用 `libver=('earliest', 'v108')`（h5py 默认），即 **v0 superblock + v1.6 sym-table 布局**。代价是文件比官网 vendor 参考大约 2 倍（73 KB vs 35 KB），但保证所有目标 CFD 工具都能读取。

### 8.3 vendor (`box_ansa_orig.cgns`) 的特殊紧凑布局

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

### 8.4 HDF5 签名速查表（用于逐字节诊断 CGNS 输出）

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

### 8.5 输出文件 sanity check

`gph2cgns.py` 生成的 CGNS 应满足：

- `data[8] == 0` （superblock v0）
- `OHDR` count = 0
- `SNOD` count > 0（与组数量相当）
- 通过 `h5py.File(path, 'r')` 可正常打开并读取 `Base` / `Base/<Zone>` / 各 `Elements_t` 子节点
- ZoneType 数据集字节序列为 `Unstructured`

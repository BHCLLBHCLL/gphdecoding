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
| 0x26B0-0x2C60 | 1456 B | LS_Nodes | R4[n,3] 顶点坐标 |
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
| LS_CvolIdOfElements | I4[n] | 每个单元的控制体 ID，n=135 |
| LS_Links | I4[] | 连接表（单元-单元或单元-面） |
| LS_Nodes | R4[n,3] | 顶点 XYZ 坐标（大端 float32） |
| LS_SurfaceRegions | 变长 | 表面区域定义 |
| Element_InformationFlag | 变长 | 单元信息标志 |

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

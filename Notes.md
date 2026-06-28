## 问题本质

GPH 里说的 **node id** 指 `LS_Links` 的 **conn 数组中的顶点索引**（0-based，指向 `LS_Nodes`），不是坐标本身。大文件里 conn 超过 ~1 GiB 会被 scFLOW 拆成多段；段与段之间用 **裸 I4 字节数** 作标记，其中 **`1073741824`（0x40000000）= 1 GiB 容量标记**。

解析时若把 **块头/字节数** 当成 conn 数据读入，这些标记就会出现在 `face_nodes`（CSR 顶点索引）里，表现为「node id = 1073741824」或类似超大值，进而导致面翘曲、体网格不封闭。

---

## 根因：conn 续接块格式与读偏移不一致

### 1. 已确认并已修复：末段双头格式（三段 conn）

三段及以上 conn 的 **最后一段** 格式为：

```
[I4=1073741824][I4=need_bytes][payload]
```

旧版 `_read_conn_continuations` 在 `bare_bc == 1GiB` 且 `need_bytes < 1GiB` 时，从 **`pos+4`** 开始读 payload，实际读到的是 **`need_bytes`**（如 `173846272`），被当成第一个顶点索引。

| 文件 | 错误值 | 含义 |
|------|--------|------|
| `box.gph` | `173846272` | 末段 payload 字节数 |
| `voxel_v9.gph` | `713679784` | 同上 |

位置：第二段 1 GiB 边界（conn entry `536870912` = 2×268435456）。

**修复**（§11.13）：若 `read_i32_be(pos+4) == need_bytes`，从 **`pos+8`** 读 payload。

### 2. 仍可能出现 `1073741824` 本身的情况

若偏移再错 4 字节（从 `pos` 而非 `pos+4` 读），第一个 I4 就是 **`1073741824`**。可能场景：

- 续接分支判断顺序不当，双头/单头格式识别错误
- 节边界 `sec_end` 截断，完整 1 GiB 块读失败后落入 fallback 分支
- 把非 conn 的大 I4 块误选为 conn 首块

### 3. 症状被静默掩盖

`parse_gph_mesh` 中有 **clamp** 逻辑：

```777:786:gph2cgns.py
        # Clamp any garbage node indices that exceed the known vertex count
        fn = link_data["face_nodes"]
        bad = fn >= n_vertices
        if bad.any():
            fn = fn.copy()
            fn[bad] = n_vertices - 1
            link_data["face_nodes"] = fn
```

超大「node id」（含 1 GiB 标记）会被改成 `n_vertices-1`，**不报错**，但对应面顶点全错 → 面翘曲/不封闭，难以定位。

---

## 续接格式速查

| 场景 | 裸续接格式 | 正确 payload 起点 |
|------|-----------|------------------|
| 完整 1 GiB 中间段 | `[I4=1GiB][payload 1GiB]` | `pos+4` |
| 两段 conn 末段（v4/v6） | `[I4=need_bytes][payload]` | `pos+4` |
| 三段+ conn 末段（box/v9） | `[I4=1GiB][I4=need_bytes][payload]` | **`pos+8`** |

首段 conn 仍是标准块 `[12, bc][payload][bc]`，由 `iter_data_blocks` 提供 payload 偏移，不含 1 GiB 标记。

---

## 可能解决方案

### A. 解析层（核心，部分已做）

1. **保留双头末段检测**（已实现）：`bare_bc==1GiB` 且 `inner_bc==need_bytes` → 从 `pos+8` 读。
2. **收紧 fallback**：分支 142–152（单头短末段）仅在 `read_i32_be(pos+4)` **不像** byte_count 时使用——若 `next4 == need_bytes` 或 `next4 == 1GiB`，不要当顶点索引。
3. **读完后校验**：对 conn 做 `max(conn) < n_vertices`；在 268435456 的整数倍边界检查是否出现 sentinel 值（`1073741824`、等于某段 `need_bytes` 等）。

### B. 失败策略（建议）

将 silent clamp 改为 **可配置 fail-fast**：

```python
SENTINEL = {1073741824, ...}
if (conn >= n_vertices).any() or np.isin(conn, list(SENTINELS)).any():
    raise ValueError("conn contains byte-count sentinel; check continuation parsing")
```

比默默 clamp 到 `n_vertices-1` 更易调试。

### C. 测试与诊断

- 已有 `tests/test_conn_continuations.py`、`tests/validate_conn.py`——对大文件跑后者，确认 `max index == n_vertices-1` 且在 chunk 边界无 sentinel。
- 转换前打印 `conn_chunks`、chunk 边界样本（`tests/diag_conn.py`）。

### D. 长期扩展

- **cvol / surface face_ids** 若未来也 >1 GiB 分块，需类似续接逻辑（当前测试样例尚未触发）。
- 文档已更新：`GPH_FORMAT_SPEC.md` §5.2、`DEV_SUMMARY.md` §11.13。

---

## LS_Nodes float32（FPH / tr03_9）

FPH 网格（如 `tests/tr03_9.fph`）的 `LS_Nodes` 描述符 **type=4**，三轴块为 **float32 BE**（`byte_count = n_verts × 4`）。

若按 float64 误读：

- 坐标幅值落在 ~**1e-13**（denormal），旧版 ``1e-30`` 阈值会误判为合理；
- 顶点数按 ``bc // 8`` 计算 → **减半**（221786 → 110893），conn 索引越界或 ParaView 挂起。

**修复**（§11.14）：`gph_model.parse_ls_nodes_xyz()` 统一 float32 / f64 BE / word-reversed；``_COORD_MIN_ABSMAX = 1e-4`` m；顶点数取自描述符 ``dim0``。

---

## fph2cgns FlowSolution / Zone 选项

| 开关 | 说明 |
|------|------|
| 默认 | FlowSolution 写 **R4 float32**（与 `LS_SPHFile` 源一致） |
| `--flow-f64` | 场变量写 R8 float64 |
| `--skip-fluid-region` | **不导出** `FluidRegion` 整个 Zone（非仅省略场数据） |
| `--clip-flow 1` | 场值 > 1e20 置 0 |

详见 `GPH_FORMAT_SPEC.md` §5.8、`DEV_SUMMARY.md` §11.15。

---

## 小结

| 项目 | 说明 |
|------|------|
| **直接原因** | conn 裸续接块头（1 GiB 标记或 `need_bytes`）被当作顶点索引读入 |
| **典型位置** | 每 268435456 个 conn entry 的 chunk 边界 |
| **已修复** | 三段 conn 末段双头 `[1GiB][need][payload]` |
| **残留风险** | 其它格式变体、clamp 静默掩盖、偶发 `1073741824`  literal |
| **优先动作** | 确认使用含 §11.13 修复的代码；对大文件跑 `validate_conn.py`；考虑 sentinel 检测 + fail-fast 替代 blind clamp |

若你看到的错误值是 **`1073741824` 本身**还是 **`173846272` 这类末段字节数**，可以进一步判断是「偏移少 4 字节」还是「双头未识别」——我可以据此给出更精确的补丁。
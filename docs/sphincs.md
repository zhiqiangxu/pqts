# SPHINCS+ (SLH-DSA) 原理详解

## 1. 概述

SPHINCS+（正式名称 SLH-DSA，Stateless Hash-Based Digital Signature Algorithm）是一种**无状态的基于哈希的签名方案**。它的安全性**仅依赖于哈希函数的安全性**，不依赖任何数学结构假设（如格、编码等），因此被认为是最保守、最安全的后量子签名方案。

**核心安全假设**：底层哈希函数（SHA-256 或 SHAKE256）的抗碰撞性、抗原像性和抗第二原像性。

**代价**：签名尺寸大（~7856 字节），签名/验签速度较慢。

## 2. 设计哲学

### 2.1 为什么基于哈希？

格基方案（Dilithium、Falcon）的安全性依赖于格问题的计算困难性。虽然目前没有已知的量子算法能破解，但**格问题的量子难度还没有被完全理解**。

哈希函数的安全性则有更强的信心：
- 几十年的密码分析研究
- 随机预言机模型下的可证明安全
- 即使发现弱点，也可以替换底层哈希函数

SPHINCS+ 是一个"保险策略"——如果未来发现格问题的量子攻击，SPHINCS+ 仍然安全。

### 2.2 有状态 vs 无状态

经典的哈希签名方案（如 XMSS、LMS）是**有状态的**——签名者必须记住已使用的密钥索引，每用一次就标记一个。如果状态丢失或重复使用索引，安全性会被破坏。

SPHINCS+ 通过巧妙的设计实现了**无状态**：用消息的哈希值伪随机选择使用哪个密钥，避免了状态管理。

## 3. 构建模块

SPHINCS+ 是一个层层堆叠的结构，从底层到顶层：

```
┌─────────────────────────────────────┐
│           Hypertree (HT)            │  ← 多层 XMSS 树的链
│  ┌──────────────┐                   │
│  │   XMSS Tree  │ × d 层            │  ← 每层一棵 Merkle 树
│  │  ┌────────┐  │                   │
│  │  │ WOTS+  │  │ × 2^h' 个叶子     │  ← 每个叶子是一个一次签名
│  │  └────────┘  │                   │
│  └──────────────┘                   │
│           +                         │
│  ┌──────────────┐                   │
│  │    FORS      │                   │  ← 少次签名（签消息摘要）
│  └──────────────┘                   │
└─────────────────────────────────────┘
```

### 3.1 WOTS+（Winternitz One-Time Signature）

WOTS+ 是最底层的一次性签名方案。

**原理**：基于哈希链。

```
私钥: sk = (sk_0, sk_1, ..., sk_{len-1})     // len 个随机值
公钥: pk_i = H^{w-1}(sk_i)                   // 每个值哈希 w-1 次
      pk = H(pk_0 || pk_1 || ... || pk_{len-1})  // 压缩为一个值
```

参数 w = 16（Winternitz 参数），表示每个链长度为 15 步。

**签名消息 M**：

```
1. 将 M 转为 base-w 表示: m_0, m_1, ..., m_{len1-1}  （每个 ∈ [0, 15]）
2. 计算校验和: C = Σ(w-1-m_i)，也转为 base-w
3. 签名: sig_i = H^{m_i}(sk_i)       // 在链上走 m_i 步
```

**验证**：

```
对每个 sig_i: pk_i' = H^{w-1-m_i}(sig_i)    // 继续走到链顶
pk' = H(pk_0' || pk_1' || ...)
检查 pk' == pk
```

**安全性**：
- 签名只暴露链的中间值，攻击者无法回退（哈希的单向性）
- 校验和防止攻击者在链上多走几步来伪造更大的消息值

**为什么是"一次性"**：两次签名可能暴露链上更多的中间值，攻击者可以组合出伪造签名。

参数计算（n=16 字节，w=16）：
```
len1 = ⌈8n / log2(w)⌉ = ⌈128/4⌉ = 32   // 消息部分
len2 = ⌊log2(len1·(w-1)) / log2(w)⌋ + 1 = 3   // 校验和部分
len = len1 + len2 = 35                    // 总链数
WOTS 签名大小 = 35 × 16 = 560 字节
```

### 3.2 XMSS（eXtended Merkle Signature Scheme）

XMSS 解决了 WOTS+ 的"一次性"限制——它用一棵 **Merkle 树**管理多个 WOTS+ 密钥。

```
            root
           /    \
         H01    H23
        / \    / \
      H0  H1  H2  H3     ← 内部节点 = H(left || right)
      |   |   |   |
     PK0 PK1 PK2 PK3     ← 叶子 = WOTS+ 公钥
```

树高 h' = 9，所以每棵 XMSS 树有 2^9 = 512 个叶子，可以签 512 条消息。

**签名 = WOTS+ 签名 + 认证路径（auth path）**

```
要证明 PK2 在树中:
  auth_path = [H3, H01]          // 兄弟节点

验证:
  H2 = WOTS_PK_from_sig(sig, msg)   // 从 WOTS 签名恢复公钥
  H23 = H(H2 || H3)                 // 用认证路径爬上去
  root' = H(H01 || H23)
  检查 root' == root
```

XMSS 签名大小 = WOTS 签名 + h'×n = 560 + 9×16 = **704 字节**

### 3.3 Hypertree（超级树）

单棵 XMSS 树只能签 512 条消息。要支持足够多的签名（2^63），SPHINCS+ 使用**多层 XMSS 树的链式结构**——Hypertree。

```
层 6 (顶层):    XMSS 树 → root = 公钥
                  |
                叶子 WOTS 签名 →

层 5:           XMSS 树 → root
                  |
                叶子 WOTS 签名 →

  ...（共 7 层）

层 0 (底层):    XMSS 树 → root
                  |
                叶子 WOTS 签名 → 签 FORS 公钥
```

每层有一棵 XMSS 树（高度 9），共 d=7 层。

- 顶层树的 root 就是**公钥**
- 每层的叶子（WOTS+ 密钥）签署下一层树的 root
- 底层的叶子签署消息的 FORS 签名

总签名容量：2^{h'×d} = 2^{9×7} = 2^{63}（远超实际需要）

Hypertree 签名大小 = d × XMSS 签名 = 7 × 704 = **4928 字节**

### 3.4 FORS（Forest of Random Subsets）

FORS 是一种**少次签名方案**（few-time signature），用于签署消息摘要。它比 WOTS+ 更高效地处理较长的消息摘要。

FORS 由 k=14 棵二叉树组成，每棵高度 a=12：

```
树 0:         root_0
             /      \
           ...      ...
          /              \
        leaf_0  ...  leaf_{4095}     ← 2^12 = 4096 个叶子

树 1:         root_1
             ...

...（共 14 棵树）

FORS 公钥 = H(root_0 || root_1 || ... || root_13)
```

**签名消息摘要 M**：

```
1. 将 M 分为 k=14 段，每段 a=12 bit → 14 个索引，各 ∈ [0, 4095]
2. 对每棵树，揭示对应索引的叶子值 + 认证路径
```

**安全性**：每棵树只揭示一个叶子。要伪造签名，需要知道未揭示叶子的原像（哈希的抗原像性）。FORS 允许有限次数的重复使用（安全性会缓慢退化），这对 SPHINCS+ 的无状态设计至关重要。

FORS 签名大小 = k × (1 + a) × n = 14 × 13 × 16 = **2912 字节**

## 4. 完整签名流程

### 4.1 密钥生成

```
KeyGen():
  1. SK.seed ← random(n)          // 用于派生所有子密钥
  2. SK.prf ← random(n)           // 用于消息随机化
  3. PK.seed ← random(n)          // 用于哈希的公开随机性

  4. PK.root = HT_Root(PK.seed, SK.seed)
     // 构建整个 Hypertree 的顶层 XMSS 树的根
     // 这需要生成 512 个 WOTS+ 公钥（每个涉及 35 次哈希链 × 15 步）
     // 所以 KeyGen 很慢：~70ms

  公钥 PK = (PK.seed, PK.root)    // 32 字节
  私钥 SK = (SK.seed, SK.prf, PK)  // 64 字节
```

公钥只有 32 字节——这是所有后量子签名中最小的！

### 4.2 签名

```
Sign(SK, msg):
  1. R = PRF_msg(SK.prf, PK.seed, msg)    // 消息随机化
  2. digest = H_msg(R, PK.seed, PK.root, msg)  // 长摘要

  3. 从 digest 提取:
     - FORS 消息: fors_msg (k×a = 168 bit)
     - 树索引: tree_idx (54 bit) → 决定用 Hypertree 的哪个位置
     - 叶索引: leaf_idx (9 bit)  → 决定用 XMSS 树的哪个叶子

  4. FORS_SIG = FORS_Sign(fors_msg)        // FORS 签名 = 2912 B
  5. FORS_PK = FORS_PK_from_sig(FORS_SIG, fors_msg)  // FORS 公钥

  6. HT_SIG = HT_Sign(FORS_PK, tree_idx, leaf_idx)   // Hypertree 签名 = 4928 B
     // 从底层到顶层逐层签名:
     //   层 0: XMSS_Sign(FORS_PK, leaf_idx, tree=tree_idx)
     //   层 1: XMSS_Sign(root_0, leaf=tree_idx[8:0], tree=tree_idx[17:9])
     //   ...
     //   层 6: XMSS_Sign(root_5, leaf=..., tree=0)

  7. return σ = (R, FORS_SIG, HT_SIG)     // 16 + 2912 + 4928 = 7856 B
```

**无状态的关键**：tree_idx 和 leaf_idx 是从**消息的哈希值**确定性派生的（伪随机）。不同的消息大概率映射到不同的叶子，不需要状态管理。只有当两条不同消息碰巧映射到同一叶子时，FORS 的安全性才会退化——但 2^63 个可能的叶子使碰撞概率极低。

### 4.3 验证

```
Verify(PK, msg, σ = (R, FORS_SIG, HT_SIG)):
  1. digest = H_msg(R, PK.seed, PK.root, msg)
  2. 提取 fors_msg, tree_idx, leaf_idx

  3. FORS_PK = FORS_PK_from_sig(FORS_SIG, fors_msg)
     // 从 FORS 签名恢复公钥（爬 14 棵树）

  4. root = HT_Verify(FORS_PK, tree_idx, leaf_idx, HT_SIG)
     // 从底到顶逐层验证:
     //   层 0: 从 XMSS 签名恢复 root_0
     //   层 1: 从 XMSS 签名恢复 root_1
     //   ...
     //   层 6: 从 XMSS 签名恢复 root_6

  5. return root == PK.root
```

## 5. 可调哈希函数（Tweakable Hash）

SPHINCS+ 不直接使用裸的 SHA-256。为了保证**域分离**（不同用途的哈希不会碰撞），使用可调哈希函数：

```
F(PK.seed, ADRS, M) = SHA-256(PK.seed || ADRS || M)
```

ADRS（Address）是一个 32 字节的结构，编码了：
- 层地址（layer）：在 Hypertree 的哪一层
- 树地址（tree）：在该层的哪棵树
- 类型（type）：WOTS 哈希 / WOTS 公钥压缩 / 树节点 / FORS
- 键对地址（keypair）：哪个叶子
- 链地址（chain）：WOTS 的哪条链
- 哈希地址（hash）：链中的第几步

这确保了 Hypertree 中任意两次哈希调用都使用不同的 ADRS，即使输入相同也产生不同的输出。

## 6. Simple vs Robust 变体

SPHINCS+ 有两种实例化方式：

- **Simple**（本实现）：F(seed, adrs, M) = H(seed || adrs || M)
  - 更快（约 2×）
  - 安全性在标准模型下需要更强的哈希函数假设

- **Robust**：F(seed, adrs, M) = H(seed || adrs || (M ⊕ H(seed || adrs)))
  - 更慢（多一次哈希）
  - 安全性归约更紧，在更弱的假设下可证明安全

NIST 标准化选择了 Simple 变体作为默认。

## 7. 安全性分析

### 7.1 各组件的安全性

| 组件 | 安全性依赖 | 破解方式 |
|------|-----------|---------|
| WOTS+ | 哈希链的单向性 | 找到哈希的原像 |
| XMSS | Merkle 树的碰撞抗性 | 找到哈希碰撞 |
| FORS | 哈希的抗原像性 | 猜测未揭示的叶子值 |
| Hypertree | 各层 XMSS 的安全性组合 | 攻破任一层 |

### 7.2 量子安全性

量子计算对哈希函数的最佳攻击：
- **Grover 算法**：将搜索加速到 O(√N)
- 对 SHA-256 的原像攻击：2^256 → 2^128（仍然安全）
- 对 SHA-256 的碰撞攻击：2^128 → 2^85（BHT 算法）

SPHINCS+-128s 的参数选择保证了 128 位量子安全。

### 7.3 与格基方案的比较

```
安全性信心:
  SPHINCS+ ████████████████████  最高（仅依赖哈希）
  Dilithium █████████████████    高（格问题研究充分）
  Falcon   ████████████████     高（NTRU 30年历史）
```

## 8. 参数总结

### SPHINCS+-128s-simple (SHA-256)

| 参数 | 值 | 含义 |
|------|-----|------|
| n | 16 | 哈希输出长度（字节）|
| h | 63 | 总树高度 |
| d | 7 | Hypertree 层数 |
| h' | 9 | 每层 XMSS 树高度 |
| k | 14 | FORS 树的数量 |
| a | 12 | 每棵 FORS 树的高度 |
| w | 16 | Winternitz 参数 |

### 尺寸

| | 大小 | 构成 |
|---|------|------|
| 公钥 | 32 B | PK.seed (16) + PK.root (16) |
| 私钥 | 64 B | SK.seed (16) + SK.prf (16) + PK (32) |
| 签名 | 7856 B | R (16) + FORS (2912) + HT (4928) |

### 性能特征

| 操作 | 时间 | 原因 |
|------|------|------|
| KeyGen | ~70 ms | 构建顶层 XMSS 树（512 个 WOTS 密钥） |
| Sign | ~500 ms | 构建完整路径上的 XMSS 树 + FORS |
| Verify | ~0.5 ms | 只需爬树验证路径，不需构建 |

验证比签名快 **1000 倍**——这对在线验证场景非常有利。

## 9. SPHINCS+ 的命名约定

```
SPHINCS+-128s-simple
         │  │  │
         │  │  └── simple: 哈希变体
         │  └───── s: small（优化签名大小）vs f: fast（优化速度）
         └──────── 128: 量子安全级别（bit）
```

- **128s**（small）：签名 7856 B，签名慢但小
- **128f**（fast）：签名 17088 B，签名快但大
- 还有 192s/192f、256s/256f 等更高安全级别的变体

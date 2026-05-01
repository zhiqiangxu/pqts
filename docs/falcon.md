# Falcon-512 原理详解

## 1. 概述

Falcon（Fast-Fourier Lattice-based Compact Signatures over NTRU）是基于 **NTRU 格**的后量子签名方案，使用 **Hash-and-Sign（哈希签名）** 范式。它的核心优势是**签名尺寸最小**（~666 字节），因此被 Solana 选择用于链上签名验证。

**核心安全假设**：NTRU 问题和 SIS（Short Integer Solution）问题在量子计算机下仍然困难。

**核心技术贡献**：高效的 FFT（快速傅里叶变换）上离散高斯采样，使得在 NTRU 格上的 Hash-and-Sign 变得实际可行。

## 2. 数学基础

### 2.1 多项式环

Falcon-512 工作在：

```
R_q = Z_q[X] / (X^512 + 1)，q = 12289
```

多项式次数 n = 512，模数 q = 12289（一个素数，且 q ≡ 1 mod 2048，支持 NTT）。

### 2.2 NTRU 格

NTRU 格是由一对短多项式 (f, g) 定义的特殊格结构：

```
公钥：h = g · f^{-1} mod q

NTRU 格：L = { (u, v) ∈ R^2 : u + v·h ≡ 0 (mod q) }
```

关键性质：
- **(f, g)** 是格 L 中的一个**短向量**（f, g 的系数很小，约 ±1）
- 知道 (f, g) → 容易在格中找短向量（签名）
- 只知道 h → 在格中找短向量是困难的（NTRU 问题）

### 2.3 NTRU 方程

完整的 NTRU 密钥需要四个多项式 (f, g, F, G) 满足：

```
f · G - g · F = q    （在 R = Z[X]/(X^n+1) 中）
```

这四个多项式构成 NTRU 格的一组短基：

```
B = | g   -f |
    | G   -F |
```

这组基的行向量都是短的（系数小），是格的"好基"（陷门）。而公钥 h 只给出格的"坏基"。

### 2.4 求解 NTRU 方程

给定短多项式 f, g，如何找到 F, G 满足 f·G - g·F = q？

Falcon 使用**递归场范数（Field Norm）**方法：

```
1. 计算 f' = N(f), g' = N(g)    // 场范数：降低多项式次数到 n/2
2. 递归求解 f'·G' - g'·F' = q   // 在更小的环中求解
3. 从 (F', G') 提升回 (F, G)    // 通过 Galois 群作用
```

递归到底（n=1）时，问题变为简单的扩展欧几里得算法。这是 Falcon 实现中最复杂的部分。

## 3. GPV 框架：Hash-and-Sign

Falcon 基于 **GPV（Gentry-Peikert-Vaikuntanathan）** 框架：

### 3.1 直觉

```
格中的一个"陪集"：所有满足 s1 + s2·h ≡ t (mod q) 的 (s1, s2) 构成一个集合

问题：在这个集合中找到一个短的 (s1, s2)
- 没有陷门 → 这是 CVP（最近向量问题），是困难的
- 有陷门 (f,g,F,G) → 可以用 Babai/Klein 算法高效求解
```

### 3.2 签名流程概览

```
Sign(sk, msg):
  1. nonce ← random
  2. t = HashToPoint(nonce || msg)           // 哈希到 R_q 中的多项式
  3. (s1, s2) = SamplePreimage(sk, t)        // 用陷门在格中采样短向量
  4. return σ = (nonce, s2)                  // s1 可从 s2 和 msg 恢复

Verify(pk, msg, σ = (nonce, s2)):
  1. t = HashToPoint(nonce || msg)
  2. s1 = t - s2 · h mod q                   // 恢复 s1
  3. return ||(s1, s2)||² < bound            // 检查短向量
```

### 3.3 为什么只存 s2？

因为 s1 + s2·h = t mod q，知道 (s2, t, h) 就能算出 s1。签名只需要存储 s2，验证时恢复 s1。这是 Falcon 签名紧凑的关键原因之一。

## 4. 离散高斯采样：Falcon 的核心

### 4.1 为什么需要高斯采样？

简单的 Babai 最近平面算法（四舍五入到最近格点）会**泄露私钥信息**。攻击者可以从多个签名的分布中统计推断出陷门。

解决方案：不是取最近格点，而是**按照离散高斯分布采样**格点。这保证签名的分布独立于私钥（统计零知识）。

### 4.2 Klein/GPV 采样

在 NTRU 格的 Gram-Schmidt 正交化基 B̃ 上逐层采样：

```
SamplePreimage(B, t, σ):
  1. 计算 B 的 Gram-Schmidt 正交化 B̃ = {b̃_1, ..., b̃_{2n}}
  2. 令 c = t（目标向量）
  3. for i = 2n down to 1:
       c_i = <c, b̃_i> / <b̃_i, b̃_i>        // 投影系数
       z_i = SampleZ(c_i, σ/||b̃_i||)        // 离散高斯采样
       c = c - z_i · b_i                     // 减去已采样的分量
  4. return t - c                             // 短向量
```

### 4.3 Falcon 的 FFT 采样（ffSampling）

上述 Klein 采样需要 O(n²) 复杂度。Falcon 的关键创新是利用 NTRU 格的**特殊代数结构**，通过 FFT 将采样加速到 O(n log n)：

```
ffSampling(t, T):    // T 是 Falcon 树（LDL 分解的 FFT 表示）
  if n == 1:
    return SampleZ(t, σ·T.value)

  (t0, t1) = splitFFT(t)           // 将 n 维分为两个 n/2 维
  (T0, T1) = T.children

  z1 = ffSampling(t1, T1)          // 递归采样右半
  t0' = t0 + (t1 - z1) · T.middle  // 调整左半的中心
  z0 = ffSampling(t0', T0)         // 递归采样左半

  return mergeFFT(z0, z1)          // 合并
```

这利用了 NTRU 格的 LDL* 分解在 FFT 域中具有递归结构（Falcon Tree）。

### 4.4 SampleZ：整数高斯采样

最底层需要从离散高斯 D_{Z,c,σ} 中采样整数：

```
Pr[z] ∝ exp(-|z - c|² / (2σ²))
```

实现方法有多种：
- **CDT（Cumulative Distribution Table）**：预计算累积分布表，二分查找
- **Knuth-Yao 算法**：基于 bit 的最优采样
- **拒绝采样**：简单但较慢

Falcon 要求采样器是**常数时间**的（防止侧信道攻击），且使用高精度浮点运算。

## 5. 签名压缩

Falcon 使用特殊的压缩编码使签名尽可能紧凑：

### 5.1 系数分布

s2 的系数近似服从离散高斯分布，大部分集中在 0 附近：

```
Pr[|s2_i| = 0] ≈ 30%
Pr[|s2_i| = 1] ≈ 25%
Pr[|s2_i| = 2] ≈ 18%
Pr[|s2_i| ≥ 5] < 3%
```

### 5.2 压缩方案

Falcon 使用**类 Golomb-Rice 编码**：
- 1 bit 符号
- 低位直接存储
- 高位用一元编码

这是一种变长编码，对接近 0 的值非常高效，平均每个系数约 1.1~1.3 bit。

总签名大小：1 byte header + 40 bytes nonce + ~625 bytes s2 ≈ **666 bytes**。

## 6. 验证流程详解

```
Verify(pk = h, msg, σ = (nonce, s2)):

  1. t = HashToPoint(nonce || msg)
     // 用 SHAKE-256 将消息哈希到 R_q 中的均匀随机多项式
     // 生成 512 个系数，每个在 [0, q) 中

  2. s1 = t - s2 · h mod q
     // 多项式乘法（用 NTT 加速）+ 多项式减法
     // 关键：这一步的计算量决定了验证速度

  3. ||(s1, s2)||² = Σ s1[i]² + Σ s2[i]²
     // 计算欧几里得范数的平方
     // s1, s2 用中心表示：系数在 (-q/2, q/2]

  4. return ||(s1, s2)||² < β²
     // β² = 34034726（Falcon-512 的范数界）
     // 短向量 → 合法签名；长向量 → 伪造
```

**验证为什么快**：只需要 1 次 NTT 乘法 + 1 次范数计算，没有复杂的采样或分解。这就是为什么 Firedancer/Solana 只实现验证——它在链上非常高效（~4 μs）。

## 7. 安全性分析

### 7.1 伪造签名的难度

要伪造签名，攻击者需要：
1. 给定 h 和 t = Hash(msg)
2. 找到短 (s1, s2) 使得 s1 + s2·h = t mod q

这等价于在 NTRU 格中解 CVP（最近向量问题），是 NP-hard 问题。

### 7.2 NTRU 问题的难度

给定 h = g·f^{-1} mod q，恢复 (f, g) 等价于在 NTRU 格中找最短向量（SVP）。最好的已知攻击是格基约化算法（BKZ），对 n=512 的安全性为 ~128 位。

### 7.3 与 Dilithium 的安全性比较

| | Falcon | Dilithium |
|---|---|---|
| 安全假设 | NTRU + SIS | Module-LWE + Module-SIS |
| 假设强度 | NTRU 有 30 年研究历史 | LWE 有更强的最坏情况归约 |
| 量子安全级别 | NIST Level 1 (128 bit) | NIST Level 2 (128 bit) |

## 8. 参数总结

| 参数 | 值 | 含义 |
|------|-----|------|
| n | 512 | 多项式次数 |
| q | 12289 | 模数 |
| σ | 165.74 | 高斯标准差 |
| β² | 34034726 | 签名范数上界 |
| 签名大小 | ~666 B | 最紧凑的后量子签名 |
| 公钥大小 | ~897 B | h 的压缩表示 |
| 私钥大小 | ~1281 B | (f, g, F, G) 的压缩表示 |

## 9. Falcon vs Dilithium：为什么 Solana 选 Falcon

```
签名大小对比：
  Falcon-512:    ~666 B    ████
  Dilithium-2:  ~2420 B    ███████████████

验证速度对比（越短越快）：
  Falcon-512:    ~4 μs     ██
  Dilithium-2:  ~40 μs     ████████████
```

Solana 每秒处理数千笔交易，每笔交易附带签名。签名越小 → 区块越小 → 带宽压力越低 → 吞吐量越高。Falcon 的签名只有 Dilithium 的 **~1/4**，在高吞吐场景下优势明显。

缺点：Falcon 的签名生成需要高精度浮点运算和复杂的 NTRU 求解器，实现难度远高于 Dilithium。但对 Solana 来说，签名生成在客户端/钱包完成（离线），链上只需验证——而验证是 Falcon 最快的环节。

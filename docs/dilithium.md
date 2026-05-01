# Dilithium-2 (ML-DSA-44) 原理详解

## 1. 概述

Dilithium 是 NIST 后量子密码标准化中的**首选数字签名方案**（正式名称 ML-DSA，Module Lattice Digital Signature Algorithm）。它基于**模格（Module Lattice）**上的困难问题，使用 **Fiat-Shamir with Aborts** 范式构造签名。

**核心安全假设**：Module-LWE（Learning With Errors）和 Module-SIS（Short Integer Solution）问题在量子计算机下仍然困难。

## 2. 数学基础

### 2.1 多项式环

Dilithium 工作在多项式环上：

```
R_q = Z_q[X] / (X^256 + 1)
```

其中 `q = 8380417 = 2^23 - 2^13 + 1`。这个环中的元素是次数 < 256 的多项式，系数在 Z_q 中，乘法在模 X^256+1 下进行。

选择 q 的原因：
- q ≡ 1 (mod 512)，保证 NTT（数论变换）可以高效执行
- q 足够大以保证安全性，又足够小以保证效率

### 2.2 模格结构

Dilithium-2 使用 **4×4 模格**（k=4 行，l=4 列）。公钥矩阵 A 是 k×l 的多项式矩阵：

```
A ∈ R_q^{4×4}
```

每个条目是 R_q 中的一个多项式。这比单个多项式的格（如 NTRU）提供了更灵活的参数选择。

### 2.3 NTT 加速

多项式乘法是最核心的运算。朴素乘法是 O(n²)，通过 NTT（Number Theoretic Transform）可以加速到 O(n log n)：

```
a * b = NTT^{-1}(NTT(a) ⊙ NTT(b))
```

其中 ⊙ 是逐点乘法。NTT 本质上是有限域上的 FFT。

Dilithium 使用 **Montgomery 乘法**进一步优化模运算，避免昂贵的除法：

```
MontgomeryReduce(a) = a * 2^{-32} mod q
```

预计算的 twiddle factor 存储为 Montgomery 形式：`zeta_mont = zeta * 2^32 mod q`。

## 3. 方案构造

### 3.1 密钥生成

```
KeyGen():
  1. 随机种子 ζ ← {0,1}^256
  2. (ρ, ρ', K) = SHAKE256(ζ)     // 展开为三个种子
  3. A = ExpandA(ρ)                // 从 ρ 确定性生成 4×4 矩阵 A（存在 NTT 域）
  4. s1 ← S_η^l, s2 ← S_η^k      // 从 ρ' 采样短向量（系数在 [-η, η]，η=2）
  5. t = A · NTT(s1) + s2          // 计算公钥向量 t
  6. (t1, t0) = Power2Round(t, d)  // 分解：t = t1 · 2^d + t0，d=13
  7. tr = H(ρ || pack(t1))         // 公钥哈希

  公钥 pk = (ρ, t1)               // ~1312 字节
  私钥 sk = (ρ, K, tr, s1, s2, t0)
```

**Power2Round** 的作用：t1 保留高位（公开），t0 保留低位（秘密）。这减小了公钥大小，同时低位 t0 在签名时用于提示（hint）计算。

### 3.2 签名

签名使用 **Fiat-Shamir with Aborts** 范式——"先承诺、后挑战、再响应"，加上拒绝采样保证安全：

```
Sign(sk, msg):
  1. A = ExpandA(ρ)
  2. μ = H(tr || msg)                        // 消息摘要
  3. ρ' = H(K || μ)                          // 确定性随机性（可选随机化）
  4. κ = 0                                   // 重试计数器

  loop:                                      // 拒绝采样循环
    5. y ← S_{γ1}^l using (ρ', κ)           // 采样掩码向量，γ1 = 2^17
    6. w = A · NTT(y)                        // 承诺
    7. (w1, w0) = Decompose(w, 2γ2)         // 分解为高位和低位
    8. c̃ = H(μ || pack(w1))                 // 挑战哈希
    9. c = SampleChallenge(c̃)               // 从 c̃ 采样稀疏多项式（τ=39 个 ±1）
    10. z = y + c · s1                       // 响应

    // 拒绝条件（保证零知识性）：
    11. if ||z||∞ ≥ γ1 - β: continue         // z 太大，泄露 s1 信息
    12. if ||w0 - c·s2||∞ ≥ γ2 - β: continue // 低位检查
    13. if ||c·t0||∞ ≥ γ2: continue          // t0 检查

    14. 计算 hint h = MakeHint(w - c·s2, c·t0)  // 提示向量
    15. if weight(h) > ω: continue           // 提示太多

    16. return σ = (c̃, z, h)                // 签名 ~2420 字节
```

**为什么需要拒绝采样？**

如果直接输出 z = y + c·s1，攻击者可以从多个签名中统计分析出 s1。拒绝采样确保 z 的分布独立于 s1——只有当 z "看起来像"随机的均匀分布时才接受，否则重试。平均需要约 4~7 次重试。

### 3.3 验证

```
Verify(pk, msg, σ = (c̃, z, h)):
  1. if ||z||∞ ≥ γ1 - β: return false      // 范数检查
  2. A = ExpandA(ρ)
  3. μ = H(H(pk) || msg)
  4. c = SampleChallenge(c̃)
  5. w' = A·NTT(z) - c·t1·2^d              // 重新计算 w'
  6. w1' = UseHint(w', h)                   // 用 hint 恢复 w1
  7. c̃' = H(μ || pack(w1'))
  8. return c̃ == c̃'                        // 哈希比较
```

**验证为什么成立？**

```
w' = A·z - c·t1·2^d
   = A·(y + c·s1) - c·t1·2^d
   = A·y + c·(A·s1) - c·t1·2^d
   = A·y + c·(t - s2) - c·t1·2^d          // 因为 t = A·s1 + s2
   = A·y + c·t - c·s2 - c·t1·2^d
   = A·y + c·(t1·2^d + t0) - c·s2 - c·t1·2^d
   = A·y + c·t0 - c·s2
   = w + c·t0 - c·s2                       // 因为 w = A·y

HighBits(w') = HighBits(w + c·t0 - c·s2)
```

Hint h 的作用是让验证者在不知道 t0, s2 的情况下从 w' 恢复出正确的 HighBits(w) = w1。

## 4. Decompose 和 Hint 机制

这是 Dilithium 最精巧的设计之一。

### 4.1 Decompose

将值 a 分解为高位 a1 和低位 a0：

```
a = a1 · 2γ2 + a0，其中 -γ2 < a0 ≤ γ2
```

γ2 = (q-1)/88 = 95232。高位 a1 的范围是 [0, 43]，只需 6 bit 表示。

### 4.2 MakeHint / UseHint

问题：验证者计算的是 `w' = w - c·s2 + c·t0`，而不是 `w`。如何从 w' 恢复 HighBits(w)？

答案：签名者计算一个 **hint 向量 h**，告诉验证者 HighBits 在哪些位置需要调整。h 是一个 0/1 向量，只在进位发生的位置为 1。

```
MakeHint(a0, a1): 如果加上 a0 会改变 a1 的 HighBits → 返回 1
UseHint(r, h):    根据 h 调整 HighBits(r)
```

## 5. 安全性分析

### 5.1 安全归约

Dilithium 的安全性归约到两个格问题：

- **Module-LWE**：给定 (A, t = A·s1 + s2)，找到短向量 s1, s2 是困难的 → 密钥不可伪造
- **Module-SIS**：给定 A，找到短向量 z 使得 A·z ≈ 0 是困难的 → 签名不可伪造

### 5.2 抗量子安全性

经典计算机上最好的格攻击算法是 BKZ（Block Korkine-Zolotarev），时间复杂度约 2^{0.292·d}，其中 d 是格维度。

量子计算机上，Grover 搜索可以将搜索问题加速到平方根，但格问题的核心（格基约化）不是搜索问题，Shor 算法也不适用。已知最好的量子格攻击仅比经典快常数因子。

Dilithium-2 的参数选择保证了 **128 位量子安全**。

### 5.3 参数总结

| 参数 | 值 | 含义 |
|------|-----|------|
| q | 8380417 | 模数 |
| n | 256 | 多项式次数 |
| k, l | 4, 4 | 模格维度 |
| η | 2 | 秘密多项式系数范围 [-2, 2] |
| γ1 | 2^17 | 掩码范围 |
| γ2 | 95232 | 分解参数 |
| τ | 39 | 挑战多项式中 ±1 的个数 |
| β | 78 | τ·η，范数界 |
| d | 13 | Power2Round 丢弃位数 |
| ω | 80 | hint 中最多 1 的个数 |

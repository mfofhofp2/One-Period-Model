# Bond, Stock, Option — No-Arbitrage Bounds via Grassmann Algebra

## Setup

Three instruments: a bond with realized return $R$, a stock with current price $s$, and a call option with strike $K$ and current price $v$.

$$x = O + \rho + s\sigma + v\kappa$$

$$X(\omega) = O + R\rho + \omega\sigma + \max(\omega-K,0)\kappa$$

Here $O$ is a fixed origin, and $\rho, \sigma, \kappa$ are the direction vectors for the bond, stock, and option respectively (i.e. $\rho = \text{bond}-O$, and so on). $x$ is the current price point; $X(\omega)$ is the terminal price point in outcome $\omega$.

We pick three reference outcomes $L < K < H$. Since $\max(\omega-K,0)=0$ whenever $\omega \le K$:

$$X(L) = O + R\rho + L\sigma$$
$$X(K) = O + R\rho + K\sigma$$
$$X(H) = O + R\rho + H\sigma + (H-K)\kappa$$

Only $X(H)$ carries an option component, since only in that outcome does the call finish in the money.

**Goal**: write $x = a·O + b·X(L) + c·X(H) + d·X(K)$. By the fundamental theorem of asset pricing, the model is arbitrage-free iff $b,c,d \ge 0$ (the coefficient on the origin $O$ is unconstrained — a cone doesn't restrict movement along its own apex direction).

## The computational tool

Let $\Pi = O·X(L)·X(H)·X(K)$. Replacing the $j$-th factor with $x$ gives $\Pi_j(x)$, and each coefficient is recovered by

$$b = \frac{\Pi_1(x)}{\Pi}, \qquad c = \frac{\Pi_2(x)}{\Pi}, \qquad d = \frac{\Pi_3(x)}{\Pi}$$

Since every factor carries a copy of $O$, and $O\wedge O = 0$, any term drawing the $O$-part from more than one factor vanishes. So once the leading $O$ is factored out, we only need the **direction parts**:

$$X(L)\to R\rho+L\sigma, \quad X(K)\to R\rho+K\sigma, \quad X(H)\to R\rho+H\sigma+(H-K)\kappa, \quad x\to \rho+s\sigma+v\kappa$$

Any wedge product where a vector repeats (e.g. $\rho\wedge\rho$, or $\rho\sigma\rho$ after reordering) is automatically zero. The only surviving 3-fold products are permutations of $\rho\sigma\kappa$, with sign given by the permutation's parity:

$$\rho\sigma\kappa=1,\quad \rho\kappa\sigma=-1,\quad \sigma\kappa\rho=+1,\quad \kappa\sigma\rho=-1,\quad \sigma\rho\kappa=-1,\quad \kappa\rho\sigma=+1$$

## Step 1 — Compute Π = O·X(L)·X(H)·X(K)

**First multiplication** — $(R\rho+L\sigma)\wedge(R\rho+H\sigma+(H-K)\kappa)$, expanded term by term (6 cross-terms from distributing $2\times3$):

| Term | Result | Kept? |
|---|---|---|
| $R\rho\wedge R\rho$ | $R^2\rho\rho=0$ | vanishes ($\rho$ repeats) |
| $R\rho\wedge H\sigma$ | $RH\rho\sigma$ | kept |
| $R\rho\wedge(H-K)\kappa$ | $R(H-K)\rho\kappa$ | kept |
| $L\sigma\wedge R\rho$ | $LR\sigma\rho=-LR\rho\sigma$ | kept |
| $L\sigma\wedge H\sigma$ | $LH\sigma\sigma=0$ | vanishes ($\sigma$ repeats) |
| $L\sigma\wedge(H-K)\kappa$ | $L(H-K)\sigma\kappa$ | kept |

Collecting the surviving terms:

$$(R\rho+L\sigma)\wedge(R\rho+H\sigma+(H-K)\kappa) = (RH-LR)\rho\sigma + R(H-K)\rho\kappa + L(H-K)\sigma\kappa$$

**Second multiplication** — wedge the above with $(R\rho+K\sigma)$ (6 more cross-terms):

| Term | Result | Kept? |
|---|---|---|
| $(RH-LR)\rho\sigma\wedge R\rho$ | $\rho\sigma\rho$ | vanishes ($\rho$ repeats) |
| $(RH-LR)\rho\sigma\wedge K\sigma$ | $\rho\sigma\sigma$ | vanishes ($\sigma$ repeats) |
| $R(H-K)\rho\kappa\wedge R\rho$ | $\rho\kappa\rho$ | vanishes ($\rho$ repeats) |
| $R(H-K)\rho\kappa\wedge K\sigma$ | $R(H-K)K\rho\kappa\sigma = -R(H-K)K\rho\sigma\kappa$ | kept |
| $L(H-K)\sigma\kappa\wedge R\rho$ | $L(H-K)R\sigma\kappa\rho = +L(H-K)R\rho\sigma\kappa$ | kept |
| $L(H-K)\sigma\kappa\wedge K\sigma$ | $\sigma\kappa\sigma$ | vanishes ($\sigma$ repeats) |

Summing the two surviving terms and factoring out $(H-K)$:

$$\Pi = \big[-R(H-K)K + L(H-K)R\big]O\rho\sigma\kappa = R(H-K)(L-K)O\rho\sigma\kappa$$

Since $L<K<H$ and $R>0$: $(H-K)>0$, $(L-K)<0$, so **$\Pi$ is negative**. This sign carries through every coefficient computed by dividing by $\Pi$ below.

## Step 2 — Compute Π1(x) = O·x·X(H)·X(K) (coefficient of X(L))

$$(\rho+s\sigma+v\kappa)\wedge(R\rho+H\sigma+(H-K)\kappa) = (H-sR)\rho\sigma + (H-K-vR)\rho\kappa + (s(H-K)-vH)\sigma\kappa$$

Wedge with $(R\rho+K\sigma)$, keeping only the cross terms:

$$(H-K-vR)\rho\kappa\wedge K\sigma = -(H-K-vR)K\rho\sigma\kappa$$
$$(s(H-K)-vH)\sigma\kappa\wedge R\rho = (s(H-K)-vH)R\rho\sigma\kappa$$

$$\Pi_1(x) = \big[-(H-K-vR)K + (s(H-K)-vH)R\big]O\rho\sigma\kappa = (H-K)\big[R(s-v)-K\big]O\rho\sigma\kappa$$

Dividing by $\Pi$:

$$b = \frac{R(s-v)-K}{R(L-K)}$$

Since $R(L-K)<0$ (denominator negative), coefficient $\ge0$ requires the numerator $\le0$:

$$R(s-v)-K\le0  \Longrightarrow  v \ge s-\frac{K}{R}$$

This is the discounted-intrinsic-value lower bound — it doesn't depend on $L$ or $H$ at all, which is a useful check: this bound must hold in *any* arbitrage-free model, not just this trinomial one.

## Step 3 — Compute Π2(x) = O·X(L)·x·X(K) (coefficient of X(H))

$$(R\rho+L\sigma)\wedge(\rho+s\sigma+v\kappa) = (Rs-L)\rho\sigma + Rv\rho\kappa + Lv\sigma\kappa$$

Wedge with $(R\rho+K\sigma)$:

$$Rv\rho\kappa\wedge K\sigma = -RvK\rho\sigma\kappa \qquad Lv\sigma\kappa\wedge R\rho = LvR\rho\sigma\kappa$$

$$\Pi_2(x) = (LvR - RvK)O\rho\sigma\kappa = vR(L-K)O\rho\sigma\kappa$$

Dividing by $\Pi = R(H-K)(L-K)O\rho\sigma\kappa$ — the $(L-K)$ factor cancels directly, no sign flip needed here:

$$c = \frac{v}{H-K} \ge 0  \Longrightarrow  v \ge 0$$

The trivial bound: a call's payoff is never negative, so its price shouldn't be either.

## Step 4 — Compute Π3(x) = O·X(L)·X(H)·x (coefficient of X(K))

Reusing the intermediate result from Step 1:

$$(R\rho+L\sigma)\wedge(R\rho+H\sigma+(H-K)\kappa) = (RH-LR)\rho\sigma + R(H-K)\rho\kappa + L(H-K)\sigma\kappa$$

Wedge with $(\rho+s\sigma+v\kappa)$:

$$(RH-LR)\rho\sigma\wedge v\kappa = (RH-LR)v\rho\sigma\kappa$$
$$R(H-K)\rho\kappa\wedge s\sigma = -R(H-K)s\rho\sigma\kappa$$
$$L(H-K)\sigma\kappa\wedge\rho = L(H-K)\rho\sigma\kappa$$

$$\Pi_3(x) = \big[R(H-L)v - (H-K)(Rs-L)\big]O\rho\sigma\kappa$$

Dividing by $\Pi$ (negative denominator):

$$d = \frac{R(H-L)v-(H-K)(Rs-L)}{R(H-K)(L-K)}$$

Coefficient $\ge0$ requires the numerator $\le0$:

$$R(H-L)v \le (H-K)(Rs-L)  \Longrightarrow  v \le \frac{(H-K)(Rs-L)}{R(H-L)}$$

## Combined result

$$\boxed{\max\left(0,\ s-\frac{K}{R}\right)  \le  v  \le  \frac{(H-K)(Rs-L)}{R(H-L)}}$$

## Interpretation

- **Lower bound**: the discounted intrinsic value $s-K/R$, floored at $0$. Model-free — holds regardless of what $L,H$ are assumed to be.
- **Upper bound**: depends on the full trinomial structure. If $s$ is at the money ($s=K$) and $R=1$, this reduces to $(H-K)(K-L)/(H-L)$ — the classic "90-100-110" bound.
- Both bounds together define the range of option prices $v$ consistent with no arbitrage in this three-outcome model. A price outside this range implies a static arbitrage trade combining the bond, stock, and option.

## Worked numeric example

Take $L=90, K=100, H=110, R=1, s=100$ (the classic textbook setup):

$$\text{Lower bound} = \max(0,\ 100-100) = 0$$
$$\text{Upper bound} = \frac{(110-100)(1\cdot100-90)}{1\cdot(110-90)} = \frac{10\times10}{20} = 5$$

So the no-arbitrage range for this option's price is $0 \le v \le 5$.

---

# Bond, Stock, Two Options — No-Arbitrage Bounds via Grassmann Algebra

## Setup

Four instruments: a bond with realized return $R$, a stock with current price $s$, and two call options with strikes $K_1<K_2$ and current prices $v_1,v_2$.

$$x = O + \rho + s\sigma + v_1\kappa_1 + v_2\kappa_2$$

$$X(\omega) = O + R\rho + \omega\sigma + \max(\omega-K_1,0)\kappa_1 + \max(\omega-K_2,0)\kappa_2$$

We pick four reference outcomes $L<K_1<K_2<H$:

$$X(L) = O+R\rho+L\sigma \qquad X(K_1) = O+R\rho+K_1\sigma$$
$$X(K_2) = O+R\rho+K_2\sigma+(K_2-K_1)\kappa_1 \qquad X(H) = O+R\rho+H\sigma+(H-K_1)\kappa_1+(H-K_2)\kappa_2$$

We write $x = a·O + b·X(L) + c·X(K_1) + d·X(K_2) + e·X(H)$. Arbitrage-free requires $b,c,d,e\ge0$ (coefficient on $O$ unconstrained). Four instruments, four reference outcomes — this is exactly the "complete market" count, one dimension up from the single-option case.

## The mechanical tool, via determinants

With five points ($O$ plus four references), $\Pi = O·X(L)·X(K_1)·X(K_2)·X(H)$ is a 4-fold wedge in the direction space $(\rho,\sigma,\kappa_1,\kappa_2)$. A wedge product of $n$ vectors, given by their coordinates in a fixed basis, equals the determinant of the matrix whose rows are those coordinates, times the top-form (here $\rho\sigma\kappa_1\kappa_2$). This is the natural extension of the reduction-identity bookkeeping used for three vectors, and keeps calculations tractable as dimension grows — row operations (which leave a determinant unchanged) do the work that manually tracking permutation signs did before.

Coordinates in the basis $(\rho,\sigma,\kappa_1,\kappa_2)$:

$$X(L)\to(R,L,0,0) \quad X(K_1)\to(R,K_1,0,0) \quad X(K_2)\to(R,K_2,K_2-K_1,0) \quad X(H)\to(R,H,H-K_1,H-K_2) \quad x\to(1,s,v_1,v_2)$$

## Step 1 — Compute Π <-> det(X(L),X(K1),X(K2),X(H))

```
| R    L      0        0     |
| R    K1     0        0     |
| R    K2     K2-K1    0     |
| R    H      H-K1     H-K2  |
```

Subtract row 1 from rows 2, 3, 4 (determinant unchanged):

```
| R    L      0        0     |
| 0    K1-L   0        0     |
| 0    K2-L   K2-K1    0     |
| 0    H-L    H-K1     H-K2  |
```

Expanding along column 1, then the resulting matrix is lower-triangular — the determinant is the product of the diagonal entries:

$$\boxed{\Pi = R(K_1-L)(K_2-K_1)(H-K_2)\ O\rho\sigma\kappa_1\kappa_2}$$

Since $L<K_1<K_2<H$ and $R>0$, every factor is positive — **$\Pi$ is positive** here (unlike the single-option case, where it came out negative). No sign flips are needed when dividing by $\Pi$ below.

## Step 2 — Compute Π1(x) <-> det(x,X(K1),X(K2),X(H)) (coefficient b of X(L))

```
| 1    s      v1       v2    |
| R    K1     0        0     |
| R    K2     K2-K1    0     |
| R    H      H-K1     H-K2  |
```

Eliminate column 1 by subtracting R·row 1 from rows 2–4:

```
| 1    s        v1          v2        |
| 0    K1-Rs    -Rv1        -Rv2      |
| 0    K2-Rs    K2-K1-Rv1   -Rv2      |
| 0    H-Rs     H-K1-Rv1    H-K2-Rv2  |
```

Expand along column 1, leaving the $3\times3$ minor. Subtract row 1 (of the minor) from rows 2, 3:

```
| K1-Rs    -Rv1     -Rv2   |
| K2-K1    K2-K1    0      |
| H-K1     H-K1     H-K2   |
```

Expand along column 3 (only two nonzero entries: $-Rv_2$ in row 1, $H-K_2$ in row 3):

$$= -(-Rv_2)\cdot M_1 + (H-K_2)\cdot M_2$$

where the two $2\times2$ minors are:

```
M1 = | K2-K1   K2-K1 |      M2 = | K1-Rs   -Rv1  |
     | H-K1    H-K1  |           | K2-K1   K2-K1 |
```

The first $2\times2$ has identical columns (both $K_2-K_1$ vs $H-K_1$, but as a column pair $(K_2-K_1,H-K_1)$ appearing twice) — its determinant is $0$. The second:

$$(K_1-Rs)(K_2-K_1) - (-Rv_1)(K_2-K_1) = (K_2-K_1)\big[(K_1-Rs)+Rv_1\big] = (K_2-K_1)\big[K_1-R(s-v_1)\big]$$

So the $3\times3$ determinant is $(H-K_2)(K_2-K_1)\big[K_1-R(s-v_1)\big]$, and:

$$\det\big(x,X(K_1),X(K_2),X(H)\big) = (K_2-K_1)(H-K_2)\big[K_1-R(s-v_1)\big]$$

Dividing by $\Pi$:

$$b = \frac{K_1-R(s-v_1)}{R(K_1-L)}$$

$$b\ge0  \Longrightarrow  v_1 \ge s-\frac{K_1}{R}$$

This is exactly the single-option discounted-intrinsic-value lower bound, reappearing unchanged — and notably, $v_2$ doesn't appear at all. A strong consistency check.

## Step 3 — Compute Π2(x) <-> det(X(L),x,X(K2),X(H)) (coefficient c of X(K1))

```
| R    L      0        0     |
| 1    s      v1       v2    |
| R    K2     K2-K1    0     |
| R    H      H-K1     H-K2  |
```

Use row 2 (leading $1$) to eliminate column 1: subtract R·row 2 from rows 1, 3, 4:

```
| 0    L-Rs     -Rv1        -Rv2      |
| 1    s        v1          v2        |
| 0    K2-Rs    K2-K1-Rv1   -Rv2      |
| 0    H-Rs     H-K1-Rv1    H-K2-Rv2  |
```

Expand along column 1 — only row 2 is nonzero, at position $(2,1)$, giving sign $(-1)^{2+1}=-1$:

```
      | L-Rs    -Rv1        -Rv2      |
 = -  | K2-Rs   K2-K1-Rv1   -Rv2      |
      | H-Rs    H-K1-Rv1    H-K2-Rv2  |
```

Subtract row 1 from rows 2, 3:

```
      | L-Rs    -Rv1     -Rv2   |
 = -  | K2-L    K2-K1    0      |
      | H-L     H-K1     H-K2   |
```

Expand along column 3 (nonzero entries: $-Rv_2$ in row 1, $H-K_2$ in row 3):

$$= -\Big[(-Rv_2)\cdot M_1 + (H-K_2)\cdot M_2\Big]$$

where:

```
M1 = | K2-L   K2-K1 |      M2 = | L-Rs   -Rv1  |
     | H-L    H-K1  |           | K2-L   K2-K1 |
```

First $2\times2$: $(K_2-L)(H-K_1)-(K_2-K_1)(H-L)$. Expanding both products and simplifying (the $HK_2, HL, LK_2$ cross-terms cancel in pairs) gives $(K_1-L)(H-K_2)$.

Second $2\times2$: $(L-Rs)(K_2-K_1)-(-Rv_1)(K_2-L) = (L-Rs)(K_2-K_1)+Rv_1(K_2-L)$.

$$= -\Big[-Rv_2(K_1-L)(H-K_2) + (H-K_2)\big[(L-Rs)(K_2-K_1)+Rv_1(K_2-L)\big]\Big]$$

$$= (H-K_2)\Big[Rv_2(K_1-L) + (Rs-L)(K_2-K_1) - Rv_1(K_2-L)\Big]$$

Dividing by $\Pi = R(K_1-L)(K_2-K_1)(H-K_2)$:

$$c = \frac{R\big[v_2(K_1-L)-v_1(K_2-L)\big] + (Rs-L)(K_2-K_1)}{R(K_1-L)(K_2-K_1)}$$

Splitting the $(Rs-L)(K_2-K_1)$ term (its $R$ cancels against the denominator's $R$) gives the cleaner form:

$$c = \frac{(K_2-K_1)(s-L/R) - \big[(K_2-L)v_1-(K_1-L)v_2\big]}{(K_1-L)(K_2-K_1)}$$

$c\ge0$ rearranges to:

$$(K_2-L)v_1 - (K_1-L)v_2  \le  (K_2-K_1)\left(s-\frac{L}{R}\right)$$

A joint bound linking $v_1$ and $v_2$ — it only exists once both strikes are in the picture.

## Step 4 — Compute Π3(x) <-> det(X(L),X(K1),x,X(H)) (coefficient d of X(K2))

```
| R    L      0        0     |
| R    K1     0        0     |
| 1    s      v1       v2    |
| R    H      H-K1     H-K2  |
```

Subtract row 1 from row 2: row 2 becomes $(0,K_1-L,0,0)$. Expand along row 2 — only entry $(2,2)=K_1-L$ is nonzero, sign $(-1)^{2+2}=+1$:

```
                | R    0        0     |
 = (K1-L)  ·    | 1    v1       v2    |
                | R    H-K1     H-K2  |
```

Expand this $3\times3$ along the top row (only the leading $R$ is nonzero):

```
                          | v1     v2   |
 = (K1-L) · R  ·          | H-K1   H-K2 |
   = (K1-L)R [v1(H-K2) - v2(H-K1)]
```

Dividing by $\Pi = R(K_1-L)(K_2-K_1)(H-K_2)$:

$$d = \frac{v_1(H-K_2)-v_2(H-K_1)}{(K_2-K_1)(H-K_2)}$$

$$d\ge0  \Longrightarrow  \boxed{v_1 \ge v_2\cdot\frac{H-K_1}{H-K_2}}$$

This is the cross-strike relationship — it falls straight out of the coefficient calculation, no separate argument required. Since $H-K_1>H-K_2>0$, the ratio exceeds $1$, so this is strictly stronger than $v_1\ge v_2$ (the lower-strike call must be worth more): it pins down exactly how much more.

## Step 5 — Compute Π4(x) <-> det(X(L),X(K1),X(K2),x) (coefficient e of X(H))

```
| R    L      0        0     |
| R    K1     0        0     |
| R    K2     K2-K1    0     |
| 1    s      v1       v2    |
```

Same first step: subtract row 1 from row 2, giving $(0,K_1-L,0,0)$. Expanding along row 2:

```
                | R    0        0   |
 = (K1-L)  ·    | R    K2-K1    0   |
                | 1    v1       v2  |
```

This is lower-triangular after noting only the diagonal-ish structure contributes — expanding along the top row (only the leading $R$ survives, since the other two entries are $0$):

```
                          | K2-K1   0   |
 = (K1-L) · R  ·          | v1      v2  |
   = (K1-L) R (K2-K1) v2
```

Dividing by $\Pi$:

$$e = \frac{v_2}{H-K_2} \ge 0  \Longrightarrow  v_2\ge0$$

The trivial bound — a call's payoff is never negative.

## Combined result

$$v_2\ge0 \qquad v_1 \ge v_2\cdot\frac{H-K_1}{H-K_2} \qquad v_1\ge s-\frac{K_1}{R} \qquad (K_2-L)v_1-(K_1-L)v_2\le(K_2-K_1)\left(s-\frac{L}{R}\right)$$

## Interpretation

- **$v_2\ge0$**: trivial, same as the single-option case.
- **Cross-strike bound**: obtained "for free" from the $d\ge0$ coefficient — no separate proof needed to see that the lower-strike option must cost more, and by how much.
- **$v_1\ge s-K_1/R$**: the discounted-intrinsic-value lower bound reappears unchanged from the single-option case, independent of $v_2$ — a strong consistency check.
- **Joint bound**: links $v_1$ and $v_2$ together; only appears once both strikes coexist.

## Worked numeric example

Take $L=80,\ K_1=100,\ K_2=110,\ H=130,\ R=1,\ s=105$:

$$v_2\ge0 \qquad v_1\ge1.5v_2 \qquad v_1\ge5 \qquad 30v_1-20v_2\le250$$

Testing $(v_1,v_2)=(10,0)$: the first three bounds pass, but $30(10)-20(0)=300>250$ — **violated**. Solving the underlying linear system directly for this $(v_1,v_2)$ confirms $c=-0.25<0$, matching exactly.

Testing $(v_1,v_2)=(8,3)$: $3\ge0$ ✓, $8\ge4.5$ ✓, $8\ge5$ ✓, $30(8)-20(3)=240\le250$ ✓ — all four pass, so this pair is arbitrage-free.

## Why this matters for extending further

The same recipe — one more reference outcome, one more row/column in the determinant — extends to any number of strikes $K_1<\cdots<K_n$. The row-reduction pattern seen throughout (subtracting row 1 to zero out a column, then repeating on the shrinking minor) stays mechanical regardless of how large the matrix gets, though the joint bound analogous to Step 3's condition will involve more cross-terms with each additional strike.

---

# Bond, Stock, Three Options — No-Arbitrage Bounds via Grassmann Algebra

## Setup

Five instruments: a bond with realized return $R$, a stock with current price $s$, and three call options with strikes $K_1<K_2<K_3$ and current prices $v_1,v_2,v_3$.

$$x = O + \rho + s\sigma + v_1\kappa_1 + v_2\kappa_2 + v_3\kappa_3$$

$$X(\omega) = O + R\rho + \omega\sigma + \max(\omega-K_1,0)\kappa_1 + \max(\omega-K_2,0)\kappa_2 + \max(\omega-K_3,0)\kappa_3$$

We pick five reference outcomes $L<K_1<K_2<K_3<H$:

$$X(L) = O+R\rho+L\sigma \qquad X(K_1) = O+R\rho+K_1\sigma$$
$$X(K_2) = O+R\rho+K_2\sigma+(K_2-K_1)\kappa_1$$
$$X(K_3) = O+R\rho+K_3\sigma+(K_3-K_1)\kappa_1+(K_3-K_2)\kappa_2$$
$$X(H) = O+R\rho+H\sigma+(H-K_1)\kappa_1+(H-K_2)\kappa_2+(H-K_3)\kappa_3$$

We write $x = a·O + b·X(L) + c·X(K_1) + d·X(K_2) + e·X(K_3) + f·X(H)$. Arbitrage-free requires $b,c,d,e,f\ge0$. Five instruments, five reference outcomes — again exactly the "complete market" count, one dimension up from the two-strike case.

## The mechanical tool

$\Pi = O·X(L)·X(K_1)·X(K_2)·X(K_3)·X(H)$ is now a 5-fold wedge in $(\rho,\sigma,\kappa_1,\kappa_2,\kappa_3)$, equal to the determinant of the $5\times5$ coordinate matrix times the top-form $\rho\sigma\kappa_1\kappa_2\kappa_3$. Coordinates:

$$X(L)\to(R,L,0,0,0) \quad X(K_1)\to(R,K_1,0,0,0) \quad X(K_2)\to(R,K_2,K_2-K_1,0,0)$$
$$X(K_3)\to(R,K_3,K_3-K_1,K_3-K_2,0) \quad X(H)\to(R,H,H-K_1,H-K_2,H-K_3) \quad x\to(1,s,v_1,v_2,v_3)$$

The same row-reduction pattern used for two strikes extends directly: subtracting row 1 (or whichever row carries the leading $1$ or $R$) to zero out a column, then repeating on the shrinking minor.

## Step 1 — Compute Π <-> det(X(L),X(K1),X(K2),X(K3),X(H))

```
| R   L    0       0       0     |
| R   K1   0       0       0     |
| R   K2   K2-K1   0       0     |
| R   K3   K3-K1   K3-K2   0     |
| R   H    H-K1    H-K2    H-K3  |
```

Subtract row 1 from rows 2–5:

```
| R   L      0       0       0     |
| 0   K1-L   0       0       0     |
| 0   K2-L   K2-K1   0       0     |
| 0   K3-L   K3-K1   K3-K2   0     |
| 0   H-L    H-K1    H-K2    H-K3  |
```

Lower-triangular after expanding along column 1 — determinant is the product of the diagonal:

$$\boxed{\Pi = R(K_1-L)(K_2-K_1)(K_3-K_2)(H-K_3)\ O\rho\sigma\kappa_1\kappa_2\kappa_3}$$

Every factor positive (as in the two-strike case) — **$\Pi>0$**, no sign flips needed below.

## Step 2 — Coefficient b of X(L)

$\det\big(x,X(K_1),X(K_2),X(K_3),X(H)\big)$: eliminate column 1 using row 1 ($x$ has leading coordinate $1$), leaving a $4\times4$ minor whose first column is $(K_1-Rs, K_2-Rs, K_3-Rs, H-Rs)$ and whose remaining columns mirror the $\kappa_1,\kappa_2,\kappa_3$ structure of $X(K_1),X(K_2),X(K_3),X(H)$ shifted by $-Rv_1,-Rv_2,-Rv_3$ in row 1. Subtracting row 1 from rows 2–4 of this minor produces zeros in column 1 for all but the top row, collapsing the calculation to exactly the same triangular pattern as Step 1, but with $K_1-Rs$ (adjusted by $-R(s-v_1)$ tracking through) replacing $K_1-L$:

$$\det\big(x,X(K_1),X(K_2),X(K_3),X(H)\big) = (K_2-K_1)(K_3-K_2)(H-K_3)\big[K_1-R(s-v_1)\big]$$

Dividing by $\Pi$:

$$b = \frac{K_1-R(s-v_1)}{R(K_1-L)}$$

$$b\ge0  \Longrightarrow  v_1 \ge s-\frac{K_1}{R}$$

Exactly the single-option discounted-intrinsic-value bound — unchanged by the presence of $v_2,v_3$. Consistent with both prior cases.

## Step 3 — Coefficient c of X(K1)

$\det\big(X(L),x,X(K_2),X(K_3),X(H)\big)$: use row 2 (leading $1$) to eliminate column 1, giving a sign $(-1)^{2+1}=-1$, then reduce the remaining $4\times4$ the same way as Step 3 in the two-strike case. Carrying the algebra through (the $v_2,v_3$ contributions again arrange themselves into a single $(K_1-L)$-weighted term, exactly mirroring the two-strike pattern):

$$c = \frac{(K_2-K_1)(s-L/R) - \big[(K_2-L)v_1-(K_1-L)v_2\big]}{(K_1-L)(K_2-K_1)}$$

$$c\ge0  \Longrightarrow  (K_2-L)v_1-(K_1-L)v_2 \le (K_2-K_1)\left(s-\frac{L}{R}\right)$$

Identical in form to the two-strike joint bound — $v_3$ does not appear. This condition only "sees" its immediate neighbors $L,K_1,K_2$.

## Step 4 — Coefficient d of X(K2) (the new, middle condition)

$\det\big(X(L),X(K_1),x,X(K_3),X(H)\big)$: rows 1 and 2 ($X(L)$ and $X(K_1)$) agree in every $\kappa$-coordinate (both zero), so subtracting row 1 from row 2 collapses to $(0,K_1-L,0,0,0)$, exactly as in the two-strike case. Working through the resulting reduced determinant:

$$d = \frac{v_1(K_3-K_2) - v_2(K_3-K_1) + v_3(K_2-K_1)}{(K_2-K_1)(K_3-K_2)}$$

$$d\ge0  \Longrightarrow  \boxed{v_2 \le \frac{(K_3-K_2)v_1+(K_2-K_1)v_3}{K_3-K_1}}$$

**This is the butterfly-spread convexity condition** — and it involves only $K_1,K_2,K_3,v_1,v_2,v_3$. Neither $L,H,s$ nor the fourth instrument's price appears. It falls out of the same mechanical coefficient-extraction as every other bound, with no separate construction.

## Step 5 — Coefficient e of X(K3)

By the mirror-image argument to Step 3: rows for $X(L),X(K_1),X(K_2)$ all agree in the $\kappa_3$ coordinate (all zero), collapsing the determinant the same way. Working through:

$$e = \frac{v_2(H-K_3)-v_3(H-K_2)}{(K_3-K_2)(H-K_3)}$$

$$e\ge0  \Longrightarrow  \boxed{v_2 \ge v_3\cdot\frac{H-K_2}{H-K_3}}$$

Same shape as the two-strike cross-strike condition, shifted one strike to the right. Only involves $K_2,K_3,H,v_2,v_3$ — $L,K_1,v_1$ drop out entirely.

## Step 6 — Coefficient f of X(H)

Rows $X(L),X(K_1),X(K_2),X(K_3)$ all agree in the last coordinate before $x$ — same collapse pattern:

$$f = \frac{v_3}{H-K_3} \ge 0  \Longrightarrow  v_3\ge0$$

The trivial bound.

## Combined result

$$v_1\ge s-\frac{K_1}{R} \qquad (K_2-L)v_1-(K_1-L)v_2\le(K_2-K_1)\left(s-\frac{L}{R}\right)$$

$$v_2 \le \frac{(K_3-K_2)v_1+(K_2-K_1)v_3}{K_3-K_1} \qquad v_2 \ge v_3\cdot\frac{H-K_2}{H-K_3} \qquad v_3\ge0$$

## Interpretation — the locality pattern

Laying all five results side by side:

| Reference point | Condition | Depends on |
|---|---|---|
| $L$ | $v_1\ge s-K_1/R$ | $s,K_1$ only |
| $K_1$ | $(K_2-L)v_1-(K_1-L)v_2\le(K_2-K_1)(s-L/R)$ | $L,s,K_1,K_2$ only |
| $K_2$ | $v_2\le\dfrac{(K_3-K_2)v_1+(K_2-K_1)v_3}{K_3-K_1}$ (butterfly) | $K_1,K_2,K_3$ only |
| $K_3$ | $v_2\ge v_3\cdot\dfrac{H-K_2}{H-K_3}$ | $K_2,K_3,H$ only |
| $H$ | $v_3\ge0$ | trivial |

Every condition depends only on its immediate neighbors in the chain $L,K_1,K_2,K_3,H$ — none reaches further than one or two steps away. The single interior strike $K_2$ produces exactly one butterfly-spread condition, matching the general pattern: $n$ strikes give $n+2$ total conditions, of which $n-2$ are interior (butterfly-type) conditions.

## Worked numeric example

Take $L=70,\ K_1=90,\ K_2=100,\ K_3=110,\ H=130,\ R=1,\ s=105$. Check $(v_1,v_2,v_3)=(20,14,9)$:

$$v_1\ge s-K_1/R:\quad 20\ge15 \checkmark$$
$$(K_2-L)v_1-(K_1-L)v_2\le(K_2-K_1)(s-L/R):\quad 30(20)-20(14)=320\le10(35)=350 \checkmark$$
$$v_2\le\frac{(K_3-K_2)v_1+(K_2-K_1)v_3}{K_3-K_1}:\quad 14\le\frac{10(20)+10(9)}{20}=14.5 \checkmark$$
$$v_2\ge v_3\cdot\frac{H-K_2}{H-K_3}:\quad 14\ge9\cdot\frac{30}{20}=13.5 \checkmark$$
$$v_3\ge0:\quad 9\ge0 \checkmark$$

All five conditions pass — $(20,14,9)$ is a consistent, arbitrage-free set of option prices for this model.

---

# General $n$ Strikes — The Full Characterization

The three worked examples above suggest a conjecture: every no-arbitrage coefficient is *local*, depending only on its immediate neighbors in the chain $L,K_1,\ldots,K_n,H$, with the interior coefficients taking the form of a discrete-convexity ("butterfly") condition. This section states and proves that conjecture for arbitrary $n$.

## Setup, general $n$

Reference outcomes $L=Y_0 < Y_1 < \cdots < Y_n < Y_{n+1}=H$, where $Y_j=K_j$ for $1\le j\le n$. Generators:

$$x = O+\rho+s\sigma+\sum_{j=1}^n v_j\kappa_j, \qquad X(Y_p) = O+R\rho+Y_p\sigma+\sum_{m=1}^n\max(Y_p-K_m,0)\,\kappa_m$$

Write $x = a\,O + \sum_{p=0}^{n+1} b_p\,X(Y_p)$ (so $b_0=b_L$, $b_{n+1}=b_H$). Arbitrage-free $\iff$ every $b_p\ge0$.

## Theorem 1. Closed-form determinant

$$\Pi := O\wedge X(Y_0)\wedge\cdots\wedge X(Y_{n+1}) = R\,(K_1-L)\left(\prod_{i=2}^n(K_i-K_{i-1})\right)(H-K_n)\ \, O\rho\sigma\kappa_1\cdots\kappa_n$$

**Proof.** By induction on $n$. The $(n{+}1)$-strike determinant factors as $(H-K_n)$ times the $n$-strike determinant with $K_n$ substituted for $H$ — the same "peel off the last strike" row-reduction carried out by hand for $n=2,3$ above, generalized to arbitrary $n$. Concretely, expanding along the row corresponding to $\kappa_n$ isolates a single nonzero cofactor (since every generator except $X(H)$ has zero $\kappa_n$-coordinate), whose value is exactly $H-K_n$, and the remaining minor is precisely the $n$-strike matrix with $H$ replaced by $K_n$. $\blacksquare$

Since $L<K_1<\cdots<K_n<H$ and $R>0$, every factor is strictly positive, so $\Pi\ne0$ for every $n$.

## Theorem 2. FTAP at the cone level

Since $\Pi\ne0$, the $n+2$ generators are linearly independent, giving a linear isomorphism $\mathbb{R}^{n+2}\to\mathbb{R}^{n+2}$ under which the non-negative orthant maps to a closed convex cone $C$. A separating-hyperplane argument then gives, for every $n$:

$$\text{no-arbitrage at } x \iff x\in C$$

This is the abstract Fundamental Theorem of Asset Pricing, valid for arbitrary $n$. It confirms *existence* of a non-negative decomposition but does not yet hand back the explicit inequalities that Theorem 3 supplies.

## Theorem 3. The full local characterization

$$\text{no-arbitrage} \iff b_p\ge0 \text{ for every } p=0,1,\ldots,n+1$$

with each condition depending only on immediate neighbors — the exact conditions, and which ones exist, depend on $n$:

$$\boxed{\begin{array}{ll}
b_L\ge0 &\iff v_1\ge s-\dfrac{K_1}{R} \\[6pt]
b_{K_1}\ge0 &\iff (K_2-L)v_1-(K_1-L)v_2 \le (K_2-K_1)\left(s-\dfrac{L}{R}\right) \qquad (n\ge2)\\[6pt]
b_{K_i}\ge0 &\iff v_i \le \dfrac{(K_{i+1}-K_i)v_{i-1}+(K_i-K_{i-1})v_{i+1}}{K_{i+1}-K_{i-1}}, \quad 2\le i\le n-1 \qquad (n\ge3)\\[6pt]
b_{K_n}\ge0 &\iff v_{n-1}\ge v_n\cdot\dfrac{H-K_{n-1}}{H-K_n} \qquad (n\ge2)\\[6pt]
b_H\ge0 &\iff v_n\ge0
\end{array}}$$

For $n=1$ there is no interior condition and $K_1=K_n$ is simultaneously the near-$L$ and near-$H$ boundary strike, so the single condition on $b_{K_1}$ is derived on its own — it is exactly the $n=1$ worked example's coefficient $d$ above. For $n=2$, both boundary formulas apply directly and there is still no interior condition (matching the two-strike worked example above exactly). The general boundary and interior formulas are stated for $n\ge3$, matching the three-strike worked example as the base case.

### Proof: a direct linear-system argument

Rather than expand a large determinant and its cofactors for each coefficient separately — the natural but unwieldy generalization of the hand-computations above — the cleanest route is to solve the defining linear system directly.

**Setting up the system.** Matching coordinates in $x = a\,O+\sum_{p=0}^{n+1}b_p X(Y_p)$ against $x=O+\rho+s\sigma+\sum_m v_m\kappa_m$ and $X(Y_p)=O+R\rho+Y_p\sigma+\sum_m\max(Y_p-K_m,0)\kappa_m$ gives, coordinate by coordinate:

$$\rho:\quad \sum_{p=0}^{n+1} b_p = \frac1R \qquad\qquad \sigma:\quad \sum_{p=0}^{n+1} b_p Y_p = s$$

$$\kappa_m\ (1\le m\le n):\quad \sum_{p:\,Y_p>K_m} b_p\,(Y_p-K_m) = v_m$$

(the $O$-coordinate equation only pins down $a$ and plays no further role). This is $n+2$ equations in the $n+2$ unknowns $b_0,\ldots,b_{n+1}$; Theorem 1 guarantees a unique solution.

**Solving by substitution.** Define the tail sum $S_m := \sum_{p=m+1}^{n+1} b_p$ for $m=0,1,\ldots,n$ (so $S_n=b_{n+1}=b_H$, and $S_0=\frac1R-b_L$). The $\kappa_m$-equation reads $v_m = T_m - K_mS_m$ where $T_m:=\sum_{p>m}b_pY_p$.

*Step 1 (rightmost strike).* For $m=n$: only $p=n+1$ survives the sum, giving $v_n = b_{n+1}(H-K_n)$, i.e.

$$b_H = \frac{v_n}{H-K_n} \implies b_H\ge0\iff v_n\ge0$$

*Step 2 (consecutive-difference telescoping).* For $1\le m\le n-1$, subtract the $\kappa_{m+1}$-equation from the $\kappa_m$-equation. Using $T_m-T_{m+1}=b_{m+1}Y_{m+1}$, $S_m-S_{m+1}=b_{m+1}$, and — the single fact doing all the work — $Y_{m+1}=K_{m+1}$ for $m+1\le n$:

$$v_m-v_{m+1} = (T_m-T_{m+1}) - K_mS_m + K_{m+1}S_{m+1} = b_{m+1}(Y_{m+1}-K_{m+1}) + (K_{m+1}-K_m)S_m = (K_{m+1}-K_m)S_m$$

(the $b_{m+1}$ term vanishes identically, since $Y_{m+1}-K_{m+1}=0$). So:

$$S_m = \frac{v_m-v_{m+1}}{K_{m+1}-K_m}, \qquad 1\le m\le n-1$$

*Step 3 (interior coefficients).* For $2\le i\le n-1$: $b_{K_i} = b_i = S_{i-1}-S_i$, and both are given by Step 2:

$$b_{K_i} = \frac{v_{i-1}-v_i}{K_i-K_{i-1}} - \frac{v_i-v_{i+1}}{K_{i+1}-K_i} = \frac{(K_{i+1}-K_i)v_{i-1} - (K_{i+1}-K_{i-1})v_i + (K_i-K_{i-1})v_{i+1}}{(K_i-K_{i-1})(K_{i+1}-K_i)}$$

Since both factors in the denominator are positive, $b_{K_i}\ge0$ iff the stated butterfly inequality holds.

*Step 4 ($b_{K_n}$).* $b_{K_n}=b_n=S_{n-1}-S_n$:

$$b_{K_n} = \frac{v_{n-1}-v_n}{K_n-K_{n-1}} - \frac{v_n}{H-K_n} = \frac{(H-K_n)v_{n-1}-(H-K_{n-1})v_n}{(K_n-K_{n-1})(H-K_n)}$$

Both factors in the denominator are positive, giving $b_{K_n}\ge0 \iff v_{n-1}\ge v_n\cdot\frac{H-K_{n-1}}{H-K_n}$.

*Step 5 ($b_L$).* Use the $\sigma$-equation together with $S_0=\frac1R-b_L$ and the $\kappa_1$-equation. Writing the $\sigma$-equation as $s=b_LL+T_0$ with $T_0=b_1K_1+T_1$, and the $\kappa_1$-equation as $T_1=v_1+K_1S_1$, and $b_1=S_0-S_1$:

$$s = b_LL + (S_0-S_1)K_1 + v_1+K_1S_1 = b_LL+S_0K_1+v_1$$

so $s-v_1=b_LL+S_0K_1$. Using $S_0=\frac1R-b_L$:

$$s-v_1 = b_L(L-K_1)+\frac{K_1}{R} \implies b_L = \frac{K_1-R(s-v_1)}{R(K_1-L)} \implies b_L\ge0\iff v_1\ge s-\frac{K_1}{R}$$

*Step 6 ($b_{K_1}$, when $n\ge2$).* $b_{K_1}=b_1=S_0-S_1 = \left(\frac1R-b_L\right) - \frac{v_1-v_2}{K_2-K_1}$. Substituting $b_L$ from Step 5 and simplifying:

$$b_{K_1} = \frac{(K_1-K_2)(Rs-L)+Rv_1(K_2-L)-Rv_2(K_1-L)}{R(K_1-L)(K_1-K_2)}$$

The denominator is negative ($R(K_1-L)>0$, $(K_1-K_2)<0$), giving

$$b_{K_1}\ge0 \iff (K_2-L)v_1-(K_1-L)v_2 \le (K_2-K_1)\left(s-\frac{L}{R}\right). \qquad\blacksquare$$

Every step above is ordinary substitution in a linear system — no determinant cofactors, no case-by-case sign bookkeeping. The entire mechanism reduces to one identity, used repeatedly: $Y_{m+1}=K_{m+1}$ makes the $b_{m+1}$-term in the difference of consecutive $\kappa$-equations vanish automatically, which is *why* the system telescopes into the local tail-sum recursion above.

### Recovering $n=1,2,3$

Setting $n=1$: only $b_L$ and $b_{K_1}$ (computed directly from $\rho,\sigma,\kappa_1$, since there is no $S_1$ intermediate) and $b_H$ exist — matching the first worked example's $b,d,c$ exactly. Setting $n=2,3$ and reading off Steps 1, 4, 5, 6 (and, for $n=3$, the single interior application of Step 3 at $i=2$) reproduces the second and third worked examples' combined results term for term.

---

## Conclusion

The central claim of this note is a **locality theorem**: in the $n$-strike single-period model, every no-arbitrage condition — on the discounted bond, the stock, and each of the $n$ options — depends on at most three consecutive points in the ordered reference chain $L,K_1,\ldots,K_n,H$, regardless of how large $n$ is. This is proved, not merely observed: the mechanism is a single algebraic fact (that $Y_{m+1}=K_{m+1}$ exactly, so that a piecewise-linear payoff's discrete second difference vanishes outside its own support), applied uniformly via a telescoping tail-sum substitution in the defining linear system.

Three things follow from this.

**First**, the abstract statement of the Fundamental Theorem of Asset Pricing — no-arbitrage iff the price lies in a certain cone — is upgraded, for this model and for any number of strikes, to a fully explicit system of $n+2$ checkable inequalities (Theorem 3), each as easy to verify as a single arithmetic comparison between neighboring option prices and strikes.

**Second**, the interior condition recovered here is exactly the classical discrete-convexity (butterfly-spread) no-arbitrage condition on option prices across strikes, familiar from the state-price-density literature (Breeden–Litzenberger). What this note adds is not the qualitative fact of convexity, which is well known, but (i) a single uniform derivation supplying that condition *together with* the two boundary families (near $L$ and near $H$) under one proof, and (ii) a demonstration that the locality is exact and general — three neighboring points and no more, for every interior strike, for every $n$ — rather than a pattern checked case by case.

**Third**, the practical upshot is computational: checking no-arbitrage for a ladder of $n$ strikes is an $O(n)$ task — one local check per reference point — rather than requiring the evaluation of an $(n+2)$-dimensional determinant or the solution of a general linear program. The determinant $\Pi$ and the cone-membership characterization (Theorems 1–2) remain useful as the underlying reason the local system is well-posed and has a unique solution, but the day-to-day content of the result is the boxed system in Theorem 3.

The derivation throughout — from the single-option warm-up to the general-$n$ theorem — follows one consistent method: represent price points as vectors in an exterior algebra over the market's basic instruments, extract the unique non-negative combination via Cramer's rule (or, in the general case, by direct substitution in the equivalent linear system), and read off no-arbitrage as non-negativity of that combination's coefficients. That the same method scales, without modification, from three instruments to arbitrarily many, and yields a closed form rather than an ever-growing case analysis, is the paper's principal technical contribution.

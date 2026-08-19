# Bond, Stock, Option — No-Arbitrage Bounds via Grassmann Algebra

## Setup

Three instruments: a bond with realized return $R$, a stock with current price $s$, and a call option with strike $K$ and current price $v$.

$$x = O + \rho + s\sigma + v\kappa$$

$$X(\omega) = O + R\rho + \omega\sigma + \max(\omega-K,0)\kappa$$

Here $O$ is a fixed origin, and $\rho, \sigma, \kappa$ are the direction vectors for the bond, stock, and option respectively (i.e. $\rho = $ bond$-O$, and so on). $x$ is the current price point; $X(\omega)$ is the terminal price point in outcome $\omega$.

We pick three reference outcomes $L < K < H$. Since $\max(\omega-K,0)=0$ whenever $\omega \le K$:

$$X(L) = O + R\rho + L\sigma$$
$$X(K) = O + R\rho + K\sigma$$
$$X(H) = O + R\rho + H\sigma + (H-K)\kappa$$

Only $X(H)$ carries an option component, since only in that outcome does the call finish in the money.

**Goal**: write $x = a\,O + b\,X(L) + c\,X(H) + d\,X(K)$. By the fundamental theorem of asset pricing, the model is arbitrage-free iff $b,c,d \ge 0$ (the coefficient on the origin $O$ is unconstrained — a cone doesn't restrict movement along its own apex direction).

## The computational tool

Let $\Pi = O\,X(L)\,X(H)\,X(K)$. Replacing the $j$-th factor with $x$ gives $\Pi_j(x)$, and each coefficient is recovered by

$$b = \frac{\Pi_1(x)}{\Pi}, \qquad c = \frac{\Pi_2(x)}{\Pi}, \qquad d = \frac{\Pi_3(x)}{\Pi}$$

Since every factor carries a copy of $O$, and $O\wedge O = 0$, any term drawing the $O$-part from more than one factor vanishes. So once the leading $O$ is factored out, we only need the **direction parts**:

$$X(L)\to R\rho+L\sigma, \quad X(K)\to R\rho+K\sigma, \quad X(H)\to R\rho+H\sigma+(H-K)\kappa, \quad x\to \rho+s\sigma+v\kappa$$

Any wedge product where a vector repeats (e.g. $\rho\wedge\rho$, or $\rho\sigma\rho$ after reordering) is automatically zero. The only surviving 3-fold products are permutations of $\rho\sigma\kappa$, with sign given by the permutation's parity:

$$\rho\sigma\kappa=1,\quad \rho\kappa\sigma=-1,\quad \sigma\kappa\rho=+1,\quad \kappa\sigma\rho=-1,\quad \sigma\rho\kappa=-1,\quad \kappa\rho\sigma=+1$$

## Step 1 — Compute $\Pi = O\,X(L)\,X(H)\,X(K)$

**First multiplication** — $(R\rho+L\sigma)\wedge(R\rho+H\sigma+(H-K)\kappa)$, expanded term by term (6 cross-terms from distributing $2\times3$):

| Term | Result | Kept? |
|---|---|---|
| $R\rho\wedge R\rho$ | $R^2\,\rho\rho=0$ | vanishes ($\rho$ repeats) |
| $R\rho\wedge H\sigma$ | $RH\,\rho\sigma$ | kept |
| $R\rho\wedge(H-K)\kappa$ | $R(H-K)\,\rho\kappa$ | kept |
| $L\sigma\wedge R\rho$ | $LR\,\sigma\rho=-LR\,\rho\sigma$ | kept |
| $L\sigma\wedge H\sigma$ | $LH\,\sigma\sigma=0$ | vanishes ($\sigma$ repeats) |
| $L\sigma\wedge(H-K)\kappa$ | $L(H-K)\,\sigma\kappa$ | kept |

Collecting the surviving terms:

$$(R\rho+L\sigma)\wedge(R\rho+H\sigma+(H-K)\kappa) = (RH-LR)\,\rho\sigma + R(H-K)\,\rho\kappa + L(H-K)\,\sigma\kappa$$

**Second multiplication** — wedge the above with $(R\rho+K\sigma)$ (6 more cross-terms):

| Term | Result | Kept? |
|---|---|---|
| $(RH-LR)\rho\sigma\wedge R\rho$ | $\rho\sigma\rho$ | vanishes ($\rho$ repeats) |
| $(RH-LR)\rho\sigma\wedge K\sigma$ | $\rho\sigma\sigma$ | vanishes ($\sigma$ repeats) |
| $R(H-K)\rho\kappa\wedge R\rho$ | $\rho\kappa\rho$ | vanishes ($\rho$ repeats) |
| $R(H-K)\rho\kappa\wedge K\sigma$ | $R(H-K)K\,\rho\kappa\sigma = -R(H-K)K\,\rho\sigma\kappa$ | kept |
| $L(H-K)\sigma\kappa\wedge R\rho$ | $L(H-K)R\,\sigma\kappa\rho = +L(H-K)R\,\rho\sigma\kappa$ | kept |
| $L(H-K)\sigma\kappa\wedge K\sigma$ | $\sigma\kappa\sigma$ | vanishes ($\sigma$ repeats) |

Summing the two surviving terms and factoring out $(H-K)$:

$$\Pi = \big[-R(H-K)K + L(H-K)R\big]\,O\rho\sigma\kappa = R(H-K)(L-K)\,O\rho\sigma\kappa$$

Since $L<K<H$ and $R>0$: $(H-K)>0$, $(L-K)<0$, so **$\Pi$ is negative**. This sign carries through every coefficient computed by dividing by $\Pi$ below.

## Step 2 — Compute $\Pi_1(x) = O\,x\,X(H)\,X(K)$ (coefficient of $X(L)$)

$$(\rho+s\sigma+v\kappa)\wedge(R\rho+H\sigma+(H-K)\kappa) = (H-sR)\rho\sigma + (H-K-vR)\rho\kappa + (s(H-K)-vH)\sigma\kappa$$

Wedge with $(R\rho+K\sigma)$, keeping only the cross terms:

$$(H-K-vR)\rho\kappa\wedge K\sigma = -(H-K-vR)K\,\rho\sigma\kappa$$
$$(s(H-K)-vH)\sigma\kappa\wedge R\rho = (s(H-K)-vH)R\,\rho\sigma\kappa$$

$$\Pi_1(x) = \big[-(H-K-vR)K + (s(H-K)-vH)R\big]\,O\rho\sigma\kappa = (H-K)\big[R(s-v)-K\big]\,O\rho\sigma\kappa$$

Dividing by $\Pi$:

$$b = \frac{R(s-v)-K}{R(L-K)}$$

Since $R(L-K)<0$ (denominator negative), coefficient $\ge0$ requires the numerator $\le0$:

$$R(s-v)-K\le0 \;\;\Longrightarrow\;\; v \ge s-\frac{K}{R}$$

This is the discounted-intrinsic-value lower bound — it doesn't depend on $L$ or $H$ at all, which is a useful check: this bound must hold in *any* arbitrage-free model, not just this trinomial one.

## Step 3 — Compute $\Pi_2(x) = O\,X(L)\,x\,X(K)$ (coefficient of $X(H)$)

$$(R\rho+L\sigma)\wedge(\rho+s\sigma+v\kappa) = (Rs-L)\rho\sigma + Rv\,\rho\kappa + Lv\,\sigma\kappa$$

Wedge with $(R\rho+K\sigma)$:

$$Rv\,\rho\kappa\wedge K\sigma = -RvK\,\rho\sigma\kappa \qquad Lv\,\sigma\kappa\wedge R\rho = LvR\,\rho\sigma\kappa$$

$$\Pi_2(x) = (LvR - RvK)\,O\rho\sigma\kappa = vR(L-K)\,O\rho\sigma\kappa$$

Dividing by $\Pi = R(H-K)(L-K)\,O\rho\sigma\kappa$ — the $(L-K)$ factor cancels directly, no sign flip needed here:

$$c = \frac{v}{H-K} \ge 0 \;\;\Longrightarrow\;\; v \ge 0$$

The trivial bound: a call's payoff is never negative, so its price shouldn't be either.

## Step 4 — Compute $\Pi_3(x) = O\,X(L)\,X(H)\,x$ (coefficient of $X(K)$)

Reusing the intermediate result from Step 1:

$$(R\rho+L\sigma)\wedge(R\rho+H\sigma+(H-K)\kappa) = (RH-LR)\rho\sigma + R(H-K)\rho\kappa + L(H-K)\sigma\kappa$$

Wedge with $(\rho+s\sigma+v\kappa)$:

$$(RH-LR)\rho\sigma\wedge v\kappa = (RH-LR)v\,\rho\sigma\kappa$$
$$R(H-K)\rho\kappa\wedge s\sigma = -R(H-K)s\,\rho\sigma\kappa$$
$$L(H-K)\sigma\kappa\wedge\rho = L(H-K)\,\rho\sigma\kappa$$

$$\Pi_3(x) = \big[R(H-L)v - (H-K)(Rs-L)\big]\,O\rho\sigma\kappa$$

Dividing by $\Pi$ (negative denominator):

$$d = \frac{R(H-L)v-(H-K)(Rs-L)}{R(H-K)(L-K)}$$

Coefficient $\ge0$ requires the numerator $\le0$:

$$R(H-L)v \le (H-K)(Rs-L) \;\;\Longrightarrow\;\; v \le \frac{(H-K)(Rs-L)}{R(H-L)}$$

## Combined result

$$\boxed{\max\left(0,\ s-\frac{K}{R}\right) \;\le\; v \;\le\; \frac{(H-K)(Rs-L)}{R(H-L)}}$$

## Interpretation

- **Lower bound**: the discounted intrinsic value $s-K/R$, floored at $0$. Model-free — holds regardless of what $L,H$ are assumed to be.
- **Upper bound**: depends on the full trinomial structure. If $s$ is at the money ($s=K$) and $R=1$, this reduces to $(H-K)(K-L)/(H-L)$ — the classic "90-100-110" bound.
- Both bounds together define the range of option prices $v$ consistent with no arbitrage in this three-outcome model. A price outside this range implies a static arbitrage trade combining the bond, stock, and option.

## Worked numeric example

Take $L=90, K=100, H=110, R=1, s=100$ (the classic textbook setup):

$$\text{Lower bound} = \max(0,\ 100-100) = 0$$
$$\text{Upper bound} = \frac{(110-100)(1\cdot100-90)}{1\cdot(110-90)} = \frac{10\times10}{20} = 5$$

So the no-arbitrage range for this option's price is $0 \le v \le 5$.


# Bond, Stock, Two Options — No-Arbitrage Bounds via Grassmann Algebra

## Setup

Four instruments: a bond with realized return $R$, a stock with current price $s$, and two call options with strikes $K_1<K_2$ and current prices $v_1,v_2$.

$$x = O + \rho + s\sigma + v_1\kappa_1 + v_2\kappa_2$$

$$X(\omega) = O + R\rho + \omega\sigma + \max(\omega-K_1,0)\kappa_1 + \max(\omega-K_2,0)\kappa_2$$

We pick four reference outcomes $L<K_1<K_2<H$:

$$X(L) = O+R\rho+L\sigma \qquad X(K_1) = O+R\rho+K_1\sigma$$
$$X(K_2) = O+R\rho+K_2\sigma+(K_2-K_1)\kappa_1 \qquad X(H) = O+R\rho+H\sigma+(H-K_1)\kappa_1+(H-K_2)\kappa_2$$

We write $x = a\,O + b\,X(L) + c\,X(K_1) + d\,X(K_2) + e\,X(H)$. Arbitrage-free requires $b,c,d,e\ge0$ (coefficient on $O$ unconstrained). Four instruments, four reference outcomes — this is exactly the "complete market" count, one dimension up from the single-option case.

## The mechanical tool, via determinants

With five points ($O$ plus four references), $\Pi = O\,X(L)\,X(K_1)\,X(K_2)\,X(H)$ is a 4-fold wedge in the direction space $(\rho,\sigma,\kappa_1,\kappa_2)$. A wedge product of $n$ vectors, given by their coordinates in a fixed basis, equals the determinant of the matrix whose rows are those coordinates, times the top-form (here $\rho\sigma\kappa_1\kappa_2$). This is the natural extension of the reduction-identity bookkeeping used for three vectors, and keeps calculations tractable as dimension grows — row operations (which leave a determinant unchanged) do the work that manually tracking permutation signs did before.

Coordinates in the basis $(\rho,\sigma,\kappa_1,\kappa_2)$:

$$X(L)\to(R,L,0,0) \quad X(K_1)\to(R,K_1,0,0) \quad X(K_2)\to(R,K_2,K_2-K_1,0) \quad X(H)\to(R,H,H-K_1,H-K_2) \quad x\to(1,s,v_1,v_2)$$

## Step 1 — Compute $\Pi \leftrightarrow \det\big(X(L),X(K_1),X(K_2),X(H)\big)$

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

$$\boxed{\Pi = R(K_1-L)(K_2-K_1)(H-K_2)\; O\rho\sigma\kappa_1\kappa_2}$$

Since $L<K_1<K_2<H$ and $R>0$, every factor is positive — **$\Pi$ is positive** here (unlike the single-option case, where it came out negative). No sign flips are needed when dividing by $\Pi$ below.

## Step 2 — Compute $\Pi_1(x) \leftrightarrow \det\big(x,X(K_1),X(K_2),X(H)\big)$ (coefficient $b$ of $X(L)$)

```
| 1    s      v1       v2    |
| R    K1     0        0     |
| R    K2     K2-K1    0     |
| R    H      H-K1     H-K2  |
```

Eliminate column 1 by subtracting $R\times$row 1 from rows 2–4:

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

$$b\ge0 \;\;\Longrightarrow\;\; v_1 \ge s-\frac{K_1}{R}$$

This is exactly the single-option discounted-intrinsic-value lower bound, reappearing unchanged — and notably, $v_2$ doesn't appear at all. A strong consistency check.

## Step 3 — Compute $\Pi_2(x) \leftrightarrow \det\big(X(L),x,X(K_2),X(H)\big)$ (coefficient $c$ of $X(K_1)$)

```
| R    L      0        0     |
| 1    s      v1       v2    |
| R    K2     K2-K1    0     |
| R    H      H-K1     H-K2  |
```

Use row 2 (leading $1$) to eliminate column 1: subtract $R\times$row 2 from rows 1, 3, 4:

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

$$(K_2-L)v_1 - (K_1-L)v_2 \;\le\; (K_2-K_1)\left(s-\frac{L}{R}\right)$$

A joint bound linking $v_1$ and $v_2$ — it only exists once both strikes are in the picture.

## Step 4 — Compute $\Pi_3(x) \leftrightarrow \det\big(X(L),X(K_1),x,X(H)\big)$ (coefficient $d$ of $X(K_2)$)

```
| R    L      0        0     |
| R    K1     0        0     |
| 1    s      v1       v2    |
| R    H      H-K1     H-K2  |
```

Subtract row 1 from row 2: row 2 becomes $(0,K_1-L,0,0)$. Expand along row 2 — only entry $(2,2)=K_1-L$ is nonzero, sign $(-1)^{2+2}=+1$:

```
                | R    0        0     |
 = (K1-L)  ×    | 1    v1       v2    |
                | R    H-K1     H-K2  |
```

Expand this $3\times3$ along the top row (only the leading $R$ is nonzero):

```
                          | v1     v2   |
 = (K1-L) · R  ×          | H-K1   H-K2 |
   = (K1-L)R [v1(H-K2) - v2(H-K1)]
```

Dividing by $\Pi = R(K_1-L)(K_2-K_1)(H-K_2)$:

$$d = \frac{v_1(H-K_2)-v_2(H-K_1)}{(K_2-K_1)(H-K_2)}$$

$$d\ge0 \;\;\Longrightarrow\;\; \boxed{v_1 \ge v_2\cdot\frac{H-K_1}{H-K_2}}$$

This is the cross-strike relationship — it falls straight out of the coefficient calculation, no separate argument required. Since $H-K_1>H-K_2>0$, the ratio exceeds $1$, so this is strictly stronger than "$v_1\ge v_2$" (the lower-strike call must be worth more): it pins down exactly how much more.

## Step 5 — Compute $\Pi_4(x) \leftrightarrow \det\big(X(L),X(K_1),X(K_2),x\big)$ (coefficient $e$ of $X(H)$)

```
| R    L      0        0     |
| R    K1     0        0     |
| R    K2     K2-K1    0     |
| 1    s      v1       v2    |
```

Same first step: subtract row 1 from row 2, giving $(0,K_1-L,0,0)$. Expanding along row 2:

```
                | R    0        0   |
 = (K1-L)  ×    | R    K2-K1    0   |
                | 1    v1       v2  |
```

This is lower-triangular after noting only the diagonal-ish structure contributes — expanding along the top row (only the leading $R$ survives, since the other two entries are $0$):

```
                          | K2-K1   0   |
 = (K1-L) · R  ×          | v1      v2  |
   = (K1-L) R (K2-K1) v2
```

Dividing by $\Pi$:

$$e = \frac{v_2}{H-K_2} \ge 0 \;\;\Longrightarrow\;\; v_2\ge0$$

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

$$v_2\ge0 \qquad v_1\ge1.5\,v_2 \qquad v_1\ge5 \qquad 30v_1-20v_2\le250$$

Testing $(v_1,v_2)=(10,0)$: the first three bounds pass, but $30(10)-20(0)=300>250$ — **violated**. Solving the underlying linear system directly for this $(v_1,v_2)$ confirms $c=-0.25<0$, matching exactly.

Testing $(v_1,v_2)=(8,3)$: $3\ge0$ ✓, $8\ge4.5$ ✓, $8\ge5$ ✓, $30(8)-20(3)=240\le250$ ✓ — all four pass, so this pair is arbitrage-free.

## Why this matters for extending further

The same recipe — one more reference outcome, one more row/column in the determinant — extends to any number of strikes $K_1<\cdots<K_n$. The row-reduction pattern seen throughout (subtracting row 1 to zero out a column, then repeating on the shrinking minor) stays mechanical regardless of how large the matrix gets, though the joint bound analogous to Step 3's condition will involve more cross-terms with each additional strike.


# Bond, Stock, Two Options — No-Arbitrage Bounds via Grassmann Algebra

## Setup

Four instruments: a bond with realized return $R$, a stock with current price $s$, and two call options with strikes $K_1<K_2$ and current prices $v_1,v_2$.

$$x = O + \rho + s\sigma + v_1\kappa_1 + v_2\kappa_2$$

$$X(\omega) = O + R\rho + \omega\sigma + \max(\omega-K_1,0)\kappa_1 + \max(\omega-K_2,0)\kappa_2$$

We pick four reference outcomes $L<K_1<K_2<H$:

$$X(L) = O+R\rho+L\sigma \qquad X(K_1) = O+R\rho+K_1\sigma$$
$$X(K_2) = O+R\rho+K_2\sigma+(K_2-K_1)\kappa_1 \qquad X(H) = O+R\rho+H\sigma+(H-K_1)\kappa_1+(H-K_2)\kappa_2$$

We write $x = a\,O + b\,X(L) + c\,X(K_1) + d\,X(K_2) + e\,X(H)$. Arbitrage-free requires $b,c,d,e\ge0$ (coefficient on $O$ unconstrained). Four instruments, four reference outcomes — this is exactly the "complete market" count, one dimension up from the single-option case.

## The mechanical tool, via determinants

With five points ($O$ plus four references), $\Pi = O\,X(L)\,X(K_1)\,X(K_2)\,X(H)$ is a 4-fold wedge in the direction space $(\rho,\sigma,\kappa_1,\kappa_2)$. A wedge product of $n$ vectors, given by their coordinates in a fixed basis, equals the determinant of the matrix whose rows are those coordinates, times the top-form (here $\rho\sigma\kappa_1\kappa_2$). This is the natural extension of the reduction-identity bookkeeping used for three vectors, and keeps calculations tractable as dimension grows — row operations (which leave a determinant unchanged) do the work that manually tracking permutation signs did before.

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

$$\boxed{\Pi = R(K_1-L)(K_2-K_1)(H-K_2)\; O\rho\sigma\kappa_1\kappa_2}$$

Since $L<K_1<K_2<H$ and $R>0$, every factor is positive — **$\Pi$ is positive** here (unlike the single-option case, where it came out negative). No sign flips are needed when dividing by $\Pi$ below.

## Step 2 — Compute Π1(x) <-> det(x,X(K1),X(K2),X(H)) (coefficient b of X(L))

```
| 1    s      v1       v2    |
| R    K1     0        0     |
| R    K2     K2-K1    0     |
| R    H      H-K1     H-K2  |
```

Eliminate column 1 by subtracting $R\times$row 1 from rows 2–4:

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

$$b\ge0 \;\;\Longrightarrow\;\; v_1 \ge s-\frac{K_1}{R}$$

This is exactly the single-option discounted-intrinsic-value lower bound, reappearing unchanged — and notably, $v_2$ doesn't appear at all. A strong consistency check.

## Step 3 — Compute Π2(x) <-> det(X(L),x,X(K2),X(H)) (coefficient c of X(K1))

```
| R    L      0        0     |
| 1    s      v1       v2    |
| R    K2     K2-K1    0     |
| R    H      H-K1     H-K2  |
```

Use row 2 (leading $1$) to eliminate column 1: subtract $R\times$row 2 from rows 1, 3, 4:

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

$$(K_2-L)v_1 - (K_1-L)v_2 \;\le\; (K_2-K_1)\left(s-\frac{L}{R}\right)$$

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
 = (K1-L)  ×    | 1    v1       v2    |
                | R    H-K1     H-K2  |
```

Expand this $3\times3$ along the top row (only the leading $R$ is nonzero):

```
                          | v1     v2   |
 = (K1-L) · R  ×          | H-K1   H-K2 |
   = (K1-L)R [v1(H-K2) - v2(H-K1)]
```

Dividing by $\Pi = R(K_1-L)(K_2-K_1)(H-K_2)$:

$$d = \frac{v_1(H-K_2)-v_2(H-K_1)}{(K_2-K_1)(H-K_2)}$$

$$d\ge0 \;\;\Longrightarrow\;\; \boxed{v_1 \ge v_2\cdot\frac{H-K_1}{H-K_2}}$$

This is the cross-strike relationship — it falls straight out of the coefficient calculation, no separate argument required. Since $H-K_1>H-K_2>0$, the ratio exceeds $1$, so this is strictly stronger than "$v_1\ge v_2$" (the lower-strike call must be worth more): it pins down exactly how much more.

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
 = (K1-L)  ×    | R    K2-K1    0   |
                | 1    v1       v2  |
```

This is lower-triangular after noting only the diagonal-ish structure contributes — expanding along the top row (only the leading $R$ survives, since the other two entries are $0$):

```
                          | K2-K1   0   |
 = (K1-L) · R  ×          | v1      v2  |
   = (K1-L) R (K2-K1) v2
```

Dividing by $\Pi$:

$$e = \frac{v_2}{H-K_2} \ge 0 \;\;\Longrightarrow\;\; v_2\ge0$$

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

$$v_2\ge0 \qquad v_1\ge1.5\,v_2 \qquad v_1\ge5 \qquad 30v_1-20v_2\le250$$

Testing $(v_1,v_2)=(10,0)$: the first three bounds pass, but $30(10)-20(0)=300>250$ — **violated**. Solving the underlying linear system directly for this $(v_1,v_2)$ confirms $c=-0.25<0$, matching exactly.

Testing $(v_1,v_2)=(8,3)$: $3\ge0$ ✓, $8\ge4.5$ ✓, $8\ge5$ ✓, $30(8)-20(3)=240\le250$ ✓ — all four pass, so this pair is arbitrage-free.

## Why this matters for extending further

The same recipe — one more reference outcome, one more row/column in the determinant — extends to any number of strikes $K_1<\cdots<K_n$. The row-reduction pattern seen throughout (subtracting row 1 to zero out a column, then repeating on the shrinking minor) stays mechanical regardless of how large the matrix gets, though the joint bound analogous to Step 3's condition will involve more cross-terms with each additional strike.

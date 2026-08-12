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

We pick four reference outcomes $L<K_1<K_2<H$. Since each $\max(\cdot,0)$ term is zero below its strike:

$$X(L) = O+R\rho+L\sigma \qquad X(K_1) = O+R\rho+K_1\sigma$$
$$X(K_2) = O+R\rho+K_2\sigma+(K_2-K_1)\kappa_1 \qquad X(H) = O+R\rho+H\sigma+(H-K_1)\kappa_1+(H-K_2)\kappa_2$$

Only $X(K_2)$ and $X(H)$ carry option components — $X(L)$ and $X(K_1)$ sit below both strikes, so neither option is in the money there.

**Goal**: write $x = a\,O + b\,X(L) + c\,X(K_1) + d\,X(K_2) + e\,X(H)$. Arbitrage-free requires $b,c,d,e\ge0$.

## The computational tool

Let $\Pi = O\,X(L)\,X(K_1)\,X(K_2)\,X(H)$. Replacing the $j$-th factor with $x$ gives $\Pi_j(x)$, and each coefficient is recovered by dividing by $\Pi$, exactly as before — just with one more factor in the product.

Once the leading $O$ is factored out (using $O\wedge O=0$), we work with direction parts only:

$$X(L)\to R\rho+L\sigma, \quad X(K_1)\to R\rho+K_1\sigma, \quad X(K_2)\to R\rho+K_2\sigma+(K_2-K_1)\kappa_1$$
$$X(H)\to R\rho+H\sigma+(H-K_1)\kappa_1+(H-K_2)\kappa_2, \quad x\to \rho+s\sigma+v_1\kappa_1+v_2\kappa_2$$

Any wedge product with a repeated vector vanishes. The only surviving 4-fold products are permutations of $\rho\sigma\kappa_1\kappa_2$, with sign given by permutation parity. Identities used below:

$$\rho\sigma\kappa_1\kappa_2=1,\quad \rho\kappa_1\sigma=-\rho\sigma\kappa_1,\quad \sigma\kappa_1\rho=+\rho\sigma\kappa_1,\quad \rho\kappa_2\sigma=-\rho\sigma\kappa_2,\quad \sigma\kappa_2\rho=+\rho\sigma\kappa_2$$
$$\rho\kappa_1\kappa_2\sigma=+\rho\sigma\kappa_1\kappa_2, \quad \sigma\kappa_1\kappa_2\rho=-\rho\sigma\kappa_1\kappa_2, \quad \rho\sigma\kappa_2\kappa_1=-\rho\sigma\kappa_1\kappa_2$$

## Step 1 — Compute $\Pi = O\,X(L)\,X(K_1)\,X(K_2)\,X(H)$

**First multiplication** — $(R\rho+L\sigma)\wedge(R\rho+K_1\sigma)$:

| Term | Result | Kept? |
|---|---|---|
| $R\rho\wedge R\rho$ | $0$ | vanishes ($\rho$ repeats) |
| $R\rho\wedge K_1\sigma$ | $RK_1\,\rho\sigma$ | kept |
| $L\sigma\wedge R\rho$ | $-LR\,\rho\sigma$ | kept |
| $L\sigma\wedge K_1\sigma$ | $0$ | vanishes ($\sigma$ repeats) |

$$(R\rho+L\sigma)\wedge(R\rho+K_1\sigma) = R(K_1-L)\,\rho\sigma$$

**Second multiplication** — wedge with $(R\rho+K_2\sigma+(K_2-K_1)\kappa_1)$. Since we already have $\rho\sigma$, only the $\kappa_1$ term avoids a repeat:

$$R(K_1-L)\,\rho\sigma \wedge (K_2-K_1)\kappa_1 = R(K_1-L)(K_2-K_1)\,\rho\sigma\kappa_1$$

**Third multiplication** — wedge with $(R\rho+H\sigma+(H-K_1)\kappa_1+(H-K_2)\kappa_2)$. Only the $\kappa_2$ term avoids a repeat:

$$R(K_1-L)(K_2-K_1)\,\rho\sigma\kappa_1 \wedge (H-K_2)\kappa_2 = R(K_1-L)(K_2-K_1)(H-K_2)\,\rho\sigma\kappa_1\kappa_2$$

$$\boxed{\Pi = R(K_1-L)(K_2-K_1)(H-K_2)\; O\rho\sigma\kappa_1\kappa_2}$$

Since $L<K_1<K_2<H$ and $R>0$, every factor is positive — **$\Pi$ is positive** this time (unlike the single-option case). No sign flips are needed when dividing by $\Pi$ below.

## Step 2 — Compute $\Pi_1(x) = O\,x\,X(K_1)\,X(K_2)\,X(H)$ (coefficient $b$ of $X(L)$)

**First multiplication** — $(\rho+s\sigma+v_1\kappa_1+v_2\kappa_2)\wedge(R\rho+K_1\sigma)$:

| Term | Result | Kept? |
|---|---|---|
| $\rho\wedge R\rho$ | $0$ | vanishes |
| $\rho\wedge K_1\sigma$ | $K_1\,\rho\sigma$ | kept |
| $s\sigma\wedge R\rho$ | $-sR\,\rho\sigma$ | kept |
| $s\sigma\wedge K_1\sigma$ | $0$ | vanishes |
| $v_1\kappa_1\wedge R\rho$ | $-v_1R\,\rho\kappa_1$ | kept |
| $v_1\kappa_1\wedge K_1\sigma$ | $-v_1K_1\,\sigma\kappa_1$ | kept |
| $v_2\kappa_2\wedge R\rho$ | $-v_2R\,\rho\kappa_2$ | kept |
| $v_2\kappa_2\wedge K_1\sigma$ | $-v_2K_1\,\sigma\kappa_2$ | kept |

$$= (K_1-sR)\rho\sigma - v_1R\,\rho\kappa_1 - v_1K_1\,\sigma\kappa_1 - v_2R\,\rho\kappa_2 - v_2K_1\,\sigma\kappa_2$$

**Second multiplication** — wedge each term with $(R\rho+K_2\sigma+(K_2-K_1)\kappa_1)$, keeping only terms that avoid repeats:

| Source term | Surviving cross-term | Result |
|---|---|---|
| $(K_1-sR)\rho\sigma$ | $\wedge(K_2-K_1)\kappa_1$ | $(K_1-sR)(K_2-K_1)\,\rho\sigma\kappa_1$ |
| $-v_1R\,\rho\kappa_1$ | $\wedge K_2\sigma$ | $+v_1RK_2\,\rho\sigma\kappa_1$ |
| $-v_1K_1\,\sigma\kappa_1$ | $\wedge R\rho$ | $-v_1K_1R\,\rho\sigma\kappa_1$ |
| $-v_2R\,\rho\kappa_2$ | $\wedge K_2\sigma,\ \wedge(K_2-K_1)\kappa_1$ | $v_2RK_2\,\rho\sigma\kappa_2 + v_2R(K_2-K_1)\,\rho\kappa_1\kappa_2$ |
| $-v_2K_1\,\sigma\kappa_2$ | $\wedge R\rho,\ \wedge(K_2-K_1)\kappa_1$ | $-v_2K_1R\,\rho\sigma\kappa_2 + v_2K_1(K_2-K_1)\,\sigma\kappa_1\kappa_2$ |

Collecting by basis 3-vector:

$$\rho\sigma\kappa_1:\ \ (K_2-K_1)\big[K_1-R(s-v_1)\big] \qquad \rho\sigma\kappa_2:\ \ v_2R(K_2-K_1)$$
$$\rho\kappa_1\kappa_2:\ \ v_2R(K_2-K_1) \qquad \sigma\kappa_1\kappa_2:\ \ v_2K_1(K_2-K_1)$$

**Third multiplication** — wedge each 3-vector term with $(R\rho+H\sigma+(H-K_1)\kappa_1+(H-K_2)\kappa_2)$, keeping only the complementary basis vector each time:

| Source 3-vector | Complementary term | Result |
|---|---|---|
| $\rho\sigma\kappa_1$ | $\wedge(H-K_2)\kappa_2$ | $(K_2-K_1)[K_1-R(s-v_1)](H-K_2)\,\rho\sigma\kappa_1\kappa_2$ |
| $\rho\sigma\kappa_2$ | $\wedge(H-K_1)\kappa_1 = -(H-K_1)\,\rho\sigma\kappa_1\kappa_2$ | $-v_2R(K_2-K_1)(H-K_1)\,\rho\sigma\kappa_1\kappa_2$ |
| $\rho\kappa_1\kappa_2$ | $\wedge H\sigma = +H\,\rho\sigma\kappa_1\kappa_2$ | $v_2R(K_2-K_1)H\,\rho\sigma\kappa_1\kappa_2$ |
| $\sigma\kappa_1\kappa_2$ | $\wedge R\rho = -R\,\rho\sigma\kappa_1\kappa_2$ | $-v_2K_1(K_2-K_1)R\,\rho\sigma\kappa_1\kappa_2$ |

The three $v_2$-carrying terms combine as $v_2R(K_2-K_1)\big[-(H-K_1)+H-K_1\big]=0$ — **they cancel completely**. What survives is:

$$\Pi_1(x) = (K_2-K_1)(H-K_2)\big[K_1-R(s-v_1)\big]\,O\rho\sigma\kappa_1\kappa_2$$

Dividing by $\Pi$:

$$b = \frac{K_1-R(s-v_1)}{R(K_1-L)}$$

$$b\ge0 \;\;\Longrightarrow\;\; v_1 \ge s-\frac{K_1}{R}$$

This is exactly the single-option discounted-intrinsic-value lower bound, reappearing unchanged — and notably, $v_2$ drops out entirely. A clean consistency check.

## Step 3 — Compute $\Pi_2(x) = O\,X(L)\,x\,X(K_2)\,X(H)$ (coefficient $c$ of $X(K_1)$)

**First multiplication** — $(R\rho+L\sigma)\wedge(\rho+s\sigma+v_1\kappa_1+v_2\kappa_2)$:

| Term | Result | Kept? |
|---|---|---|
| $R\rho\wedge\rho$ | $0$ | vanishes |
| $R\rho\wedge s\sigma$ | $Rs\,\rho\sigma$ | kept |
| $R\rho\wedge v_1\kappa_1$ | $Rv_1\,\rho\kappa_1$ | kept |
| $R\rho\wedge v_2\kappa_2$ | $Rv_2\,\rho\kappa_2$ | kept |
| $L\sigma\wedge\rho$ | $-L\,\rho\sigma$ | kept |
| $L\sigma\wedge s\sigma$ | $0$ | vanishes |
| $L\sigma\wedge v_1\kappa_1$ | $Lv_1\,\sigma\kappa_1$ | kept |
| $L\sigma\wedge v_2\kappa_2$ | $Lv_2\,\sigma\kappa_2$ | kept |

$$= (Rs-L)\rho\sigma + Rv_1\,\rho\kappa_1 + Rv_2\,\rho\kappa_2 + Lv_1\,\sigma\kappa_1 + Lv_2\,\sigma\kappa_2$$

**Second multiplication** — wedge with $(R\rho+K_2\sigma+(K_2-K_1)\kappa_1)$:

| Source term | Surviving cross-term(s) | Result |
|---|---|---|
| $(Rs-L)\rho\sigma$ | $\wedge(K_2-K_1)\kappa_1$ | $(Rs-L)(K_2-K_1)\,\rho\sigma\kappa_1$ |
| $Rv_1\,\rho\kappa_1$ | $\wedge K_2\sigma$ | $-Rv_1K_2\,\rho\sigma\kappa_1$ |
| $Rv_2\,\rho\kappa_2$ | $\wedge K_2\sigma,\ \wedge(K_2-K_1)\kappa_1$ | $-Rv_2K_2\,\rho\sigma\kappa_2 - Rv_2(K_2-K_1)\,\rho\kappa_1\kappa_2$ |
| $Lv_1\,\sigma\kappa_1$ | $\wedge R\rho$ | $Lv_1R\,\rho\sigma\kappa_1$ |
| $Lv_2\,\sigma\kappa_2$ | $\wedge R\rho,\ \wedge(K_2-K_1)\kappa_1$ | $Lv_2R\,\rho\sigma\kappa_2 - Lv_2(K_2-K_1)\,\sigma\kappa_1\kappa_2$ |

Collecting:

$$\rho\sigma\kappa_1:\ \ (K_2-K_1)(Rs-L) - Rv_1(K_2-K_1) \qquad \rho\sigma\kappa_2:\ \ -Rv_2(K_2-K_1)$$
$$\rho\kappa_1\kappa_2:\ \ -Rv_2(K_2-K_1) \qquad \sigma\kappa_1\kappa_2:\ \ -Lv_2(K_2-K_1)$$

**Third multiplication** — wedge with $(R\rho+H\sigma+(H-K_1)\kappa_1+(H-K_2)\kappa_2)$, same complementary-vector rule as Step 2. Collecting all $v_2$ contributions:

$$v_2\text{-terms}: \ -Rv_2(K_2-K_1)(H-K_1) + \big(-Rv_2(K_2-K_1)\big)\cdot H + \big(-Lv_2(K_2-K_1)\big)\cdot(-R)$$
$$= Rv_2(K_2-K_1)\big[-(H-K_1)-H+L\big] = -Rv_2(K_2-K_1)(K_1-L)$$

and the $v_1,s,L$ term:

$$(K_2-K_1)\big[(Rs-L)-Rv_1\big](H-K_2)$$

$$\Pi_2(x) = (K_2-K_1)\Big[(H-K_2)(Rs-L-Rv_1) - R(K_1-L)v_2\Big]\,O\rho\sigma\kappa_1\kappa_2$$

Dividing by $\Pi = R(K_1-L)(K_2-K_1)(H-K_2)\,O\rho\sigma\kappa_1\kappa_2$:

$$c = \frac{(H-K_2)(Rs-L-Rv_1) - R(K_1-L)v_2}{R(K_1-L)(H-K_2)}$$

$c\ge0$ rearranges to:

$$(K_2-L)v_1 - (K_1-L)v_2 \;\le\; (K_2-K_1)\left(s-\frac{L}{R}\right)$$

A joint bound linking $v_1$ and $v_2$ — it only exists once both strikes are in the picture.

## Step 4 — Compute $\Pi_3(x) = O\,X(L)\,X(K_1)\,x\,X(H)$ (coefficient $d$ of $X(K_2)$)

**First multiplication** is the same as Step 1's: $(R\rho+L\sigma)\wedge(R\rho+K_1\sigma) = R(K_1-L)\,\rho\sigma$.

**Second multiplication** — wedge with $x=(\rho+s\sigma+v_1\kappa_1+v_2\kappa_2)$:

| Term | Result | Kept? |
|---|---|---|
| $\rho\sigma\wedge\rho$ | $0$ | vanishes |
| $\rho\sigma\wedge s\sigma$ | $0$ | vanishes |
| $\rho\sigma\wedge v_1\kappa_1$ | $v_1\,\rho\sigma\kappa_1$ | kept |
| $\rho\sigma\wedge v_2\kappa_2$ | $v_2\,\rho\sigma\kappa_2$ | kept |

$$R(K_1-L)\,\rho\sigma \wedge x = R(K_1-L)\big[v_1\,\rho\sigma\kappa_1 + v_2\,\rho\sigma\kappa_2\big]$$

**Third multiplication** — wedge with $(R\rho+H\sigma+(H-K_1)\kappa_1+(H-K_2)\kappa_2)$:

$$v_1\,\rho\sigma\kappa_1\wedge(H-K_2)\kappa_2 = v_1(H-K_2)\,\rho\sigma\kappa_1\kappa_2$$
$$v_2\,\rho\sigma\kappa_2\wedge(H-K_1)\kappa_1 = v_2(H-K_1)\,\rho\sigma\kappa_2\kappa_1 = -v_2(H-K_1)\,\rho\sigma\kappa_1\kappa_2$$

$$\Pi_3(x) = R(K_1-L)\big[v_1(H-K_2)-v_2(H-K_1)\big]\,O\rho\sigma\kappa_1\kappa_2$$

Dividing by $\Pi$:

$$d = \frac{v_1(H-K_2)-v_2(H-K_1)}{(K_2-K_1)(H-K_2)}$$

$$d\ge0 \;\;\Longrightarrow\;\; \boxed{v_1 \ge v_2\cdot\frac{H-K_1}{H-K_2}}$$

This is the cross-strike relationship — it falls straight out of the coefficient calculation, no separate argument required. Since $H-K_1>H-K_2>0$, the ratio exceeds $1$, so this is strictly stronger than "$v_1\ge v_2$" (the lower-strike call must be worth more): it pins down exactly how much more.

## Step 5 — Compute $\Pi_4(x) = O\,X(L)\,X(K_1)\,X(K_2)\,x$ (coefficient $e$ of $X(H)$)

**First two multiplications** match Step 1 exactly, through $R(K_1-L)(K_2-K_1)\,\rho\sigma\kappa_1$.

**Third multiplication** — wedge with $x=(\rho+s\sigma+v_1\kappa_1+v_2\kappa_2)$. Only the $v_2\kappa_2$ term avoids a repeat:

$$R(K_1-L)(K_2-K_1)\,\rho\sigma\kappa_1 \wedge v_2\kappa_2 = R(K_1-L)(K_2-K_1)v_2\,\rho\sigma\kappa_1\kappa_2$$

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

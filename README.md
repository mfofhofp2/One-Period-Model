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

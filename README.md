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

# Bond, Stock, Two Options — Extending to a Second Strike

## Setup

Four instruments: bond (return $R$), stock (price $s$), call option with strike $K_1$ (price $v_1$), call option with strike $K_2$ (price $v_2$), with $K_1<K_2$.

$$x = O + \rho + s\sigma + v_1\kappa_1 + v_2\kappa_2$$
$$X(\omega) = O + R\rho + \omega\sigma + \max(\omega-K_1,0)\kappa_1 + \max(\omega-K_2,0)\kappa_2$$

We take four reference outcomes $L<K_1<K_2<H$:

$$X(L) = O+R\rho+L\sigma \qquad (\text{neither option ITM})$$
$$X(K_1) = O+R\rho+K_1\sigma \qquad (\text{neither option ITM})$$
$$X(K_2) = O+R\rho+K_2\sigma+(K_2-K_1)\kappa_1 \qquad (\text{only the } K_1\text{ call is ITM})$$
$$X(H) = O+R\rho+H\sigma+(H-K_1)\kappa_1+(H-K_2)\kappa_2$$

We write $x = a\,O + b\,X(L) + c\,X(K_1) + d\,X(K_2) + e\,X(H)$. Arbitrage-free requires $b,c,d,e\ge0$ (coefficient on $O$ unconstrained). Four instruments, four reference outcomes — this is exactly the "complete market" count from the single-option case, now one dimension up.

## The mechanical tool, extended

With five points ($O$ plus four references), $\Pi = O\,X(L)\,X(K_1)\,X(K_2)\,X(H)$ is a 4-fold wedge in the direction space $(\rho,\sigma,\kappa_1,\kappa_2)$. As before, factor out $O$ (using $O\wedge O=0$) and work with direction parts only. The only surviving top-form is $\rho\sigma\kappa_1\kappa_2$ (any repeated vector kills a term), with sign given by permutation parity — same rule as before, just one more slot to track.

A convenient shortcut for a 4-fold wedge of vectors, each given by its coordinates in $(\rho,\sigma,\kappa_1,\kappa_2)$: the wedge product equals the determinant of the matrix whose rows are those coordinates, times $\rho\sigma\kappa_1\kappa_2$. This is the natural extension of the $\rho\sigma\kappa$ reduction identities used for three vectors, and keeps the bookkeeping manageable as the dimension grows.

Coordinates (in the basis $\rho,\sigma,\kappa_1,\kappa_2$):

$$X(L)\to(R,L,0,0) \quad X(K_1)\to(R,K_1,0,0) \quad X(K_2)\to(R,K_2,K_2-K_1,0) \quad X(H)\to(R,H,H-K_1,H-K_2) \quad x\to(1,s,v_1,v_2)$$

## Step 1 — Compute $\Pi \leftrightarrow \det\big(X(L),X(K_1),X(K_2),X(H)\big)$

$$\begin{vmatrix} R&L&0&0\\ R&K_1&0&0\\ R&K_2&K_2-K_1&0\\ R&H&H-K_1&H-K_2 \end{vmatrix}$$

Subtract row 1 from rows 2, 3, 4 (determinant unchanged):

$$\begin{vmatrix} R&L&0&0\\ 0&K_1-L&0&0\\ 0&K_2-L&K_2-K_1&0\\ 0&H-L&H-K_1&H-K_2 \end{vmatrix}$$

This is block lower-triangular after expanding along column 1 — the determinant is the product of the diagonal entries:

$$\boxed{\Pi = R(K_1-L)(K_2-K_1)(H-K_2)\; O\rho\sigma\kappa_1\kappa_2}$$

Unlike the single-option case (where $\Pi$ came out negative), here every factor is positive since $L<K_1<K_2<H$ and $R>0$ — so $\Pi>0$, and no sign-flip bookkeeping is needed below.

## Step 2 — Compute $\Pi_1(x) \leftrightarrow \det\big(x,X(K_1),X(K_2),X(H)\big)$ (coefficient $b$ of $X(L)$)

$$\begin{vmatrix} 1&s&v_1&v_2\\ R&K_1&0&0\\ R&K_2&K_2-K_1&0\\ R&H&H-K_1&H-K_2 \end{vmatrix}$$

Eliminate column 1 using row 1 (subtract $R\times$row 1 from rows 2–4):

$$\begin{vmatrix} 1&s&v_1&v_2\\ 0&K_1-Rs&-Rv_1&-Rv_2\\ 0&K_2-Rs&K_2-K_1-Rv_1&-Rv_2\\ 0&H-Rs&H-K_1-Rv_1&H-K_2-Rv_2 \end{vmatrix}$$

Expand along column 1 — reduces to the $3\times3$ minor. Subtract row 2 from rows 3, 4 of that minor:

$$\begin{vmatrix} K_1-Rs&-Rv_1&-Rv_2\\ K_2-K_1&K_2-K_1&0\\ H-K_1&H-K_1&H-K_2-Rv_2+Rv_2 \end{vmatrix} = \begin{vmatrix} K_1-Rs&-Rv_1&-Rv_2\\ K_2-K_1&K_2-K_1&0\\ H-K_1&H-K_1&H-K_2 \end{vmatrix}$$

Expand along the third column (only two nonzero entries):

$$-(-Rv_2)\cdot\big[(K_2-K_1)^2-(K_2-K_1)(H-K_1)\big]\cdots$$

This is getting heavy to track by hand across every term — collecting all contributions (row/column expansion of the $3\times3$, then simplifying) gives:

$$\det\big(x,X(K_1),X(K_2),X(H)\big) = (K_2-K_1)(H-K_2)\big[K_1-R(s-v_1)\big]$$

Dividing by $\Pi$:

$$b = \frac{K_1-R(s-v_1)}{R(K_1-L)}$$

**$b\ge0$** (denominator positive):
$$K_1-R(s-v_1)\ge0 \;\Longrightarrow\; v_1 \ge s-\frac{K_1}{R}$$

This is exactly the single-option discounted-intrinsic-value lower bound, reappearing unchanged with a second option in the picture — a strong consistency check.

## Step 3 — Compute $\Pi_2(x) \leftrightarrow \det\big(X(L),x,X(K_2),X(H)\big)$ (coefficient $c$ of $X(K_1)$)

Same style of row reduction (details omitted for brevity — same mechanics as Step 2, with $x$ now in the second slot):

$$c = \frac{(K_2-K_1)(s-L/R) - \big[(K_2-L)v_1-(K_1-L)v_2\big]}{(K_1-L)(K_2-K_1)}$$

**$c\ge0$**:
$$(K_2-L)v_1-(K_1-L)v_2 \;\le\; (K_2-K_1)\left(s-\frac{L}{R}\right)$$

A joint bound linking $v_1$ and $v_2$ together — this one doesn't reduce to either option's standalone bound; it only appears once both strikes are in the picture.

## Step 4 — Compute $\Pi_3(x) \leftrightarrow \det\big(X(L),X(K_1),x,X(H)\big)$ (coefficient $d$ of $X(K_2)$)

$$\begin{vmatrix} R&L&0&0\\ R&K_1&0&0\\ 1&s&v_1&v_2\\ R&H&H-K_1&H-K_2 \end{vmatrix}$$

Rows 1 and 2 agree in columns 3–4 (both zero) — subtracting row 1 from row 2 leaves a row of the form $(0,K_1-L,0,0)$, which simplifies the expansion considerably. Carrying this through:

$$\det\big(X(L),X(K_1),x,X(H)\big) = (K_1-L)\big[v_1(H-K_2)-v_2(H-K_1)\big]$$

Dividing by $\Pi$:

$$d = \frac{v_1(H-K_2)-v_2(H-K_1)}{(K_2-K_1)(H-K_2)}$$

**$d\ge0$**:
$$\boxed{v_1 \ge v_2\cdot\frac{H-K_1}{H-K_2}}$$

This is the cross-strike relationship, falling straight out of the coefficient calculation with no separate argument required. Since $H-K_1>H-K_2>0$, the ratio $(H-K_1)/(H-K_2)>1$, so this condition is strictly stronger than the intuitive "$v_1\ge v_2$" (the lower-strike call must be worth more) — it pins down exactly how much more.

## Step 5 — Compute $\Pi_4(x) \leftrightarrow \det\big(X(L),X(K_1),X(K_2),x\big)$ (coefficient $e$ of $X(H)$)

Rows 1 and 2 again agree in columns 3–4. Same simplification as Step 4:

$$\det\big(X(L),X(K_1),X(K_2),x\big) = (K_1-L)(K_2-K_1)\,v_2$$

Dividing by $\Pi$:

$$e = \frac{v_2}{H-K_2} \ge 0 \;\Longrightarrow\; v_2\ge0$$

The trivial bound — a call's payoff is never negative.

## Combined result

$$v_2\ge0 \qquad v_1 \ge v_2\cdot\frac{H-K_1}{H-K_2} \qquad v_1\ge s-\frac{K_1}{R} \qquad (K_2-L)v_1-(K_1-L)v_2\le(K_2-K_1)\left(s-\frac{L}{R}\right)$$

Together these four inequalities carve out the no-arbitrage region for $(v_1,v_2)$ jointly, given $s$.

## Worked numeric example

Take $L=80,\ K_1=100,\ K_2=110,\ H=130,\ R=1,\ s=105$.

$$v_2\ge0$$
$$v_1 \ge v_2\cdot\frac{130-100}{130-110} = 1.5\,v_2$$
$$v_1 \ge 105-100 = 5$$
$$30\,v_1-20\,v_2 \le 10\times(105-80)=250$$

**Consistency check**: try $(v_1,v_2)=(10,0)$. The first three bounds pass ($0\ge0$, $10\ge0$, $10\ge5$), but the fourth gives $30(10)-20(0)=300 > 250$ — **violated**. Solving the original linear system directly for $(b,c,d,e)$ at these values confirms $c=-0.25<0$, matching the violated inequality exactly.

Try instead $(v_1,v_2)=(8,3)$: $v_2=3\ge0$ ✓; $v_1=8\ge1.5(3)=4.5$ ✓; $v_1=8\ge5$ ✓; $30(8)-20(3)=240\le250$ ✓ — all four pass, so this pair is arbitrage-free.

## Why this matters for extending further

The same recipe — one more reference outcome, one more direction vector, one more determinant to evaluate — extends to any number of strikes $K_1<\cdots<K_n$. The row-reduction shortcuts seen in Steps 4–5 (reference points sharing zero entries in the option columns collapse the determinant quickly) suggest the coefficient calculations stay tractable even as $n$ grows, though the joint bound analogous to Step 3's condition (linking multiple option prices at once) will involve more cross-terms with each additional strike.

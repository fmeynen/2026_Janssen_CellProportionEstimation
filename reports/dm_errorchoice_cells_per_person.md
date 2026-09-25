# Why success rate barely increases with cells per person

Dirichlet-multinomial errorchoice simulation (`scripts/simulations/simulation_dm_errorchoice.R`)

## Observation

In the *success rate vs cells per person* plots, the lines are practically flat for both AE and ARE. Settings:

- concentration $c = 50$
- $K = 10$ cell types
- alpha $\in \{2, 3, 4, 5\}$
- $n_{\text{people}} \in \{1, 2, 3, 5, 10\}$
- $n_{\text{per person}}$ from $10^4$ to $10^8$
- $\tau_{AE} = 0.02$, $\tau_{ARE} = 0.5$

Going from $10^4$ to $10^8$ cells per person hardly changes the success rate. The success rate does change clearly with the number of people.

## Success rule

For each replicate and cell type $j$, the estimated proportions are averaged over the $N$ people in the replicate and compared with the population-level true proportion $p_j$:

$$\bar p_j = \frac{1}{N}\sum_{i=1}^{N} \hat p_{ij}, \qquad AE_j = |\bar p_j - p_j|, \qquad ARE_j = \frac{|\bar p_j - p_j|}{p_j}.$$

A replicate succeeds if $\max_j AE_j \le \tau_{AE}$ (and separately for ARE).

## Two sources of error

The data are generated in two steps:

1. Person $i$ has a true composition $\theta_i \sim \text{Dirichlet}(c \cdot p)$.
2. From that person, $n$ cells are counted: $\text{counts}_i \sim \text{Multinomial}(n, \theta_i)$, so $\hat p_{ij} = \text{count}_{ij} / n$.

The error of $\bar p_j$ therefore has two parts:

- **Differences between people (biological variation).** Each $\theta_{ij}$ scatters around $p_j$ with variance $\dfrac{p_j(1-p_j)}{c+1}$.
- **Counting (sampling variation).** Each $\hat p_{ij}$ scatters around $\theta_{ij}$. On average this adds a variance of $\dfrac{p_j(1-p_j)\,c}{(c+1)\,n}$.

Averaged over $N$ independent people:

$$\operatorname{Var}(\bar p_j - p_j) = p_j(1-p_j)\left[\frac{1}{N(c+1)} + \frac{c}{(c+1)\,N n}\right].$$

- The first term depends only on the **number of people** $N$ and the concentration $c$.
- The second term depends only on the **total number of cells counted**, $N \cdot n$.

## How large each part is at $c = 50$

Per person (the factor $p_j(1-p_j)/N$ is left out):

| $n$ per person | People term $1/(c+1)$ | Counting term $c/((c+1)n)$ | Counting share of variance |
|---|---|---|---|
| 100 | 0.0196 | 0.0098 | 33% |
| 1,000 | 0.0196 | 0.00098 | 5% |
| 10,000 | 0.0196 | 0.000098 | 0.5% |
| $10^8$ | 0.0196 | $\approx 0$ | $\approx 0\%$ |

The two terms are equal at $n = c = 50$.

- At $10^4$ cells per person, counting is already only 0.5% of the variance.
- Increasing $n$ from $10^4$ to $10^8$ (or even to infinity) lowers the standard deviation by only about 0.25%, which is invisible in the plots.
- Even with infinitely many cells per person, $\bar p_j$ converges to the average composition of *those particular* people, not to the population proportion $p_j$. The people term remains.

## Why the success rates are where they are

**AE.** Take the most abundant cell type, $p \approx 0.4$ (alpha = 5). With $N = 10$ people:

$$\text{SD}(\bar p_j - p_j) \approx \sqrt{\frac{0.4 \cdot 0.6}{51 \cdot 10}} \approx 0.022.$$

That is already larger than $\tau_{AE} = 0.02$, and the rule takes the worst of the 10 cell types. Simulated success rates of roughly 25–37% for 10 people are therefore expected.

**ARE.** Dividing by $p_j$ makes the rare cell types dominate. For $p_j = 0.01$ and $N = 10$:

$$\text{SD}(ARE_j) \approx \sqrt{\frac{1-p_j}{p_j\,(c+1)\,N}} = \sqrt{\frac{0.99}{0.01 \cdot 51 \cdot 10}} \approx 0.44,$$

which is close to $\tau_{ARE} = 0.5$. Larger alpha makes the smallest $p_j$ smaller, so ARE success drops to 0 for alpha = 4 and 5.

**The jump at ARE = 1.** For rare cell types $c \cdot p_j$ is very small. Many people then have essentially none of that cell type in their true composition. If none of the people in a replicate have it, $\bar p_j = 0$ and $ARE_j = |0 - p_j| / p_j = 1$ exactly. This produces the jump at $\tau = 1$ in the ARE *success rate vs threshold* plots; pooled over all scenarios at 100,000 cells per person, the median ARE is exactly 1.

## Conclusions

1. **At concentration 50, more people is the only effective lever once a person has more than a few hundred cells.**
   - Around $n \approx 100$, counting error is still about a third of the variance, so going to about 1,000 cells per person helps somewhat.
   - Beyond about 1,000 cells per person, counting is 5% of the variance or less, and more cells per person barely helps.
2. **For a fixed total number of cells, spreading them over more people is always better.** The counting term depends only on the total $N \cdot n$. The people term decreases with $N$. So 100,000 cells from 50 people give a much more accurate estimate than 100,000 cells from 5 people. In practice the limit is the cost of including extra people, not the cost of counting cells.
3. **This depends on the estimation target.** The conclusions above hold when the goal is the *population* composition $p$. If the goal is each individual's own composition $\theta_i$, the people term disappears and cells per person is the lever.
4. **It depends strongly on the concentration.** $c = 50$ means large differences between people. With a higher concentration (more similar people) the people term shrinks, and cells per person keeps mattering up to roughly $n \approx c$. A realistic value of $c$, estimated from data, is therefore essential for interpreting these results.

## Possible next step

Show success rate against the number of people (e.g. 1 to 200) at a moderate fixed depth such as 1,000 cells per person. That shows how many people are needed to reach 95% success for AE and ARE at each alpha.

## MCMC for fitting the wp

## auto corr (gal-gal)
$$
w_p = w_{p(1h)} + w_{p(2h)} = w_{p(1h, cs)} + w_{p(1h, sc)} + w_{p(1h, ss)} + w_{p(2h,cc)} + w_{p(2h,ss)} + w_{p(2h,sc)} + w_{p(2h,cs)}
$$
### 1 halo
$$
cs = sc = \sum_i \frac{N_h(M_i)}{N_g^2} N_c N_s w_{cs}(M_i)\\
ss = \sum_i \frac{N_h(M_i)}{N_g^2}N_s(N_s-1) w_{ss}(M_i)
$$
Note that $w_{cs}, w_{ss}$ is a projected probability distribution function, which means it is normalized.
*I have to note that in code, the $w_{cs}$ is called $w_p$, but it is totally **different** thing with $w_p$.*

### 2 halo
$$
cc = \sum_{i,j} \frac{N_h(M_i) N_h (M_j)}{N_g^2} N_c(M_i)N_c(M_j) w_{p(2h, cc)}(M_i, M_j)\\
cs = \sum_{i,j} \frac{N_h(M_i) N_h (M_j)}{N_g^2} N_c(M_i) N_s(M_j) w_{p(2h, cs)}(M_i, M_j)\\
sc = \sum_{i,j} \frac{N_h(M_i) N_h (M_j)}{N_g^2} N_s(M_i) N_c(M_j) w_{p(2h, sc)}(M_i, M_j)\\
ss = \sum_{i,j} \frac{N_h(M_i) N_h (M_j)}{N_g^2} N_s(M_i)N_s(M_j) w_{p(2h, ss)}(M_i, M_j)
$$

Here, $N_c, N_s$ is the average central galaxy number and satellite galaxy number per halo, so it is a function of $M$. $N_h(M_i)$ is halo mass function, which is **number density** of halos in mass bin $\log M_h \pm \Delta \log M$.
$$
N_c(M) = \frac 12 \left( 1+erf\left(\frac{\log M - \log M_{min}}{\sigma_{\log M}}\right)\right)\\
N_s(M) = \begin{cases}\left(\frac{M-M_0}{M_1'}\right)^\alpha N_c(M),& M>M_0\\ 0,&M\le M_0\end{cases}\\
$$
Therefore, the number density of all galaxies and satellite galaxies are 
$$
N_g = \sum_i N_h(M_i) (N_c(M_i) + N_s(M_i))\\
N_g^{sat} = \sum_i N_h(M_i) N_s(M_i)
$$

## crosss corr (gal-matter)
$$
w_p = w_{p(c, m)} + w_{p(s, m)}
$$
where
$$
w_{p(c,m)} = \sum_i \frac{N_h(M_i)}{N_g} N_c w_{p(c, m)}(M_i)\\
w_{p(s,m)} = \sum_i \frac{N_h(M_i)}{N_g} N_s w_{p(s, m)}(M_i)
$$

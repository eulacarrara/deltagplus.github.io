---  
title:  
---  
  
<br>  
  
### Solving mixed model equations with unknown parent groups (UPG)  
  
##### 1) Direct solving (R)  
##### 2) Iterative solving (BLUPF90+)  
  
---  
  
**Author:** Eula Carrara (eulacarrara@gmail.com)  
**Created on:** 22-Feb-2025  
**Modified on:** 15-Mar-2025  
**Theory:** lecture notes Dr. Daniela Lourenco, UGA, 2025  
**Numerical example taken from** [tutorial_blupf90.pdf](https://nce.ads.uga.edu/wiki/lib/exe/fetch.php?media=tutorial_blupf90.pdf)  
**Required R packages:** igraph, MASS, optiSel, tidyverse  
```r  
install.packages(c("igraph", "MASS", "optiSel", "tidyverse"))  
```  
```{r, echo=FALSE, results='hide', message=FALSE, warning=FALSE}  
options(width = 200)  
```  
---  
  
Considering an individual **a**, progeny of **s** and **d**,  
  
\begin{array}{c}  
s_i \quad \quad d_i \\  
\backslash \quad \, / \\  
\, a_i  
\end{array}  
  
its expected breeding value will be:  
  
- **if both parents known:** $a_i = \frac{a_{s_i} + a_{d_i}}{2}$  
  
- **if only sire known:** $a_i = \frac{a_{s_i} + 0}{2}$  
  
- **if only dam known:** $a_i = \frac{0 + a_{d_i}}{2}$  
  
- **if both parents unknown:** $a_i = \frac{0 + 0}{2}$  
  
That is, the breeding value of an unknown sire/dam is zero. This does not correspond to reality.  
One way to fill these gaps in the pedigree is to build Unknown Parent Groups (UPG; or phantom parents, or genetic group):  
  
- **only known sire:** $a_i = \frac{a_{s_i} + \boldsymbol{\mathit{UPG_i}}}{2}$  
  
- **only known dam:** $a_i = \frac{\boldsymbol{\mathit{UPG_i}} + a_{d_i}}{2}$  
  
- **both unknown parents:** $a_i = \frac{\boldsymbol{\mathit{UPG_i}} + \boldsymbol{\mathit{UPG_i}}}{2}$  
  
<br>  
  
**How ​​to create UPGs?**  
  
Using information about breed, year of birth, sex, country, etc...  
UPGs are not animals! They are most often treated as a fixed effect in the model.  
  
<br>  
    
##### Animal mixed model  
The mixed animal model can be described as:  
$y = Xb + Za + e$  
  
in which:  
- `y` is the vector of phenotypes,  
- `b` is the vector of fixed effects,  
- `a` is the vector of direct additive genetic effects (breeding values),  
- `X` and `Z` relate `b` and `a` to `y`,  
- `e` is the vector of residuals.  
  
<br>  
  
##### Animal mixed model, but with **UPG**:  
$y = Xb + \mathbf{\mathit{ZQg}} + Za + e$  
  
in which:  
- `Q` is the matrix that relates animals to UPGs,  
- `g` is the effect of UPGs.  
  
So we need to build `Q`. How?  
The `Q` matrix is ​​constructed with the <u>contribution fractions</u> of each UPG in the expected breeding values ​​of the individuals related to the UPG.  
  
<br>  
  
##### Example of constructing the Q matrix:  
  
**Pedigree:**  
```{r, echo=FALSE, results='hide', message=FALSE, warning=FALSE}  
library(igraph)  
library(tidyverse)  
pedigree_ex <- data.frame(  
ID = c('a1', 'a2', 'a3', 'a4', 'a5', 'a6'),  
Sire = c('g2', 'g3', 'g3', 'g1', 'a3', 'g3'),  
Dam = c('a4', 'a1', 'a1', 'g1', 'a6', 'g3')  
)  
  
edges <- pedigree_ex %>%  
pivot_longer(cols = c(Sire, Dam), values_drop_na = TRUE) %>%  
select(from = value, to = ID)  
  
gr <- graph_from_data_frame(edges, directed = TRUE)  
  
generation_levels <- c("g1" = 1,   
 "a4" = 2,   
 "g2" = 3,   
 "a1" = 4,   
 "g3" = 5,   
 "a2" = 6, "a3" = 6, "a6" = 6,   
 "a5" = 7)  

V(gr)$y <- generation_levels[V(gr)$name]  
  
layout <- layout_as_tree(gr)  
  
unique_generations <- sort(unique(V(gr)$y))  
gen_map <- setNames(seq(0, -length(unique_generations) + 1, length.out = length(unique_generations)), unique_generations)  
layout[,2] <- sapply(V(gr)$y, function(y) gen_map[y])  
  
height_scale <- 2  
layout[, 2] <- layout[, 2] * height_scale  
  
width_scale <- 2   
for (gen in unique_generations) {  
nodes_in_gen <- which(V(gr)$y == gen)  
if (length(nodes_in_gen) > 1) {  
layout[nodes_in_gen, 1] <- seq(-width_scale, width_scale, length.out = length(nodes_in_gen))  
}  
}  
```  
```{r, echo=FALSE}  
print(pedigree_ex)  
```  
```{r, echo=FALSE, results='hide', message=FALSE, warning=FALSE}  
par(mar = c(1, 1, 1, 1))   
plot(gr, layout = layout,  
 vertex.size = 30,   
 vertex.color = "lightblue",  
 vertex.label.color = "black",   
 vertex.label.font = 2,   
 vertex.label.cex = 1,  
 edge.width = 1,   
 edge.arrow.size = 0.5,# Uniform arrow size  
 edge.color = "black")  
```  
  
First, we need to put all the expected breeding values as a function of the UPGs:  
  
**Defining the expected breeding values depending on the UPG**  
\begin{aligned}  
E(a_4) &= \frac{1}{2} \cdot g_1 + \frac{1}{2} \cdot g_1 = g_1 \\  
\\[0.1cm]  
E(a_1) &= \frac{1}{2} \cdot g_2 + \frac{1}{2} \cdot E(a_4) = \frac{1}{2} g_1 + \frac{1}{2} g_2 \\  
\\[0.1cm]  
E(a_3) &= \frac{1}{2} \cdot g_3 + \frac{1}{2} \cdot E(a_1) = \frac{1}{4} g_1 + \frac{1}{4} g_2 + \frac{1}{2} g_3 \\  
\\[0.1cm]  
E(a_6) &= \frac{1}{2} \cdot g_3 + \frac{1}{2} \cdot g_3 = g_3 \\  
\\[0.1cm]  
E(a_5) &= \frac{1}{2} \cdot E(a_3) + \frac{1}{2} \cdot E(a_6) = \frac{1}{8} g_1 + \frac{1}{8} g_2 + \frac{3}{4} g_3 \\  
\\[0.1cm]  
E(a_2) &= \frac{1}{2} \cdot g_3 + \frac{1}{2} \cdot E(a_1) = \frac{1}{4} g_1 + \frac{1}{4} g_2 + \frac{1}{2} g_3  
\end{aligned}  
  
**Rearranging from a1 to a6** (just for better visualization)  
\begin{aligned}  
E(a_1) &= \frac{1}{2} \cdot g_1 + \frac{1}{2} \cdot g_2 + 0 \cdot g_3 \\  
\\[0.1cm]  
E(a_2) &= \frac{1}{4} \cdot g_1 + \frac{1}{4} \cdot g_2 + \frac{1}{2} \cdot g_3 \\  
\\[0.1cm]  
E(a_3) &= \frac{1}{4} \cdot g_1 + \frac{1}{4} \cdot g_2 + \frac{1}{2} \cdot g_3 \\  
\\[0.1cm]  
E(a_4) &= 1 \cdot g_1 + 0 \cdot g_2 + 0 \cdot g_3 \\  
\\[0.1cm]  
E(a_5) &= \frac{1}{8} \cdot g_1 + \frac{1}{8} \cdot g_2 + \frac{3}{4} \cdot g_3 \\  
\\[0.1cm]  
E(a_6) &= 0 \cdot g_1 + 0 \cdot g_2 + 1 \cdot g_3  
\end{aligned}  
  
**Q matrix with the fractions of each UPG (cols) for each animal (rows)**  
\[  
\begin{bmatrix}  
1/2 & 1/2 & 0 \\  
1/4 & 1/4 & 1/2 \\  
1/4 & 1/4 & 1/2 \\  
1 & 0 & 0 \\  
1/8 & 1/8 & 3/4 \\  
0 & 0 & 1  
\end{bmatrix}  
\]  
  
This is the `Q` matrix that will be incorporated into the mixed model equations.  
  
<br>  
  
#### Example of solving mixed model equations considering UPG  
##### Solving via R program, step by step  
  
**Data**  
```{r, echo=TRUE}  
data1 <- data.frame(  
ID = c('ID006', 'ID009', 'ID012', 'ID007', 'ID010', 'ID013', 'ID008', 'ID011', 'ID014', 'ID015'),  
A = c('A', 'A', 'A', 'B', 'B', 'B', 'C', 'C', 'C', 'C'),  
S = c(1, 2, 1, 2, 1, 2, 1, 2, 1, 2),  
cov = c(1.0, 1.0, 2.0, 2.0, 1.0, 2.0, 2.0, 1.0, 1.0, 2.0),  
obs = c(3.0, 2.0, 4.0, 6.0, 3.0, 6.0, 6.0, 6.0, 8.0, 4.0))  
print(data1)  
```  
  
<br>  
  
**Pedigree without UPG**  
```{r, echo=TRUE}  
pedigree <- data.frame(  
ID = c('ID001', 'ID002', 'ID003', 'ID004', 'ID005', 'ID006', 'ID007', 'ID008', 'ID009', 'ID010', 'ID011', 'ID012', 'ID013', 'ID014', 'ID015'),  
Sire = c(NA, NA, NA, NA, NA, NA, 'ID002', 'ID001', 'ID002', 'ID007', 'ID007', 'ID011', 'ID011', 'ID009', 'ID011'),  
Dam = c(NA, NA, NA, NA, NA, NA, 'ID005', 'ID004', 'ID003', 'ID006', 'ID004', 'ID008', 'ID010', 'ID013', 'ID010')  
)  
print(pedigree)  
```  
  
<br>  
  
**Pedigree with UPG**  
```{r, echo=TRUE}  
pedigree_upg <- data.frame(  
ID = c('ID001', 'ID002', 'ID003', 'ID004', 'ID005', 'ID006', 'ID007', 'ID008', 'ID009', 'ID010', 'ID011', 'ID012', 'ID013', 'ID014', 'ID015'),  
Sire = c('g1', 'g2', 'g1', 'g2', 'g2', 'g1', 'ID002', 'ID001', 'ID002', 'ID007', 'ID007', 'ID011', 'ID011', 'ID009', 'ID011'),  
Dam = c('g4', 'g3', 'g3', 'g3', 'g4', 'g3', 'ID005', 'ID004', 'ID003', 'ID006', 'ID004', 'ID008', 'ID010', 'ID013', 'ID010')  
)  
print(pedigree_upg)  
  
# Ordering pedigree without UPG by ID  
pedigree <- pedigree[order(pedigree$ID), ]  
  
# Ordering data by ID  
data1 <- data1[order(data1$ID), ]  
```  
  
<br>  
  
##### Model  
$y = Xb + ZQg + Za + e$  
  
<br>  
  
##### Mixed Model Equations  
\[  
\begin{bmatrix}  
X'X & 0 & X'Z \\  
0 & Q'A^{-1}Q\alpha & -Q'A^{-1}\alpha \\  
Z'X & -A^{-1}Q\alpha & Z'Z + A^{-1}\alpha  
\end{bmatrix}  
*  
\begin{bmatrix}  
b\\  
0\\  
Qg + a  
\end{bmatrix}  
=  
\begin{bmatrix}  
X'y \\  
0 \\  
Z'y  
\end{bmatrix}  
\]  
  
<br>  
  
##### Elements of mixed model equations  
**y vector**  
```{r, echo=TRUE}  
y <- as.matrix(data1$obs)  
print(y)  
```  
  
<br>  
  
**X matrix**  
```{r, echo=TRUE}  
library(dplyr)  
X <- data1 %>%  
mutate(A_A = ifelse(A == "A", 1, 0),  
 A_B = ifelse(A == "B", 1, 0),  
 A_C = ifelse(A == "C", 1, 0),  
 S_1 = ifelse(S == 1, 1, 0),  
 S_2 = ifelse(S == 2, 1, 0)) %>%  
select(A_A, A_B, A_C, S_1, S_2, cov)  
  
X <- as.matrix(X)  
print(X)  
```  
  
<br>  
  
**Z matrix**  
(expanding Z for all animals to sum with the inverse of numerator relationship matrix, \( A^{-1} \))  
```{r, echo=TRUE}  
animals <- pedigree$ID  
N <- nrow(data1)  
Np <- length(animals)  
  
Z <- matrix(0, nrow = N, ncol = Np)  
for (i in 1:N) {  
animal_index <- match(data1$ID[i], animals)  
Z[i, animal_index] <- 1  
}  
print(Z)  
```  
  
<br>  
  
**Constructing matrices A and inverse A**  
(note that the inverse A is constructed using the pedigree without UPG)  
```{r, echo=TRUE, results='hide', message=FALSE, warning=FALSE}  
library(optiSel)  
# A matrix  
A <- makeA(pedigree)  
```  
```{r, echo=TRUE}  
# Show first 5 col/row  
print(A[1:5, 1:5], digits = 4)  
```  
```{r, echo=TRUE,results='hide', message=FALSE, warning=FALSE}  
# Inverse matrix A  
Ainv <- solve(A)  
```  
```{r, echo=TRUE}  
# Show first 5 col/row  
print(Ainv[1:5, 1:5], digits = 4)  
```  
  
<br>  
  
**Q matrix**  
- Function to compute Q matrix:  
```{r, echo=TRUE}  
n_UPG = 4  
  
compute_expectation <- function(pedigree, N) {  
expectations <- list()  
  
base_groups <- unique(c(pedigree$Sire, pedigree$Dam))  
base_groups <- base_groups[grepl("^g[0-9]+$", base_groups)]  
  
for (group in base_groups) {  
expectations[[group]] <- setNames(as.list(rep(0, N)), paste0("g", 1:N))  
expectations[[group]][[group]] <- 1  
}  
  
get_expectation <- function(id) {  
if (!is.null(expectations[[id]])) {  
return(expectations[[id]])  
}  
  
row <- pedigree[pedigree$ID == id, ]  
if (nrow(row) == 0) {  
return(setNames(as.list(rep(0, N)), paste0("g", 1:N)))  
}  
  
sire_exp <- get_expectation(row$Sire)  
dam_exp <- get_expectation(row$Dam)  
  
expectations[[id]] <- setNames(lapply(1:N, function(i) {  
(1/2) * sire_exp[[paste0("g", i)]] + (1/2) * dam_exp[[paste0("g", i)]]  
}), paste0("g", 1:N))  
  
return(expectations[[id]])  
}  
  
animals <- pedigree$ID[grepl("^ID[0-9]+$", pedigree$ID)]  
for (animal in animals) {  
expectations[[animal]] <- get_expectation(animal)  
}  
  
return(expectations)  
}  
```  
  
<br>  
  
**Applying the function in pedigree with UPG**  
(format: ID, Sire, Dam, Sire/Dam with UPG codes)  
```{r, echo=TRUE}  
expectations <- compute_expectation(pedigree_upg, n_UPG)  
  
coefficients <- data.frame(  
ID = names(expectations),  
g1 = sapply(expectations, function(x) x$g1),  
g2 = sapply(expectations, function(x) x$g2),  
g3 = sapply(expectations, function(x) x$g3),  
g4 = sapply(expectations, function(x) x$g4))  
```  
  
**Final Q matrix:**  
```{r, echo=TRUE}  
Q <- as.matrix(coefficients[-(1:4), -1])  
print(Q)  
```  
  
<br>  
  
Let's assume that the variance components are known:  
  
Direct additive variance: $0.5$  
```{r, echo=TRUE}  
vara <- 0.5  
```  
  
Residual variance: $2.0$  
```{r, echo=TRUE}  
vare <- 2.0  
```  
  
Therefore, alpha is equal to:  
```{r, echo=TRUE}  
alpha <- vare / vara  
print(alpha)  
```  
  
<br>  
  
We have the individual elements of the mixed model equations, namely \( X \), \( Z \), \( Q \), \( A^{-1} \), \( \alpha \), and \( y \).  
Now, we need to calculate each "block":  
  
\[  
\begin{bmatrix}  
X'X & 0 & X'Z \\  
0 & Q'A^{-1}Q\alpha & -Q'A^{-1}\alpha \\  
Z'X & -A^{-1}Q\alpha & Z'Z + A^{-1}\alpha  
\end{bmatrix}  
*  
\begin{bmatrix}  
b \\  
0 \\  
Qg + a  
\end{bmatrix}  
=  
\begin{bmatrix}  
X'y \\  
0 \\  
Z'y  
\end{bmatrix}  
\]  
  
\[  
\begin{bmatrix}  
\text{block11} & \text{block12} & \text{block13} \\  
\text{block21} & \text{block22} & \text{block23} \\  
\text{block31} & \text{block32} & \text{block33}  
\end{bmatrix}  
*  
\begin{bmatrix}  
b \\  
0 \\  
Qg + a  
\end{bmatrix}  
=  
\begin{bmatrix}  
\text{block1} \\  
\text{block2} \\  
\text{block3}  
\end{bmatrix}  
\]  
  
<br>  
  
Remember: $LHS * sol = RHS$  
  
<br>  
  
**LHS - left-hand side of mixed model equations**  
```{r, echo=TRUE, results='hide', message=FALSE, warning=FALSE}  
block_11 <- t(X) %*% X  
block_12 <- matrix(0, nrow = nrow(block_11), ncol = ncol(Q))  
block_13 <- t(X) %*% Z  
  
block_21 <- matrix(0, nrow = ncol(Q), ncol = nrow(block_11))  
block_22 <- t(Q) %*% Ainv %*% Q * alpha  
block_23 <- -t(Q) %*% Ainv * alpha  
  
block_31 <- t(Z) %*% X  
block_32 <- -Ainv %*% Q * alpha  
block_33 <- t(Z) %*% Z + Ainv * alpha  
  
LHS <- rbind(  
cbind(block_11, block_12, block_13),  
cbind(block_21, block_22, block_23),  
cbind(block_31, block_32, block_33))  
dimnames(LHS) <- NULL  
```  
```{r, echo=TRUE}  
# Show first 5 col/row 
print(LHS[1:5, 1:5], digits = 4)  
```  
  
**RHS - right-hand side of mixed model equations**  
```{r, echo=TRUE, results='hide', message=FALSE, warning=FALSE}  
block_1 <- t(X) %*% y  
block_2 <- matrix(0, nrow = ncol(Q), ncol = ncol(y))  
block_3 <- t(Z) %*% y  
  
RHS <- rbind(block_1, block_2, block_3)  
dimnames(RHS) <- NULL  
```  
```{r, echo=TRUE}  
# Show first 5 col/row  
print(RHS[1:5, 1], digits = 4)  
```  
  
<br>  
  
We need to <u>invert the LHS</u>, because:  
$\text{sol} = \text{LHS}^{-1} \times \text{RHS}$  
```r  
solve(LHS)  
```  
<span style="color:red;">`! system is computationally singular: reciprocal condition number = 1.58711e-18`</span>  
  
<br>  
  
Notice that the LHS **is not a positive definite matrix** (`system is computationally singular`), so it does not have an inverse.  
We need to calculate the generalized inverse...  
  
We will use the `ginv()` function from the `MASS` package, which calculates the generalized Moore-Penrose inverse (one of the most common ginv)  
```{r, echo=TRUE, results='hide', message=FALSE, warning=FALSE}  
library(MASS)  
LHSinv <- ginv(LHS)  
```  
```{r, echo=TRUE}  
# Show first 5 col/row  
print(LHSinv[1:5, 1:5], digits = 4)  
```  
  
We multiply the inverse of the LHS by the RHS and we will have our vector of solutions.  
  
**Final solutions (R):**  
```{r, echo=TRUE}  
solutions_R <- LHSinv %*% RHS  
print(solutions_R)  
```  
  
The solutions order is:  
$$  
\begin{bmatrix}  
b \\  
g \\  
Qg + a  
\end{bmatrix}  
$$  
in which:  
- `b` is the vector of solutions for the fixed effects,  
- `g` is the vector of solutions for the UPG effects (also fixed here),  
- `Qg + a` is the vector of breeding values.  
Note that the **breeding value of the animal** will be the **effect of the UPG (Qg) + breeding value (a)**.  
  
that is:  
  
<br>  
  
**Vector of fixed effects (b)**  
```{r, echo=TRUE}  
beta_hat <- solutions_R[1:ncol(X)]  
b <- as.matrix(beta_hat)  
print(b)  
```  
  
<br>  
  
**Vector of UPG effects (g)**  
```{r, echo=TRUE}  
upg_hat <- solutions_R[(ncol(X) + 1):(ncol(X) + ncol(Q))]  
g <- as.matrix(upg_hat)  
print(g)  
```  
  
<br>  
  
**Vector of direct additive genetic effects - breeding value (Qg+a)**  
```{r, echo=TRUE}  
animal_hat <- solutions_R[(ncol(X) + ncol(Q) + 1):(ncol(X) + ncol(Q) + ncol(Z))]  
Qg_a <- as.matrix(animal_hat)  
print(Qg_a)  
```  
And those were our final solutions!  
  
<br>  
  
##### y predicted  
Now, let's calculate the vector `y` (predicted y) from the solutions we found.  
To do this, we will use the first and third rows of the mixed model equations, since the second row of the RHS is equal to zero.  
  
**From first row:**  
$$  
X'y = X'Xb + X'Z(Qg+a)  
$$  
```{r, echo=TRUE}  
rhs1 <- t(X) %*% X %*% b + t(X) %*% Z %*% Qg_a  
```  
  
**From third row:**  
$$  
Z'y = Z'Xb - A^{-1}Q\alpha g + (Z'Z + A^{-1} \alpha)(Qg + a)  
$$  
```{r, echo=TRUE}  
rhs3 <- t(Z) %*% X %*% b - Ainv %*% Q %*% (alpha * g) + (t(Z) %*% Z + Ainv * alpha) %*% Qg_a  
```  
  
**Let's create a system with these two equations, to solve for `y`**  
\[  
\left\{  
\begin{aligned}  
X'y &= X'bX + X'Z(Qg+a) \\  
Z'y &= Z'Xb - A^{-1}Q\alpha g + (Z'Z + A^{-1} \alpha)(Qg + a)  
\end{aligned}  
\right.  
\]  
  
Whatever the animal model, the part that multiplies `y` will always be an `mxn` matrix. Let's call it `M`.  
Similarly, the whole part after the equal sign will be an `nx1` vector. Let's call it `c`.  
  
So, we have that $My=c$.  
  
<br>  
  
**M matrix**  
```{r, echo=TRUE}  
M <- rbind(t(X), t(Z))  
```  
  
<br>  
  
**c vector**  
```{r, echo=TRUE}  
c <- rbind(rhs1, rhs3)  
```  
  
<br>  
  
Solving for y...  
If we try to invert the matrix `M`, we see that it is also not positive definite, so we need a generalized inverse: `ginv()`.  
  
```{r, echo=TRUE}  
y_hat <- ginv(M) %*% c  
  
final <- data.frame(y_predicted = y_hat, y_original = y)  
print(final)  
```  
The `predicted y` is identical to the `original y`.  
  
<br>  
  
This is easy for 15 animals. But if we have 1 million of them?  
We need a **more efficient** program.  
  
<br>  
  
##### Solving via program [`BLUPF90+`](https://nce.ads.uga.edu/wiki/doku.php?id=readme.blupf90plus)  
Let's solve this same example using the BLUPF90 family of programs  
*(Misztal I., Tsuruta S., Lourenco D.A.L., Aguilar I., Legarra A., and Vitezica Z. 2014. Manual for BLUPF90 family of programs.)*.  
  
<br>  
  
**Data**  
`ID006 A 1 1.0 3.0`  
`ID009 A 2 1.0 2.0`  
`ID012 A 1 2.0 4.0`  
`ID007 B 2 2.0 6.0`  
`ID010 B 1 1.0 3.0`  
`ID013 B 2 2.0 6.0`  
`ID008 C 1 2.0 6.0`  
`ID011 C 2 1.0 6.0`  
`ID014 C 1 1.0 8.0`  
`ID015 C 2 2.0 4.0`  
  
**Pedigree with UPG**  
`ID001 0 0 -1 -4`  
`ID002 0 0 -2 -3`  
`ID003 0 0 -1 -3`  
`ID004 0 0 -2 -3`  
`ID005 0 0 -2 -4`  
`ID006 0 0 -1 -3`  
`ID007 ID002 ID005 ID002 ID005`  
`ID008 ID001 ID004 ID001 ID004`  
`ID009 ID002 ID003 ID002 ID003`  
`ID010 ID007 ID006 ID007 ID006`  
`ID011 ID007 ID004 ID007 ID004`  
`ID012 ID011 ID008 ID011 ID008`  
`ID013 ID011 ID010 ID011 ID010`  
`ID014 ID009 ID013 ID009 ID013`  
`ID015 ID011 ID010 ID011 ID010`  
  
Here we encode the UPGs with <u>negative numbers</u> instead of letters. This is a requirement of the BLUPF90+ program!  
  
  
**Parameter file**  
`DATAFILE`  
`data1.txt`  
`TRAITS`  
`5`  
`FIELDS_PASSED TO OUTPUT`  
` `  
`WEIGHT(S)`  
` `  
`RESIDUAL_VARIANCE`  
`2.0`  
`EFFECT`  
`2 cross alpha`  
`EFFECT`  
`3 cross alpha`  
`EFFECT`  
`4 cov`  
`EFFECT`  
`1 cross alpha`  
`RANDOM`  
`animal`  
`FILE`  
`ped2.txt`  
`FILE_POS`  
`1 4 5 0 0   # id, sire, dam - with upg code`  
`UPG_TYPE`  
`in_ped`  
`INBREEDING`  
`no-inbreeding   # no inbreeding`  
`(CO)VARIANCES`  
`0.5`  
  
<br>  
  
We run the program [`RENUMF90`](https://nce.ads.uga.edu/wiki/doku.php?id=readme.RENUMF90) and then the program [`BLUPF90+`](https://nce.ads.uga.edu/wiki/doku.php?id=readme.blupf90plus).  
The solutions (`solutions` file) will be:  
```{r, echo=TRUE}  
sol_blup <- data.frame(  
trait = rep(1, 25),  
effect = c(rep(1, 3), rep(2, 2), 3, rep(4, 19)),  
level = c(1,2,3, 1,2, 1, 1:19),  
solution = c(  
-0.09929111, 1.82312339, 2.56875029, 2.43588849, 0.65341111,   
1.01898764, 0.31554213, -1.15906490, 0.54803039, -1.31841690,   
-0.53693713, 0.78918582, 1.56194108, 0.06661379, 0.68188768,   
0.87340751, -1.52607734, 2.63594054, 2.60934675, -5.25998297,  
-1.10486630, -4.38182541, 3.11359067, 2.12479542, -6.15194294  
))  
print(sol_blup)  
```  
  
<br>  
  
The solutions of the `BLUPF90+` program have the order:  
`[fixed effects]`, `[breeding values]`, `[UPG effects]`.  
Let's put them in the same order as the solutions we obtained in the R program, for comparison.  
You can find the original animal ID in the `renadd04.ped` file and match it with the `solutions` file to identify the order.  
  
```{r, echo=TRUE}  
sol_blup$order <- c("A", "B", "C", "S1", "S2", "cov",   
"ID015", "ID006", "ID007", "ID008", "ID010", "ID009", "ID011", "ID012",  
"ID013", "ID014", "ID005", "ID002", "ID004", "ID001", "ID003",   
"g1", "g2", "g3", "g4")  
  
R_blup_order <- c("A", "B", "C", "S1", "S2", "cov",   
 "g1", "g2", "g3", "g4",   
 "ID001", "ID002", "ID003", "ID004", "ID005", "ID006", "ID007", "ID008",  
 "ID009", "ID010", "ID011", "ID012", "ID013", "ID014", "ID015")  
```  
  
Reordering  
```{r, echo=TRUE}  
sol_blup$order <- factor(sol_blup$order, levels = R_blup_order)  
sol_blup2 <- sol_blup %>% arrange(order)  
```  
  
Final solutions (BLUPF90+)  
```{r, echo=TRUE}  
solutions_blup <- as.matrix(sol_blup2$solution)  
print(solutions_blup)  
```  
  
Split the effects  
```{r, echo=TRUE}  
# Fixed effects  
beta_hat <- solutions_blup[1:ncol(X)]; b <- as.matrix(beta_hat)  
print(beta_hat)  
# UPG effects  
upg_hat <- solutions_blup[(ncol(X) + 1):(ncol(X) + ncol(Q))]; g <- as.matrix(upg_hat)  
print(upg_hat)  
# Breeding values  
animal_hat <- solutions_blup[(ncol(X) + ncol(Q) + 1):(ncol(X) + ncol(Q) + ncol(Z))]; Qg_a <- as.matrix(animal_hat)  
print(animal_hat)  
```  
  
<br>  
  
The solutions are **different** from those obtained directly, since we used the **generalized inverse**.  
In other words, we will have different solutions, depending on which inverse we use.  
However, our `predicted y` must be the same in both predictions, since <u>estimable functions are invariant to the choice of the generalized inverse</u> (this is a topic for another post).  
  
Let's get the `predicted y` vector, just like we did for the R solutions.  
```{r, echo=TRUE}  
rhs1 <- t(X) %*% X %*% b + t(X) %*% Z %*% Qg_a  
rhs3 <- t(Z) %*% X %*% b - Ainv %*% Q %*% (alpha * g) + (t(Z) %*% Z + Ainv * alpha) %*% Qg_a  
```  
  
$My=c$  
  
M matrix  
```{r, echo=TRUE}  
M <- rbind(t(X), t(Z))  
```  
  
c vector  
```{r, echo=TRUE}  
c <- rbind(rhs1, rhs3)  
```  
  
Solving y...  
```{r, echo=TRUE}  
y_hat <- ginv(M) %*% c  
  
final_blup <- data.frame(y_predicted = y_hat, y_original = y)  
print(final_blup)  
```  
  
#### Let's compare the results via R and via blupf90+!  
```{r, echo=TRUE}  
final_solutions <- data.frame(  
 ID = data1$ID,  
 y_original = data1$obs,  
 y_predicted_R = final$y_predicted,  
 y_predicted_blupf90 = round(final_blup$y_predicted,3),  
 residual_blupf90 = round((final_blup$y_predicted - data1$obs),4))  
print(final_solutions)  
```  
  
Note that the predicted `y` obtained either directly (`R`) or iteratively (`BLUPF90+`) matches the `original y`.  
The solutions obtained via `BLUPF90+` are not exact, as they were generated by an iterative process. However, the residual is minimal.  
R is limited in many aspects, while BLUPF90+ is much more efficient.  
  
<br><br>  

#a  install.packages("expm")

library(expm) # 用於矩陣次方運算 %^%

G <- matrix(c(3, 0, 2, 0, 1, 1, 2, 1, 9), nrow=3, byrow=TRUE)
b <- c(2, 3, 1)

ans_i_val <- eigen(G)$values     # (i) 特徵值
ans_i_vec <- eigen(G)$vectors    # (i) 特徵向量
ans_ii <- det(G)                 # (ii) 行列式
ans_iii <- solve(G, b)           # (iii) 求解 Gx=b
# (iv) 計算 7*G^6 + 2*G^-5 + 3*G'

# 先計算 G^5 的反矩陣 (即 G^-5)
G_inv5 <- solve(G %^% 5)

# 重新組合方程式
ans_iv_package <- 7 * (G %^% 6) + 2 * G_inv5 + 3 * t(G)

print("第 (iv) 題的答案：")
print(ans_iv_package)
ans_v_colVar <- apply(G, 2, var)    # (v) 各行變異數
ans_v_rowMed <- apply(G, 1, median) # (v) 各列中位數


#b
x <- seq(0.1, 3, by=0.1)
fx <- exp(x) + 2*sin(x)
gx <- 2*cos(x) - x^2 + 2
hx <- log(x) + 2

plot(x, fx, type="o", pch=1, lty=1, col="red", 
     xlab="x", ylab="Mathematical Functions", 
     main="The Plot of Different Functions", 
     ylim=range(c(fx, gx, hx)))
lines(x, gx, type="o", pch=2, lty=2, col="blue")
lines(x, hx, type="o", pch=3, lty=3, col="green")
legend("topleft", legend=c("f(x)", "g(x)", "h(x)"), 
       pch=1:3, lty=1:3, col=c("red", "blue", "green"))

#c
set.seed(2026) # 固定亂數種子以便重現結果
reps <- 2000
n <- 50
true_param <- (3^2) / 9  # 真實值為 1

t1_est <- numeric(reps)
t2_est <- numeric(reps)

for(i in 1:reps) {
  samp <- rnorm(n, mean=3, sd=3) # sd = sqrt(variance)
  x_bar <- mean(samp)
  s2 <- var(samp)
  
  t1_est[i] <- (50 * x_bar^2) / (49 * s2)
  t2_est[i] <- (x_bar^2 - 1/50) / s2
}

# (i) 兩個估計量的平均值
mean_t1 <- mean(t1_est)
mean_t2 <- mean(t2_est)

# (ii) 絕對誤差平均值
mae_t1 <- mean(abs(t1_est - true_param))
mae_t2 <- mean(abs(t2_est - true_param))


#d
initial <- c(2, 0.25)
tol <- 1e-6

repeat {
  x_val <- initial[1]
  y_val <- initial[2]
  
  # 計算函數值 F(x)
  f1 <- x_val^2 - 2*x_val - y_val + 0.5
  f2 <- x_val^2 + 4*y_val^2 - 4
  F_mat <- c(f1, f2)
  
  # 計算 Jacobian 矩陣
  J <- matrix(c(2*x_val - 2, -1, 
                2*x_val,     8*y_val), nrow=2, byrow=TRUE)
  
  # 牛頓法迭代公式
  x_new <- initial - solve(J) %*% F_mat
  
  # 判斷是否收斂 (誤差小於 1e-6)
  if (sqrt(sum((x_new - initial)^2)) < tol) break
  
  initial <- x_new
}
ans_d <- round(x_new, 6) # 四捨五入至小數點後6位


#e
check_clt <- function(n, lambda, reps=1000) {
  z_scores <- numeric(reps)
  for(i in 1:reps) {
    samp <- rpois(n, lambda)
    x_bar <- mean(samp)
    var_xbar <- lambda / n
    z_scores[i] <- (x_bar - lambda) / sqrt(var_xbar)
  }
  # 畫 QQ 圖檢視常態性
  qqnorm(z_scores, main=paste("n=", n, ", lambda=", lambda))
  qqline(z_scores, col="red")
}

par(mfrow=c(2,2)) # 切割畫布為 2x2
check_clt(10, 1)   # (a)
check_clt(50, 0.1) # (b) 
check_clt(100, 0.1)# (c)
check_clt(1000, 1) # (d) 結論：n=1000 時常態近似最完美


#2
#a
check_clt <- function(n, lambda, reps=1000) {
  z_scores <- numeric(reps)
  for(i in 1:reps) {
    samp <- rpois(n, lambda)
    x_bar <- mean(samp)
    var_xbar <- lambda / n
    z_scores[i] <- (x_bar - lambda) / sqrt(var_xbar)
  }
  # 畫 QQ 圖檢視常態性
  qqnorm(z_scores, main=paste("n=", n, ", lambda=", lambda))
  qqline(z_scores, col="red")
}

par(mfrow=c(2,2)) # 切割畫布為 2x2
check_clt(10, 1)   # (a)
check_clt(50, 0.1) # (b) 
check_clt(100, 0.1)# (c)
check_clt(1000, 1) # (d) 結論：n=1000 時常態近似最完美


#b
sample1 <- c(37, 35, 35, 31, 37)
sample2 <- c(39, 38, 39, 41, 43)
sample3 <- c(33, 36, 35, 36, 40)
data_b <- data.frame(
  val = c(sample1, sample2, sample3),
  group = factor(rep(1:3, each=5))
)
summary(aov(val ~ group, data=data_b)) # 讀取 p-value 判斷是否拒絕 H0

#c
s1 <- c(11, 12, 8, 7, 7, 9)
s2 <- c(5, 8, 6, 7, 4, 8)
# (i) 單樣本 H0: mu >= 10 (左尾檢定)
t.test(s1, mu=10, alternative="less")
# (ii) 獨立雙樣本 H0: mu1 - mu2 <= 2 (右尾檢定)
t.test(s1, s2, mu=2, alternative="greater", var.equal=TRUE)
# (iii) 成對雙樣本 85% 信賴區間
t.test(s1, s2, paired=TRUE, conf.level=0.85)$conf.int

#d
# (i) QQ plots
d_data <- c(1.86, 0.51, -1.14, 0.33, 0.17, 0.67, -1.69, 1.02, 3.32, -0.20, 
            -0.27, 1.00, 0.39, -1.80, 0.68, 0.40, -2.50, 0.62, 0.49, 1.85)
par(mfrow=c(1,2))
qqplot(qt(ppoints(length(d_data)), df=10), d_data, main="QQ t-plot (df=10)")
qqnorm(d_data, main="QQ Normal N(0, 2)") # 預設標準常態，僅檢視形狀

# (ii) Poisson Goodness of Fit
freq_data <- c(3, 2, 3, 6, 2, 3, 0, 2, 1, 5, 5, 3, 6, 3, 5)
obs_table <- table(factor(freq_data, levels=0:6))
# 依照 Poisson(0.9) 理論機率進行檢定
chisq.test(obs_table, p=dpois(0:6, 0.9), rescale.p=TRUE, simulate.p.value=TRUE)


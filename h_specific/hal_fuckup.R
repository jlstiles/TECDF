library(hal9001)
library(hal)
set.seed(10)

# defining a propensity score
g0 = function(W1) plogis(-.1-.5*sin(W1) - .4*(abs(W1)>1)*W1^2)

# generating pscores and outcomes for sample of 1000
n=1000
W = rnorm(n)
ps = g0(W)
A = rbinom(n, 1, ps)

# going to fit over 1/2 the data
df_train = data.frame(W = W[1:(n/2)])

# defining validation over other half, one df including A
df_valA = data.frame(A = A[(n/2+1):n], W = W[(n/2+1):n])
df_val = data.frame(W = W[(n/2+1):n])

fitg_hal = fit_hal(X = df_train, Y = A[1:(n/2)], degrees = NULL, fit_type = "glmnet",
                   n_folds = 10, use_min = TRUE, family = "binomial",
                   return_lasso = FALSE, yolo = TRUE)

# predicting with dataframe including A should give an error instead of wrong answer
pred_halA = predict(fitg_hal, new_data = df_valA, type = 'response')
# prediction with dataframe not including A
pred_hal = predict(fitg_hal, new_data = df_val, type = 'response')

# big difference
pred_hal - pred_halA

# ridiculous risk for the wrong preds
risk_hal = mean(with(df_valA, -A*log(pred_hal)-(1-A)*log(1-pred_hal)))
risk_hal

risk_hal_wrong = mean(with(df_valA, -A*log(pred_halA)-(1-A)*log(1-pred_halA)))
risk_hal_wrong

# glm fcn gives correct preds either way
fit_glm = glm(A[1:(n/2)] ~ W, data = df_train, family = binomial)
pred_glmA = predict(fit_glm, newdata = df_valA, type = 'response')
pred_glm = predict(fit_glm, newdata = df_val, type = 'response')
pred_glmA - pred_glm

fitg_hal1 = hal(X = df, Y = A[1:(n/2)], family = binomial(), parallel = TRUE)

# gives error
pred_hal1 = predict(fitg_hal1, newdata = newX_withA, type = 'response')

# correct
pred_hal1 = predict(fitg_hal1, newdata = newX, type = 'response')
pred_hal1 - pred_hal

risk_hal1 = mean(with(newX_withA, -A*log(pred_hal1)-(1-A)*log(1-pred_hal1)))

risk_glm
risk_hal
risk_hal1

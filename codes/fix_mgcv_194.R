# fix_mgcv_194.R
# 临时修复 mgcv 1.9-4 中 bs = "re"（随机效应光滑）的 bug：
#   gam(y ~ s(g, bs = "re")) 报错 "data must be a data.frame"
#
# 原因：mgcv 1.9-x 用 ExtractData() 把平滑项数据以纯 list 传给
#   smooth.construct.re.smooth.spec / Predict.matrix.random.effect，
#   而这两个函数调用 model.matrix() 时要求 data.frame。
#
# 用法：在使用 gam() 之前 source 本文件即可，仅对当前 R 会话生效，
#   重启 R 后需要重新 source。
#
# 注：bs = "fs" 可能也有同样的问题，如遇到可照此扩展。

if (packageVersion("mgcv") == "1.9.4") {

  ns <- asNamespace("mgcv")

  # 原函数及其依赖放在独立环境里，避免污染全局环境
  patch_env <- new.env(parent = ns)
  patch_env$re_orig <- get("smooth.construct.re.smooth.spec", envir = ns)
  patch_env$pm_orig <- get("Predict.matrix.random.effect", envir = ns)

  # 构造阶段：把 list 数据强制转回 data.frame 再交给原函数
  re_fixed <- function(object, data, knots) {
    if (is.list(data) && !is.data.frame(data)) data <- as.data.frame(data)
    re_orig(object, data, knots)
  }
  environment(re_fixed) <- patch_env

  # 预测阶段：同样的处理
  pm_fixed <- function(object, data) {
    if (is.list(data) && !is.data.frame(data)) data <- as.data.frame(data)
    pm_orig(object, data)
  }
  environment(pm_fixed) <- patch_env

  unlockBinding("smooth.construct.re.smooth.spec", ns)
  unlockBinding("Predict.matrix.random.effect", ns)
  assignInNamespace("smooth.construct.re.smooth.spec", re_fixed, ns = "mgcv")
  assignInNamespace("Predict.matrix.random.effect", pm_fixed, ns = "mgcv")
  lockBinding("smooth.construct.re.smooth.spec", ns)
  lockBinding("Predict.matrix.random.effect", ns)

  message("mgcv 1.9-4 patch applied (bs='re' fixed for this session).")

} else {
  message("mgcv version is ", as.character(packageVersion("mgcv")),
          " -- patch not applied (only targets 1.9-4).")
}

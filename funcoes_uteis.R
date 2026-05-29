# --------------------------- FUNÇÃO DE ESCALA Y -----------------------------
addUnits <- function(n) {
  ifelse(n < 1e3, n,
         ifelse(n < 1e6, paste0(round(n / 1e3), "k"),
                ifelse(n < 1e9, paste0(round(n / 1e6), "M"),
                       ifelse(n < 1e12, paste0(round(n / 1e9), "B"),
                              ifelse(n < 1e15, paste0(round(n / 1e12), "T"), "too big!")))))
}

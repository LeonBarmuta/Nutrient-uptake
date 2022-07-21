NutUpTS <- function (times, y, parm){
  # Advection-dispersion with transient storage and first order uptake in main channel
  require(ReacTran)
  Cs <- y[1:200]
  Cts <- y[201:400]
  AFDW <- fiadeiro(v = parm[3], D = parm[2], grid = grid)
  trans <- tran.1D(C = Cs,
                   C.up = parm[1],
                   C.down = parm[1],
                   D = parm[2],
                   v = parm[3],
                   dx = grid,
                   AFDW = AFDW)
  ts.exchange <- parm[5] * (Cs-Cts)
  uptake <- parm[4] * Cs
  input <- parm[4] * parm[1]
  dCs <- trans$dC - ts.exchange - uptake + input
  dCts <- parm[6] * ts.exchange
  list(c(dCs, dCts))
}

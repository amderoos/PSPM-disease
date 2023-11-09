ReqFunc <- function(R, pars) {
  with(as.list(c(pars)), {
    z <- Sb / Sm
    funcresp <- M * R / (H + R)
    nuj <- SIGMA * (2 - Q) * funcresp - TS
    nua <- SIGMA *      Q  * funcresp - TS
    LRO <- (nua / MUS) * z^(MUS / nuj - 1)
    
    return(LRO-1)
  })
}

BeqFunc <- function(R, pars) {
  with(as.list(c(pars)), {
    z <- Sb / Sm
    funcresp <- M * R / (H + R)
    nuj <- SIGMA * (2 - Q) * funcresp - TS
    nua <- SIGMA *      Q  * funcresp - TS
    Fm  <- z^(MUS / nuj)
    G <- funcresp * ((2 - Q) * (Sm * Fm - Sb) / (nuj - MUS) + Q * Fm * Sm / MUS)
    Beq <- Rho * (Rmax - R) / G
    return(Beq)
  })
}

JnrFunc <- function(R, B, pars) {
  with(as.list(c(pars)), {
    z <- Sb / Sm
    funcresp <- M * R / (H + R)
    nuj <- SIGMA * (2 - Q) * funcresp - TS
    Fm  <- z^(MUS / nuj)
    return(B * (1 - Fm) / MUS)
  })
}

AnrFunc <- function(R, B, pars) {
  with(as.list(c(pars)), {
    z <- Sb / Sm
    funcresp <- M * R / (H + R)
    nuj <- SIGMA * (2 - Q) * funcresp - TS
    Fm  <- z^(MUS / nuj)
    return(B * Fm / MUS)
  })
}

JbioFunc <- function(R, B, pars) {
  with(as.list(c(pars)), {
    z <- Sb / Sm
    funcresp <- M * R / (H + R)
    nuj <- SIGMA * (2 - Q) * funcresp - TS
    Fm  <- z^(MUS / nuj)
    return(B * (Sm * Fm - Sb)/(nuj - MUS))
  })
}

AbioFunc <- function(R, B, pars) {
  with(as.list(c(pars)), {
    z <- Sb / Sm
    funcresp <- M * R / (H + R)
    nuj <- SIGMA * (2 - Q) * funcresp - TS
    Fm  <- z^(MUS / nuj)
    return(B * Fm * Sm / MUS)
  })
}

NtotFunc <- function(R, B, pars) {
  with(as.list(c(pars)), {
    z <- Sb / Sm
    funcresp <- M * R / (H + R)
    nuj <- SIGMA * (2 - Q) * funcresp - TS
    nua <- SIGMA *      Q  * funcresp - TS
    Fm <- z^(MUS / nuj)
    Ntot <- B / MUS
    return(Ntot)
  })
}

R0Func <- function(R, B, pars) {
  with(as.list(c(pars)), {
    z <- Sb / Sm
    funcresp <- M * R / (H + R)
    nuj  <- SIGMA * (2 - Q) * funcresp - TS
    nua  <- SIGMA *      Q  * funcresp - TS
    nuij <- SIGMA * (2 - Q) * funcresp - (TS + TI + TIJ)
    nuia <- SIGMA *      Q  * funcresp - (TS + TI + TIA)
    nuijplus <- max(nuij, 0)
    nuijmin  <- min(nuij, 0)
    nuiaplus <- max(nuia, 0)
    nuiamin  <- min(nuia, 0)
    dij <- (MUS + MUI + MUIJ - nuijmin)
    dia <- (MUS + MUI + MUIA - nuiamin)
    Fm  <- z^(MUS / nuj)

    R0 <- (1 - Fm) / (MUS * dij)
    R0 <- R0 + Fm  / (MUS * dia)

    if (nuijplus > 0) {
      Fmi  <- z^((MUS + MUI + MUIJ) / nuijplus)
      denom <- 1 / (nuj * dij - MUS * nuijplus)
      if (is.nan(denom) || is.infinite(denom)) {
        denom <- -log(z) * z^(MUS / nuj) / (nuijplus * nuj)
      }
      
      R0 <- R0 + nuijplus *  (1/dia - 1/dij) * (Fm - Fmi) * denom 
    }
    return(BETA * B * R0)
  })
}


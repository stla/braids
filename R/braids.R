isPositiveInteger <- function(n) {
  length(n) == 1L && is.numeric(n) && !is.na(n) && n != 0 && floor(n) == n
}

areIntegers <- function(x) {
  is.numeric(x) && !any(is.na(x)) && all(floor(x) == x)
}

Sigma <- function(i) {
  stopifnot(isPositiveInteger(i))
  c(as.integer(i), 1L)
}

SigmaInv <- function(i) {
  stopifnot(isPositiveInteger(i))
  c(as.integer(i), -1L)
}

#' @title Make a braid
#' @description Make a braid.
#'
#' @param n number of strands, an integer, at least 2
#' @param artingens Artin generators given by a vector of non-zero
#'   integers; a positive integer \code{i} corresponds to the
#'   standard positive Artin generator of a braid which
#'   represents twisting the neighbour strands \code{i} and \code{i+1},
#'   such that strand \code{i} goes \emph{under} strand \code{i+1}; a
#'   negative integer \code{-i} corresponds to the inverse.
#'
#' @return A \code{braid} object.
#' @export
#'
#' @examples
#' mkBraid(n = 4, c(2, -3))
mkBraid <- function(n, artingens) {
  stopifnot(isPositiveInteger(n), n >= 2)
  stopifnot(areIntegers(artingens), all(artingens != 0))
  if(any(abs(artingens) >= n)) {
    stop("Found a generator with a too large index.")
  }
  gens <- lapply(as.integer(artingens), function(i) {
    if(i > 0L) {
      Sigma(i)
    } else {
      SigmaInv(-i)
    }
  })
  mkBraid0(n, gens)
}

mkBraid0 <- function(n, gens) {
  out <- gens
  attr(out, "n") <- as.integer(n)
  class(out) <- "braid"
  out
}

#' @exportS3Method print braid
print.braid <- function(x, ...) {
  idx <- vapply(x, `[[`, integer(1L), 1L)
  signs <- vapply(x, `[[`, integer(1L), 2L)
  print(paste0(vapply(seq_along(idx), function(i) {
    if(signs[i] == 1L) {
      sprintf("sigma_%d", idx[i])
    } else {
      sprintf("sigmaInv_%d", idx[i])
    }
  }, character(1L)), collapse = "*"))
  invisible()
}

#' @title Number of strands
#' @description The number of strands of a braid.
#'
#' @param braid a \code{braid} object (e.g. created with \code{\link{mkBraid}})
#'
#' @return An integer.
#' @export
numberOfStrands <- function(braid) {
  attr(braid, "n")
}

#' @title Free reduction of a braid
#' @description Applies free reduction to a braid, i.e. removes pairs of
#'   consecutive generators inverse of each other.
#'
#' @param braid a \code{braid} object (e.g. created with \code{\link{mkBraid}})
#'
#' @return A \code{braid} object.
#' @export
#' @importFrom maybe just nothing is_nothing from_just
#'
#' @examples
#' braid <- mkBraid(4, c(2, -3, 3))
#' freeReduceBraidWord(braid)
freeReduceBraidWord <- function(braid) {
  reduceStep <- function(gens) {
    go <- function(changed, w) {
      if(length(w) >= 2L) {
        w1 <- w[[1L]]
        w2 <- w[[2L]]
        x <- w1[1L]
        y <- w2[1L]
        s1 <- w1[2L]
        s2 <- w2[2L]
        if(x == y && ((s1 == 1L && s2 == -1L) || (s1 == -1L && s2 == 1L))) {
          go(TRUE, w[-c(1L, 2L)])
        } else {
          gg <- go(changed, w[-1L])
          if(is_nothing(gg)) {
            nothing()
          } else {
            just(c(list(w1), from_just(gg)))
          }
        }
      } else if(length(w) == 1L) {
        if(changed) {
          just(c(list(w[[1L]]), w))
        } else {
          nothing()
        }
      } else {
        if(changed) {
          just(w)
        } else {
          nothing()
        }
      }
    }
    go(FALSE, gens)
  }
  loop <- function(w) {
    x <- reduceStep(w)
    if(is_nothing(x)) {
      w
    } else {
      loop(from_just(x))
    }
  }
  mkBraid0(numberOfStrands(braid), loop(braid))
}

#' @title Braid permutation
#' @description Returns the left-to-right permutation associated to a braid.
#'
#' @param braid a \code{braid} object (e.g. created with \code{\link{mkBraid}})
#'
#' @return A permutation.
#' @export
#'
#' @examples
#' braid <- mkBraid(4, c(2, -3, 3))
#' braidPermutation(braid)
braidPermutation <- function(braid) {
  n <- numberOfStrands(braid)
  idxs <- vapply(braid, `[[`, integer(1L), 1L)
  worker <- function(arr, idxs) {
    if(length(idxs) == 0L) {
      arr
    } else {
      i <- idxs[1L]
      is <- idxs[-1L]
      a <- arr[i]
      b <- arr[i + 1L]
      arr[i] <- b
      arr[i + 1L] <- a
      worker(arr, is)
    }
  }
  worker(1L:n, idxs)
}

#' @title Double generator
#' @description Generator \eqn{\sigma_{s,t}} in the Birman-Ko-Lee new
#'   presentation. It twists the strands \code{s} and \code{t} while going over
#'   all other strands (for \code{t=s+1}, this is \eqn{\sigma_s}).
#'
#' @param n number of strands, integer \code{>=2}
#' @param s,t indices of two strands, \code{s < t}
#'
#' @return A \code{braid} object.
#' @export
#'
#' @examples
#' doubleSigma(5, 1, 3)
doubleSigma <- function(n, s, t) {
  stopifnot(isPositiveInteger(n), isPositiveInteger(s), isPositiveInteger(t))
  stopifnot(n >= 2)
  if(s > n) {
    stop("The `s` index is out of range.")
  }
  if(t > n) {
    stop("The `t` index is out of range.")
  }
  if(s >= t) {
    stop("The `s` index must be strictly smaller than the `t` index.")
  }
  if(t - 1L <= s) {
    gens <- lapply((t-1L):s, Sigma)
  } else {
    gens <- lapply((s+1L):(t-1L), SigmaInv)
  }
  mkBraid0(n, gens)
}

#' @title Half-twist
#' @description The (positive) half-twist of all the braid strands, usually
#'   denoted by \eqn{\Delta}.
#'
#' @param n number of strands, integer \code{>=2}
#'
#' @return A \code{braid} object.
#' @export
#'
#' @examples
#' halfTwist(4)
halfTwist <- function(n) {
  stopifnot(isPositiveInteger(n), n >= 2)
  # if(n == 1L) {
  #   gens <- list()
  # } else {
  subs <- lapply(1L:(n-1L), function(k) {
    (n-1L):k
  })
  gens <- lapply(do.call(c, subs), Sigma)
  # }
  mkBraid0(n, gens)
}

#' @title Inner automorphism
#' @description The inner automorphism defined by
#'   \eqn{\tau X = \Delta^{-1} X \Delta}, where \eqn{\Delta} is the
#'   positive half-twist; it sends each generator \eqn{\sigma_j} to
#'   \eqn{\sigma_{n-j}}.
#'
#' @param braid a \code{braid} object
#'
#' @return A \code{braid} object.
#' @export
#'
#' @examples
#' braid <- mkBraid(4, c(2, -3, 3))
#' tau(braid)
tau <- function(braid) {
  n <- numberOfStrands(braid)
  gens <- lapply(braid, function(gen) {
    i <- gen[2L]
    if(i == 1L) {
      Sigma(n - i)
    } else {
      SigmaInv(n - i)
    }
  })
  mkBraid0(n, gens)
}

#' @title Inverse braid
#' @description The inverse of a braid (without performing reduction).
#'
#' @param braid a \code{braid} object
#'
#' @return A \code{braid} object.
#' @export
#'
#' @examples
#' braid <- mkBraid(4, c(2, -3, 3))
#' ibraid <- inverseBraid(braid)
#' composeTwoBraids(braid, ibraid)
inverseBraid <- function(braid) {
  n <- numberOfStrands(braid)
  invgens <- lapply(braid, function(gen) {
    c(gen[1L], -gen[2L])
  })
  mkBraid0(n, rev(invgens))
}

#' @title Composition of two braids
#' @description Composes two braids, doing free reduction on the result.
#'
#' @param braid1,braid2 \code{braid} objects with the same number of strands
#'
#' @return A \code{braid} object.
#' @export
#'
#' @examples
#' braid <- mkBraid(4, c(2, -3, 3))
#' composeTwoBraids(braid, braid)
composeTwoBraids <- function(braid1, braid2) {
  n <- numberOfStrands(braid1)
  if(n != numberOfStrands(braid2)) {
    stop("Unequal numbers of strands.")
  }
  freeReduceBraidWord(mkBraid0(n, c(braid1, braid2)))
}

#' @title Composition of many braids.
#' @description Composes many braids, doing free reduction on the result.
#'
#' @param braids list of \code{braid} objects with the same number of strands
#'
#' @return A \code{braid} object.
#' @export
#'
#' @examples
#' braid <- mkBraid(4, c(2, -3, 3))
#' composeManyBraids(list(braid, braid, braid))
composeManyBraids <- function(braids) {
  ns <- vapply(braids, numberOfStrands, integer(1L))
  n <- ns[1L]
  if(any(ns != n)) {
    stop("Unequal numbers of strands.")
  }
  freeReduceBraidWord(mkBraid0(n, do.call(c, braids)))
}

#' @title Whether a braid is pure
#' @description Checks whether a braid is pure, i.e. its permutation is trivial.
#'
#' @param braid a \code{braid} object
#'
#' @return A Boolean value.
#' @export
#'
#' @examples
#' braid <- mkBraid(4, c(2, -3, 3))
#' isPureBraid(braid)
isPureBraid <- function(braid) {
  identical(braidPermutation(braid), seq_len(numberOfStrands(braid)))
}

#' @title Whether a braid is positive
#' @description Checks whether a braid has only positive Artin generators.
#'
#' @param braid a \code{braid} object
#'
#' @return A Boolean value.
#' @export
#'
#' @examples
#' braid <- mkBraid(4, c(2, -3, 3))
#' isPositiveBraidWord(braid)
isPositiveBraidWord <- function(braid) {
  signs <- vapply(braid, `[[`, integer(1L), 2L)
  all(signs == 1L)
}

#' @title Linking matrix
#' @description Linking numbers between all pairs of strands of a braid.
#'
#' @param braid a \code{braid} object
#'
#' @return A matrix.
#' @export
#' @seealso \code{\link{strandLinking}} to get the linking number between
#'   two strands of the braid.
#'
#' @examples
#' braid <- mkBraid(4, c(2, -3, 3))
#' linkingMatrix(braid)
linkingMatrix <- function(braid) {
  n <- numberOfStrands(braid)
  perm <- 1L:n
  doSwap <- function(i) {
    a <- perm[i]
    b <- perm[i+1L]
    perm[i]    <<- b
    perm[i+1L] <<- a
    invisible()
  }
  mat <- matrix(0L, nrow = n, ncol = n)
  doAdd <- function(i, j, pm1) {
    x <- mat[i, j]
    mat[i, j] <<- mat[j, i] <<- x + pm1
    invisible()
  }
  for(gen in braid) {
    k <- gen[1L]
    u <- perm[k]
    v <- perm[k+1L]
    doAdd(u, v, gen[2L])
    doSwap(k)
  }
  mat
}

#' @title Linking number between two strands
#' @description The linking number between two strands of a braid.
#'
#' @param braid a \code{braid} object
#' @param i,j indices of two strands
#'
#' @return An integer.
#' @export
#' @seealso \code{\link{linkingMatrix}} to get the linking numbers between
#'   all pairs of strands of the braid.
#'
#' @examples
#' braid <- mkBraid(4, c(2, -3, 3))
#' strandLinking(braid, 1, 3)
strandLinking <- function(braid, i, j) {
  n <- numberOfStrands(braid)
  stopifnot(isPositiveInteger(i), isPositiveInteger(j), i <= n, j <= n)
  if(i == j) {
    0L
  } else {
    go <- function(s, t, gens) {
      if(length(gens) == 0L) {
        0L
      } else {
        g <- gens[[1L]]
        gs <- gens[-1L]
        k <- g[1L]
        sgn <- g[2L]
        if(s == k && t == k+1L) {
          sgn + go(s + 1L, t - 1L, gs)
        } else if(t == k && s == k+1L) {
          sgn + go(s - 1L, t + 1L, gs)
        } else if(s == k) {
          go(s + 1L, t, gs)
        } else if(s == k+1L) {
          go(s - 1L, t, gs)
        } else if(t == k) {
          go(s, t + 1L, gs)
        } else if(t == k + 1L) {
          go(s, t - 1L, gs)
        } else {
          go(s, t, gs)
        }
      }
    }
    go(i, j, braid)
  }
}

#' @title Whether a braid is a permutation braid
#' @description Checks whether a braid is a permutation braid, that is,
#'   a positive braid where any two strands cross at most one, and positively.
#'
#' @param braid a \code{braid} object
#'
#' @return A Boolean value.
#' @export
#'
#' @examples
#' braid <- mkBraid(4, c(2, -3, 3))
#' isPermutationBraid(braid)
isPermutationBraid <- function(braid) {
  if(isPositiveBraidWord(braid)) {
    lkMatrix <- linkingMatrix(braid)
    all(lkMatrix[upper.tri(lkMatrix)] %in% c(0L, 1L))
  } else {
    FALSE
  }
}

isPermutation <- function(x) {
  setequal(x, seq_along(x))
}

.permutationBraid <- function(perm) {
  n <- length(perm)
  cfwd <- cinv <- 1L:n
  doSwap <- function(i) {
    a <- cinv[i]
    b <- cinv[i+1L]
    cinv[i]    <<- b
    cinv[i+1L] <<- a
    u <- cfwd[a]
    v <- cfwd[b]
    cfwd[a] <<- v
    cfwd[b] <<- u
    invisible()
  }
  worker <- function(phase) {
    if(phase >= n) {
      list()
    } else {
      tgt <- perm[phase]
      src <- cfwd[tgt]
      this <- (src-1L):phase
      lapply(this, doSwap)
      rest <- worker(phase + 1L)
      c(list(this), rest)
    }
  }
  worker(1L)
}

#' @title Permutation braid
#' @description Makes a permutation braid from a permutation.
#'
#' @param perm a permutation
#'
#' @return A \code{braid} object.
#' @export
#'
#' @examples
#' perm <- c(3, 1, 4, 2)
#' braid <- permutationBraid(perm)
#' isPermutationBraid(braid)
#' braidPermutation(braid)
permutationBraid <- function(perm) {
  stopifnot(isPermutation(perm), length(perm) >= 2)
  gens <- lapply(do.call(c, .permutationBraid(perm)), Sigma)
  mkBraid0(length(perm), gens)
}

.allPositiveBraidWords <- function(n, l) {
  go <- function(k) {
    if(k == 0L) {
      list(list())
    } else {
      do.call(c, lapply(1L:(n-1L), function(i) {
        lapply(go(k - 1L), function(rest) {
          c(list(Sigma(i)), rest)
        })
      }))
    }
  }
  go(l)
}

#' @title Positive braid words of given length
#' @description All positive braid words of the given length.
#'
#' @param n number of strands, integer \code{>=2}
#' @param l length of the words
#'
#' @return A list of \code{braid} objects.
#' @export
#'
#' @examples
#' allPositiveBraidWords(3, 4)
allPositiveBraidWords <- function(n, l) {
  stopifnot(isPositiveInteger(n), n >= 2)
  lapply(.allPositiveBraidWords(n, l), function(gens) {
    mkBraid0(n, gens)
  })
}

.allBraidWords <- function(n, l) {
  go <- function(k) {
    if(k == 0L) {
      list(list())
    } else {
      gens <- do.call(c, lapply(1L:(n-1L), function(i) {
        c(list(Sigma(i)), list(SigmaInv(i)))
      }))
      do.call(c, lapply(go(k - 1L), function(rest) {
        lapply(gens, function(gen) {
          c(list(gen), rest)
        })
      }))
    }
  }
  go(l)
}

#' @title Braid words of given length
#' @description All braid words of the given length.
#'
#' @param n number of strands, integer \code{>=2}
#' @param l length of the words
#'
#' @return A list of \code{braid} objects.
#' @export
#'
#' @examples
#' allBraidWords(3, 2)
allBraidWords <- function(n, l) {
  stopifnot(isPositiveInteger(n), n >= 2)
  lapply(.allBraidWords(n, l), function(gens) {
    mkBraid0(n, gens)
  })
}

longZipWith <- function(a0, b0, f, a_, b_) {
  go <- function(xxs, yys) {
    if(length(xxs) == 0L) {
      vapply(yys, function(y) {
        f(a0, y)
      }, integer(1L))
    } else if(length(yys) == 0L) {
      vapply(xxs, function(x) {
        f(x, b0)
      }, integer(1L))
    } else {
      c(f(xxs[[1L]], yys[[1L]]), go(xxs[-1L], yys[-1L]))
    }
  }
  go(a_, b_)
}

addSeries <- function(xs, ys) {
  longZipWith(0L, 0L, `+`, xs, ys)
}

sumSeries <- function(xs) {
  if(length(xs) == 0L) {
    0L
  } else {
    Reduce(addSeries, xs)
  }
}

#' @title Bronfman polynomials
#' @description The Bronfman polynomial of a braid group is the reciprocal of
#'   the growth function of the positive braids. This function computes the
#'   Bronfman polynomial of the braid group on \code{n} strands for \code{n}
#'   going to \code{1} to \code{N}.
#'
#' @param N maximum number of strands
#'
#' @return A list of integer vectors representing the Bronfman polynomials;
#'   each vector represents the polynomial coefficients in increasing order.
#' @export
#'
#' @examples
#' bronfmanPolynomials(3) # 1, 1 - X, 1 - 2X + X^3
bronfmanPolynomials <- function(N) {
  sgn <- function(i, x) {
    if(i %% 2L == 0L) {
      -x
    } else {
      x
    }
  }
  choose2 <- function(k) {
    (k * (k - 1L)) %/% 2L
  }
  go <- function(n) {
    if(n == 0L) {
      1L
    } else {
      sumSeries(lapply(1L:n, function(i) {
        e <- go(n - i)
        v <- c(rep(0L, choose2(i)), e)
        sgn(i, v)
      }))
    }
  }
  lapply(1L:N, go)
}

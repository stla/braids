isPositiveInteger <- function(n) {
  length(n) == 1L && is.numeric(n) && !is.na(n) && n != 0 && floor(n) == n
}

#' @title Artin generators
#' @description A standard Artin generator of a braid: \code{Sigma(i)}
#'   represents twisting the neighbour strands \code{i} and \code{i+1},
#'   such that strand \code{i} goes \emph{under} strand \code{i+1}.
#'
#' @param i index of the strand, a positive integer
#'
#' @returns A vector of two integers.
#' @export
#' @name ArtinGenerators
#' @rdname ArtinGenerators
Sigma <- function(i) {
  stopifnot(isPositiveInteger(i))
  c(as.integer(i), 1L)
}

#' @export
#' @rdname ArtinGenerators
SigmaInv <- function(i) {
  stopifnot(isPositiveInteger(i))
  c(as.integer(i), -1L)
}

#' @title Make a braid
#' @description Make a braid.
#'
#' @param n number of strands, an integer, at least 2
#' @param brgens list of generators obtained with \code{\link{Sigma}} or
#'   \code{\link{SigmaInv}}
#'
#' @return A \code{braid} object.
#' @export
#'
#' @examples
#' mkBraid(4, list(Sigma(2), SigmaInv(3)))
mkBraid <- function(n, brgens) {
  stopifnot(isPositiveInteger(n), n >= 2)
  out <- brgens
  idx <- vapply(brgens, `[[`, integer(1L), 1L)
  if(any(idx) >= n) {
    stop("Found a generator with a too large index.")
  }
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
  }, character(1L)), collapse = " "))
  invisible()
}

#' @title Number of strands
#' @description The number of strands of a braid.
#'
#' @param braid a \code{braid} object created with \code{\link{mkBraid}}
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
#' @param braid a \code{braid} object created with \code{\link{mkBraid}}
#'
#' @return A \code{braid} object.
#' @export
#' @importFrom maybe just nothing is_nothing from_just
#'
#' @examples
#' braid <- mkBraid(4, list(Sigma(2), SigmaInv(3), Sigma(3)))
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
  mkBraid(numberOfStrands(braid), loop(braid))
}

#' @title Braid permutation
#' @description Returns the left-to-right permutation associated to a braid.
#'
#' @param braid a \code{braid} object created with \code{\link{mkBraid}}
#'
#' @return A permutation.
#' @export
#'
#' @examples
#' braid <- mkBraid(4, list(Sigma(2), SigmaInv(3), Sigma(3)))
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
#' @description Generator \code{sigma_{s,t}} in the Birman-Ko-Lee new
#'   presentation. It twistes the strands \code{s} and \code{t} whie going over
#'   all other strands (for \code{t=s+1}, this is \code{sigma_s}).
#'
#' @param n number of strands, integer \code{>= 2}
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
    stop("`s` must be strictly smaller than `t`.")
  }
  if(t - 1L <= s) {
    gens <- lapply((t-1L):s, Sigma)
  } else {
    gens <- lapply((s+1L):(t-1L), SigmaInv)
  }
  mkBraid(n, gens)
}

#' @title Half-twist
#' @description The (positive) half-twist of all the braid strands, usually
#'   denoted by \eqn{\Delta}.
#'
#' @param n number of strands, integer \code{>= 2}
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
  mkBraid(n, gens)
}

#' @title Inner automorphism
#' @description The inner automorphism defined by
#'   \eqn{\tau X = \Delta^{-1} X \Delta}, where \eqn{\Delta} is the
#'   positive half-twist; it send each generator \eqn{\sigma_j} to
#'   \eqn{\sigma_{n-j}}.
#'
#' @param braid a \code{braid} object created with \code{\link{mkBraid}}
#'
#' @return A \code{braid} object.
#' @export
#'
#' @examples
#' braid <- mkBraid(4, list(Sigma(2), SigmaInv(3), Sigma(3)))
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
  mkBraid(n, gens)
}

#' @title Inverse braid
#' @description The inverse of a braid (without performing reduction).
#'
#' @param braid a \code{braid} object created with \code{\link{mkBraid}}
#'
#' @return A \code{braid} object.
#' @export
#'
#' @examples
#' braid <- mkBraid(4, list(Sigma(2), SigmaInv(3), Sigma(3)))
#' ibraid <- inverseBraid(braid)
#' composeTwoBraids(braid, ibraid)
inverseBraid <- function(braid) {
  n <- numberOfStrands(braid)
  invgens <- lapply(braid, function(gen) {
    c(gen[1L], -gen[2L])
  })
  mkBraid(n, rev(invgens))
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
#' braid <- mkBraid(4, list(Sigma(2), SigmaInv(3), Sigma(3)))
#' composeTwoBraids(braid, braid)
composeTwoBraids <- function(braid1, braid2) {
  n <- numberOfStrands(braid1)
  if(n != numberOfStrands(braid2)) {
    stop("Unequal numbers of strands.")
  }
  freeReduceBraidWord(mkBraid(n, c(braid1, braid2)))
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
#' braid <- mkBraid(4, list(Sigma(2), SigmaInv(3), Sigma(3)))
#' composeManyBraids(list(braid, braid, braid))
composeManyBraids <- function(braids) {
  ns <- vapply(braids, numberOfStrands, integer(1L))
  n <- ns[1L]
  if(any(ns != n)) {
    stop("Unequal numbers of strands.")
  }
  freeReduceBraidWord(mkBraid(n, do.call(c, braids)))
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
#' braid <- mkBraid(4, list(Sigma(2), SigmaInv(3), Sigma(3)))
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
#' braid <- mkBraid(4, list(Sigma(2), SigmaInv(3), Sigma(3)))
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
#'
#' @examples
#' braid <- mkBraid(4, list(Sigma(2), SigmaInv(3), Sigma(3)))
#' linkingMatrix(braid)
linkingMatrix <- function(braid) {
  n <- numberOfStrands(braid)
  perm <- 1L:n
  doSwap <- function(i) {
    a <- perm[i]
    b <- perm[i+1L]
    perm[i] <<- b
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
#' braid <- mkBraid(4, list(Sigma(2), SigmaInv(3), Sigma(3)))
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
    cinv[i] <<- b
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
  mkBraid(length(perm), gens)
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
#' @param n number of strands, positive integer \code{>= 2}
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
    mkBraid(n, gens)
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
#' @param n number of strands, positive integer \code{>= 2}
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
    mkBraid(n, gens)
  })
}
hSepString <- function(hsep) {
  sep <- attr(hsep, "sep")
  switch(
    sep,
    "empty"  = "",
    "spaces" = paste0(rep(" ", hsep), collapse = ""),
    "string" = hsep
  )
}

vSepString <- function(vsep) {
  sep <- attr(vsep, "sep")
  switch(
    sep,
    "empty"  = character(0L),
    "spaces" = paste0(rep(" ", vsep), collapse = ""),
    "string" = vsep
  )
}

vSepSpaces <- function(k) {
  out <- k
  attr(out, "sep") <- "spaces"
  out
}

vSepSize <- function(vsep) {
  nchar(vSepString(vsep))
}

HSepEmpty <- function() {
  out <- "."
  attr(out, "sep") <- "empty"
  out
}

ASCII <- function(x, y, lines) {
  list("x" = x, "y" = y, "lines" = lines)
}

asciiLines <- function(ascii) {
  ascii[["lines"]]
}

vExtendTo <- function(valign, n0, rect) {
  y <- rect[["y"]]
  vExtendWith(valign, max(n0, y) - y, rect)
}

vExtendWith <- function(valign, d, rect) {
  x <- rect[["x"]]
  y <- rect[["y"]]
  lines <- rect[["lines"]]
  emptyLine <- paste0(rep(" ", x), collapse = "")
  f <- function(ls) {
    a <- d %/% 2L
    switch(
      valign,
      "Vtop"    = c(ls, rep(emptyLine, d)),
      "Vbottom" = c(rep(emptyLine, d), ls),
      "Vcenter" = c(rep(emptyLine, a), ls, rep(emptyLine, d - a))
    )
  }
  ASCII(x, y + d, f(lines))
}

hExtendTo <- function(halign, n0, rect) {
  x <- rect[["x"]]
  hExtendWith(halign, max(n0, x) - x, rect)
}

hExtendWith <- function(halign, d, rect) {
  x <- rect[["x"]]
  y <- rect[["y"]]
  lines <- rect[["lines"]]
  f <- function(l) {
    a <- d %/% 2L
    switch(
      halign,
      "Hleft"   = paste0(c(l, rep(" ", d)), collapse = ""),
      "Hright"  = paste0(c(rep(" ", d), l), collapse = ""),
      "Hcenter" = paste0(c(rep(" ", a), l, rep(" ", d - a)), collapse = "")
    )
  }
  ASCII(x + d, y, vapply(lines, f, character(1L)))
}

asciiFromLines <- function(ls) {
  y <- length(ls)
  x <- max(vapply(ls, nchar, integer(1L)))
  f <- function(l) {
    paste0(c(l, rep(" ", x - nchar(l))), collapse = "")
  }
  ASCII(x, y, vapply(ls, f, character(1L)))
}

hCatTop <- function(asciis) {
  hCatWith("Vtop", HSepEmpty(), asciis)
}

hCatWith <- function(valign, hsep, rects) {
  n <- length(rects)
  maxy <- max(vapply(rects, `[[`, integer(1L), "y"))
  xsz <- vapply(rects, `[[`, integer(1L), "x")
  sep <- hSepString(hsep)
  sepx <- length(sep)
  rects1 <- lapply(rects, function(rect) {
    vExtendTo(valign, maxy, rect)
  })
  x2 <- sum(xsz) + (n - 1L) * sepx
  M <- do.call(rbind, lapply(rects1, asciiLines))
  final <- apply(M, 2L, paste0, collapse = sep)
  ASCII(x2, maxy, final)
}

intercalate <- function(sep, x) {
  if(length(x) == 0L) {
    list()
  } else if(length(x) == 1L) {
    x[1L]
  } else {
    unlist(c(x[1L], list(sep), intercalate(sep, x[-1L])))
  }
}

vCatWith <- function(halign, vsep, rects) {
  n <- length(rects)
  maxx <- max(vapply(rects, `[[`, integer(1L), "x"))
  ysz <- vapply(rects, `[[`, integer(1L), "y")
  sepy <- vSepSize(vsep)
  vsepstring <- vSepString(vsep)
  fullsep <- apply(
    rbind(
      vapply(strsplit(vsepstring, "")[[1L]], rep, character(maxx), times = maxx)
    ),
    2L, paste0, collapse = ""
  )
  rects1 <- lapply(rects, function(rect) {
    hExtendTo(halign, maxx, rect)
  })
  y2 <- sum(ysz) + (n - 1L) * sepy
  final <- intercalate(fullsep, lapply(rects1, asciiLines))
  ASCII(maxx, y2, final)
}

filledBox <- function(c, x0, y0) {
  x <- max(0L, x0)
  y <- max(0L, y0)
  asciiFromLines(rep(paste0(rep(c, x), collapse = ""), y))
}

transparentBox <- function(x0, y0) {
  filledBox(" ", x0, y0)
}

asciiShow <- function(x) {
  asciiFromLines(x)
}

horizBraidASCII <- function(flipped, braid) {
  under  <- c("\\ /" , " / "  , "/ \\")
  over   <- c("\\ /" , " \\ " , "/ \\")
  horiz  <- c("   "  , "   "  , "___")
  space3 <- c("   "  , "   "  , "   ")
  n <- numberOfStrands(braid)
  block2 <- function(i, middle) {
    if(flipped) {
      a <- n - i - 1L
      b <- i - 1L
    } else {
      a <- i - 1L
      b <- n - i - 1L
    }
    x <- c(rep(horiz, a), c(space3, middle), rep(horiz, b))
    asciiFromLines(x[-c(1L, 2L)])
  }
  block <- function(g) {
    i <- g[1L]
    if(g[2L] == 1L) {
      block2(i, if(flipped) over else under)
    } else {
      block2(i, if(flipped) under else over)
    }
  }
  mkNumbers <- function(x) {
    if(flipped) x <- rev(x)
    vCatWith(
      "Hright",
      vSepSpaces(2L),
      lapply(x, asciiShow)
    )
  }
  spaceBlock    <- transparentBox(1L, 3L*n - 2L)
  beginEndBlock <- asciiFromLines(rep(horiz, n)[-c(1L, 2L)])
  numberBlock   <- mkNumbers(1L:n)
  numberBlock2  <- mkNumbers(braidPermutation(braid))
  prelude       <- list(numberBlock, spaceBlock, beginEndBlock)
  epilogue      <- list(beginEndBlock, spaceBlock, numberBlock2)
  middleBlocks  <- lapply(braid, block)
  allBlocks     <- c(prelude, middleBlocks, epilogue)
  vExtendWith("Vtop", 1L, hCatTop(allBlocks))
}

#' @title ASCII braid
#' @description Prints an ASCII figure of a braid.
#'
#' @param braid a \code{braid} object
#'
#' @return No value is returned, just prints the ASCII figure.
#' @export
#'
#' @examples
#' braid <- mkBraid(4, c(1, -2))
#' braidASCII(braid)
braidASCII <- function(braid) {
  ascii <- horizBraidASCII(FALSE, braid)
  cat(ascii[["lines"]], sep = "\n")
  invisible()
}

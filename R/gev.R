is.Numeric <-
function (x, length.arg = Inf, integer.valued = FALSE, positive = FALSE)
if (all(is.numeric(x)) && all(is.finite(x)) && (if (is.finite(length.arg)) length(x) ==
    length.arg else TRUE) && (if (integer.valued) all(x == round(x)) else TRUE) &&
    (if (positive) all(x > 0) else TRUE)) TRUE else FALSE

dgev <-
function (x, location = 0, scale = 1, shape = 0, log = FALSE,
    tolshape0 = sqrt(.Machine$double.eps))
{
    oobounds.log <- -Inf
    if (!is.logical(log.arg <- log) || length(log) != 1)
        stop("bad input for argument 'log'")
    rm(log)
    if (!is.Numeric(tolshape0, length.arg = 1, positive = TRUE))
        stop("bad input for argument 'tolshape0'")
    use.n <- max(length(x), length(location), length(scale),
        length(shape))
    if (length(shape) != use.n)
        shape <- rep_len(shape, use.n)
    if (length(location) != use.n)
        location <- rep_len(location, use.n)
    if (length(scale) != use.n)
        scale <- rep_len(scale, use.n)
    if (length(x) != use.n)
        x <- rep_len(x, use.n)
    logdensity <- rep_len(log(0), use.n)
    scase <- (abs(shape) < tolshape0)
    nscase <- sum(scase)
    if (use.n - nscase) {
        zedd <- 1 + shape * (x - location)/scale
        xok <- (!scase) & (zedd > 0)
        logdensity[xok] <- -log(scale[xok]) - zedd[xok]^(-1/shape[xok]) -
            (1 + 1/shape[xok]) * log(zedd[xok])
        outofbounds <- (!scase) & (zedd <= 0)
        if (any(outofbounds)) {
            logdensity[outofbounds] <- oobounds.log
            no.oob <- sum(outofbounds)
        }
    }
    if (nscase) {
        logdensity[scase] <- dgumbel(x[scase], location = location[scase],
            scale = scale[scase], log = TRUE)
    }
    logdensity[scale <= 0] <- NaN
    logdensity[is.infinite(x)] <- log(0)
    if (log.arg)
        logdensity
    else exp(logdensity)
}

pgev <-
function (q, location = 0, scale = 1, shape = 0, lower.tail = TRUE,
    log.p = FALSE)
{
    if (!is.logical(lower.tail) || length(lower.tail) != 1)
        stop("bad input for argument 'lower.tail'")
    if (!is.logical(log.arg <- log.p) || length(log.p) != 1)
        stop("bad input for argument 'log.p'")
    use.n <- max(length(q), length(location), length(scale),
        length(shape))
    if (length(shape) != use.n)
        shape <- rep_len(shape, use.n)
    if (length(location) != use.n)
        location <- rep_len(location, use.n)
    if (length(scale) != use.n)
        scale <- rep_len(scale, use.n)
    if (length(q) != use.n)
        q <- rep_len(q, use.n)
    scase0 <- abs(shape) < sqrt(.Machine$double.eps)
    zedd <- (q - location)/scale
    use.zedd <- pmax(0, 1 + shape * zedd)
    if (lower.tail) {
        if (log.p) {
            ans <- -use.zedd^(-1/shape)
        }
        else {
            ans <- exp(-use.zedd^(-1/shape))
        }
    }
    else {
        if (log.p) {
            ans <- log(-expm1(-use.zedd^(-1/shape)))
        }
        else {
            ans <- -expm1(-use.zedd^(-1/shape))
        }
    }
    if (any(scase0)) {
        ans[scase0] <- pgumbel(q[scase0], location = location[scase0],
            scale = scale[scase0], lower.tail = lower.tail, log.p = log.p)
    }
    ans[scale <= 0] <- NaN
    ans
}

qgev <-
function (p, location = 0, scale = 1, shape = 0, lower.tail = TRUE,
    log.p = FALSE)
{
    if (!is.logical(log.p) || length(log.p) != 1)
        stop("bad input for argument 'log.p'")
    use.n <- max(length(p), length(location), length(scale),
        length(shape))
    if (length(shape) != use.n)
        shape <- rep_len(shape, use.n)
    if (length(location) != use.n)
        location <- rep_len(location, use.n)
    if (length(scale) != use.n)
        scale <- rep_len(scale, use.n)
    if (length(p) != use.n)
        p <- rep_len(p, use.n)
    scase0 <- abs(shape) < sqrt(.Machine$double.eps)
    if (lower.tail) {
        if (log.p) {
            ln.p <- p
            ans <- location + scale * ((-ln.p)^(-shape) - 1)/shape
            ans[ln.p > 0] <- NaN
        }
        else {
            ans <- location + scale * ((-log(p))^(-shape) - 1)/shape
            ans[p == 1] <- Inf
            ans[p > 1] <- NaN
        }
    }
    else {
        if (log.p) {
            ln.p <- p
            ans <- location + scale * ((-log1p(-exp(ln.p)))^(-shape) -
                1)/shape
            ans[ln.p > 0] <- NaN
        }
        else {
            ans <- location + scale * ((-log1p(-p))^(-shape) -
                1)/shape
            ans[p == 1] <- Inf
            ans[p > 1] <- NaN
            ans[p < 0] <- NaN
        }
    }
    if (any(scase0))
        ans[scase0] <- qgumbel(p[scase0], location = location[scase0],
            scale = scale[scase0], lower.tail = lower.tail, log.p = log.p)
    ans[scale <= 0] <- NaN
    ans
}

rgev <-
function (n, location = 0, scale = 1, shape = 0)
{
    use.n <- if ((length.n <- length(n)) > 1)
        length.n
    else if (!is.Numeric(n, integer.valued = TRUE, length.arg = 1,
        positive = TRUE))
        stop("bad input for argument 'n'")
    else n
    if (!is.Numeric(location))
        stop("bad input for argument argument 'location'")
    if (!is.Numeric(shape))
        stop("bad input for argument argument 'shape'")
    ans <- numeric(use.n)
    if (length(shape) != use.n)
        shape <- rep_len(shape, use.n)
    if (length(location) != use.n)
        location <- rep_len(location, use.n)
    if (length(scale) != use.n)
        scale <- rep_len(scale, use.n)
    scase <- abs(shape) < sqrt(.Machine$double.eps)
    nscase <- sum(scase)
    if (use.n - nscase)
        ans[!scase] <- location[!scase] + scale[!scase] * ((-log(runif(use.n -
            nscase)))^(-shape[!scase]) - 1)/shape[!scase]
    if (nscase)
        ans[scase] <- rgumbel(nscase, location = location[scase],
            scale = scale[scase])
    ans[scale <= 0] <- NaN
    ans
}

dgumbel <-
function (x, location = 0, scale = 1, log = FALSE)
{
    if (!is.logical(log.arg <- log) || length(log) != 1)
        stop("bad input for argument 'log'")
    rm(log)
    zedd <- (x - location)/scale
    logdensity <- -zedd - exp(-zedd) - log(scale)
    logdensity[is.infinite(x)] <- log(0)
    if (log.arg)
        logdensity
    else exp(logdensity)
}

pgumbel <-
function (q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
{
    if (!is.logical(lower.tail) || length(lower.tail) != 1)
        stop("bad input for argument 'lower.tail'")
    if (!is.logical(log.p) || length(log.p) != 1)
        stop("bad input for argument 'log.p'")
    if (lower.tail) {
        if (log.p) {
            ans <- -exp(-(q - location)/scale)
            ans[q <= -Inf] <- -Inf
            ans[q == Inf] <- 0
        }
        else {
            ans <- exp(-exp(-(q - location)/scale))
            ans[q <= -Inf] <- 0
            ans[q == Inf] <- 1
        }
    }
    else {
        if (log.p) {
            ans <- log(-expm1(-exp(-(q - location)/scale)))
            ans[q <= -Inf] <- 0
            ans[q == Inf] <- -Inf
        }
        else {
            ans <- -expm1(-exp(-(q - location)/scale))
            ans[q <= -Inf] <- 1
            ans[q == Inf] <- 0
        }
    }
    ans[scale <= 0] <- NaN
    ans
}

qgumbel <-
function (p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
{
    if (!is.logical(lower.tail) || length(lower.tail) != 1)
        stop("bad input for argument 'lower.tail'")
    if (!is.logical(log.p) || length(log.p) != 1)
        stop("bad input for argument 'log.p'")
    if (lower.tail) {
        if (log.p) {
            ln.p <- p
            ans <- location - scale * log(-ln.p)
        }
        else {
            ans <- location - scale * log(-log(p))
            ans[p == 0] <- -Inf
            ans[p == 1] <- Inf
        }
    }
    else {
        if (log.p) {
            ln.p <- p
            ans <- location - scale * log(-log(-expm1(ln.p)))
            ans[ln.p > 0] <- NaN
        }
        else {
            ans <- location - scale * log(-log1p(-p))
            ans[p == 0] <- Inf
            ans[p == 1] <- -Inf
        }
    }
    ans[scale <= 0] <- NaN
    ans
}

rgumbel <-
function (n, location = 0, scale = 1)
{
    answer <- location - scale * log(-log(runif(n)))
    answer[scale <= 0] <- NaN
    answer
}

pointinTR <-
function (dpoint, tregion) 
{
    delta = tregion[, 1:(ncol(tregion) - 1)] %*% dpoint + tregion[, 
        ncol(tregion)]
    if (max(delta) > 0) 
        FALSE
    else TRUE
}


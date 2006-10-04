plotsingle <- function (cm, labels, method, ndigits, ...)
{
	plot.new()

	oldmar = par ("mar")

	par (mar = rep(1, 4))
	p = ncol (cm)
	lim = c(-1, p + 1)

	plot.window(xlim = lim, ylim = lim, xaxs = "i", yaxs = "i")

	for (i in 1:p)
	{
		text (i - 0.5, p + 0.5, labels[i], srt=90)
		text (-0.5, p - i + 0.5, labels[i])
	}

	for (i in 2:p)
	{
		for (j in 1:(i-1))
		{
			doEllipses (cm[c(i,j), c(i,j)], pos = c(i - 1.5,p - j - 0.5), lwd = 2, ...)
			text (j - 0.5 ,p - i + 0.5, round (cm[i,j], ndigits), ...)
		}
	}

	lines (c(0.5, p-0.5), c(p - 0.5, 0.5), lwd = 3)

	lines (lim[2] - 2 + c(-0.5, 0.3), c(-0.5, -0.5), col = col, ...)
	text (lim[2] - 2 - 0.7, -0.475, method, pos = 2)

	par (mar = oldmar)
}

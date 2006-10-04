plotcomp <- function (cor1, cor2, labels1, labels2, method1, method2, ndigits = 4, ...)	## internal function
{
	plot.new()
# 	plot (0,0, pch = "", axes = F, xlab = "", ylab = "")
	oldmar = par ("mar")
	par (mar = rep(1, 4))
	p = ncol (cor1)
	lim = c(-1, p + 1)
	plot.window(xlim=lim, ylim= lim, xaxs="i", yaxs="i")

	for (i in 1:p)
	{
		text (i - 0.5, p + 0.5, labels1[i], srt=90)
		text (-0.5, p - i + 0.5, labels2[i])
	}

	for (i in 2:p)
	{
		for (j in 1:(i-1))
		{
			doEllipses (cor1[c(i,j), c(i,j)], pos = c(i - 1.5,p - j - 0.5), lwd = 2)
			doEllipses (cor2[c(i,j), c(i,j)], pos = c(i - 1.5, p - j - 0.5), lwd = 2, ...)

			text (j - 0.5,p - i + 0.5, round (cor1[i,j], ndigits), adj = c(0.5,-0.1))
			text (j - 0.5,p - i + 0.5, round (cor2[i,j], ndigits), adj = c(0.5,1.1), ...)
		}
	}

	lines (c(0.5, p-0.5), c(p - 0.5, 0.5), lwd = 3)

#	lines (mean (lim) + 

	lines (lim[2] - 2 + c(-0.5, 0.3), c(-0.3, -0.3))
	lines (lim[2] - 2 + c(-0.5, 0.3), c(-0.7, -0.7), ...)
	text (lim[2] - 2 - 0.7, -0.275, method1, pos = 2)
	text (lim[2] - 2 - 0.7, -0.675, method2, pos = 2)

	par (mar = oldmar)
}

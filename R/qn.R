"qn" <-
function (x)
{
	x = as.double (x)
	.C("Qn", x, as.integer (length(x)), qn = double (1), PACKAGE = "pcaPP")$qn
}


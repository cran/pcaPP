qn = function (x)
{
	x = as.double (x)
	.C("Qn", PACKAGE="pcaPP",x, as.integer (length(x)), qn = double (1))$qn
}

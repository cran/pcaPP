EvalScaleFunction <- function (method)
{
	if (is.function (method))
		return (method) 
	return  (eval (parse (text = method)))
}

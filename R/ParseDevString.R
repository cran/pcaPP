"ParseDevString" <-
function (method)
{
	if (method == "mad")  return (0)
	if (method == "sd")  return (1)
	if (method == "Qn" | method == "qn" )  return (2)
	return (1)
}


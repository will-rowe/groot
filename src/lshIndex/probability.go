// copy of https://github.com/ekzhu/lshensemble/blob/0322dae1f4d960f6fb3f9e6e2870786b9f4239ed/probability.go
package lshIndex

import "math"

// Compute the integral of function f, lower limit a, upper limit l, and
// precision defined as the quantize step
func integral(f func(float64) float64, a, b, precision float64) float64 {
	var area float64
	for x := a; x < b; x += precision {
		area += f(x+0.5*precision) * precision
	}
	return area
}

/*
  The following are using Jaccard similarity
*/
// Probability density function for false positive
func falsePositive(l, k int) func(float64) float64 {
	return func(j float64) float64 {
		return 1.0 - math.Pow(1.0-math.Pow(j, float64(k)), float64(l))
	}
}

// Probability density function for false negative
func falseNegative(l, k int) func(float64) float64 {
	return func(j float64) float64 {
		return 1.0 - (1.0 - math.Pow(1.0-math.Pow(j, float64(k)), float64(l)))
	}
}

// Compute the cummulative probability of false negative given threshold t
func probFalseNegative(l, k int, t, precision float64) float64 {
	return integral(falseNegative(l, k), t, 1.0, precision)
}

// Compute the cummulative probability of false positive given threshold t
func probFalsePositive(l, k int, t, precision float64) float64 {
	return integral(falsePositive(l, k), 0, t, precision)
}

/*
  The following are using Jaccard containment TODO: consolidate these functions with the above
*/
// Probability density function for false positive
func falsePositiveC(x, q, l, k int) func(float64) float64 {
	return func(t float64) float64 {
		return 1.0 - math.Pow(1.0-math.Pow(t/(1.0+float64(x)/float64(q)-t), float64(k)), float64(l))
	}
}

// Probability density function for false negative
func falseNegativeC(x, q, l, k int) func(float64) float64 {
	return func(t float64) float64 {
		return 1.0 - (1.0 - math.Pow(1.0-math.Pow(t/(1.0+float64(x)/float64(q)-t), float64(k)), float64(l)))
	}
}

// Compute the cummulative probability of false negative
func probFalseNegativeC(x, q, l, k int, t, precision float64) float64 {
	fn := falseNegativeC(x, q, l, k)
	xq := float64(x) / float64(q)
	if xq >= 1.0 {
		return integral(fn, t, 1.0, precision)
	}
	if xq >= t {
		return integral(fn, t, xq, precision)
	} else {
		return 0.0
	}
}

// Compute the cummulative probability of false positive
func probFalsePositiveC(x, q, l, k int, t, precision float64) float64 {
	fp := falsePositiveC(x, q, l, k)
	xq := float64(x) / float64(q)
	if xq >= 1.0 {
		return integral(fp, 0.0, t, precision)
	}
	if xq >= t {
		return integral(fp, 0.0, t, precision)
	} else {
		return integral(fp, 0.0, xq, precision)
	}
}

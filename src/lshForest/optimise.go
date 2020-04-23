package lshforest

import "math"

//  the following funcs are taken from https://github.com/ekzhu/minhash-lsh
// optimise returns the optimal number of hash functions and the optimal number of buckets for Jaccard similarity search, as well as  the false positive and negative probabilities.
func optimise(sketchSize int, jsThresh float64) (int, int, float64, float64) {
	optimumK, optimumL := 0, 0
	fp, fn := 0.0, 0.0
	minError := math.MaxFloat64
	for l := 1; l <= sketchSize; l++ {
		for k := 1; k <= sketchSize; k++ {
			if l*k > sketchSize {
				break
			}
			currFp := probFalsePositive(l, k, jsThresh, 0.01)
			currFn := probFalseNegative(l, k, jsThresh, 0.01)
			currErr := currFn + currFp
			if minError > currErr {
				minError = currErr
				optimumK = k
				optimumL = l
				fp = currFp
				fn = currFn
			}
		}
	}
	return optimumK, optimumL, fp, fn
}

// integral of function f, lower limit a, upper limit l, and precision defined as the quantize step
func integral(f func(float64) float64, a, b, precision float64) float64 {
	var area float64
	for x := a; x < b; x += precision {
		area += f(x+0.5*precision) * precision
	}
	return area
}

// falsePositive is the probability density function for false positive
func falsePositive(l, k int) func(float64) float64 {
	return func(j float64) float64 {
		return 1.0 - math.Pow(1.0-math.Pow(j, float64(k)), float64(l))
	}
}

// falseNegative is the probability density function for false negative
func falseNegative(l, k int) func(float64) float64 {
	return func(j float64) float64 {
		return 1.0 - (1.0 - math.Pow(1.0-math.Pow(j, float64(k)), float64(l)))
	}
}

// probFalseNegative to compute the cummulative probability of false negative given threshold t
func probFalseNegative(l, k int, t, precision float64) float64 {
	return integral(falseNegative(l, k), t, 1.0, precision)
}

// probFalsePositive to compute the cummulative probability of false positive given threshold t
func probFalsePositive(l, k int, t, precision float64) float64 {
	return integral(falsePositive(l, k), 0, t, precision)
}

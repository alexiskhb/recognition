// var mapExp = {}; // memo
var dftExponent = function(k, N) {
	let x = -2 * Math.PI * (k / N);
	// mapExp[N] = mapExp[N] || {};
	// mapExp[N][k] = mapExp[N][k] || [Math.cos(x), Math.sin(x)];// [Re, Im]
	return [Math.cos(x), Math.sin(x)];
};

// complexNumber = [Re, Im]
var cAdd = function (a, b) {
	return [a[0] + b[0], a[1] + b[1]];
};

var cSub = function (a, b) {
	return [a[0] - b[0], a[1] - b[1]];
};

var cMultiply = function (a, b) {
	return [(a[0] * b[0] - a[1] * b[1]), 
	        (a[0] * b[1] + a[1] * b[0])];
};

var cMagnitude = function (c) {
	return Math.sqrt(c[0]*c[0] + c[1]*c[1]); 
};

var DFT = function(vector) {
	let X = [], N = vector.length;

	for (let k = 0; k < N; k++) {
		X[k] = [0, 0];
		for (let n = 0; n < N; n++) {
			let exp = dftExponent(k * n, N);
			let term = cMultiply([vector[n], 0], exp);
			X[k] = cAdd(X[k], term);
		}
	}

	return X;
};

var IDFT = function(vector) {
	let X = [], N = vector.length;

	for (let n = 0; n < N; n++) {
		X[n] = [0, 0];
		for (let k = 0; k < N; k++) {
			let exp = dftExponent(-n * k, N);
			let term = cMultiply(vector[k], exp);
			X[n] = cAdd(X[n], term);
		}
		X[n][0] /= N;
		X[n][1] /= N;
	}

	return X;
};


$(document).ready(function() {
	var precision = 4;
	var pfc = 96/2.54; // approx pixels per centimiter
	var displayedSignals = [];

	function splitNumsByWhitespace(s) {
		return s.split(/[,\s]+/).filter( function(e) { return e.trim().length > 0; } ).map(Number);
	}

	$('#drawSignal').on('click', function() {
		displaySignal(false);
	});
	
	$('#addSignal').on('click', function() {
		displaySignal(true);
	});

	function displaySignal(add) {
		function screenLinspace(width) {
			let sampleRate = Number($('#signalParam_sampleRate').val());
			sampleRate /= pfc;
			let delta = 1/pfc/sampleRate;
			let r = [];
			for (let x = 0; x < width/pfc && x >= 0 && delta != 0; x += delta) {
				r.push(x);
			}
			return r;
		} 

		let signals = [
			function (n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma) { //1
				let x = screenLinspace(n);
				let y = x.map(function(t){return {x:t, y:0};});
				for (let i = 0; i < y.length; i++) {
					if (x[i] >= n0) {
						y[i].y = 1;
						break;
					}
				}
				return y;
			},
			function (n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma) { //2
				let x = screenLinspace(n);
				let y = x.map(function(t){return {x:t, y:0};});
				for (let i = 0; i < y.length; i++) {
					if (x[i] >= n0) {
						y[i].y = 1;
					}
				}
				return y;
			},
			function (n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma) { //3
				let x = screenLinspace(n);
				let y = x.map(function(t){return {x:t, y:Math.pow(a, t)};});
				return y;
			},
			function (n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma) { //4
				let x = screenLinspace(n);
				let y = x.map(function(t){return {x:t, y:a*Math.sin(2*Math.PI*omega*t + phi)};});
				return y;
			},
			function (n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma) { //5
				L = Math.trunc(L);
				let x = screenLinspace(n);
				let y = x.map(function(t) {
					let tr = Math.trunc(t);
					return {x:t, y:tr%L < L/2 ? 1 : -1};
				});
				return y;
			},
			function (n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma) { //6 
				L = Math.trunc(L);
				let x = screenLinspace(n);
				let y = x.map(function(t) {
					let tr = Math.trunc(t);
					return {x:t, y:(tr%L)/L};
				});
				return y;
			},
			function (n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma) { //7 
				let x = screenLinspace(n);
				let y = x.map(function(t) {
					return {x:t, y:a*Math.pow(Math.E, -t/tau)*Math.cos(2*Math.PI*omega*t + phi)};
				});
				return y;
			},
			function (n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma) { //8
				let x = screenLinspace(n);
				let y = x.map(function(t) {
					return {x:t, y:a*Math.cos(u*t)*Math.cos(2*Math.PI*omega*t + phi)};
				});
				return y;
			},
			function (n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma) { //9
				let x = screenLinspace(n);
				let y = x.map(function(t) {
					return {x:t, y:a*(1 + m*Math.cos(u*t))*Math.cos(2*Math.PI*omega*t + phi)};
				});
				return y;
			},
			function (n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma) { //10
				let x = screenLinspace(n);
				let y = x.map(function(t) {
					return {x:t, y: A + (B - A)*Math.random()};
				});
				return y;
			},
			function (n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma) { //11
				let x = screenLinspace(n);
				let y = x.map(function(t) {
					return {x:t, y: jStat.normal.sample(a, sigma)};
				});
				return y;
			},
			function (n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma) { //12
				while (q.length < p.length) {
					q.push(0);
				}
				while (p.length < q.length) {
					p.push(0);
				}
				let len = p.length;
				let x = screenLinspace(n);
				let rds = [];
				let y = [];
				for (let i = 0; i < x.length; i++) {
					let t = x[i];
					let r = jStat.normal.sample(0, 1);
					rds.push(r);
					for (let j = 1; j <= len && i - j >= 0; j++) {
						r += q[j - 1]*rds[i - j] - p[j - 1]*y[i - j].y;
					}
					y.push({x:t, y:r});
				};
				return y;
			}
		];

		let n = Number($('#signalParam_N').val());
		let n0 = Number($('#signalParam_n0').val());
		let a = Number($('#signalParam_a').val());
		let omega = Number($('#signalParam_omega').val());
		let phi = Number($('#signalParam_phi').val());
		let L = Number($('#signalParam_L').val());
		let tau = Number($('#signalParam_tau').val());
		let u = Number($('#signalParam_u').val());
		let m = Number($('#signalParam_m').val());
		let A = Number($('#signalParam_A').val());
		let B = Number($('#signalParam_B').val());
		let sigma = Number($('#signalParam_sigma').val());
		let showLine = $('#showLine').is(":checked");
		let p = splitNumsByWhitespace($('#signalParam_p').val());
		let q = splitNumsByWhitespace($('#signalParam_q').val());
		let color = $('#signalParam_color').val();

		let selectId = Number($("select#signal").val());
		let selectText = $('#signal>option:selected').text();
		let points = signals[selectId](n, n0, a, omega, phi, L, p, q, tau, u, m, A, B, sigma);
		let axisY = {};
		if ($('#signalParam_yMin').val() != "" && !isNaN(Number($('#signalParam_yMin').val()))) {
			axisY.minimum = Number($('#signalParam_yMin').val());
		}
		if ($('#signalParam_yMax').val() != "" && !isNaN(Number($('#signalParam_yMax').val()))) {
			axisY.maximum = Number($('#signalParam_yMax').val());
		}
		let newData = {
			color: color,
			type: showLine ? "line" : "column",
			dataPoints: points,
			markerSize: 0
		};
		if (add) {
			displayedSignals.push(newData);
		} else {
			displayedSignals = [newData];
		}
		(new CanvasJS.Chart("signalPlotContainer", {
			zoomEnabled: true,
			zoomType: "x",
			title: {text: selectText},
			data: displayedSignals,
			axisX:{minimum: 0},
			axisY:axisY
		})).render();
	}


	$('#generateNthNormal').on('click', function(event) {
		function getCovariationMatrix(vars) {
			let size = vars.length;
			function validate(mat) {
				function minor(mat, n) {
					if (n > mat.length) {
						return mat;
					}
					let r = [];
					for (let i = 0; i < n; i++) {
						r.push([]);
						for (let j = 0; j < n; j++) {
							r[i].push(mat[i][j]);
						}
					}
					return r;
				}

				if (!jStat(mat).symmetric()) {
					alert('Матрица не симметричная');
					return false;
				}
				for (let i = 2; i < mat.length; i++) {
					let m = minor(mat, i);
					if (jStat.det(m) <= 0) {
						alert('Матрица не положительно определена');
						return false;
					}
				}
				return true;
			}

			let corr = splitNumsByWhitespace($('#covariationMatrix').val());
			if (size*size != corr.length && size*(size + 1)/2 != corr.length) {
				alert('Размеры вектора дисперсий и матрицы ковариации несовместимы');
				return null;
			}
			let mat = [];
			for (let i = 0; i < size; i++) {
				mat.push([]);
			}
			if (size*size == corr.length) {
				for (let i = 0; i < size; i++) {
					for (let j = 0; j < size; j++) {
						mat[i].push(corr[i*size + j]);
					}
				}
			} else {
				let k = 0;
				for (let i = 0; i < size; i++) {
					for (let j = 0; j < size; j++) {
						mat[i].push(j <= i ? corr[k++] : 0);
					}
				}
				for (let i = 0; i < size; i++) {
					for (let j = i + 1; j < size; j++) {
						mat[i][j] = mat[j][i];
					}
				}
			}
			for (let i = 0; i < size; i++) {
				for (let j = 0; j < size; j++) {
					mat[i][j] *= Math.sqrt(vars[i])*Math.sqrt(vars[j]);
				}
			}
			if (!validate(mat)) {
				return null;
			}
			return mat;
		}

		function cholesky(corr) {
			let B = [];
			for (let i = 0; i < corr.length; i++) {
				B.push([]);
				for (let j = 0; j < corr.length; j++) {
					B[i].push(0);
				}
			}
			for (let i = 0; i < corr.length; i++) {
				for (let j = 0; j <= i; j++) {
					let nom = corr[i][j];
					for (let k = 0; k < j - 1; k++) {
						nom -= B[i][k]*B[j][k];
					}
					let denom = corr[i][j];
					for (let k = 0; k < j - 1; k++) {
						denom -= B[j][k]*B[j][k];
					}
					B[i][j] = nom/Math.sqrt(denom);
					if (isNaN(B[i][j])) {
						continue;
					}
				}
			}
			return B;
		}

		function getSample(N, means, B) {
			let X = [];
			for (let i = 0; i < N; i++) {
				X.push([jStat.normal.sample(0, 1)]);
			}
			let prod = jStat.transpose(jStat.multiply(B, X));
			let r = jStat.add(prod, means);
			return r;
		}

		function getSamples(nOfVectors, N, means, B) {
			let samples = [];
			for (let i = 0; i < nOfVectors; i++) {
				samples.push(getSample(N, means, B))
			}
			return samples;
		}

		function f2D(x, y, sigmaX, sigmaY, meanX, meanY, corr) {
			let expPow = -1/(1-Math.pow(corr, 2))/2*(Math.pow((x-meanX)/sigmaX, 2) - 2*corr*(x-meanX)*(y-meanY)/sigmaX/sigmaY + Math.pow((y-meanY)/sigmaY, 2));
			return Math.pow(Math.E, expPow)/Math.PI/sigmaX/sigmaY/Math.sqrt(1-Math.pow(corr, 2))/2;
		}

		function getNormal2DHeatmap(x0, y0, d, sigmaX, sigmaY, meanX, meanY, corr) {
			let delta = d/($('#normalPlotContainer').width()/2);
			let xs = [], ys = [], zs = [];
			for (let t = 0; t < d; t += delta) {
				xs.push(x0 + t);
				ys.push(y0 + t);
			}
			for (let x = x0; x < x0 + d; x += delta) {
				zs.push([]);
				for (let y = y0; y < y0 + d; y += delta) {
					zs[zs.length - 1].push(f2D(x, y, sigmaX, sigmaY, meanX, meanY, corr));
				}
			}
			return {
				x: xs, 
				y: ys, 
				z: jStat.transpose(zs),
				type: 'heatmap',
				// colorscale: 'Greys'
				colorscale: 'YIGnBu'
			};
		}

		let means = splitNumsByWhitespace($('#meanVector').val());
		if (!means) {
			return null;
		}
		let vars = splitNumsByWhitespace($('#varVector').val());
		if (!vars) {
			return null;
		}
		if (vars.length != means.length) {
			alert('Вектора дисперсий и средних несовместимы');
			return null;
		}
		let corr = getCovariationMatrix(vars);
		if (!corr) {
			return null;
		}
		let radioX = Number($('input[name=radioX]:checked').val());
		let radioY = Number($('input[name=radioY]:checked').val());
		let B = cholesky(corr);
		let B2 = jStat.cholesky(corr); 
		let samples = getSamples(Number($('#nOfVectors').val()), means.length, means, B2);
		let samplesTr = jStat.transpose(samples);
		let sampleMeans = samplesTr.map(jStat.mean).map(function(k){return Number(k.toFixed(precision))});
		let sampleVars = samplesTr.map(jStat.variance).map(function(k){return Number(k.toFixed(precision))});
		let sampleMeansDelta = jStat.subtract(means, sampleMeans).map(function(k){return Number(k.toFixed(precision))});
		let sampleVarsDelta = jStat.subtract(vars, sampleVars).map(function(k){return Number(k.toFixed(precision))});
		let sampleCorrcoeff = jStat.corrcoeff(samplesTr[radioX], samplesTr[radioY]);
		$('#sampleMean').html(JSON.stringify(sampleMeans));
		$('#sampleVar').html(JSON.stringify(sampleVars));
		$('#sampleMeanDelta').html(JSON.stringify(sampleVarsDelta));
		$('#sampleVarDelta').html(JSON.stringify(sampleMeansDelta));
		$('#sampleCorrcoeff').html(sampleCorrcoeff.toFixed(precision));
		$('#sampleCorrcoeffDelta').html((corr[radioX][radioY]/Math.sqrt(vars[radioX]*vars[radioY]) - sampleCorrcoeff).toFixed(precision));
		let minX = jStat.min(samplesTr[radioX]);
		let maxX = jStat.max(samplesTr[radioX]);
		let minY = jStat.min(samplesTr[radioY]);
		let maxY = jStat.max(samplesTr[radioY]);
		let d = Math.max(maxX - minX, maxY - minY);
		let heatmap = getNormal2DHeatmap(means[radioX] - vars[radioX], means[radioY] - vars[radioY], 2*Math.max(vars[radioX], vars[radioY]), 
		                                 Math.sqrt(vars[radioX]), Math.sqrt(vars[radioY]), 
		                                 means[radioX], means[radioY], 
		                                 corr[radioX][radioY]/Math.sqrt(vars[radioX])/Math.sqrt(vars[radioY]));

		let estimatedHeatmap = getNormal2DHeatmap(means[radioX] - vars[radioX], means[radioY] - vars[radioY], 2*Math.max(vars[radioX], vars[radioY]), 
		                                 Math.sqrt(sampleVars[radioX]), Math.sqrt(sampleVars[radioY]), 
		                                 sampleMeans[radioX], sampleMeans[radioY], 
		                                 sampleCorrcoeff);
		let data = [heatmap, {
			x: samplesTr[radioX],
			y: samplesTr[radioY],
			mode: 'markers',
			type: 'scatter',
			marker: {size:5, color: 'orange'}
		}];
		let layout = {
			xaxis: {range: [means[radioX] - vars[radioX], means[radioX] - vars[radioX] + 2*Math.max(vars[radioX], vars[radioY])]},
			yaxis: {range: [means[radioY] - vars[radioY], means[radioY] - vars[radioY] + 2*Math.max(vars[radioX], vars[radioY])]}
		};
		Plotly.newPlot('normalPlotContainer', data, layout);
		data[0] = estimatedHeatmap;
		Plotly.newPlot('estimatedNormalPlotContainer', data, layout);
	});




	var seriesN = 50;

	function distrSeries() {
		let dataPoints = [];
		for (let i = 1; i <= seriesN; i++) {
			dataPoints.push({x:i, y:Math.pow(2/3, i)/2});
		}
		return dataPoints;
	}

	function init() {
		displayedSignals = [];
		updateSample();
		recalculateDPos();
		recalculateContPos(); 
		redrawCont();
		updateSample();
		updateContSample();
	}

	var ds = distrSeries();

	let chartDSPolygon = new CanvasJS.Chart("chartContainerDSPolygon",
	{
		zoomEnabled: true,
        zoomType: "x",
		title: {
			text: "Многоугольник распределения"
		},
		data: [
		{
			type: "line",
			dataPoints: ds
		}],
		axisY:{
			minimum: 0,
		},
		axisX:{
			minimum: 0,
		}
	});
	chartDSPolygon.render();

	function getFunc(ds) {
		let dataPoints = [{x:0, y:0}];
		for (let i = 0; i < Math.min(ds.length, ds.length); i++) {
			dataPoints.push({x:i + 1, y: dataPoints[i].y + ds[i].y});
		}
		return dataPoints;
	}

	let f = getFunc(ds);
	let chartF = new CanvasJS.Chart("chartContainerF",
	{
		zoomEnabled: true,
        zoomType: "x",
		title: {
			text: "Функция распределения"
		},
		data: [
		{
			color: "blue",
			type: "column",
			dataPoints: f
		}],
		axisY:{
			minimum: 0,
		},
		axisX:{
			minimum: 0,
		}
	});
	chartF.render();

	function myLog(x) {
		return Math.max(Math.log(x) / Math.log(2.0/3), 0.1);
	}

	function getSampleDataPoints(N) {
		let dataPoints = [];
		for (let i = 0; i < N; i++) {
			let r = 1 - Math.random();
			let l = Math.ceil(myLog(r));
			if (l >= dataPoints.length) {
				while (l >= dataPoints.length) {
					dataPoints.push({x:dataPoints.length, y:0});
				}
			}
			dataPoints[l].y += 1;
		}
		return dataPoints;
	}

	function getStats(data, N) {
		let mean = 0;
		for (let i = 0; i < data.length; i++) {
			mean += data[i].x*data[i].y;
		}
		
		let variance = 0;
		for (let i = 0; i < data.length; i++) {
			variance += Math.pow(data[i].x - mean, 2);
		}
		variance /= data.length;

		let m2 = 0, m3 = 0, m4 = 0;
		for (let i = 0; i < data.length; i++) {
			m2 += Math.pow(data[i].x - mean, 2);
			m3 += Math.pow(data[i].x - mean, 3);
			m4 += Math.pow(data[i].x - mean, 4);
		}
		m2 = m2/(N - 1);
		m3 = m3/(N - 1)*N/(N - 2);
		let asym = m3/Math.pow(Math.sqrt(m3), 3);
		let kurtosis = N*(N+1)/(N-1)/(N-2)*m4/(N-3)/Math.pow(variance, 2) - 3*(N-1)*(N-1)/(N-2)/(N-3);

		return {
			mean: mean,
			variance: variance,
			asym: asym,
			kurtosis: kurtosis
		};
	}

	function updateSample() {
		let N = $('#N').val();
		let sample = getSampleDataPoints(N);
		
		for (let i = 0 ; i < sample.length; i++) {
			sample[i].y /= N;
		}

		let chart = new CanvasJS.Chart("chartContainer",
		{
			zoomEnabled: true,
			zoomType: "x",
			title: {
				text: "Выборка"
			},
			data: [
			{
				type: "column",
				dataPoints: sample
			}
			],
			axisY:{
				minimum: 0,
			},
			axisX:{
				minimum: 0,
			}
		});
		chart.render();
	
		let st = getStats(sample, N);
		let variance = 6;
		let mean = 3;
		$('#resultsTbl').empty();
		$('#resultsTbl').
			append($('<tr>').append($('<td>').html('.')).append($('<td>').html('Вычисл.')).append($('<td>').html('Аналит.')).append($('<td>').html('Ошибка'))).
			append($('<tr>').append($('<td>').html('Среднее')).append($('<td>').html(st.mean.toFixed(precision))).append($('<td>').html(mean.toFixed(precision))).append($('<td>').html((st.mean-mean).toFixed(precision)))).
			append($('<tr>').append($('<td>').html('Дисперсия')).append($('<td>').html(Math.sqrt(st.variance).toFixed(precision))).append($('<td>').html(variance.toFixed(precision))).append($('<td>').html((Math.sqrt(st.variance)-variance).toFixed(precision))));

		$('#results').html(
			'<br>Коэф. асимметрии ' + st.asym.toFixed(precision) + 
			'<br>Коэф. эксцесса ' + st.kurtosis.toFixed(precision) +'<br>');

		let stF = [{x:0, y:0}];
		for (let i = 1; i < sample.length; i++) {
			stF.push({x:i, y: stF[i-1].y + sample[i].y});
		}
		let chartStF = new CanvasJS.Chart("chartContainerStF",
		{
			zoomEnabled: true,
			zoomType: "x",
			title: {
				text: "Статистическая функция распределения"
			},
			data: [
			{
				type: "column",
				dataPoints: stF
			}
			],
			axisY:{
				minimum: 0,
			},
			axisX:{
				minimum: 0,
			}
		});
		chartStF.render();
	}

	$('#N').on('change', function(event) {
		updateSample();
	});

	$('#refresh').on('click', function(event) {
		updateSample();
	});

	function recalculateDPos() {
		let A = Number($('#A').val());
		let B = Number($('#B').val());
		let result = 0;
		if (A <= B) {
			for (let n = A; n <= B; n++) {
				result += Math.pow(2/3, n);
			}
		}
		result /= 2;
		result = Math.max(Math.min(result, 1), 0);
		$('#d_pos').html(result.toFixed(6));
	}

	$('#A').on('change', function(event) { 
		recalculateDPos(); 
	});

	$('#B').on('change', function(event) { 
		recalculateDPos(); 
	});

	function F(x, sigma) {
		return 1 - Math.pow(Math.E, -x*x/2/sigma/sigma);
	}

	function p(x, sigma) {
		return x*Math.pow(Math.E, -x*x/2/sigma/sigma)/sigma/sigma;
	}

	function recalculateContPos() {
		let A = Number($('#cA').val());
		let B = Number($('#cB').val());
		let sigma = Number($('#sigma').val());
		let result = F(B, sigma) - F(A, sigma);
		$('#cont_pos').html(result.toFixed(precision));
	}

	function redrawCont() {
		let sigma = Number($('#sigma').val());
		let L = Math.sqrt(2*sigma*sigma*Math.log(1000/999))-1;
		let R = Math.sqrt(2*sigma*sigma*Math.log(1000))+1;
		let delta = (R - L)/3000;
		let f = [];
		let contMode = sigma;
		let contMean = sigma*Math.sqrt(Math.PI/2);
		let contMedian = Math.sqrt(2*sigma*sigma*Math.log(2));

		for (let x = L; x <= R; x += delta) {
			f.push({x:x, y:F(x, sigma)});
		}

		let chartF = new CanvasJS.Chart("chartContainerContF",
		{
			zoomEnabled: true,
			zoomType: "x",
			title: {
				text: "F(x)"
			},
			data: [
			{
				type: "line",
				dataPoints: f
			}
			],
			axisY:{
				minimum: 0,
				maximum: 1.1
			},
			axisX:{
				minimum: 0
			}
		});
		chartF.render();

		let f2 = [];
		for (let x = L; x <= R; x += delta) {
			f2.push({x:x, y:p(x, sigma)});
		}
		
		$('#contMode').html(contMode.toFixed(precision));
		$('#contMean').html(contMean.toFixed(precision));
		$('#contMedian').html(contMedian.toFixed(precision));

		let chartP = new CanvasJS.Chart("chartContainerContP",
		{
			zoomEnabled: true,
			zoomType: "x",
			title: {
				text: "f(x)"
			},
			data: [
			{
				type: "line",
				dataPoints: f2
			}
			],
			axisY:{
				minimum: 0,
				maximum: 1.1/sigma/Math.sqrt(Math.E)
			},
			axisX:{
				minimum: 0,
			    stripLines:[
				{
					value: contMode,
					color:"blue",
					label : "m",
					labelFontColor: "blue",
					labelPlacement: "outside"
				},
				{
					value: contMedian,
					color:"blue",
					label : "d",
					labelFontColor: "blue",
					labelPlacement: "outside"
				},
				{
					value: contMean,
					color:"blue",
					label : "n",
					labelFontColor: "blue",
					labelPlacement: "outside"
				}
				]
			}
		});
		chartP.render();
	}

	$('#cA').on('change', function(event) { 
		recalculateContPos(); 
		redrawCont();
	});

	$('#cB').on('change', function(event) { 
		recalculateContPos(); 
		redrawCont();
	});

	$('#sigma').on('change', function(event) { 
		recalculateContPos(); 
		redrawCont();
		updateContSample();
	});

	function randFx(sigma) {
		return Math.sqrt(-2*sigma*sigma*Math.log(1 - Math.random()));
	}

	function getContSampleDataPoints() {
		let sigma = Number($('#sigma').val());
		let columns = Number($('#contNBin').val());
		let L = Math.sqrt(2*sigma*sigma*Math.log(1000/999));
		let R = Math.sqrt(2*sigma*sigma*Math.log(1000));
		let contN = Number($('#contN').val());
		let delta = (R - L)/columns;
		let dataPoints = [];
		for (let i = 0; i < contN; i++) {
			let x = randFx(sigma);
			x = Math.trunc((Math.max(x - L, 0))/delta) + 1;
			if (x >= dataPoints.length) {
				while (x >= dataPoints.length) {
					dataPoints.push({x:dataPoints.length, y:0});
				}
			}
			dataPoints[x].y += 1;
		}
		for (let i = 0; i < dataPoints.length; i++) {
			dataPoints[i].x *= delta;
		}
		return dataPoints;
	}

	function updateContSample() {
		let N = $('#contN').val();
		let sigma = Number($('#sigma').val());
		let sample = getContSampleDataPoints(N);
		for (let i = 0 ; i < sample.length; i++) {
			sample[i].y /= N;
		}
		let st = getStats(sample, N);
		let chart = new CanvasJS.Chart("chartContainerContSample",
		{
			zoomEnabled: true,
			zoomType: "x",
			title: {
				text: "Выборка"
			},
			data: [
			{
				color: "blue",
				type: "column",
				dataPoints: sample
			}
			],
			axisY:{
				minimum: 0,
			},
			axisX:{
				minimum: 0,
			}
		});
		chart.render();

		let contMean = sigma*Math.sqrt(Math.PI/2);
		let contVar = sigma*sigma*(2-Math.PI/2);
		$('#contResultsTbl').empty();
		$('#contResultsTbl').
			append($('<tr>').append($('<td>').html('.')).append($('<td>').html('Вычисл.')).append($('<td>').html('Аналит.')).append($('<td>').html('Ошибка'))).
			append($('<tr>').append($('<td>').html('Среднее')).append($('<td>').html(st.mean.toFixed(precision))).append($('<td>').html(contMean.toFixed(precision))).append($('<td>').html((st.mean-contMean).toFixed(precision)))).
			append($('<tr>').append($('<td>').html('Дисперсия')).append($('<td>').html(Math.sqrt(st.variance).toFixed(precision))).append($('<td>').html(contVar.toFixed(precision))).append($('<td>').html((Math.sqrt(st.variance)-contVar).toFixed(precision))));

		$('#contResults').html(
			'<br>Коэф. асимметрии ' + st.asym.toFixed(precision) +
			'<br>Коэф. эксцесса ' + st.kurtosis.toFixed(precision) +'<br>');

		let stF = [{x:0, y:0}];
		for (let i = 1; i < sample.length; i++) {
			stF.push({x:sample[i].x, y: stF[i-1].y + sample[i].y});
		}

		let chartStF = new CanvasJS.Chart("chartContainerContStF",
		{
			zoomEnabled: true,
			zoomType: "x",
			title: {
				text: "Статистическая функция распределения"
			},
			data: [
			{
				color: "blue",
				type: "column",
				dataPoints: stF
			}
			],
			axisY:{
				minimum: 0,
			},
			axisX:{
				minimum: 0,
			}
		});
		chartStF.render();
	}

	$('#contRefresh').on('click', function(event) {
		updateContSample();
	});
	
	$('#contN').on('change', function(event) {
		updateContSample();
	});
	
	$('#contNBin').on('change', function(event) {
		updateContSample();
	});

	init();
});

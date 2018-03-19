$(document).ready(function() {
	var seriesN = 50;
	var precision = 5;
	updateSample();

	function distrSeries() {
		let dataPoints = [];
		for (let i = 1; i <= seriesN; i++) {
			dataPoints.push({x:i, y:Math.pow(2/3, i)/2});
		}
		return dataPoints;
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
			variance += (data[i].x - mean)*(data[i].x - mean);
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
		$('#results').html(
			'Среднее: ' + st.mean.toFixed(precision) + 
			'<br>Дисперсия ' + st.variance.toFixed(precision) +
			'<br>Срквадр ' + Math.sqrt(st.variance).toFixed(precision) + 
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
	$('#A').on('change', function(event) { recalculateDPos(); });
	$('#B').on('change', function(event) { recalculateDPos(); });
	recalculateDPos();

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
		$('#cont_pos').html(result.toFixed(4));
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
				minimum: 0,
				stripLines:[
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
		chartF.render();

		let f2 = [];
		for (let x = L; x <= R; x += delta) {
			f2.push({x:x, y:p(x, sigma)});
		}
		
		$('#contMode').html(contMode.toFixed(4));
		$('#contMean').html(contMean.toFixed(4));
		$('#contMedian').html(contMedian.toFixed(4));

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
					label : "mode",
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
	recalculateContPos(); 
	redrawCont();

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
		$('#contResults').html(
			'Среднее: ' + st.mean.toFixed(precision) + 
			'<br>Дисперсия ' + st.variance.toFixed(precision) +
			'<br>Срквадр ' + Math.sqrt(st.variance).toFixed(precision) + 
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

	$('#contN').on('change', function(event) {
		updateContSample();
	});
	$('#contNBin').on('change', function(event) {
		updateContSample();
	});

	updateSample();
	updateContSample();
});

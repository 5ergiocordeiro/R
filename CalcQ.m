function q = CalcQ(V, iout, jout, jd, epsilon)
% Calcula a carga por meio da lei de Gauss a partir do potencial 'V'
	% Parâmetros para o cálculo
	% ... intervalos a considerar
	interval = [
		1, 1, jout + 1, jout + 1;
		1, 1, jout, jout;
		2, iout, jout + 1, jout + 1;
		2, iout, jout, jout;
		iout + 1, iout + 1, jd + 1, jout;
		iout, iout, jd + 1, jout;
		iout + 1, iout + 1, 2, jd - 1;
		iout, iout, 2, jd - 1;
		iout + 1, iout + 1, 1, jd - 1;
		iout, iout, 1, jd - 1;
		iout + 1, iout + 1, jd, jd;
		iout, iout, jd, jd;
		iout + 1, iout + 1, jd, jd;
		iout, iout, jd, jd;		
		];
	% ... pesos e permissividades em cada intervalo
	mult = epsilon([1 1 1 1 1 1 2 2 2 2 1 1 2 2]) .* [-0.5 0.5 -1 1 -1 1 -1 1 -0.5 0.5 -0.5 0.5 -0.5 0.5];
	% Cálculo pela lei de Gauss
	q = 0; 
	for k = 1:14
		q = q + sum(V(interval(k,1):interval(k,2), interval(k,3):interval(k,4))) * mult(k);
	end
end

function geramattri(n)
% Gera dois sistemas triangulares 'n' x 'n', um superior e outro inferior, a partir de valores aleatórios na faixa ]0 - 1[, grava-os em disco, resolve-os e calcula a norma 2 do resultado.
	A = B = zeros(n, n);
	for m = 1:n
		nvals = 1 + n - m;
		A(m, m:n) = rand(1, nvals);
		B(nvals, 1:nvals) = rand(1, nvals);
	end
	bA =  rand(1, n)';
	bB =  rand(1, n)';
	TS = [A bA];
	TI = [B bB];	
	save('TS', 'TS');
	save('TI', 'TI');
	disp(sprintf("Sistemas gerados e gravados."));
	INVA = A ^ (-1);
	INVB = B ^ (-1);
	xA = INVA * bA;
	xB = INVB * bB;
	disp(sprintf("Norma 2 dos resultados = %f e %f", norm(xA, 2), norm(xB, 2)));
end
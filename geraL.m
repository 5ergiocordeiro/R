function geraL(n)
%Gera uma matriz de estados para o posicionador de antena de satélite e calcula seus autovalores.
	C = 7.72e-3;
	a = 1.14e-2;
	for i = 1:5
		K = 10^(i-1);
		for j = 1:5
			alpha = 10^(-j);
			L = [0, 1, 0; 0, 0, 1; - K * C * alpha, - K * C, - a];
			v = eig(L);
			disp(sprintf("K = %f, alpha = %f => v = ", K, alpha));
			disp(v);
		end
	end
end
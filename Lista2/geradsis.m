function geradsis(n)
% Gera um sistema 'n' x 'n' diagonalmente dominante a partir de valores aleatórios, grava-a em disco, resolve-o e calcula a norma 2 do resultado.
	A = rand(n);
	for i = 1:n
		A(i,i) = 10 + A(i,i) * 10;
	end
	b = rand(n,1);
	S = [A, b];
	save('D', 'S');
	disp(sprintf("Sistema gerado e gravado."));
	INVA = A^(-1);
	X = INVA * b;
	disp(sprintf("Norma 2 da solução = %f", norm(X, 2)));
end
function geramat(n)
% Gera duas matrizes 'n' x 'n' a partir de valores aleatórios na faixa ]0 - 1[, grava-as em disco, calcula o produto e a norma 2 do resultado.
	A = rand(n);
	save('MatA', 'A');
	B = rand(n);
	save('MatB', 'B');
	disp(sprintf("Matrizes geradas e gravadas."));
	C = A * B;
	disp(sprintf("Norma 2 do produto = %f", norm(C, 2)));
end
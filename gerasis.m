function gerasis(n)
%Gera um sistema a partir de valores aleatórios na faixa ]0 - 1[, grava-o em disco, resolve-o e calcula o determinante e a norma 2 do resultado.
	A = rand(n);
	bA =  rand(1, n)';
	S = [A bA];
	save('S', 'S');
	disp(sprintf("Sistema gerado e gravado."));
	disp(sprintf("Determinante = %f", det(A)));
	INVA = A ^ (-1);
	xA = INVA * bA;
	disp(sprintf("Norma 2 do resultado = %f", norm(xA, 2)));
end
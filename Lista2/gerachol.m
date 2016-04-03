function gerachol(n)
%Gera um sistema a partir de valores aleatórios na faixa ]0 - 1[, grava-o em disco, resolve-o e calcula o determinante e a norma 2 do resultado.
%A matriz gerada é simétrica e definida positiva.
	% Gera matriz não singular
	test = 0;
	while test == 0
		A = rand(n);
		test = det(A);
	end	
	% Gera matriz simétrica definida positiva
	S = A * A';
	% Gera o sistema
	b = rand(n,1);
	C = [S, b];
	save('C', 'C');
	disp(sprintf("Sistema gerado e gravado."));
	disp(sprintf("Determinante = %f", test^2));
	% Resolve o sistema
	INVS = S ^ (-1);
	x = INVS * b;
	disp(sprintf("Norma 2 do resultado = %f", norm(x, 2)));
end
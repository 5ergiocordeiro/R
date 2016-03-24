function gerachol(n)
%Gera um sistema a partir de valores aleat�rios na faixa ]0 - 1[, grava-o em disco, resolve-o e calcula o determinante e a norma 2 do resultado.
%A matriz gerada � sim�trica e definida positiva.
	% Gera matriz n�o singular
	test = 0;
	while test == 0
		A = rand(n);
		test = det(A);
	end	
	% Gera matriz sim�trica definida positiva
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
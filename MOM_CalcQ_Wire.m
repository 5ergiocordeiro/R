% Calcula a distribuição de carga 'rho' em um condutor linear de raio 'aa' e comprimento 'l', mantido em um potencial 'Vo', pelo método dos momentos, com 'n' subdivisões da região.
function [q, rho, y] = MOM_CalcQ_Wire(V0, er, aa, l, n)
	% Constantes
	epsilon0 = 8.8541e-12;
	h = l / n;
	doisvezesloghsobreaa = 2 * log(h / aa);
	% Calcula a matrix de coeficientes 'A'
	y = h * ((1:n) - 0.5);
	A = zeros(n);
	for i = 1:n
		for j = 1:n
			if (i ~= j)
				A(i,j) = h / abs(y(i) - y(j));
			else
				A(i,j) = doisvezesloghsobreaa;
			end
		end
	end
	% Calcula a matrix de momentos 'b'
	b = 4.0 * pi * epsilon0 * er * V0 * ones(n, 1);
	% Calcula a distribuição de carga
	rho = A \ b;
	% Calcula a carga total
	q = h * sum(rho);
end
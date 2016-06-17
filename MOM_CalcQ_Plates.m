% Calcula a distribuição de carga 'rho' nas placas de um capacitor plano de arestas 'aa' e 'bb', separadas por uma distância 'd', submetidas a uma diferença de potencial 'Vo', pelo método dos momentos, com 'n' subdivisões em cada dimensão. O vcalor de 'n' precisa ser um quadrado perfeito.
function [C, q, rho, x, y] = MOM_CalcQ_Plates(V0, er, aa, bb, d, n)
	% Testa as entradas
	m = sqrt(n);
	if (m > floor(m))
		fprintf('O valor de "n" deve ser um quadrado perfeito! \n');
		C = NaN; q = NaN; rho = NaN; x = NaN; y = NaN;
		return;
	end
	% Constantes
	epsilon0 = 8.8541e-12;
	h = [aa bb] / m;
	nt = 2 * n;
	dl2 = h(1) ^2;
	V0sobre2 = V0 / 2;
	quatropiepsilon = 4 * pi * epsilon0 * er;
	mandrake = h(1) * log(1 + sqrt(2)) / (pi  * epsilon0 * er);
	% Calcula a matrix de coeficientes 'A'
	A = zeros(nt);
	p = 0;
	x = zeros(nt, 1);
	y = x; z = x;
	for k = 1:2
		for i = 1:m
			for j = 1:m
				p = p + 1;
				x(p) = h(1) * (i - 0.5);
				y(p) = h(2) * (j - 0.5);
			end
		end
	end
	z((n+1):nt) = d;
	for i = 1:nt
		for j = 1:nt
			if (i == j)
				A(i,j) = mandrake;
			else
				r = sqrt(sum([x(i) - x(j), y(i) - y(j), z(i) - z(j)] .^ 2));
				A(i,j) = dl2 / (quatropiepsilon * r);
			end
		end
	end	
	% Calcula a matrix de momentos 'b'
	b = ones(nt, 1);
	b((n+1):nt) = -1;
	b = b * V0sobre2;
	% Calcula a distribuição de carga
	rho = A \ b;
	% Calcula a carga total
	q = dl2 * sum(rho(1:n));
	% Calcula a capacitância
	C = abs(q/V0);
end
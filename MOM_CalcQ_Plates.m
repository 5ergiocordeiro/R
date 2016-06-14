% Calcula a distribui��o de carga 'rho' nas placas de um capacitor plano de arestas 'aa' e 'bb', separadas por uma dist�ncia 'd', submetidas a uma diferen�a de potencial 'Vo', pelo m�todo dos momentos, com 'n' subdivis�es da regi�o.
function [C, q, rho, x, y] = MOM_CalcQ_Plates(V0, er, aa, bb, d, n)
	% Constantes
	epsilon0 = 8.8541e-12;
	m = sqrt(n);
	h = [aa bb] / m;
	nt = 2 * n;
	dl2 = h(1) ^2;
	V0sobre2 = V0 / 2;
	quatropiepsilon = 4 * pi * epsilon0 * er;
	mandrake = h(1) * log(1 + sqrt(2)) / (pi  * epsilon0 * er);
	% Calcula a matrix de coeficientes 'A'
	A = zeros(nt);
	p = 0;
	x = zeros(2 * n, 1);
	y = x;
	for k = 1:2
		for i = 1:m
			for j = 1:m
				p = p + 1;
				x(p) = h(1) * (i - 0.5);
				y(p) = h(2) * (j - 0.5);
			end
		end
	end
	z = zeros(nt, 1);
	z(n+1:nt) = d;
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
	b(n+1:nt) = -1;
	b = b * V0sobre2;
	% Calcula a distribui��o de carga
	rho = A \ b;
	% Calcula a carga total
	q = dl2 * sum(rho);
	% Calcula a capacit�ncia
	C = abs(q/V0);
end
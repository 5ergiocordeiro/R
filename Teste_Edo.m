%{
Funções para obtenção de raízes de polinômios e funções trasncendentes.
Coeficientes de polinômios devem começar pelo de maior ordem.
%}

function Teste_Edo(testes)
	%feature('DefaultCharacterSet', 'UTF8');
	fprintf('Testando as funções de obtenção de solução de equações diferenciais. \n');	
	% Funções de teste
	% ... ordem, derivada, valor inicial, intervalo, número de subintervalos, métodos, tolerâncias
	% ... especial:
	% ....... ordem para os algoritmos de Runge-Kutta ou
	% ....... índice para os algoritmos de Adams
	global f1 f2 f3 f4 f5 f6 metodos metAdams;
	metodos = {'EMIRDA', {'de Euler'; 'de Euler modificado'; 'de Euler melhorado';'de Runge-Kutta'; 'de Dormand-Price'; 'de Adams'}};
	metAdams = {'-Bashforth de segunda ordem'; '-Bashforth de terceira ordem'; '-Bashforth-Moulton (quarta ordem)'; '-Moulton de segunda ordem'; '-Moulton de segunda ordem'};
	xf1 = @(x,y) x - 2*y + 1;
	f1 = {1, xf1, 1, [0 1], [10 100 10 10 10], ['E';'E';'R';'D';'A'], {[], []}, [0 0 4 0 5]};
	xf2 = @(x,y) y - x^2 + 1;
	f2 = {1, xf2, 0.5, [0 2], [10], ['E'], {[], []}};
	xf3 = @(x,y) - 2 * x * y;
	f3 = {1, xf3, 0.5, [0 1], [10 10 10], ['E';'M';'I'], {[], []}};

	len = size(testes,1);
	for i = 1:len
		Edo_Teste(testes(i));
	end
end

function Edo_Teste(teste)
% Despacha a função de execução do 'teste' pedido
	global f1 f2 f3 f4 f5 f6;
	fprintf('Teste %d: \n', teste);
	switch teste
		case 1
			Edo_Teste1(f1);
			%Edo_Teste1(f2);
			%Edo_Teste1(f3);
		case 2
		case 3
		case 4
	end
end
	
function Edo_Teste1(f)
	global metodos metAdams;
	[ordem, fxy, y0, limx, n, meth, tol, rk] = f{1:8};
	fstr = func2str(fxy);
	fprintf('\t Testando a solução da equação diferencial y = %s \n', fstr(7:end));
	ntestes = size(meth, 1);
	methstr = metodos{2};
	for i = 1:ntestes
		methnum = findstr(metodos{1}, meth(i));
		strmeth = methstr(methnum);
		name = strmeth{1};
		if meth(i) == 'R'
			name = sprintf('%s de ordem %d', name, rk(i))
		end
		if meth(i) == 'A'
			parmAdams = metAdams(rk(i));
			strAdams = parmAdams{1};
			name = sprintf('%s%s', name, strAdams);
		end
		fprintf('\t \t Método %s no intervalo (%d,%d) com %d subintervalos: \n', name, limx(1), limx(2), n(i));
		[x, y, res, niter, err] = Edo_Sol(ordem, fxy, y0, limx, n(i), meth(i), 0, rk(i));
		if res == 0
			size(err)
			if size(err, 2) > 1
				fprintf('Resultado: \n \t \t x \t \t y \t \t erro \n');
				table = [x' y' err'];
			else
				fprintf('Resultado: \n \t \t x \t \t y \n');
				table = [x' y'];
			end
			disp(table);
		end
	end
end
	

function [x, y, res, niter, err] = Edo_Sol(ordem, fxy, y0, limx, n, meth, tol, rk)
	err = 0;
	switch meth
		case 'E'
			[x, y, res, niter] = Edo_Sol_Euler(ordem, fxy, y0, limx, n);
		case 'M'
			[x, y, res, niter] = Edo_Sol_Euler_Mod(ordem, fxy, y0, limx, n);
		case 'I'
			[x, y, res, niter] = Edo_Sol_Euler_Melh(ordem, fxy, y0, limx, n);
		case 'R'
			[x, y, res, niter] = Edo_Sol_RK(ordem, fxy, y0, limx, n, rk);
		case 'D'
			[x, y, res, niter, err] = Edo_Sol_DP(ordem, fxy, y0, limx, n);
		case 'A'
			[x, y, res, niter] = Edo_Sol_Adams(ordem, fxy, y0, limx, n, rk);
	end
end

function [x, y, res, niter] = Edo_Sol_Euler(ordem, fxy, y0, limx, n)
% Resolve a equação diferencial pelo método de Euler
	h = (limx(2) - limx(1)) / n;
	x = linspace(limx(1), limx(2), n + 1);
	y = zeros(1, n + 1);
	y(1) = y0;
	for i = 1:n
		y(i+1) = y(i) + h * fxy(x(i), y(i));
	end
	res = 0;
	niter = 0;
end

function [x, y, res, niter] = Edo_Sol_Euler_Mod(ordem, fxy, y0, limx, n)
% Resolve a equação diferencial pelo método de Euler modificado
	h = (limx(2) - limx(1)) / n;
	hsobre2 = h / 2;
	x = linspace(limx(1), limx(2), n + 1);
	y = zeros(1, n + 1);
	y(1) = y0;
	for i = 1:n
		yaux = y(i) + hsobre2 * fxy(x(i), y(i));
		y(i+1) = y(i) + h * fxy(x(i) + hsobre2, yaux);
	end
	res = 0;
	niter = 0;
end

function [x, y, res, niter] = Edo_Sol_Euler_Melh(ordem, fxy, y0, limx, n)
% Resolve a equação diferencial pelo método de Euler melhorado
	h = (limx(2) - limx(1)) / n;
	hsobre2 = h / 2;
	x = linspace(limx(1), limx(2), n + 1);
	y = zeros(1, n + 1);
	y(1) = y0;
	for i = 1:n
		faux = fxy(x(i), y(i));
		yaux = y(i) + h * faux;
		y(i+1) = y(i) + hsobre2 * (faux + fxy(x(i+1), yaux));
	end
	res = 0;
	niter = 0;
end

function [x, y, res, niter] = Edo_Sol_RK(ordem, fxy, y0, limx, n, rk)
% Resolve a equação diferencial pelo método de Runge-Kutta clássico
	RKParms = {
	%		b					c			
		{[1/2 1/2],			[0 1],			[0 0; 1 0]},								% segunda ordem (Euler melhorado)
		{[1/6 2/3 1/6],		[0 1/2 1],		[0 0 0; 1/2 0 0; -1 2 0]},					% terceira ordem
		{[1/6 1/3 1/3 1/6],	[0 1/2 1/2 1],	[0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0]}	% quarta ordem
		};
	m = zeros(1, rk + 1);
	parms = RKParms{rk - 1};
	b = parms{1};
	c = parms{2};
	a = parms{3};
	%drawButcher(a, b, c);
	h = (limx(2) - limx(1)) / n;
	x = linspace(limx(1), limx(2), n + 1);
	y = zeros(1, n + 1);
	y(1) = y0;
	for i = 1:n
		for j = 1:rk
			xj = x(i) + c(j) * h;
			yaux = m(1:j) * a(j,1:j)';
			yj = y(i) + h * yaux;
			m(j+1) = fxy(xj, yj);
		end
		y(i+1) = y(i) + h * m(2:rk+1) * b';
	end
	res = 0;
	niter = 0;
end

function drawButcher(a, b, c)
	n = size(a,1);
	fprintf('Matriz de Butcher: \n');
	for i = 1:n
		fprintf('%f |', c(i));
		for j = 1:(i-1)
			fprintf(' %f', a(i,j));
		end
		fprintf('\n');
	end
	fprintf('%s \n \t \t | ', repmat('-', 1, n * 11));
	for j = 1:n
		fprintf('%f ', b(j));
	end
	fprintf('\n');
end

function [x, y, res, niter, err] = Edo_Sol_DP(ordem, fxy, y0, limx, n)
% Resolve a equação diferencial pelo método de Dormand-Price
	a = [0 0 0 0 0 0 0; 1/5 0 0 0 0 0 0; 3/40 9/40 0 0 0 0 0; 44/45 -56/15 32/9 0 0 0 0; 
		19372/6561 -25360/2187 64448/6561 -212/729 0 0 0; 9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0;
		35/384 0 500/1113 125/192 -2187/6784 11/84 0];
	b = a(7,:);
	c = [0 1/5 3/10 4/5 8/9 1 1];
	e = [71/57600 0 -71/16695 71/1920 -17253/339200 22/525 -1/40];
	m = zeros(1,8);
	%drawButcher(a, b, c);
	h = (limx(2) - limx(1)) / n;
	x = linspace(limx(1), limx(2), n + 1);
	y = zeros(1, n + 1);
	err = zeros(1, n + 1);
	y(1) = y0;
	err(1) = 0;
	for i = 1:n
		for j = 1:7
			xj = x(i) + c(j) * h;
			yaux = m(1:j) * a(j,1:j)';
			yj = y(i) + h * yaux;
			m(j+1) = fxy(xj, yj);
		end
		y(i+1) = y(i) + h * m(2:8) * b';
		err(i+1) = m(2:8) * e';
	end
	res = 0;
	niter = 0;
end

function [x, y, res, niter] = Edo_Sol_Adams(ordem, fxy, y0, limx, n, idx)
% Resolve a equação diferencial pelo método de Adams
	AdamsParms = [
		%	b							k
		% Adams-Bashforth
		{[0 3/2 -1/2 0 0], 				2};
		{[0 23/12 -16/12 5/12 0], 		3};
		{[0 55/24 -59/24 37/24 -9/24], 	4};
		% Adams-Moulton
		{[5/12 8/12 -1/12 0 0], 		2};
		{[9/24 19/24 -5/24 1/24 0], 	3};
		];
	parms = AdamsParms(idx,1:2);
	[b, k] = parms{1:2};
	h = (limx(2) - limx(1)) / n;
	x = linspace(limx(1), limx(2), n + 1);
	y = zeros(1, n + 1);
	f = zeros(1, 5);
	y(1) = y0;
	[lixo1, est, lixo2, lixo3, lixo4] = Edo_Sol_DP(0, fxy, y(1), [x(1) x(k)], k - 1);
	hest = est;
	for i = 2:k
		f(k + 1) = fxy(x(i), est(i));
	end
	for i = 1:n
		if b(1) ~= 0
			[lixo1, est, lixo2, lixo3, lixo4] = Edo_Sol_DP(0, fxy, y(i), [x(i) x(i+1)], 2);
			hest(i+1) = est(2);
			f(1) = fxy(x(i+1), est(2));
		end
		f(2) = fxy(x(i), y(i));
		y(i+1) = y(i) + h * f * b';
		f(2:5) = f(1:4);
	end
	res = 0;
	niter = 0;
end

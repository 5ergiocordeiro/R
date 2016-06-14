%{
Funções para obtenção de raízes de polinômios e funções trasncendentes.
Coeficientes de polinômios devem começar pelo de maior ordem.
%}

function Teste_Roots(testes)
	%feature('DefaultCharacterSet', 'UTF8');
	fprintf('Testando as funções de obtenção de raízes de funções. Exige o toolbox de processamento simbólico. \n');	
	% Funções de teste
	% ... polinômios:
	% ...... coeficientes, limites para as raízes, número de raízes,
	% ...... raiz real, número de interações, tolerâncias, métodos usados, intervalo de busca,
	% ...... derivada primeira, derivada segunda e fator de aceleração
	global f1 f2 f3 f4 f5 f6 f7 f8;
	f1 = {[1 2 -13 -14 24], [4.74 0.63 -14 -0.58], {[2 0], [2 0], [0]}, [-4 3], [8 6], [1e-3 1e-6], ['V';'N'], [2 4], [4 6 -26 -14], [12 12 -26]};
	f2 = {[1 -3 -6 8], [7, 0.57, -3.83, -0.62], {[2 0], [1], [0]}, [-1.9993], [7], [0.05], ['B'], [-3.83 -0.62] [], []};	
	f3 = {[1 -5 7 19 -98 -104 0], [105 0.7 -5.61 -0.51], {[3 1], [2 0], [1]}};
	f8 = {[1 2 -12 14 -5], [], {}, [1 1], [27 4], [1e-3 1e-3], ['N';'N'], [0 2], [4 6 -24 14], [12 12 -24], [1 3]};
	% ... funções transcendentais
	% ...... expressão, raiz real, número de interações, tolerâncias e métodos usados	
	% sym x;
	xf4 = @(x) 2 * x.^3 - cos(x + 1) - 3;
	f4 = {xf4, [1.08008 1.07912 1.07903 1.0791 1.07912], [9 7 15 7 5], [5e-3 1e-3 1e-3 1e-3 1e-3], ['B';'S';'R'; 'P'; 'M'], [-1 2]};
	xf5 = @(x) 0.05 * x.^3 - 0.4 * x.^2 + 3 * sin(x) .* x;
	f5 = {xf5, [11.74390 11.74393], [13 7], [5e-3 1e-3], ['B'; 'V'], [10 12]};	
	xf6 = @(x) 3 * x.^2 + sqrt(x + 1) .* cos(x).^3 - 2;
	f6 = {xf6, [0.6812], [6], [1e-5], ['P'], [0 1]};
	xf7 = @(x) sin(x) - exp(x) - 2 * x.^2 + 10;
	xf7linha = @(x) cos(x) - exp(x) - 4 * x;
	xf72linha = @(x) cos(x) - exp(x) - 4 * x;
	f7 = {xf7, [1.67866], [4], [1e-5], ['N'], [1 2], xf7linha, xf72linha, [1]};
	

	len = size(testes,1);
	for i = 1:len
		Roots_Teste(testes(i));
	end
end

function Roots_Teste(teste)
% Despacha a função de execução do 'teste' pedido
	global f1 f2 f3 f4 f5 f6 f7 f8;
	fprintf('Teste %d: \n', teste);
	switch teste
		case 1
			Roots_Teste1(f1);
			Roots_Teste1(f2);
			Roots_Teste1(f3);
		case 2
		case 3
			%Roots_Teste3(f1);
			%Roots_Teste3(f2);
			Roots_Teste3(f8);
		case 4
			%Roots_Teste4(f4);
			%Roots_Teste4(f5);
			%Roots_Teste4(f6);
			Roots_Teste4(f7);
	end
end
	
function Roots_Teste1(f)
	fprintf('\t Testando a análise de raízes do polinômio com coeficientes ');
	% Obtém os coeficientes, os limites esperados e o número de raízes de cada tipo
	[coef, lim, nr] = f{1:3};
	fprintf('%f ', coef);
	% Executa a análise
	[limits, nroots] = RootsPolyAn(coef);
	% Verifica o acerto da análise feita
	error = sqrt(sum((lim - limits) .^ 2));
	fprintf('\n \t \t Função f1: \n');
	fprintf('\t \t \t ENCONTRADO: Raízes positivas no intervalo (%f,%f). Raízes negativas no intervalo (%f,%f). \n', limits(2), limits(1), limits(3), limits(4));
	fprintf('\t \t \t PREVISTO: Raízes positivas no intervalo (%f,%f). Raízes negativas no intervalo (%f,%f). \n', lim(2), lim(1), lim(3), lim(4));
	fprintf('\t \t \t Erro: %f. \n', error);
	[p, n, z] = nroots{:};
	fprintf('\t \t \t ENCONTRADO: Número de raízes: Positivas: ');
	fprintf('%d ', p); fprintf('\t Negativas:'); fprintf('%d ', n); fprintf('\t Nulas:'); fprintf('%d ', z);
	[p, n, z] = nr{:};
	fprintf('\n \t \t \t PREVISTO: Número de raízes: Positivas: ');
	fprintf('%d ', p); fprintf('\t Negativas:'); fprintf('%d ', n); fprintf('\t Nulas:'); fprintf('%d ', z);
	fprintf('\n');
end
	
function Roots_Teste2(f)
end

function Roots_Teste3(f)
	fprintf('\t Testando a obtenção de raiz do polinômio [');
	fprintf('%f ', f{1});
	[coef, nu1, nu2, raiz, iter, tol, meth, region, flinha, f2linha, m] = f{1:11};
	nmeth = size(meth, 1);
	for i = 1:nmeth
		[root, niter] = RootsCalc('P', meth(i), {coef, region, tol(i), flinha, f2linha, m(i)});
		fprintf('] \n \t \t Método %s, intervalo (%f, %f), com tolerância %f: ', meth(i), region(1), region(2), tol(i));
		if niter <= 0
			niter = - niter;
			fprintf('divergiu');
		else
			fprintf('%f', root);
		end
		fprintf(', com %d iterações. PREVISTO = %f. com %d iterações\n', niter, raiz(i), iter(i));
	end
end

function Roots_Teste4(f)
	[fun, raiz, iter, tol, meth, region, flinha, f2linha] = f{1:8};
	fstr = func2str(fun);
	fprintf('\t Testando a obtenção de raiz da função f = %s \n', fstr(5:end));
	nmeth = size(meth, 1);
	for i = 1:nmeth
		[root, niter] = RootsCalc('F', meth(i), {fun, region, tol(i), flinha, f2linha});
		fprintf('\n \t \t Método %s, intervalo (%f, %f), com tolerância %f: ', meth(i), region(1), region(2), tol(i));
		if niter <= 0
			niter = - niter;
			fprintf('divergiu');
		else
			fprintf('%f', root);
		end
		fprintf(', com %d iterações. PREVISTO = %f. com %d iterações\n', niter, raiz(i), iter(i));
	end
end

	
function [L, n] = RootsPolyAn(coef)
% Analisa as raízes do polinômio cujos coeficientes 'coef' são dados.
	% Deflaciona o polinômio, se possível
	len = size(coef, 2);
	last = find(coef, 1, 'last');
	poly = coef(1:last);
	nz = len - last;
	% Calcula os intervalos onde pode haver raízes
	L = RootsPolyLim(poly);
	[np, nn] = RootsPolyNum(poly);
	n = {np, nn, nz};
end

function [p, n] = RootsPolyNum(coef)
% Calcula o número de raízes positivas e negativas que um polinômio pode ter
	poly = coef(coef ~= 0);
	len = size(poly, 2);
	pos = (poly > 0);
	comp = pos(1:(len-1)) == pos(2:len);
	np = sum(comp == 0);
	p = np:-2:0;
	nn = len - np - 1;
	n = nn:-2:0;
end

function L = RootsPolyLim(coef)
% Calcula os intervalos	onde pode haver raízes do polinômio cujos coeficientes 'coef' são dados.
	len = size(coef, 2);
	n = len - 1;
	% Calcula L0 e L1
	L0 = RootsPolyCalcL(coef);
	L1 = RootsPolyCalcL(coef(len:-1:1));
	% Inverte o sinal dos coeficientes das potências ímpares
	coef = coef .* (1 - 2 * mod(n:-1:0, 2));
	% Calcula L2 e L3
	L2 = RootsPolyCalcL(coef);
	L3 = RootsPolyCalcL(coef(len:-1:1));
	% Monta a resposta
	L = [L0, 1/L1, -L2, -1/L3];
end

function [L] = RootsPolyCalcL(coef)
% Calcula um dos limites para o arranjo dado 'coef' dos coeficientes do polinômio
	% Acerta os sinais, se preciso
	if coef(1) < 0
		coef = - coef;
	end
	% Identifica os coeficientes negativos
	idxneg = find(coef(:) < 0);
	% Calcula os valores
	B = max(abs(coef(idxneg)));
	L = 1 + (B / coef(1))^(1/(idxneg(1) - 1));
end

function [root, iter] = RootsCalc(typ, meth, fun)
	switch meth
		case 'B'
			[root, iter] = RootsCalcBis(typ, fun{1}, fun{2}, fun{3});
		case 'S'
			[root, iter] = RootsCalcSec(typ, fun{1}, fun{2}, fun{3});
		case 'R'
			[root, iter] = RootsCalcReg(typ, fun{1}, fun{2}, fun{3});
		case 'P'
			[root, iter] = RootsCalcPeg(typ, fun{1}, fun{2}, fun{3});
		case 'M'
			[root, iter] = RootsCalcMil(typ, fun{1}, fun{2}, fun{3});
		case 'V'
			[root, iter] = RootsCalcvWDB(typ, fun{1}, fun{2}, fun{3});
		case 'N'
			[root, iter] = RootsCalcNew(typ, fun{1}, fun{2}, fun{3}, fun{4}, fun{5}, fun{6});
	end
end

function [root, iter] = RootsCalcBis(typ, parm, lim, tol)
% Calcula as raízes do polinômio ou função dados pelo método da bisseção
	if typ == 'P'
		flim = polyval(parm, lim);
	else
		flim = parm(lim);
	end
	slim = sign(flim);
	m = lim(2) - lim(1);
	lastm = 1e3;
	iter = 0;
	while m > tol
		if m > lastm
			root = 0;
			iter = - iter;
			return;
		end
		lastm = m;
		iter = iter + 1;
		med = mean(lim);
		if typ == 'P'
			fmed = polyval(parm, med);
		else
			fmed = parm(med);
		end
		smed = sign(fmed);
		if smed == slim(1)
			lim(1) = med;
		else
			lim(2) = med;
		end
		m = lim(2) - lim(1);
		%fprintf('\t \t \t %d %f %f %f %f %f \n', iter, lim(1), lim(2), flim(1), flim(2), m);
	end
	root = med;
end

function [root, iter] = RootsCalcSec(typ, parm, lim, tol)
% Calcula as raízes do polinômio ou função dados pelo método da secante
	old = [0 0];
	if typ == 'P'
		flim = polyval(parm, lim);
	else
		flim = parm(lim);
	end
	slim = sign(flim);
	m = lim(2) - lim(1);
	lastm = 1e3;
	iter = 0;
	while m > tol
		if slim(1) == slim(2) || m > lastm
			root = 0;
			iter = - iter;
			return;
		end
		lastm = m;
		iter = iter + 1;
		next = (lim(1)*flim(2) - lim(2)*flim(1)) / (flim(2) - flim(1));
		if typ == 'P'
			fnext = polyval(parm, next);
		else
			fnext = parm(next);
		end
		snext = sign(fnext);
		if old(1) == old(2)
			idxold = (2 - (snext == slim(1)));
			old(idxold) = 1; 
		end
		idxold = (2 - (old(1) == 1));
		lim(idxold) = next;
		flim(idxold) = fnext;
		slim(idxold) = snext;
		old(idxold) = 0;
		m = lim(2) - lim(1);
	end
	root = next;
end

function [root, iter] = RootsCalcReg(typ, parm, lim, tol)
% Calcula as raízes do polinômio polinômio ou função dados pelo método Regula Falsi
	if typ == 'P'
		flim = polyval(parm, lim);
	else
		flim = parm(lim);
	end
	slim = sign(flim);
	m = lim(2) - lim(1);
	lastm = 1e3;
	iter = 0;
	while m > tol
		if slim(1) == slim(2) || m > lastm
			root = 0;
			iter = - iter;
			return;
		end
		lastm = m;
		iter = iter + 1;
		next = (lim(1)*flim(2) - lim(2)*flim(1)) / (flim(2) - flim(1));
		if typ == 'P'
			fnext = polyval(parm, next);
		else
			fnext = parm(next);
		end
		snext = sign(fnext);
		idxchange = (2 - (slim(1) == snext));
		lim(idxchange) = next;
		flim(idxchange) = fnext;
		slim(idxchange) = snext;
		m = lim(2) - lim(1);
	end
	root = next;
end

function [root, iter] = RootsCalcPeg(typ, parm, lim, tol)
% Calcula as raízes do polinômio polinômio ou função dados pelo método Pégaso
	if typ == 'P'
		flim = polyval(parm, lim);
	else
		flim = parm(lim);
	end
	slim = sign(flim);
	m = lim(2) - lim(1);
	lastm = 1e3;
	fold = [0 0 0];
	iter = 0;
	while m > tol
		if slim(1) == slim(2) || m > lastm
			root = 0;
			iter = - iter;
			return;
		end
		lastm = m;
		iter = iter + 1;
		if iter > 3
			idxold = 3 - idxchange;
			p = fold(1) * fold(2) / (fold(2) + fold(3));
			if sign(p) ~= slim(idxchange)
				flim(idxold) = p
			end
		end
		next = (lim(1)*flim(2) - lim(2)*flim(1)) / (flim(2) - flim(1));			
		if typ == 'P'
			fnext = polyval(parm, next);
		else
			fnext = parm(next);
		end
		snext = sign(fnext);
		idxchange = (2 - (slim(1) == snext));
		lim(idxchange) = next;
		flim(idxchange) = fnext;
		slim(idxchange) = snext;
		fold = [fold(2), fold(3), fnext];
		m = lim(2) - lim(1);
	end
	root = next;
end

function [root, iter] = RootsCalcMil(typ, parm, lim, tol)
% Calcula as raízes do polinômio polinômio ou função dados pelo método de Miller
	x = [lim(1), mean(lim), lim(2)];
	if typ == 'P'
		y = polyval(parm, x);
	else
		y = parm(x);
	end
	m = lim(2) - lim(1);
	lastm = 1e3;
	iter = 0;
	while m > tol
		if m > lastm
			root = 0;
			iter = - iter;
			return;
		end		
		lastm = m;
		iter = iter + 1;
		h = x([3 2]) - x([2 1]);
		r = h(1) / h(2);
		a = (y * [r, -(r + 1), 1]') / (h(1) * sum(h));
		b = (y(3) - y(2)) / h(1) - a * h(1);
		delta = b^2 - 4 * a * y(2);
		z = (- b + sign(b) * sqrt(delta) ) / (2 * a);
		next = x(2) + z;
		if typ == 'P'
			fnext = polyval(parm, next);
		else
			fnext = parm(next);
		end
		dist = abs(next - x([1 3]));
		if dist(2) > dist(1)
			newx = [x(1), x(2), next];
			newy = [y(1), y(2), fnext];
			m = dist(1);
		else
			newx = [x(2), x(3), next];
			newy = [y(2), y(3), fnext];
			m = dist(2);
		end
		[x, idx] = sort(newx);
		y = newy(idx); 
		%fprintf('\t \t \t %d %f %f %f \n', iter, next, fnext, m);
	end
	root = next;
end

function [root, iter] = RootsCalcvWDB(typ, parm, lim, tol)
% Calcula as raízes do polinômio polinômio ou função dados pelo método de van Wijngaarden-Dekker-Brent
	x = [lim(1), mean(lim), lim(2)];
	if typ == 'P'
		y = polyval(parm, x);
	else
		y = parm(x);
	end
	m = lim(2) - lim(1);
	lastm = 1e3;
	iter = 0;
	while m > tol
		if m > lastm
			root = 0;
			iter = - iter;
			return;
		end				
		lastm = m;
		iter = iter + 1;
		next = LagPol(y', x', 0);
		if typ == 'P'
			fnext = polyval(parm, next);
		else
			fnext = parm(next);
		end
		snext = sign(fnext);
		if snext == sign(y(2))
			if snext == sign(y(1))
				idxkeep = 3;
			else
				idxkeep = 1;
			end
		else
			dist = abs(next - x([1 3]));
			if dist(2) > dist(1)
				idxkeep = 1;
			else				
				idxkeep = 3;
			end
		end
		newx = [x(2), next, x(idxkeep)];
		newy = [y(2), fnext, y(idxkeep)];
		m = abs(next - x(2));	
		[x, idx] = sort(newx); 
		y = newy(idx);
		%fprintf('\t \t \t %d %f %f %f %f %f %f %f \n', iter, x(1), x(2), x(3), y(1), y(2), y(3), m);
	end
	root = next;
end

function [root, iter] = RootsCalcNew(typ, p, lim, tol, plinha, p2linha, mult)
% Calcula as raízes do polinômio ou função dados pelo método de Newton
	if typ == 'P'
		flim = polyval(p, lim);
		flim2linha = polyval(p2linha, lim);
	else
		flim = p(lim);
		flim2linha = p2linha(lim);
	end
	slim = sign(flim);
	slim2linha = sign(flim2linha);
	sel = slim .* slim2linha
	x = lim(sel == 1)
	x = x(2)
	m = lim(2) - lim(1);
	lastm = 1e3;
	iter = 0;
	while m > tol
		if m > lastm
			root = 0;
			iter = - iter;
			return;
		end
		lastm = m;
		iter = iter + 1;
		if typ == 'P'
			fnext = polyval(p, x);
			fnextlinha = polyval(plinha, x);
		else
			fnext = p(x);
			fnextlinha = plinha(x);			
		end
		next = x - mult * fnext / fnextlinha;
		m = abs(next - x);
		fprintf('\t \t \t %d %f %f %f %f \n', iter, x, next, fnext, m);
		x = next;
	end
	root = next;
end

function [fpoint] = LagPol(x, y, point)
% Interpola o ponto 'point' pelo método de Lagrange
	n = size(x,1);
	coef = zeros(1, n);
	for i = 1:n
		sel = x([1:n] ~= i);
		coef(i) = prod(point - sel) / prod(x(i) - sel);
	end
	fpoint = coef * y;
end
function [erro,V,E] = calcevk(Po,Q)
%{
Retorna o campo elétrico e o potencial elétrico, na posição dada por Po, para uma distribuição de cargas dada por Q1, Q2 ... Qm.
Uso:
	[erro,V,E] = calcevk(Po,Q)
	onde
		Po é a posição do observador, em coordenadas cartesianas (dimensão [1,3]), com unidade m,
		Q é uma matriz de m linhas, cada linha correspondendo a uma carga.
			Qi é um vetor de dimensão [1,4] onde o primeiro elemento é o valor da carga, em nC, e as demais são a sua posição, em coordenadas cartesianas, com unidade m,
		erro é 0 quando o cálculo pode ser feito e 1 em caso contrário,
		V é o valor do potencial elétrico em Po, com unidade V,
		E é o vetor campo elétrico na posição Po, em coordenadas cartesianas (dimensão [1x3]), com unidade V/m.
Se as dimensões de Po e Q não forem consistentes, retorna NaN em V e PE, e 1 em erro.
Limitações:
	1) Funciona apenas em coordenadas cartesianas.
	2) Funciona apenas no espaço livre.
	3) Não deve haver movimento relativo entre as cargas e o observador.
	4) O campo magnético na região deve ser nulo.
	
Testado com Octave 4.0 em Windows 8.1 (ver função testcalcevk).
Sérgio Cordeiro = 05/09/2015
%}
	% Testa se o número de argumentos e as dimensões de Po e Q são consistentes
	if (nargin != 2) || (nargout != 3) || (size(Q)(2) != 4) || (prod(size(Po) == [1 3]) != 1)
		V = NaN;
		E = [V V V];
		erro = 1;
		return
	end
	% Calcula o potencial e o campo devidos a cada carga individual e depois soma para obter os resultantes
	k = 8.988;				% para carga em nC
	kq = k * Q(:,1);
	p = Q(:,2:4);
	r = Po - p;
	% Testa se alguma carga se encontra em Po
	if sum(sum(r,2) == 0)
		V = NaN;
		E = [V V V];
		erro = 1;
		return	
	end
	% Calcula V e E
	d = sqrt(sum(r.^2,2));
	d3 = d.^3;
	V = sum(kq./d);
	E = sum(r.*(kq./d3),1);
	erro = 0;
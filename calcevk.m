function [erro,V,E] = calcevk(Po,Q)
%{
Retorna o campo el�trico e o potencial el�trico, na posi��o dada por Po, para uma distribui��o de cargas dada por Q1, Q2 ... Qm.
Uso:
	[erro,V,E] = calcevk(Po,Q)
	onde
		Po � a posi��o do observador, em coordenadas cartesianas (dimens�o [1,3]), com unidade m,
		Q � uma matriz de m linhas, cada linha correspondendo a uma carga.
			Qi � um vetor de dimens�o [1,4] onde o primeiro elemento � o valor da carga, em nC, e as demais s�o a sua posi��o, em coordenadas cartesianas, com unidade m,
		erro � 0 quando o c�lculo pode ser feito e 1 em caso contr�rio,
		V � o valor do potencial el�trico em Po, com unidade V,
		E � o vetor campo el�trico na posi��o Po, em coordenadas cartesianas (dimens�o [1x3]), com unidade V/m.
Se as dimens�es de Po e Q n�o forem consistentes, retorna NaN em V e PE, e 1 em erro.
Limita��es:
	1) Funciona apenas em coordenadas cartesianas.
	2) Funciona apenas no espa�o livre.
	3) N�o deve haver movimento relativo entre as cargas e o observador.
	4) O campo magn�tico na regi�o deve ser nulo.
	
Testado com Octave 4.0 em Windows 8.1 (ver fun��o testcalcevk).
S�rgio Cordeiro = 05/09/2015
%}
	% Testa se o n�mero de argumentos e as dimens�es de Po e Q s�o consistentes
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
function [V,calc]=initv(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,RhoNint,RhoNext,ZNint,ZNext,Nrho,Nz);
% Inicializa o potencial 'V' para o c�lculo iterativo.
% Retorna a grade inicializada e a posi��o dos cilindros nela.
% Chamada pela fun��o 'calcv'.
	% Cria a grade
	V = zeros(Nrho,Nz);
	% Calcula a aproxima��o inicial.
	% Teoricamente, a aproxima��o inicial pode ser qualquer; no etantom para tentar auxiliar a converg�ncia, foram escolhidas as condi��es seguintes:
	% ... Vint dentro do cilindro interno
	V(1:RhoNint,ZNint(1):ZNint(2)) = Vint;
	% ... queda exponencial at� Vext no cilindro externo
	a = (Vint - Vext) / (RhoNint - RhoNext);
	for i = RhoNint:RhoNext
		V(i,ZNext(1):ZNext(2)) = Vint * log(RhoNext/i);
	end
	% Onde n�o calcular o Laplaciano
	calc = ones(Nrho,Nz);
	calc(RhoNext,ZNext(1):ZNext(2)) = 0 ;		% no cilindro externo
	calc(RhoNint,ZNint(1):ZNint(2)) = 0;		% no cilindro interno
end

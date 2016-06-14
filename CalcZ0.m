function Z0 = CalcZ0(h, nt)
% Resolve o exemplo 3.6 do livro "Numerical Techniques in Electromagnetic with MATLAB", de Matthew N. O. Sadiku, terceira edi��o, pp. 132 a 159}.
% 'h' � o tamanho da malha usada para o c�lculo do potencial pelo m�todo das diferen�as finitas, e 'nt' � o n�mero de itera��es para o m�todo iterativo. Quando 'nt' = 0, usa-se o m�todo de solu��o de um sistema linear.
	% Inicializa��o
	% ... constantes f�sicas:
	global epsilon0 c;
	[epsilon0, c] = deal(8.81e-12, 3e8);
	% ... par�metros do problema:
	% ...... dimens�es em cm
	[a, b, d, w, epsilonr] = deal(2.5, 2.5, 0.5, 1.0, 2.35);
	% ... posi��es fixas
	pos = num2cell(round([a b d w]/h));
	[nx, ny, jd, iw] = deal(pos{:});
	% C�lculos
	% ... calcula o potencial duas vezes, uma com e a outra sem o diel�trico
	if nt > 0
		V0 = CalcVI(nx, ny, jd, iw, 1, nt);
		Vd = CalcVI(nx, ny, jd, iw, epsilonr, nt);
	else
		V0 = CalcVS(nx, ny, jd, iw, 1);
		Vd = CalcVS(nx, ny, jd, iw, epsilonr);	
	end
	% ... calcula o campo el�trico para cada distribui��o de potencial
	E0 = Grad(V0, h * 0.01);
	Ed = Grad(Vd, h * 0.01);
	% ... plota os gr�ficos de V e E
	plotVE(V0, E0, 2);
	plotVE(Vd, Ed, 4);	
	% ... calcula a carga el�trica e a capacit�ncia, considerando que a diferen�a de potencial aplicada foi 1 Volt
	% ...... a superf�cie Gaussiana pode estar em qualquer lugar
	out = round(([iw jd] + [nx ny]) / 2);
	%out = [iw jd] + 3;
	%out = [100 100];
	C0 = 4 * CalcQ(V0, out(1), out(2), jd, [1 1] * epsilon0);
	Cd = 4 * CalcQ(Vd, out(1), out(2), jd, [1 epsilonr] * epsilon0);
	% ... calcula a imped�ncia da linha
	Z0 = 1 / (c * sqrt(C0 * Cd));
end

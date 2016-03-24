function plotheatmap(x,y,vals,texts,paleta,tipo,cbar)
%{
%}
	if ! strcmp(paleta,"")
		colormap(paleta);
	end
	switch tipo
		case "blur"
			pcolor(x,y,vals);
			shading interp;
		case "contour"
			contour(x,y,vals);			
		otherwise
			imagesc(x,y,vals);
	end
	xlabel(strtrim(texts(1,:))); ylabel(strtrim(texts(2,:)));
	last = size(x)(1);
	xlim([x(1) x(last)]);
	last = size(y)(1);
	ylim([y(1) y(last)]);
	if ! strcmp(cbar,"no")
		ab = colorbar();
		set(ab,'title',strtrim(texts(4,:)));
	end
	title(strtrim(texts(3,:)));

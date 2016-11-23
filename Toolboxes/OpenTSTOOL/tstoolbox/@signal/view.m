function view(s, fontsize, fhandle)

%tstoolbox/@signal/view
%   Syntax:
%     * view(signal) (fontsize=12)
%     * view(signal, fontsize)
%     * view(signal, fontsize, figurehandle)
%
%   Signal viewer that decides from the signal's attributes which kind of
%   plot to produce, using the signal's plothint entry to get a hint which
%   kind of plot to produce
%   Possible plothints are:
%     * 'graph'
%     * 'bar'
%     * 'surrbar'
%     * 'surrerrorbar'
%     * 'points'
%     * 'xyplot'
%     * 'xypoints'
%     * 'scatter'
%     * '3dcurve'
%     * '3dpoints'
%     * 'spectrogram'
%     * 'image'
%     * 'multigraph'
%     * 'multipoints'
%     * 'subplotgraph'
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(1,3);

if nargin < 2
	fontsize = 13;
end

nd = ndim(s);		% Anzahl der Dimensionen mit Laenge groesser eins
dlen = dlens(s);

% check which type of plot we shall produce
pmode = plothint(s); % take this as hint

switch nd			% but check if the desired way of plotting is possible
	case 1
		default = 'graph';
		switch pmode
			case {'graph', 'bar', 'points','surrbar'}
			otherwise
				pmode = default;
		end
    case 2
		default = 'image';
		switch pmode
			case {'xyplot','xypoints', 'scatter'}
				if dlen(2)~=2, pmode = default; end	
			case {'3dcurve', '3dpoints'}
				if dlen(2)~=3, pmode = default; end		
			case {'spectrogram', 'image', 'multigraph', 'multipoints', 'subplotgraph','surrerrorbar'}
			otherwise
				pmode = default;
		end	
    otherwise
	    error('Plotting for this kind of data not yet supported');
		return
end    

% set labels, titles etc.
titel = label(s);

% do the plotting
switch pmode	
	case 'surrerrorbar'      
		x = spacing(s, 1);
		dat = data(s);
		plot(x,dat(:,1));
		hold on
		errorbar(x,dat(:,2),dat(:,3),dat(:,3),'r');
		hold off
		set(gca, 'FontSize', fontsize)
		xlabel(lab(getaxis(s,1)));
		yn = yname(s);
		if ~isempty(yn)
		  ylabel(yn);
		else
			ylabel(label(yunit(s)));
		end
		zoom on
		axis tight;
		
	  
	case	'graph'
		x = spacing(s, 1);
		dat = data(s);
		if ~isreal(dat)
			phandle = plot(x, real(dat), x, imag(dat));
		else
	    	phandle = plot(x, dat);
		end
	    set(gca, 'FontSize', fontsize)
	    xlabel(lab(getaxis(s,1)));
		yn = yname(s);
		if ~isempty(yn)
	    	ylabel(yn);
	    else
			ylabel(label(yunit(s)));
		end
		zoom on
		axis tight;
	case	'points'
		x = spacing(s, 1);
		dat = data(s);
		if ~isreal(dat)
			phandle = plot(x, real(dat), x, imag(dat), '.');
		else
	    	phandle = plot(x, dat, '.');
		end
	    set(gca, 'FontSize', fontsize)
	    %title(titel);
	    xlabel(lab(getaxis(s,1)));
		yn = yname(s);
		if ~isempty(yn)
	    	ylabel(yn);
	    else
			ylabel(label(yunit(s)));
		end
		zoom on
		axis tight;
	case	'multigraph'
		x = spacing(s, 1);
		dat = data(s);
		phandle = plot(x, dat);
	    set(gca, 'FontSize', fontsize)
	    %title(titel);
	    xlabel(lab(getaxis(s,1)));
		yn = yname(s);
		if ~isempty(yn)
	    	ylabel(yn);
	    else
			ylabel(label(yunit(s)));
		end
		zoom on
		axis tight;
	case	'subplotgraph'
		x = spacing(s, 1);
		N = dlens(s,2); 	% number of graphs = number of subplots
		dat = data(s);
		subplot(N,1,1);
		for i=1:N
			subplot(N,1,i);
			phandle = plot(x, dat(:,i));
	    	set(gca, 'FontSize', ceil(fontsize/sqrt(N)));
	    	xlabel(lab(getaxis(s,1)));
			yn = yname(s);
			if ~isempty(yn)
	    		ylabel(yn);
	    	else
				ylabel(label(yunit(s)));
			end
			zoom on
			axis tight;		
		end
	case	'multipoints'
		x = spacing(s, 1);
		dat = data(s);
		phandle = plot(x, dat, '.');
	    set(gca, 'FontSize', fontsize)
	    %title(titel);
	    xlabel(lab(getaxis(s,1)));
		yn = yname(s);
		if ~isempty(yn)
	    	ylabel(yn);
	    else
			ylabel(label(yunit(s)));
		end
		zoom on
		axis tight;
	      
	case 'surrbar'
	        s1=histo(cut(s,1,2,length(data(s))),100);
		x = spacing(s1, 1);
		dat = data(s1);
		bar(x,dat);
		hold on
		plot([data(cut(s,1,1,1)) data(cut(s,1,1,1))],[0 max(dat)],'r','LineWidth',3);
	%		bar(data(cut(s,1,1,1)),1,(max(x)-min(x))/100, 'r');
		hold off
	    set(gca, 'FontSize', fontsize)
	    %title(titel);
	    ylabel(lab(getaxis(s,1)));
		yn = yname(s);
		if ~isempty(yn)
	    	xlabel(yn);
	    else
			ylabel(label(yunit(s)));
		end	
		axis tight;
	case 'bar'
		x = spacing(s, 1);
		dat = data(s);
		bar(x,dat);
	    set(gca, 'FontSize', fontsize)
	    %title(titel);
	    xlabel(lab(getaxis(s,1)));
		yn = yname(s);
		if ~isempty(yn)
	    	ylabel(yn);
	    else
			ylabel(label(yunit(s)));
		end	
		axis tight;
	case	'xyplot'
		set(gcf, 'Renderer', 'zbuffer')
		tmp = data(s);
		phandle = plot(tmp(:,1),tmp(:,2));
		set(gca, 'FontSize', fontsize);
		%title(titel);
		ylabel(yname(s));
		xlabel(yname(s));
		zoom on
	case	'xypoints'
		set(gcf, 'Renderer', 'zbuffer')
		tmp = data(s);
		phandle = plot(tmp(:,1),tmp(:,2), '.');
		set(gca, 'FontSize', fontsize);
		%title(titel);
		ylabel(yname(s));
		xlabel(yname(s));
		zoom on		
	case	'scatter'
		set(gcf, 'Renderer', 'zbuffer')
		tmp = data(s);
		phandle = scatter(tmp(:,1),tmp(:,2));
		set(gca, 'FontSize', fontsize);
		%title(titel);
		ylabel(yname(s));
		xlabel(yname(s));
		zoom on
	case	'3dcurve'
		set(gcf, 'Renderer', 'zbuffer')
		tmp = data(s);
		phandle = plot3(tmp(:,1),tmp(:,2),tmp(:,3));
		set(gca, 'FontSize', fontsize);
		rotate3d on
		set(gca, 'Projection' , 'perspective', 'LineWidth' , [1]);
		grid on
	case	'3dpoints'
		set(gcf, 'Renderer', 'zbuffer')
		tmp = data(s);
		phandle = plot3(tmp(:,1),tmp(:,2),tmp(:,3), '.');
		set(gca, 'FontSize', fontsize);
		rotate3d on
		set(gca, 'Projection' , 'perspective');
		grid on
	case    'image'
		set(gcf, 'Renderer', 'zbuffer')
		v1 = spacing(s,1);
		v2 = spacing(s,2);
		phandle = imagesc(v2, v1, data(s));
		colormap('default')
		set(gca, 'FontSize', fontsize);
		%title(titel);
		set(gca, 'YDir', 'normal');
		ylabel(lab(getaxis(s,1)));
		xlabel(lab(getaxis(s,2)));
		axis tight
		zoom on
	case 	'spectrogram'
		set(gcf, 'Renderer', 'zbuffer')
	    v1 = spacing(s,1);
		v2 = spacing(s,2);
		phandle = pcolor(v2, v1, data(s));
		shading interp
		colormap('default')
		set(gca, 'FontSize', fontsize);
		%title(titel);
		ylabel(lab(getaxis(s,1)));
		xlabel(lab(getaxis(s,2)));
		axis tight
		zoom on
end 					% switch	

title(titel);
set(gcf, 'PaperPositionMode', 'auto');	% Plotten immer im Aspect Ratio des wirklichen Bildes

% end of main routine

% helper functions

function r = lab(ax)
if ~isempty(label(ax))
	r = [name(ax) ' (' label(ax)  ')'];
else
	r = name(ax);
end

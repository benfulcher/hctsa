function out=EN_entropies(y,n,p)
	switch n
		case 'shannon' % Shannon entropy
			out=wentropy(y,'shannon')/length(y); % scales with N for large N
		case 'logenergy' % Log Energy entropy
			out=wentropy(y,'log energy')/length(y);% scales with N for large N
        case 'threshold' % magnitude of the signal greater than some value
            out=wentropy(y,'threshold',p)/length(y);
        case 'sure'
            out=wentropy(y,'sure',p)/length(y);
    end
end
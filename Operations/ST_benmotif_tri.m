function out = ST_benmotif_tri(y,trit)
% Time series y
% looks for motifs in a course-graining of the time series to an alphabet
% of three letters
% Ben Fulcher 2009

switch trit
	case 'quantile'
		yt = SUB_coursegrain(y,3,'quantile');
	case 'diffquant'
		yt = SUB_coursegrain(diff(y),3,'quantile');
end

% yt=y;
N = length(yt); % time series length

% so we have a vector yt with entries \in {1,2,3}

%% 1
r1 = cell(3,1); % stores ranges as vectors
out1 = zeros(3,1); % stores probabilities as doubles

for i = 1:3
	r1{i} = find(yt == i);
	out1(i) = length(r1{i})/N;
end


% ------ Record these -------
out.a = out1(1); % proportion of a
out.b = out1(2); % proportion of b
out.c = out1(3); % proportion of c
out.h = -sum(out1(out1>0).*log(out1(out1>0))); % entropy of this result

%% 2
% Make sure ranges are valid for looking at the next one
for i = 1:3
	if (~isempty(r1{i}) && r1{i}(end) == N)
        r1{i} = r1{i}(1:end-1);
    end
end

r2 = cell(3,3);
out2 = zeros(3,3);
for i = 1:3
	for j = 1:3
		r2{i,j} = r1{i}(yt(r1{i}+1) == j);
		out2(i,j) = length(r2{i,j})/(N-1);
	end
end


% ------ Record these -------
out.aa = out2(1,1); out.ab = out2(1,2); out.ac = out2(1,3);
out.ba = out2(2,1); out.bb = out2(2,2); out.bc = out2(2,3);
out.ca = out2(3,1); out.cb = out2(3,2); out.cc = out2(3,3);

out.hh = -sum(out2(out2>0).*log(out2(out2>0)))/2; % entropy of this result


%% 3
% Make sure ranges are valid for looking at the next one
for i=1:3
	for j=1:3
		if ~isempty(r2{i,j}) && r2{i,j}(end) == N-1; r2{i,j}=r2{i,j}(1:end-1); end
	end
end

% Do the calculation
r3 = cell(3,3,3); out3 = zeros(3,3,3);
for i = 1:3
	for j = 1:3
		for k = 1:3
			r3{i,j,k} = r2{i,j}(yt(r2{i,j}+2) == k);
			out3(i,j,k) = length(r3{i,j,k})/(N-2);
		end
	end
end

% ------ Record these -------
out.aaa = out3(1,1,1); out.aab = out3(1,1,2); out.aac = out3(1,1,3);
out.aba = out3(1,2,1); out.abb = out3(1,2,2); out.abc = out3(1,2,3);
out.aca = out3(1,3,1); out.acb = out3(1,3,2); out.acc = out3(1,3,3);

out.baa = out3(2,1,1); out.bab = out3(2,1,2); out.bac = out3(2,1,3);
out.bba = out3(2,2,1); out.bbb = out3(2,2,2); out.bbc = out3(2,2,3);
out.bca = out3(2,3,1); out.bcb = out3(2,3,2); out.bcc = out3(2,3,3);

out.caa = out3(3,1,1); out.cab = out3(3,1,2); out.cac = out3(3,1,3);
out.cba = out3(3,2,1); out.cbb = out3(3,2,2); out.cbc = out3(3,2,3);
out.cca = out3(3,3,1); out.ccb = out3(3,3,2); out.ccc = out3(3,3,3);

out.hhh = -sum(out3(out3>0).*log(out3(out3>0)))/3; % entropy of this result

%% 4
% Make sure ranges are valid for looking at the next one
for i = 1:3
	for j = 1:3
		for k = 1:3
			if (~isempty(r3{i,j,k}) && r3{i,j,k}(end) == N-2)
                r3{i,j,k}=r3{i,j,k}(1:end-1);
            end
		end
	end
end

% Do the calculation
r4 = cell(3,3,3,3); out4 = zeros(3,3,3,3);
for i=1:3
	for j=1:3
		for k=1:3
			for l=1:3
				r4{i,j,k,l} = r3{i,j,k}(yt(r3{i,j,k}+3) == l);
				out4(i,j,k,l) = length(r4{i,j,k,l})/(N-3);
			end
		end
	end
end

% ------ Record these -------
out.aaaa = out4(1,1,1,1); out.aaab = out4(1,1,1,2); out.aaac = out4(1,1,1,3);
out.aaba = out4(1,1,2,1); out.aabb = out4(1,1,2,2); out.aabc = out4(1,1,2,3);
out.aaca = out4(1,1,3,1); out.aacb = out4(1,1,3,2); out.aacc = out4(1,1,3,3);

out.abaa = out4(1,2,1,1); out.abab = out4(1,2,1,2); out.abac = out4(1,2,1,3);
out.abba = out4(1,2,2,1); out.abbb = out4(1,2,2,2); out.abbc = out4(1,2,2,3);
out.abca = out4(1,2,3,1); out.abcb = out4(1,2,3,2); out.abcc = out4(1,2,3,3);

out.acaa = out4(1,3,1,1); out.acab = out4(1,3,1,2); out.acac = out4(1,3,1,3);
out.acba = out4(1,3,2,1); out.acbb = out4(1,3,2,2); out.acbc = out4(1,3,2,3);
out.acca = out4(1,3,3,1); out.accb = out4(1,3,3,2); out.accc = out4(1,3,3,3);


out.baaa = out4(2,1,1,1); out.baab = out4(2,1,1,2); out.baac = out4(2,1,1,3);
out.baba = out4(2,1,2,1); out.babb = out4(2,1,2,2); out.babc = out4(2,1,2,3);
out.baca = out4(2,1,3,1); out.bacb = out4(2,1,3,2); out.bacc = out4(2,1,3,3);

out.bbaa = out4(2,2,1,1); out.bbab = out4(2,2,1,2); out.bbac = out4(2,2,1,3);
out.bbba = out4(2,2,2,1); out.bbbb = out4(2,2,2,2); out.bbbc = out4(2,2,2,3);
out.bbca = out4(2,2,3,1); out.bbcb = out4(2,2,3,2); out.bbcc = out4(2,2,3,3);

out.bcaa = out4(2,3,1,1); out.bcab = out4(2,3,1,2); out.bcac = out4(2,3,1,3);
out.bcba = out4(2,3,2,1); out.bcbb = out4(2,3,2,2); out.bcbc = out4(2,3,2,3);
out.bcca = out4(2,3,3,1); out.bccb = out4(2,3,3,2); out.bccc = out4(2,3,3,3);


out.caaa = out4(3,1,1,1); out.caab = out4(3,1,1,2); out.caac = out4(3,1,1,3);
out.caba = out4(3,1,2,1); out.cabb = out4(3,1,2,2); out.cabc = out4(3,1,2,3);
out.caca = out4(3,1,3,1); out.cacb = out4(3,1,3,2); out.cacc = out4(3,1,3,3);

out.cbaa = out4(3,2,1,1); out.cbab = out4(3,2,1,2); out.cbac = out4(3,2,1,3);
out.cbba = out4(3,2,2,1); out.cbbb = out4(3,2,2,2); out.cbbc = out4(3,2,2,3);
out.cbca = out4(3,2,3,1); out.cbcb = out4(3,2,3,2); out.cbcc = out4(3,2,3,3);

out.ccaa = out4(3,3,1,1); out.ccab = out4(3,3,1,2); out.ccac = out4(3,3,1,3);
out.ccba = out4(3,3,2,1); out.ccbb = out4(3,3,2,2); out.ccbc = out4(3,3,2,3);
out.ccca = out4(3,3,3,1); out.cccb = out4(3,3,3,2); out.cccc = out4(3,3,3,3);

out.hhhh = -sum(out4(out4>0).*log(out4(out4>0)))/4; % entropy of this result



% %% 4
% % Make sure ranges are valid for looking at the next one
% if ~isempty(r111) && r111(end) == N-2; r111=r111(2:end-2); end
% if ~isempty(r112) && r112(end) == N-2; r112=r112(2:end-2); end
% if ~isempty(r121) && r121(end) == N-2; r121=r121(2:end-2); end
% if ~isempty(r122) && r122(end) == N-2; r122=r122(2:end-2); end
% if ~isempty(r211) && r211(end) == N-2; r211=r211(2:end-2); end
% if ~isempty(r212) && r212(end) == N-2; r212=r212(2:end-2); end
% if ~isempty(r221) && r221(end) == N-2; r221=r221(2:end-2); end
% if ~isempty(r222) && r222(end) == N-2; r222=r222(2:end-2); end
% 
% 
% 
% % 1111
% r1111 = r111(yt(r111+3) == 1);
% % 1112
% r1112 = r111(yt(r111+3) == 2);
% % 1121
% r1121 = r112(yt(r112+3) == 1);
% % 1122
% r1122 = r112(yt(r112+3) == 2);
% % 1211
% r1211 = r121(yt(r121+3) == 1);
% % 1212
% r1212 = r121(yt(r121+3) == 2);
% % 1221
% r1221 = r122(yt(r122+3) == 1);
% % 1222
% r1222 = r122(yt(r122+3) == 2);
% % 2111
% r2111 = r211(yt(r211+3) == 1);
% % 2112
% r2112 = r211(yt(r211+3) == 2);
% % 2121
% r2121 = r212(yt(r212+3) == 1);
% % 2122
% r2122 = r212(yt(r212+3) == 2);
% % 2211
% r2211 = r221(yt(r221+3) == 1);
% % 2212
% r2212 = r221(yt(r221+3) == 2);
% % 2221
% r2221 = r222(yt(r222+3) == 1);
% % 2222
% r2222 = r222(yt(r222+3) == 2);
% 
% 
% % ----- Record these -----
% out.dddd = length(r1111)/(N-3);
% out.dddu = length(r1112)/(N-3);
% out.ddud = length(r1121)/(N-3);
% out.dduu = length(r1122)/(N-3);
% out.dudd = length(r1211)/(N-3);
% out.dudu = length(r1212)/(N-3);
% out.duud = length(r1221)/(N-3);
% out.duuu = length(r1222)/(N-3);
% out.uddd = length(r2111)/(N-3);
% out.uddu = length(r2112)/(N-3);
% out.udud = length(r2121)/(N-3);
% out.uduu = length(r2122)/(N-3);
% out.uudd = length(r2211)/(N-3);
% out.uudu = length(r2212)/(N-3);
% out.uuud = length(r2221)/(N-3);
% out.uuuu = length(r2222)/(N-3);
% out.hhhh = -(out.dddd*log(out.dddd) + out.dddu*log(out.dddu) + out.ddud*log(out.ddud) + out.dduu*log(out.dduu) + ...
% 			out.dudd*log(out.dudd) + out.dudu*log(out.dudu) + out.duud*log(out.duud) + out.duuu*log(out.duuu) +  ...
% 			out.uddd*log(out.uddd) + out.uddu*log(out.uddu) + out.udud*log(out.udud) + out.uduu*log(out.uduu) + ...
% 			out.uudd*log(out.uudd) + out.uudu*log(out.uudu) + out.uuud*log(out.uuud) + out.uuuu*log(out.uuuu) )/4;
% 


end
function [x_confinated]=FcnUsefullData(varargin) %->counts_spectrum,thr_sx,thr_dx

[~,args,nargs] = axescheck(varargin{:});

counts_spettro = args{1};

switch nargs
    case 0
        disp('ERROR: Missing Data')
        return
    case 1
        thr_sx = 0.02;
        thr_dx = 0.98;
    case 2
        thr_sx = args{2};
        thr_dx = 1-thr_sx;
    case 3
        thr_sx = args{2};
        thr_dx = args{3};
end

if isrow(counts_spettro)
    counts_spettro = counts_spettro';
end

CSpettro = cumsum(counts_spettro);
CSpettro = CSpettro / CSpettro(end);

x_confinated = dsearchn(CSpettro, thr_sx) : dsearchn(CSpettro, thr_dx);

end


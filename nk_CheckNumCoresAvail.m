function ncores = nk_CheckNumCoresAvail
global PARALLELENABLE

if isempty(PARALLELENABLE), PARALLELENABLE = false; end

ncores = 1;

a = str2double(getenv('OMP_NUM_THREADS'));
if ~isnan(a) && PARALLELENABLE, ncores = a; end
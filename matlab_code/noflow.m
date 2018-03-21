function initial = noflow(phiu, phiq, param)
%  Define the initial condition for the model with zero flow and specified
%  temperatures
%
%    Efficient Model version,  August 2008

%  Setup the Q grids with zero flow:   q = f 
initial.Qu = 2*param.f0/sqrt(2)*sin(phiq);
initial.Ql = initial.Qu;
%  Set the initial mean global temperature
initial.Tave = 11;
initial.Ts = 25-30*2/pi*phiq;
initial.Ta = initial.Tave*ones(size(phiq));
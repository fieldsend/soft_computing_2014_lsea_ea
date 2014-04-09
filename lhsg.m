function X = lhsg( n, params )

% default latin hypercube sampling -- requires matlab stats toolbox

max_v = params.maximum_values;
min_v = params.minimum_values;
d =length(min_v);
X =lhsdesign(n,d);

mx=repmat(max_v,n,1);
mn=repmat(min_v,n,1);
X=X.*(mx-mn)+mn;   % rescale samples to input range of objective function
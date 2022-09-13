function [ u, info ] = schss68_bamphi( u, k, t, A, g, opts, info )
% schss68_bamphi:
%
% The strategy we adopt to run this exponential integrator in combination with
% bamphi consists in running the (IOM)-Arnoldi procedure at the first call of
% bamphi and then recycle the information gathered in info throughout all the
% exponential integration steps. As an alternative one could delete info after
% each call to bamphi and compute (IOM)-Arnoldi from scratch at each call. 
%
  mu = g( NaN, NaN );

  u = exp( ( - 1i * mu * k / 2 ) * ( u .* conj( u ) ) ) .* u;
  [ u, info ] = bamphi( 1i * k, A, [], u, opts, info );
  u = exp( ( - 1i * mu * k / 2 ) * ( u .* conj( u ) ) ) .* u;

end

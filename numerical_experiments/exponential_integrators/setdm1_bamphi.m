function [ u, info ] = setdm1_bamphi( u, k, t, A, nonlin, opts, info, B )
% setdm1_bamphi
%
% The strategy we adopt to run this exponential integrator in combination with
% bamphi consists in running the (IOM)-Arnoldi procedure at the first call of
% bamphi and then recycle the information gathered in info throughout all the
% exponential integration steps. As an alternative one could delete info after
% each call to bamphi and compute (IOM)-Arnoldi from scratch at each call.
%
% The real trick is applied externally: one could loop over the L = 100 Monte
% Carlo run or rather form a huge matrix as described in the paper. The second
% approach is rewarding even though for growing L it might finally cause memory
% limitations to manifest. To keep this from happening, with bamphi (and bamphi
% only) one can pass B to be the repeated block in the larger Kronecker matrix.
% In doing so we run (IOM)-Arnoldi on a smaller matrix and the allocated memory
% is independent on the number of Monte Carlo launches L.
%
  persistent store_info

  if length( store_info ) > 0
    info = store_info{ 1 };
  else
    if ( nargin == 8 )
      % aux = bamphi_fov( size( u,1 ), B ); % B is a matrix we rather use for computing A's FoV
      % aux = bamphi_fov( size( u,1 ), B, [], [], u(:,1) ); % B is a matrix we rather use for computing A's FoV
      aux = bamphi_fov( size( u,1 ), B, [], [], sum( u,2 ) ); % B is a matrix we rather use for computing A's FoV
      info = aux;
    end
  end

  step = round( t / k ) + 1;
  dW = memo_stoch( step );

  [ du, info ]  = bamphi( k, A, [], [ zeros( numel( u ),1 ),  A( u(:) + nonlin.w( t,u(:) ) .* dW(:) ) + nonlin.g( t, u(:) ) ], opts, info );
  store_info{ 1 } = info;
  u = u + reshape( du + nonlin.w( t,u(:) ) .* dW(:), size( u ) );

end

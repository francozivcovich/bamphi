function [ u, info ] = exptg2s1_bamphi( u, k, t, A, g, opts, info, B, iRtfun )
% exptg2s1_bamphi
%
% The strategy we could adopt to run this exponential integrator in combination
% with bamphi consists in running the (IOM)-Arnoldi procedure at the first call
% of bamphi and then recycle the information gathered in info throughout all the
% exponential integration steps. As an alternative one could delete info after
% each call to bamphi and compute (IOM)-Arnoldi from scratch at each call.
%
% Instead, we give in input the matrix B (in the paper is Lhat), which is
% antisymmetric (thus normal) even though much less sparse thatn A. With B and
% iRtfun which is such that iRtfun(x) = R'\x we compute information on A, then
% we run computations using A.
%
% The user is ionvited to try and change test.mode in main.m inside the
% "if strcmp( which_tests{ test_idx }, 'skg' )"" statement) to see the efficiency 
% gap with a normal launch of this integrator.
%
  persistent store_info

  if length( store_info ) > 0
    info = store_info{ 1 };
  else
    if ( nargin >= 8 )
      % aux = bamphi_fov( length( u ), B ); % B is a matrix we rather use for computing A's FoV <- this is valid too
      v = [ u( 1 : length(u)/2 ); iRtfun( u(length(u)/2+1 : end ) ) ];
      aux = bamphi_fov( length( u ), B, [],[], v ); % B and v are a matrix and vector we rather use for computing A's FoV
      info = aux;
      info.A.pts.rho(:) = 1i * imag( info.A.pts.rho(:) );
    end
  end

  [ u, info ] = bamphi( k, A, [], u + ( k / 2 ) * g( t, u ), opts, info );
  store_info{ 1 } = info;

  u = u + ( k / 2 ) * g( t, u );


end

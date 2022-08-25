function out = memo_stoch( in )

  persistent dW

  if ( nargout == 0 )
    dW = in;
  end
  if ( nargout == 1 )
    out = dW{ in };
  end

end

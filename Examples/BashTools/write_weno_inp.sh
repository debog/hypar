FILE="weno.inp"

/bin/cat <<EOM > $FILE
begin
  mapped      $p_mapped
  borges      $p_borges
  yc          $p_yc
  no_limiting $p_nl
  epsilon     $p_eps
  p           $p_p
  rc          $p_rc
  xi          $p_xi
  tol         $p_tol
end
EOM

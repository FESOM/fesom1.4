subroutine dist_on_earth(lon1, lat1, lon2, lat2, dist)
  ! distance on the earth between two points
  ! input: lon1 lat2 and lon2 lat2 in radian
  ! output: dist in m
  use o_param
  implicit none
  real(kind=8)  :: lon1, lat1, lon2, lat2, alpha, dist

  alpha=acos(cos(lat1)*cos(lat2)*cos(lon1-lon2)+sin(lat1)*sin(lat2))
  dist=r_earth*abs(alpha)

end subroutine dist_on_earth

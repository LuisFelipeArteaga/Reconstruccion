function q = Matrix_to_quaternion(mT)
tr = trace(mT);
m00 = mT(1,1); m01 = mT(1,2); m02 = mT(1,3);
m10 = mT(2,1); m11 = mT(2,2); m12 = mT(2,3);
m20 = mT(3,1); m21 = mT(3,2); m22 = mT(3,3);
if tr > 0 
  S = sqrt(tr+1)*2;
  qw = 0.25*S;
  qx = (m21 - m12)/S;
  qy = (m02 - m20)/S; 
  qz = (m10 - m01)/S; 
elseif ((m00 > m11) && (m00 > m22)) 
  S = sqrt(1.0 + m00 - m11 - m22)*2; 
  qw = (m21 - m12)/S;
  qx = 0.25*S;
  qy = (m01 + m10)/S; 
  qz = (m02 + m20)/S; 
elseif (m11 > m22)  
  S = sqrt(1.0+m11-m00-m22)*2; 
  qw = (m02 - m20)/S;
  qx = (m01 + m10)/S; 
  qy = 0.25 * S;
  qz = (m12 + m21)/S; 
else 
  S = sqrt(1.0+m22-m00-m11)*2; 
  qw = (m10 - m01)/S;
  qx = (m02 + m20)/S;
  qy = (m12 + m21)/S;
  qz = 0.25*S;
end
q = [qw,qx,qy,qz];
end
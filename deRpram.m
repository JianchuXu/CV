function R=deRpram(omega)
theta=norm(omega);
R=cos(theta)*eye(3)+csinc(theta)*skew(omega)+...
    (1-cos(theta))/theta^2*(omega*omega');
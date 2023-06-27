function val=q_tilde(x)
q0=evalin('base', 'q0');
val=q(x)-q0;
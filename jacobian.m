    syms phi theta si
    syms wx wy wz real
    syms bgx bgy bgz
    syms Ts
    
    x_sym = [phi; theta; si; bgx; bgy; bgz];
    Wbi_b_sym = [wx; wy; wz];

    dx1 = state_derivative(x_sym, Wbi_b_sym) * Ts;
    
    x_next = x_sym + dx1;
    
    % state Jacobian
    Fk_sym = jacobian(x_next, x_sym);
    Lk_sym = jacobian(x_next, Wbi_b_sym);

    % Create function file for numeric Jacobian evaluation
    matlabFunction(Fk_sym, 'File', 'Fk', ...
        'Vars', {x_sym, Wbi_b_sym, Ts}); 
    matlabFunction(Lk_sym, 'File', 'Lk', ...
        'Vars', {x_sym, Wbi_b_sym, Ts}); 

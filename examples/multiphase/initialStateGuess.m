function states = initialStateGuess(obj)
    
    N = obj.ctrl.N;
    nx = obj.nx;
    
    x0 = obj.auxdata.x0;
    xf = obj.auxdata.xf;
    xf_idx = obj.auxdata.xf_idx;
    md = obj.bnds.x_min(1);

    states = zeros(nx, N);
    states(:,1) = x0;
    states(xf_idx,N) = xf;
    states(1,N) = md;
    for id = 2:N-1
        states(xf_idx,id) = interp1([0 obj.time(end)], [x0(xf_idx),xf].', obj.time(id)).';
        states(1,id) = interp1([0 obj.time(end)], [x0(1);md].', obj.time(id)).';
    end
    
end
function err = score(Iref,Vm,Para,States,Time)
    err = 0;
    Iks = mm(Vm,Para,States,Time);
    %calculate Imax for each voltage protocol
    Imax = max(abs(Iref));
    for n = 1:length(Time);
        d = ((Iks(n)-Iref(n))/Imax)^2;
        err = err + sqrt(d);
    end
end
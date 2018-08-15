function Vj_extremum = BTiterative_LP_for_inj_env (LP_prob,Vj_min,LP_opt)
%returns max(Vj_imag) when Vi_imag=0. Iteratively solves LP with
%constraint l<Vj_real<u

Vj_max=LP_prob.bu;
for iter=1:LP_opt.max_iter
    if (LP_opt.use_mex)
        [Vj_extremum, info]=greedy_LP_solution_mex( LP_prob.c, LP_prob.a, LP_prob.xl, LP_prob.xu, LP_prob.bl, LP_prob.bu );
    else
        [Vj_extremum, info]=greedy_LP_solution( LP_prob.c, LP_prob.a, LP_prob.xl, LP_prob.xu, LP_prob.bl, LP_prob.bu );
    end

    if (Vj_extremum>=0)
        b_temp=Vj_min^2-Vj_extremum^2;
        if (info~=2 && b_temp>=0 && sqrt(b_temp)-LP_prob.bl>LP_opt.eps_tol)
            LP_prob.bl=sqrt(b_temp);
        else
            break;
        end
    else
        b_temp=Vj_max^2-Vj_extremum^2;
        if (info~=1 && b_temp>=0 && LP_prob.bu-sqrt(b_temp)>LP_opt.eps_tol)
            LP_prob.bu=sqrt(b_temp);
        else
            break;
        end
    end
end

end
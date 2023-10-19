function C=copula_distribution(u,v,COPULA)

if strcmpi(COPULA.name,'Joe')
    C=copula_distribution_Joe(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'FGM')
    C=copula_distribution_FGM(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Frank')
    C=copula_distribution_Frank(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'MEV_symlog')
    C=copula_distribution_MEV_symlog(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'MEV_asymlog')
    C=copula_distribution_MEV_asymlog(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Clayton')
    C=copula_distribution_Clayton(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Nelsen16')
    C=copula_distribution_Nelsen16(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Nelsen17')
    C=copula_distribution_Nelsen17(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'HuslerReiss')
    C=copula_distribution_MEV_HuslerReiss(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Dirichlet')
    C=copula_distribution_MEV_dirichlet(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'BB9')
    C=copula_distribution_BB9(u,v,COPULA.params);
end

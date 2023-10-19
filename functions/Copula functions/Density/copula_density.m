function c=copula_density(u,v,COPULA)

if isfield(COPULA,'corner')
    u0=COPULA.corner(1);
    v0=COPULA.corner(2);
    u = u0 + (-1)^u0 * u;
    v = v0 + (-1)^v0 * v;
end

if strcmpi(COPULA.name,'indep')
    c=1;

    % Extreme value copulas    
elseif strcmpi(COPULA.name,'MEV_symlog')
    c=copula_density_MEV_symlog(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'MEV_asymlog')
    c=copula_density_MEV_asymlog(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'MEV_polynomial_pickands')
    c=copula_density_MEV_polynomial_pickands(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'HuslerReiss')
    c=copula_density_MEV_HuslerReiss(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Dirichlet')
    c=copula_density_MEV_dirichlet(u,v,COPULA.params);

    % ellipitical copulas
elseif strcmpi(COPULA.name,'Gaussian')
    c=copula_density_Gaussian(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Student')
    c=copula_density_student(u,v,COPULA.params);

    % Nelsen Archimedean Copulas
elseif strcmpi(COPULA.name,'Nelsen12')
    c=copula_density_Nelsen12(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Nelsen15')
    c=copula_density_Nelsen15(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Nelsen16')
    c=copula_density_Nelsen16(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Nelsen17')
    c=copula_density_Nelsen17(u,v,COPULA.params);

    % Others
elseif strcmpi(COPULA.name,'Joe')
    c=copula_density_Joe(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'FGM')
    c=copula_density_FGM(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Frank')
    c=copula_density_Frank(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Clayton')
    c=copula_density_Clayton(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'AMH')
    c=copula_density_AMH(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'Gumbel')
    c=copula_density_Gumbel(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'bivariate_exponential')
    c=copula_density_bivariate_exponential(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'plackett')
    c=copula_density_Plackett(u,v,COPULA.params);
elseif strcmpi(COPULA.name,'BB9')
    c=copula_density_BB9(u,v,COPULA.params);
end

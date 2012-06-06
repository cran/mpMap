\name{mpMap-internal}
\Rdversion{1.1}
\alias{check_qtl}
\alias{check_ped}
\alias{combine_chr}
\alias{combine_ibd}
\alias{combine_rf}
\alias{compute_bnrf}
\alias{convertped}
\alias{cleanrf}
\alias{CR_cross}
\alias{CR_gen_geno}
\alias{CR_estrf}
\alias{CR_calcLD}
\alias{fill}
\alias{generate_error}
\alias{generate_obs}
\alias{generate_pheno}
\alias{getpval}
\alias{reorgRIgenoprob}
\alias{wald.test.asreml}
\alias{cintern}
\alias{scan3pt}

\title{ Internal mpMap functions} 
\description{
Internal MPmap functions. These are generally not to be called by the user.}
\usage{
check_qtl(QTL, map, n.founders)
check_ped(pedigree)
combine_chr(map, map.function)
combine_ibd(obsgeno)
combine_rf(mpcross1, mpcross2, r)
compute_bnrf(mpcross1, mpcross2, r, grid)
convertped(pedigree)
CR_cross(mpcross)
CR_gen_geno(map, map.function, ped, seed, transpos, transval)
CR_estrf(finals, founders, pedigree, pairs, r)
CR_calcLD(finals, founders, pedigree, pairs, rmat)
cleanrf(mpcross)
fill(mat, missfx=2, ...)
generate_error(geno, error.prob)
generate_obs(geno, map, full.prob, fg, transpos, transval, founderld)
generate_pheno(n.founders, qtlgeno, qtleffects, vare, n.ind)
getpval(model)
reorgRIgenoprob <- function(cross)
scan3pt <- function(mpcross, chr, step=0, mapfx=c("haldane", "kosambi"))
wald.test.asreml(object, cc)
cintern(cc, tau, vrb, sigma2)
}
\keyword{internal}

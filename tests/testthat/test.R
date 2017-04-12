
sim.norm <- function(){
    n <- 100
    z <- rep(c(0,1),each=(n/2))
    y <- stats::rnorm(n=n,mean=2,sd=1)
    z[(n/4):n] <- NA
    list(y=y,z=z)
}

sim.nbinom <- function(){
    n <- 100
    z <- rep(c(0,1),each=(n/2))
    y <- stats::rnbinom(n=n,mu=1+2,size=1/0.05)
    y <- y + stats::runif(n=n,min=0,max=0.001)
    gamma <- stats::runif(n=n,min=0,max=2)
    z[(n/4):n] <- NA
    list(y=y,z=z,gamma=gamma)
}

sim.zinb <- function(){
    n <- 100
    z <- rep(c(0,1),each=(n/2))
    y <- stats::rnbinom(n=n,mu=1+2,size=1/0.05)
    y[sample(1:n,size=0.2*n)] <- 0
    y <- y + stats::runif(n=n,min=0,max=0.001)
    gamma <- stats::runif(n=n,min=0,max=2)
    z[(n/4):n] <- NA
    list(y=y,z=z,gamma=gamma)
}

testthat::test_that("wrapping \"norm\"",{
    set.seed(1)
    sim <- sim.norm()
    y <- sim$y; z <- sim$z
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z,dist="norm")
    set.seed(1)
    b <- semisup::fit.wrap(y=y,z=z,dist="norm",phi=NULL,pi=NULL,gamma=NULL)
    set.seed(1)
    c <- semisup::fit.norm(y=y,z=z,it.em=100,epsilon=1e-04)
    testthat::expect_identical(a,b)
    testthat::expect_identical(a,c)
})

testthat::test_that("wrapping \"nbinom\"",{
    sim <- sim.nbinom()
    y <- sim$y; z <- sim$z; gamma <- sim$gamma
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z,dist="nbinom",phi=0.05,gamma=gamma)
    set.seed(1)
    b <- semisup::fit.wrap(y=y,z=z,dist="nbinom",phi=0.05,pi=NULL,gamma=gamma)
    set.seed(1)
    c <- semisup::fit.nbinom(y=y,z=z,phi=0.05,gamma=gamma,it.em=100,epsilon=1e-04)
    testthat::expect_identical(a,b)
    testthat::expect_identical(a,c)
})

testthat::test_that("wrapping \"zinb\"",{
    sim <- sim.zinb()
    y <- sim$y; z <- sim$z; gamma <- sim$gamma
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z,dist="zinb",phi=0.05,pi=0.2,gamma=gamma)
    set.seed(1)
    b <- semisup::fit.wrap(y=y,z=z,dist="zinb",phi=0.05,pi=0.2,gamma=gamma)
    set.seed(1)
    c <- semisup::fit.zinb(y=y,z=z,phi=0.05,pi=0.2,gamma=gamma,it.em=100,epsilon=1e-04)
    testthat::expect_identical(a,b)
    testthat::expect_identical(a,c)
})

testthat::test_that("mean estimation",{
     y <- stats::runif(n=100,min=0,max=100)
     z <- rep(c(0,NA),each=50)
     a <- semisup::mixtura(y=y,z=z,dist="nbinom")
     b <- semisup::mixtura(y=y,z=z,dist="norm")
     testthat::expect_equal(a$estim0$mu0,b$estim0$mean0,mean(y))
})

testthat::test_that("default offset",{
    sim <- sim.nbinom()
    y <- sim$y; z <- sim$z; gamma <- sim$gamma
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z,dist="nbinom")
    set.seed(1)
    gamma <- rep(1,times=length(y))
    b <- semisup::mixtura(y=y,z=z,dist="nbinom",gamma=gamma)
    testthat::expect_identical(a,b)
})

testthat::test_that("default dispersion",{
    sim <- sim.nbinom()
    y <- sim$y; z <- sim$z; gamma <- sim$gamma
    a <- semisup::mixtura(y=y,z=z,gamma=gamma,dist="nbinom")$estim0$phi
    b <- semisup::estim.nbinom(y=y,z=z,gamma=gamma)$phi
    testthat::expect_identical(a,b)
})

testthat::test_that("default inflation",{
    sim <- sim.nbinom()
    y <- sim$y; z <- sim$z; gamma <- sim$gamma
    temp <- semisup::mixtura(y=y,z=z,gamma=gamma,dist="zinb")
    a <- list(phi=temp$estim0$phi,pi=temp$estim0$pi)
    temp <- semisup::estim.zinb(y=round(y),z=z,gamma=gamma)
    b <- list(phi=temp$phi,pi=temp$pi)
    testthat::expect_identical(a,b)
})


testthat::test_that("interrupt \"norm\"",{
    sim <- sim.norm()
    y <- sim$y; z <- sim$z
    kind <- 0.9
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z,kind=1,test="boot")$p.value
    set.seed(1)
    b <- semisup::mixtura(y=y,z=z,kind=kind,test="boot")$p.value
    if(a > kind){
        b <- a
    }
    testthat::expect_identical(a,b)
})

testthat::test_that("interrupt \"nbinom\"",{
    sim <- sim.nbinom()
    y <- sim$y; z <- sim$z; gamma <- sim$gamma
    kind <- 0.9
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z,dist="nbinom",kind=1,test="perm")$p.value
    set.seed(1)
    b <- semisup::mixtura(y=y,z=z,dist="nbinom",kind=kind,test="perm")$p.value
    if(a > kind){
        b <- a
    }
    testthat::expect_identical(a,b)
})

testthat::test_that("interrupt \"zinb\"",{
    sim <- sim.zinb()
    y <- sim$y; z <- sim$z; gamma <- sim$gamma
    kind <- 0.9
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z,dist="zinb",kind=1,test="perm")$p.value
    set.seed(1)
    b <- semisup::mixtura(y=y,z=z,dist="zinb",kind=kind,test="perm")$p.value
    if(a > kind){
        b <- a
    }
    testthat::expect_identical(a,b)
})


testthat::test_that("multiple starts \"norm\"",{
    sim <- sim.norm()
    y <- sim$y; z <- sim$z
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z)$lrts
    set.seed(1)
    b <- semisup::mixtura(y=y,z=z,starts=4)$lrts
    testthat::expect_gte(object=b,expected=a)
})

testthat::test_that("multiple starts \"nbinom\"",{
    sim <- sim.nbinom()
    y <- sim$y; z <- sim$z
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z,dist="nbinom")$lrts
    set.seed(1)
    b <- semisup::mixtura(y=y,z=z,starts=4,dist="nbinom")$lrts
    testthat::expect_gte(object=b,expected=a)
})

testthat::test_that("multiple starts \"zinb\"",{
    sim <- sim.zinb()
    y <- sim$y; z <- sim$z
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z,dist="zinb")$lrts
    set.seed(1)
    b <- semisup::mixtura(y=y,z=z,starts=4,dist="zinb")$lrts
    testthat::expect_gte(object=b,expected=a)
})

testthat::test_that("compare \"norm\"",{
    sim <- sim.norm()
    y <- sim$y; z <- sim$z
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z)$lrts
    set.seed(1)
    b <- semisup::scrutor(Y=y,Z=z)$lrts
    testthat::expect_identical(object=b,expected=a)
})

testthat::test_that("compare \"nbinom\"",{
    sim <- sim.nbinom()
    y <- sim$y; z <- sim$z; gamma <- sim$gamma
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z,dist="nbinom",gamma=gamma)$lrts
    set.seed(1)
    b <- semisup::scrutor(Y=y,Z=z,dist="nbinom",gamma=gamma)$lrts
    testthat::expect_identical(object=b,expected=a)
})

testthat::test_that("compare \"zinb\"",{
    sim <- sim.zinb()
    y <- sim$y; z <- sim$z; gamma <- sim$gamma
    set.seed(1)
    a <- semisup::mixtura(y=y,z=z,dist="zinb",gamma=gamma)$lrts
    set.seed(1)
    b <- semisup::scrutor(Y=y,Z=z,dist="zinb",gamma=gamma)$lrts
    testthat::expect_identical(object=b,expected=a)
})



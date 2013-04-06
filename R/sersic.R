sersic = function(mag, re, n, e=0, r=re){
    bn = qgamma(0.5,2*n)
    lumtot = 1*(re^2)*2*pi*n*((exp(bn))/(bn^(2*n)))*gamma(2*n)*(1-e)
    magtot = -2.5*log10(lumtot)
    Ie = 1/(10^(0.4*(mag-magtot)))
    x = bn*(r/re)^(1/n)
    lumr = Ie*lumtot*pgamma(x,2*n)
    a = r
    b = a*(1-e)
    intenr = Ie*exp(-bn*(((r/re)^(1/n))-1))
    lumtot = Ie*lumtot
    magtot = -2.5*log10(lumtot)
    magr = -2.5*log10(lumr)
    mur = -2.5*log10(intenr)
    muavgr = -2.5*log10(lumr/(pi*a*b))
    return(list(mag=magr, magdiff=magtot-magr, mu=mur, muavg=muavgr, inten=intenr, lum=lumr, lumtot=lumtot, lumfrac=lumr/lumtot))
}

convrad = function(n, f = 0.9, r = 1, fr = 0.5){
    # converts sersic radii containing differing amounts of total luminosity
    bn1 = qgamma(fr,2*n)
    bn2 = qgamma(f,2*n)
    r2 = r*((bn2/bn1)^n)
    return(r2)
}

convmu = function(r, mu, n, re = 1, rmu = re){
    # convert between surface brightnesses
    bn = qgamma(0.5,2*n)
    mu2 = mu + ((2.5*bn)/log(10))*(((r/re)^(1/n))-((rmu/re)^(1/n)))
    return(mu2)
}

re2h = function(n, re = 1){
    # convert re to scalelength
    bn = qgamma(0.5,2*n)
    h = re/(bn^n)
    return(h)
}

h2re = function(n, h = 1){
    # convert scalelength to re
    bn = qgamma(0.5,2*n)
    re = h*(bn^n)
    return(re)
}

igamma = function(x, s){
    # TRUE incomplete gamma function
    i = pgamma(x,s)*gamma(s)
    return(i)
}

concen = function(a, n){
    # central concentration
    bn = qgamma(0.5,2*n)
    c = (igamma(bn*(a^(1/n)),2*n))/(igamma(bn,2*n))
    return(c)
}

petro = function(mag, n, e = 0, rp = 3, i = 0.5){
    # petrosian
    temp = function(n, i, re, res=10, loop=5){
        bn = qgamma(0.5,2*n)
        r = seq(0,10*re,len=res)
        lower = 1
        upper = res
        for(j in 1:loop){
            r = seq(r[lower],r[upper],len=res)
            x = bn*((r/re)^(1/n))
            p = suppressWarnings((2*n*igamma(x,2*n))/((exp(1)^(-x))*(x^(2*n))))
            pmin = which.min(abs((1/p)-i))
            lower = max((pmin-1),1)
            upper = min((pmin+1),length(r))
        }
        return=c(r[pmin])
    }
    prad = as.numeric(Vectorize(temp)(n,i,re=1,res=10,loop=10))
    petro = sersic(mag=mag, re=1, n=n, e=e, r=rp*prad)
    return(petro)
}

petroindex = function(r, n, re = 1){
    # petrosian index
    bn = qgamma(0.5,2*n)
    x = bn*((r/re)^(1/n))
    i = (2*n*igamma(x,2*n))/((exp(1)^(-x))*(x^(2*n)))
    return(i)
}

petrorad = function(n, i = 0.5, re = 1){
    # calculate petrosian radius as multiples of re
    temp = function(n, i, re, res=10, loop=5){
        bn = qgamma(0.5,2*n)
        r = seq(0,10*re,len=res)
        lower = 1
        upper = res
        for(j in 1:loop){
            r = seq(r[lower],r[upper],len=res)
            x = bn*((r/re)^(1/n))
            p = suppressWarnings((2*n*igamma(x,2*n))/((exp(1)^(-x))*(x^(2*n))))
            pmin = which.min(abs((1/p)-i))
            lower = max((pmin-1),1)
            upper = min((pmin+1),length(r))
        }
        return=c(r[pmin])
    }
    prad = as.numeric(Vectorize(temp)(n,i,re=1,res=10,loop=10))
    return(prad)
}

kron = function(mag, n, e = 0, rk = 2.5, r=1e10, re = 1){
    # kron
    bn = qgamma(0.5,2*n)
    x = bn*((r/re)^(1/n))
    krad = suppressWarnings((re/(bn^n))*(igamma(x,3*n)/igamma(x,2*n)))
    kron = sersic(mag=mag, re=re, n=n, e=e, r=rk*krad)
    return(kron)
}

kronrad = function(n, r = 1e10, re = 1){
    # kron radii
    bn = qgamma(0.5,2*n)
    x = bn*((r/re)^(1/n))
    krad = suppressWarnings((re/(bn^n))*(igamma(x,3*n)/igamma(x,2*n)))
    return(krad)
}





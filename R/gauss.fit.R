gauss.fit = function(data, mean = 0, sigma = 1, norm = 1, mean.alt = NA, sigma.alt = NA, norm.alt = NA, meanlo = -Inf, sigmalo = 1e-9, normlo = 0, meanhi = Inf, sigmahi = Inf, normhi = Inf, fixmean = FALSE, fixsigma = FALSE, fixnorm = FALSE, range = NA, lim1 = NA, lim2 = NA, numlim = 1, method = "nlminb", volume = 1, bw = 0.1, ftype = "lin", null = 1E-9, error = "jack", subvol = 10, sampnum = subvol, ...){
    
    #data=rnorm(1e4, mean=2, sd=4); mean = 0; sigma = 1; norm = 1; mean.alt = 10; sigma.alt = 10; norm.alt = 1; meanlo = -Inf; sigmalo = 1e-9; normlo = 0; meanhi = Inf; sigmahi = Inf; normhi = Inf; fixmean = FALSE; fixsigma = FALSE; fixnorm = FALSE; range = range(data); lim1 = NA; lim2 = NA; numlim = 1; method = "nlminb"; volume = 1; bw = 0.1; ftype = "lin"; null = 1E-9; error = "jack"; subvol = 10; sampnum = subvol
    
    # setup
    if(!is.na(range[1])){
        range = c(min(range),max(range))
    }else{
        range = range(data)
    }
    
#    # export
#    .gauss.fit.bin = astro:::.gauss.fit.bin
#    .gauss.fit.dat = astro:::.gauss.fit.dat
#    .gauss.fit.chi = astro:::.gauss.fit.chi
    
    # gauss bin calculations (split data into bins)
    bindat = .gauss.fit.bin(data=data, range=range, lim1=lim1, lim2=lim2, numlim=numlim, volume=volume, bw=bw, null=null)
    
    # only continue if any data points remain
    if(length(bindat$fitden)>0){
        
        fitdat = .gauss.fit.dat(bindat=bindat, mean=mean, sigma=sigma, norm=norm, mean.alt=mean.alt, sigma.alt=sigma.alt, norm.alt=norm.alt, meanlo=meanlo, sigmalo=sigmalo, normlo=normlo, meanhi=meanhi, sigmahi=sigmahi, normhi=normhi, method=method, bw=bw, ftype=ftype, fixmean=fixmean, fixsigma=fixsigma, fixnorm=fixnorm)
        
        # errors?
        if(error=="jack"){
            
            # jack setup
            subvols = round(seq(1,(length(data)+1),len=subvol+1))
            jacklist = 1:subvol
            
            for(l in jacklist){
                
                # setup
                cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",l,"/",(length(subvols)-1),"")
                jackdat = data[-(subvols[l]:subvols[l+1])]
                
                # jackknifed binned data
                jackbin = .gauss.fit.bin(data=jackdat, range=range, lim1=lim1, lim2=lim2, numlim=numlim, volume=volume, bw=bw, null=null)
                
                # jackknifed fitted data
                jackfit = .gauss.fit.dat(bindat=jackbin, mean=mean, sigma=sigma, norm=norm, mean.alt=mean.alt, sigma.alt=sigma.alt, norm.alt=norm.alt, meanlo=meanlo, sigmalo=sigmalo, normlo=normlo, meanhi=meanhi, sigmahi=sigmahi, normhi=normhi, method=method, bw=bw, ftype=ftype, fixmean=fixmean, fixsigma=fixsigma, fixnorm=fixnorm)
                
                # add par to results
                if(l==jacklist[1]){
                    pars = jackfit$par
                    chi2s = jackfit$chi2
                    rchi2s = jackfit$rchi2
                }else{
                    pars = rbind(pars, jackfit$par)
                    chi2s = c(chi2s, jackfit$chi2)
                    rchi2s = c(rchi2s, jackfit$rchi2)
                }
                
            }
            
            # column definitions
            mcols = grep("m",colnames(pars))
            scols = grep("s",colnames(pars))
            ncols = grep("n",colnames(pars))
            
            # phi* rescaling
            pars[,ncols] = pars[,ncols]/((subvol-1)/subvol)
            
            # calculate errors
            sigma = function(N, x2, x){
                sigma = sqrt(((N-1)/N)*sum((x2-x)^2))
                return(sigma)
            }
            meanlo1 = -1*sigma(subvol, pars[,mcols[1]], fitdat$par["m1"])
            meanhi1 = +1*sigma(subvol, pars[,mcols[1]], fitdat$par["m1"])
            sigmalo1 = -1*sigma(subvol, pars[,scols[1]], fitdat$par["s1"])
            sigmahi1 = +1*sigma(subvol, pars[,scols[1]], fitdat$par["s1"])
            normlo1 = -1*sigma(subvol, pars[,ncols[1]], fitdat$par["n1"])
            normhi1 = +1*sigma(subvol, pars[,ncols[1]], fitdat$par["n1"])
            if("m2"%in%names(fitdat$par)){
                meanlo2 = -1*sigma(subvol, pars[,mcols[2]], fitdat$par["m2"])
                meanhi2 = +1*sigma(subvol, pars[,mcols[2]], fitdat$par["m2"])
            }else{
                meanlo2 = meanhi2 = NA
            }
            if("s2"%in%names(fitdat$par)){
                sigmalo2 = -1*sigma(subvol, pars[,scols[2]], fitdat$par["s2"])
                sigmahi2 = +1*sigma(subvol, pars[,scols[2]], fitdat$par["s2"])
            }else{
                sigmalo2 = sigmahi2 = NA
            }
            if("n2"%in%names(fitdat$par)){
                normlo2 = -1*sigma(subvol, pars[,ncols[2]], fitdat$par["n2"])
                normhi2 = +1*sigma(subvol, pars[,ncols[2]], fitdat$par["n2"])
            }else{
                normlo2 = normhi2 = NA
            }
            parlo = c(m1lo=meanlo1, s1lo=sigmalo1, n1lo=normlo1, m2lo=meanlo2, s2lo=sigmalo2, n2lo=normlo2)
            parhi = c(m1hi=meanhi1, s1hi=sigmahi1, n1hi=normhi1, m2hi=meanhi2, s2hi=sigmahi2, n2hi=normhi2)
            if(any(is.na(parlo))){parlo = parlo[-which(is.na(parlo))]}
            if(any(is.na(parhi))){parhi = parhi[-which(is.na(parhi))]}
            
            cat("\n")
            
        }else if(error=="boot"){
            
            NULL
            
        }else{
            
            parlo = parhi = NA
            
        }
        
        # build data list
        dat = list(binmid=bindat$binmid, num=bindat$num, den=bindat$den, err=bindat$err, errlo=bindat$errlo, errhi=bindat$errhi, par=fitdat$par, parlo=parlo, parhi=parhi, chi2=fitdat$chi2, dof=fitdat$dof, rchi2=fitdat$rchi2, denlim=(numlim/volume)/bw, hessian=fitdat$hessian)
        
    }else{
        
        warning("No data available for fitting, check limits")
        dat = list(binmid=NA, num=NA, den=NA, err=NA, errlo=NA, errhi=NA, par=NA, parlo=NA, parhi=NA, chi2=NA, dof=NA, rchi2=NA, denlim=NA, hessian=NA)
        
    }
    
    # return fitting data
    return(dat)
    
}

# prepare and fit gauss function
.gauss.fit.dat = function(bindat, mean, sigma, norm, mean.alt, sigma.alt, norm.alt, meanlo, sigmalo, normlo, meanhi, sigmahi, normhi, method, bw, ftype=ftype, fixmean, fixsigma, fixnorm){
    
    # lengths
    maxlength = max(length(mean), length(sigma), length(norm), length(mean.alt), length(sigma.alt), length(norm.alt))
    
    # fixed vectors
    fixmeans = rep(TRUE,2)
    fixsigmas = rep(TRUE,2)
    fixnorms = rep(TRUE,2)
    fixmeans[1:length(fixmean)] = fixmean
    fixsigmas[1:length(fixsigma)] = fixsigma
    fixnorms[1:length(fixnorm)] = fixnorm
    
    # grids
    blank = NA
    pri = c(m = mean[1:maxlength], s = sigma[1:maxlength], n = norm[1:maxlength])
    sec = c(m = mean.alt[1:maxlength], s = sigma.alt[1:maxlength], n = norm.alt[1:maxlength])
    grid = rbind(blank,pri,sec)
    if( sum(!is.na(sec))==0 ){
        grid = rbind(blank,pri)
    }
    fixvars = grid[1,]
    
    # remove NA columns from grid
    bad = {}
    for(i in 1:length(grid[1,])){
        if(sum(is.na(grid[,i])) == length(grid[,1])){
            bad = c(bad,i)
        }
    }
    if(length(bad)>0){
        grid = grid[,-bad]
    }
    
    k = 2
    
    # optimise
    chis = {}
    pars = {}
    hess = {}
    for(k in 2:length(grid[,1])){
        
        # create input for fit
        vars = grid[k,]
        for(i in 1:length(vars)){
            if(names(vars)[i] %in% names(fixvars)){
                bad = which(names(fixvars)==names(vars)[i])
                fixvars = fixvars[-bad]
            }
        }
        
        # single names fix
        if("m" %in% names(vars)){names(vars)[which(names(vars)=="m")]="m1"}
        if("s" %in% names(vars)){names(vars)[which(names(vars)=="s")]="s1"}
        if("n" %in% names(vars)){names(vars)[which(names(vars)=="n")]="n1"}
        
        # fixed parameters
        if(fixmeans[1]){if("m1"%in%names(vars)){fixvars=c(fixvars,vars["m1"]); vars=vars[-which(names(vars)=="m1")]}}
        if(fixsigmas[1]){if("s1"%in%names(vars)){fixvars=c(fixvars,vars["s1"]); vars=vars[-which(names(vars)=="s1")]}}
        if(fixnorms[1]){if("n1"%in%names(vars)){fixvars=c(fixvars,vars["n1"]); vars=vars[-which(names(vars)=="n1")]}}
        if(fixmeans[2]){if("m2"%in%names(vars)){fixvars=c(fixvars,vars["m2"]); vars=vars[-which(names(vars)=="m2")]}}
        if(fixsigmas[2]){if("s2"%in%names(vars)){fixvars=c(fixvars,vars["s2"]); vars=vars[-which(names(vars)=="s2")]}}
        if(fixnorms[2]){if("n2"%in%names(vars)){fixvars=c(fixvars,vars["n2"]); vars=vars[-which(names(vars)=="n2")]}}
        
        # fitting limits
        lolim = {}
        hilim = {}
        if("m1"%in%names(vars)){lolim = c(lolim, meanlo); hilim=c(hilim, meanhi)}
        if("s1"%in%names(vars)){lolim = c(lolim, sigmalo); hilim=c(hilim, sigmahi)}
        if("n1"%in%names(vars)){lolim = c(lolim, normlo); hilim=c(hilim, normhi)}
        if("m2"%in%names(vars)){lolim = c(lolim, meanlo); hilim=c(hilim, meanhi)}
        if("s2"%in%names(vars)){lolim = c(lolim, sigmalo); hilim=c(hilim, sigmahi)}
        if("n2"%in%names(vars)){lolim = c(lolim, normlo); hilim=c(hilim, normhi)}
        
        # fit the Gaussian function
        if(method!="nlminb"){
            fit = optim(vars, .gauss.fit.chi, data=bindat$fitden, binmid=bindat$fitbinmid, stdev=bindat$fiterr, bw=bw, ftype=ftype, fixvars=fixvars, hessian=TRUE, method=method, lower=lolim, upper=hilim)
        }else{
            fit = nlminb(vars, .gauss.fit.chi, data=bindat$fitden, binmid=bindat$fitbinmid, stdev=bindat$fiterr, bw=bw, ftype=ftype, fixvars=fixvars, lower=lolim, upper=hilim)#, control = list(trace=TRUE))
            fit$value=fit$objective
        }
        
        # add on any fixed parameters to $par
        if(fixmeans[1]){fit$par = c(fit$par,fixvars["k1"])}
        if(fixsigmas[1]){fit$par = c(fit$par,fixvars["s1"])}
        if(fixnorms[1]){fit$par = c(fit$par,fixvars["n1"])}
        if(fixmeans[2]){fit$par = c(fit$par,fixvars["k2"])}
        if(fixsigmas[2]){fit$par = c(fit$par,fixvars["s2"])}
        if(fixnorms[2]){fit$par = c(fit$par,fixvars["n2"])}
        
        # collect results
        chis = c(chis, fit$value)
        pars = c(pars, list(fit$par))
        hess = c(hess, list(fit$hessian))
        
    }
    
    # choose best fit
    chi = chis[[which.min(chis)]]
    par = pars[[which.min(chis)]]
    hessian = hess[[which.min(chis)]]
    chi2 = chi
    dof = (length(bindat$fitden)-length(par))
    rchi2 = chi2/dof
    
    # return results
    return(list(par=par, chi2=chi2, dof=dof, rchi2=rchi2, hessian=hessian))
    
}

# calculate number (total), number density (phi) and poissonian errors
.gauss.fit.bin = function(data, range, lim1, lim2, numlim, volume, bw, null){
    
    # data vectors
    bins = seq(range[1], range[2], by=bw)
    binmid = bins[-1]-(bw/2)
    binlo = bins[1:(length(bins)-1)]
    binhi = bins[2:length(bins)]
    
    # calculate bin sample statistics
    num = {}
    den = {}
    err = {}
    for(i in 1:length(binmid)){
        
        # sub-sample
        samp = data[data>binlo[i] & data<binhi[i]]
        
        # calculate bin values
        sampnum = length(samp)
        sampden = (sampnum/volume)/bw
        samperr = (sqrt(sampnum)/volume)/bw
        
        # add to data vectors
        num = c(num,sampnum)
        den = c(den,sampden)
        err = c(err,samperr)
        
    }

    # calculate upper and lower errors
    errlo = den-err
    errhi = den+err
    if(any(errlo==0)){
        errlo[which(errlo==0)] = null
    }
    if(any(errhi==0)){
        errhi[which(errhi==0)] = null
    }

    # gaussian (to be) fit values
    if(any(num<=numlim)){
        bad = which(num<=numlim)
        fitbinmid = binmid[-bad]
        fitbinlo = binlo[-bad]
        fitbinhi = binhi[-bad]
        fitnum = num[-bad]
        fitden = den[-bad]
        fiterr = err[-bad]
        fiterrlo = errlo[-bad]
        fiterrhi = errhi[-bad]
    }else{
        fitbinmid = binmid
        fitbinlo = binlo
        fitbinhi = binhi
        fitnum = num
        fitden = den
        fiterr = err
        fiterrlo = errlo
        fiterrhi = errhi
    }

    # impose any upper and lower limits
    bad = {}
    if(!is.na(lim1)){
        if(any(fitbinlo<lim1)){
            bad = c(bad,which(fitbinlo<lim1))
        }
    }
    if(!is.na(lim2)){
        if(any(fitbinhi>lim2)){
            bad = c(bad,which(fitbinhi>lim2))
        }
    }
    if(length(bad)>0){
        fitbinmid = fitbinmid[-bad]
        fitbinlo = fitbinlo[-bad]
        fitbinhi = fitbinhi[-bad]
        fitnum = fitnum[-bad]
        fitden = fitden[-bad]
        fiterr = fiterr[-bad]
        fiterrlo = fiterrlo[-bad]
        fiterrhi = fiterrhi[-bad]
    }
    
    # return results
    return(list(bins=bins, binmid=binmid, binlo=binlo, binhi=binhi, num=num, den=den, err=err, errlo=errlo, errhi=errhi, fitbinmid=fitbinmid, fitbinlo=fitbinlo, fitbinhi=fitbinhi, fitnum=fitnum, fitden=fitden, fiterr=fiterr, fiterrlo=fiterrlo, fiterrhi=fiterrhi))
    
}

# chi2 function
.gauss.fit.chi = function(vars, data, binmid, stdev, bw, ftype, fixvars, ...){
    
    # split input
    if("m1"%in%names(vars)){m1 = vars["m1"]}else{m1 = fixvars["m1"]}
    if("s1"%in%names(vars)){s1 = vars["s1"]}else{s1 = fixvars["s1"]}
    if("n1"%in%names(vars)){n1 = vars["n1"]}else{n1 = fixvars["n1"]}
    if("m2"%in%names(vars)){m2 = vars["m2"]}else{m2 = fixvars["m2"]}
    if("s2"%in%names(vars)){s2 = vars["s2"]}else{s2 = fixvars["s2"]}
    if("n2"%in%names(vars)){n2 = vars["n2"]}else{n2 = fixvars["n2"]}
    
    # join input
    m = c(m1,m2)
    s = c(s1,s2)
    n = c(n1,n2)
    
    # gauss model
    model = gauss(x=binmid, mean=m, sigma=s, norm=n, bw=bw, ftype=ftype, ...)
    
    # calculate chi2
    chi2 = sum(((data-model)^2)/(stdev^2))
    return(chi2)
    
}




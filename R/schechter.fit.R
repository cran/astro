schechter.fit = function(data, knee, slope, norm, knee.alt = NA, slope.alt = NA, norm.alt = NA, kneelo = -Inf, slopelo = -Inf, normlo = 0, kneehi = Inf, slopehi = Inf, normhi = Inf, fixk1 = FALSE, fixs1 = FALSE, fixn1 = FALSE, fixk2 = FALSE, fixs2 = FALSE, fixn2 = FALSE, range = range(data), lim1 = NA, lim2 = NA, numlim = 1, method = "nlminb", volume = 1, bw = 0.1, mag = FALSE, log = FALSE, null = 1E-9, error = "jack", subvol = 10, sampnum = subvol, msun = solar("r"), ...){
    
    # setup
    range = sort(range)
    
#    # export
#    .schechter.fit.bin = astro:::.schechter.fit.bin
#    .schechter.fit.dat = astro:::.schechter.fit.dat
#    .schechter.fit.chi = astro:::.schechter.fit.chi
    
    # schechter bin calculations
    bindat = .schechter.fit.bin(data=data, range=range, lim1=lim1, lim2=lim2, numlim=numlim, volume=volume, bw=bw, null=null)
    
    # only continue if any data points remain
    if(length(bindat$fitden)>0){
        
        fitdat = .schechter.fit.dat(bindat=bindat, knee=knee, slope=slope, norm=norm, knee.alt=knee.alt, slope.alt=slope.alt, norm.alt=norm.alt, kneelo=kneelo, slopelo=slopelo, normlo=normlo, kneehi=kneehi, slopehi=slopehi, normhi=normhi, method=method, bw=bw, mag=mag, log=log, fixk1=fixk1, fixs1=fixs1, fixn1=fixn1, fixk2=fixk2, fixs2=fixs2, fixn2=fixn2, msun=msun)
        
        # errors?
        if(error=="jack"){
            
            # jack setup
            subvols = round(seq(1,(length(data)+1),len=subvol+1))
            jacklist = 1:subvol
            
            cat("\n")
            
            for(l in jacklist){
                
                # setup
                cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",l,"/",(length(subvols)-1),"")
                jackdat = data[-(subvols[l]:subvols[l+1])]
                
                # jackknifed binned data
                jackbin = .schechter.fit.bin(data=jackdat, range=range, lim1=lim1, lim2=lim2, numlim=numlim, volume=volume, bw=bw, null=null)
                
                # jackknifed fitted data
                jackfit = .schechter.fit.dat(bindat=jackbin, knee=knee, slope=slope, norm=norm, knee.alt=knee.alt, slope.alt=slope.alt, norm.alt=norm.alt, kneelo=kneelo, slopelo=slopelo, normlo=normlo, kneehi=kneehi, slopehi=slopehi, normhi=normhi, method=method, bw=bw, mag=mag, log=log, fixk1=fixk1, fixs1=fixs1, fixn1=fixn1, fixk2=fixk2, fixs2=fixs2, fixn2=fixn2, msun=msun)
                
                # add par to results
                if(l==jacklist[1]){
                    pars = jackfit$par
                    js = jackfit$j
                    chi2s = jackfit$chi2
                    rchi2s = jackfit$rchi2
                }else{
                    pars = rbind(pars, jackfit$par)
                    js = c(js, jackfit$j)
                    chi2s = c(chi2s, jackfit$chi2)
                    rchi2s = c(rchi2s, jackfit$rchi2)
                }
                
            }
            
            # column definitions
            kcols = grep("k",colnames(pars))
            scols = grep("s",colnames(pars))
            ncols = grep("n",colnames(pars))
            
            # phi* rescaling
            pars[,ncols] = pars[,ncols]/((subvol-1)/subvol)
            js = js/((subvol-1)/subvol)
            
            # calculate errors
            sigma = function(N, x2, x){
                sigma = sqrt(((N-1)/N)*sum((x2-x)^2))
                return(sigma)
            }
            kneelo1 = -1*sigma(subvol, pars[,kcols[1]], fitdat$par["k1"])
            kneehi1 = +1*sigma(subvol, pars[,kcols[1]], fitdat$par["k1"])
            slopelo1 = -1*sigma(subvol, pars[,scols[1]], fitdat$par["s1"])
            slopehi1 = +1*sigma(subvol, pars[,scols[1]], fitdat$par["s1"])
            normlo1 = -1*sigma(subvol, pars[,ncols[1]], fitdat$par["n1"])
            normhi1 = +1*sigma(subvol, pars[,ncols[1]], fitdat$par["n1"])
            jlo = -1*sigma(subvol, js, fitdat$j)
            jhi = +1*sigma(subvol, js, fitdat$j)
            if("k2"%in%names(fitdat$par)){
                kneelo2 = -1*sigma(subvol, pars[,kcols[2]], fitdat$par["k2"])
                kneehi2 = +1*sigma(subvol, pars[,kcols[2]], fitdat$par["k2"])
            }else{
                kneelo2 = kneehi2 = NA
            }
            if("s2"%in%names(fitdat$par)){
                slopelo2 = -1*sigma(subvol, pars[,scols[2]], fitdat$par["s2"])
                slopehi2 = +1*sigma(subvol, pars[,scols[2]], fitdat$par["s2"])
            }else{
                slopelo2 = slopehi2 = NA
            }
            if("n2"%in%names(fitdat$par)){
                normlo2 = -1*sigma(subvol, pars[,ncols[2]], fitdat$par["n2"])
                normhi2 = +1*sigma(subvol, pars[,ncols[2]], fitdat$par["n2"])
            }else{
                normlo2 = normhi2 = NA
            }
            parlo = c(k1lo=kneelo1, s1lo=slopelo1, n1lo=normlo1, k2lo=kneelo2, s2lo=slopelo2, n2lo=normlo2)
            parhi = c(k1hi=kneehi1, s1hi=slopehi1, n1hi=normhi1, k2hi=kneehi2, s2hi=slopehi2, n2hi=normhi2)
            if(any(is.na(parlo))){parlo = parlo[-which(is.na(parlo))]}
            if(any(is.na(parhi))){parhi = parhi[-which(is.na(parhi))]}
            
        }else if(error=="boot"){
            
            NULL
            
        }else{
            
            parlo = parhi = jlo = jhi = NA
            
        }
        
        # build data list
        dat = list(binmid=bindat$binmid, num=bindat$num, den=bindat$den, err=bindat$err, errlo=bindat$errlo, errhi=bindat$errhi, par=fitdat$par, parlo=parlo, parhi=parhi, j=fitdat$j, jlo=jlo, jhi=jhi, chi2=fitdat$chi2, dof=fitdat$dof, rchi2=fitdat$rchi2, denlim=(numlim/volume)/bw, hessian=fitdat$hessian)
        
    }else{
        
        warning("No data available for fitting, check limits")
        dat = list(binmid=NA, num=NA, den=NA, err=NA, errlo=NA, errhi=NA, par=NA, parlo=NA, parhi=NA, j=NA, jlo=NA, jhi=NA, chi2=NA, dof=NA, rchi2=NA, denlim=NA, hessian=NA)
        
    }
    
    # return fitting data
    return(dat)
    
}

# prepare and fit schechter function
.schechter.fit.dat = function(bindat, knee, slope, norm, knee.alt, slope.alt, norm.alt, kneelo, slopelo, normlo, kneehi, slopehi, normhi, method, bw, mag, log, fixk1, fixs1, fixn1, fixk2, fixs2, fixn2, msun){
    
    # grids
    blank = rep(NA,6)
    pri = c(k1=knee[1],s1=slope[1],n1=norm[1],k2=knee[2],s2=slope[2],n2=norm[2])
    sec = c(k1=knee.alt[1],s1=slope.alt[1],n1=norm.alt[1],k2=knee.alt[2],s2=slope.alt[2],n2=norm.alt[2])
    grid = rbind(blank,pri,sec)
    if(sum(is.na(sec))==6){
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
    lden = {}
    for(k in 2:length(grid[,1])){
        
        # create input for fit
        vars = grid[k,]
        for(i in 1:length(vars)){
            if(names(vars)[i] %in% names(fixvars)){
                bad = which(names(fixvars)==names(vars)[i])
                fixvars = fixvars[-bad]
            }
        }
        if(fixk1){if("k1"%in%names(vars)){fixvars=c(fixvars,vars["k1"]); vars=vars[-which(names(vars)=="k1")]}}
        if(fixs1){if("s1"%in%names(vars)){fixvars=c(fixvars,vars["s1"]); vars=vars[-which(names(vars)=="s1")]}}
        if(fixn1){if("n1"%in%names(vars)){fixvars=c(fixvars,vars["n1"]); vars=vars[-which(names(vars)=="n1")]}}
        if(fixk2){if("k2"%in%names(vars)){fixvars=c(fixvars,vars["k2"]); vars=vars[-which(names(vars)=="k2")]}}
        if(fixs2){if("s2"%in%names(vars)){fixvars=c(fixvars,vars["s2"]); vars=vars[-which(names(vars)=="s2")]}}
        if(fixn2){if("n2"%in%names(vars)){fixvars=c(fixvars,vars["n2"]); vars=vars[-which(names(vars)=="n2")]}}
        
        # fitting limits
        lolim = {}
        hilim = {}
        if("k1"%in%names(vars)){lolim = c(lolim, kneelo); hilim=c(hilim, kneehi)}
        if("s1"%in%names(vars)){lolim = c(lolim, slopelo); hilim=c(hilim, slopehi)}
        if("n1"%in%names(vars)){lolim = c(lolim, normlo); hilim=c(hilim, normhi)}
        if("k2"%in%names(vars)){lolim = c(lolim, kneelo); hilim=c(hilim, kneehi)}
        if("s2"%in%names(vars)){lolim = c(lolim, slopelo); hilim=c(hilim, slopehi)}
        if("n2"%in%names(vars)){lolim = c(lolim, normlo); hilim=c(hilim, normhi)}
        
        # fit the Schechter function
        if(method!="nlminb"){
            fit = optim(vars, .schechter.fit.chi, data=bindat$fitden, binmid=bindat$fitbinmid, stdev=bindat$fiterr, bw=bw, mag=mag, log=log, fixvars=fixvars, hessian=TRUE, method=method, lower=lolim, upper=hilim)
        }else{
            fit = nlminb(vars, .schechter.fit.chi, data=bindat$fitden, binmid=bindat$fitbinmid, stdev=bindat$fiterr, bw=bw, mag=mag, log=log, fixvars=fixvars, lower=lolim, upper=hilim)
            fit$value=fit$objective
        }
        
        # add on any fixed parameters to $par
        if(fixk1){fit$par = c(fit$par,fixvars["k1"])}
        if(fixs1){fit$par = c(fit$par,fixvars["s1"])}
        if(fixn1){fit$par = c(fit$par,fixvars["n1"])}
        if(fixk2){fit$par = c(fit$par,fixvars["k2"])}
        if(fixs2){fit$par = c(fit$par,fixvars["s2"])}
        if(fixn2){fit$par = c(fit$par,fixvars["n2"])}
        
        # calculate luminosity density
        j = lumdens(knee=as.numeric(fit$par[grep("k",names(fit$par))]), slope=as.numeric(fit$par[grep("s",names(fit$par))]), norm=as.numeric(fit$par[grep("n",names(fit$par))]), msun=msun, mag=mag, log=log)
        
        # collect results
        chis = c(chis, fit$value)
        pars = c(pars, list(fit$par))
        hess = c(hess, list(fit$hessian))
        lden = c(lden, j)
        
    }
    
    # choose best fit
    chi = chis[[which.min(chis)]]
    par = pars[[which.min(chis)]]
    hessian = hess[[which.min(chis)]]
    j = lden[[which.min(chis)]]
    chi2 = chi
    dof = (length(bindat$fitden)-length(par))
    rchi2 = chi2/dof
    
    # return results
    return(list(par=par, j=j, chi2=chi2, dof=dof, rchi2=rchi2, hessian=hessian))
    
}

# calculate number (total), number density (phi) and poissonian errors
.schechter.fit.bin = function(data, range, lim1, lim2, numlim, volume, bw, null){
    
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

    # schechter (to be) fit values
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
.schechter.fit.chi = function(vars, data, binmid, stdev, bw, mag, log, fixvars, ...){
    
    # split input
    if("k1"%in%names(vars)){k1 = vars["k1"]}else{k1 = fixvars["k1"]}
    if("s1"%in%names(vars)){s1 = vars["s1"]}else{s1 = fixvars["s1"]}
    if("n1"%in%names(vars)){n1 = vars["n1"]}else{n1 = fixvars["n1"]}
    if("k2"%in%names(vars)){k2 = vars["k2"]}else{k2 = fixvars["k2"]}
    if("s2"%in%names(vars)){s2 = vars["s2"]}else{s2 = fixvars["s2"]}
    if("n2"%in%names(vars)){n2 = vars["n2"]}else{n2 = fixvars["n2"]}
    
    # join input
    k = c(k1,k2)
    s = c(s1,s2)
    n = c(n1,n2)
    
    # schechter model
    model = schechter(x=binmid, knee=k, slope=s, norm=n, bw=bw, mag=mag, log=log, ...)
    
    # calculate chi2
    chi2 = sum(((data-model)^2)/(stdev^2))
    return(chi2)
    
}




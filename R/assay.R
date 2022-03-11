# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'


# also:
# http://bioconductor.org/developers/how-to/buildingPackagesForBioc/#the-r-dir
# need to pass the BiocCheck function.

# currently there are errors in the example, which need to be remedied somehow.
# needs msaccess to do the examples, so maybe I need to change the interface
# to use mzR, despite the problems encountered previously.



# Error where it crashes if no peaks detected on first metabolite. --- done ? (actually, might have already been done... double done?)
# tested okay
#
# In isotope columns, replace blank with 0. --- done?
#
# 'Interactive' with Yes (capital Y) is not recognized. --- done?
#
# Handle question mark in metabolite name. --- done? --- handle slashes!
#
# Change 'samples' to datapoints - perhaps convert to RT width of peak??  --- done?
# done, now requires seconds instead of samples.  tested okay
#
# RTmax less than RTmin freeze & crash --- done?
#
# RT range smaller than hat error. --- done?
#
# Don't include standard if different intensity - can there be an 'omit' function during peak picking? e.g. omit file from peak picking if 'standard' is in the filename.
#  --- done, will omit names containing "std" or "standard".  Currently no configurability.
#


## CODE LOCATIONS MARKED (with "THIS IS WHERE")

# Within a channel - separate peak picking parameters (fructose and glucose)  --- location marked in code???
#
# save changes to config incremently --- location marked in code
#
# Minima error (GAA plot) - only detecting start of peak? --- location marked in code
#

## NOT DONE

# Peak identity (i.e. which line is which sample/isotope on the XIC)
# - no! Instead, give option for integrated output (stacked bar) during peak picking.
#       Would help if baseline is plotted too - pick up duff samples.
# give option in interactive mode for barplot for a chosen peak.
#
# better.table not giving bar plots (when no isotopes chosen)
#
# Instructions: Clear instructions and video (Andy).
#

# QUERY:
#
# Also, give RTs as real values (rather than hat start/end) - allow more accurate peak RT measurement.
### you mean in report or config?  config can be for many peaks in one line and therefore needs a range
#
#


# a function to help with reshaping... I keep figuring this out and not writing it down!

# the steps in reconfiguring a table...
# 1) potentially set a column from the rownames
# 2) melt
# 3) split var column
# 4) delete var column
# 5) optionally delete one of the columns from the split
# 6) optionally combine one of the column from split with the original row name column, and delete those
# 7) dcast, potentially with function to combine data (e.g. mean, sum, etc)
## parameters required...
### has.row.name.column .. can guess a matrix would not, and a data.frame that converts to a numerix matrix does not
### column.name.split.pattern ... will require one or two captures ... number of captures should determine number of columns
###                               this takes care of optional deletion of one of the columns
### combine.row.names.with.match ... 0 means none, 1 means the first capture, 2 the second... the combined columns are removed
### FUN = function for combining
# this approach simplifies the parameters to a single regex, a boolean and an integer.  Pretty good if it works.

min.RT.fin <- 0
max.RT.fin <- 1
med.diff.fin <- .01

reconfigure.table <-
  function(X,
           column.name.split.pattern = "^([,]+),([,]+)",
           combine.row.names.with.match = 1,
           has.row.name.column = false,
           FUN = sum,
           ...) {
    m <- melt(X)
    t <-
      str_match(as.character(m$variable), "(.*)(0|\\.\\d+[CNH]\\.\\d+)")  # matrix of matches
    # join them up.
    
    # TODO... implement the rest of it!
  }

## actually, the above doesn't even cover what I really need, which is to remove the replicate number from the middle of the
## column names and put the thing back together again.
## NB: the following will only work for labelled expts
combine.rep.cols <- function(x, FUN = mean, ...) {
  x[is.na(x)] = 0
  x$X <- rownames(x)
  m <- melt(x)
  m$variable <-
    sub("n\\d+(.*)-\\d+(0|-\\d+[CNH]-\\d+)",
        "n99\\1-9\\2",
        as.character(m$variable))
  dcast(m, X ~ variable, fun.aggregate = mean, na.rm = TRUE) -> d
  rownames(d) <- d$X
  d$X <- NULL
  return(d)
}


assay.plotter = function(p,
                         condition.pattern = "(^.*)-[[:digit:]]+$",
                         repeat.inj.pattern = "n[[:digit:]]+[\\-_]*.*") {
  p = labels.with.component(p)
  
  write.csv(p, "p2-results.csv")
  
  pp <- read.csv("p2-results.csv")
  pp[pp == 0] <- NA
  
  melt(pp) -> mm
  gsub("-0$|-[[:digit:]]+[NCH]-[[:digit:]]+$", "", mm$X) -> mm$X
  
  dcast(mm, X ~ variable, sum, na.rm = TRUE) -> dd
  rownames(dd) <- dd$X
  dd$X <- NULL
  meds <- apply(dd, 2, median, na.rm = TRUE)
  meds <- meds / median(meds, na.rm = TRUE)
  p <- as.data.frame(t(apply(p, 1, function(r) {
    return(r / meds)
  })))
  
  #p = p[complete.cases(p),]
  
  #x = p
  #x[x == 0] = 1
  #heatmap(scale(log(x)), margins=c(9,7))
  
  # in this case, exlcuding n1 5min rep 1, as it's clearly an outlier...
  #q = p[,grep('minutepulse',colnames(p))]
  #q = q[,-grep('_n1.*5minutepulse1',colnames(q))]
  #q = q[,-grep('_n2',colnames(q))]
  
  rns = rownames(p)
  rns = gsub('-0$', '', rns[grep("-0$", rns)])
  
  for (rn in rns) {
    rnp = unlist(paste('^', rn, '-(.*)', sep = ''))
    rn <- gsub("[/\\^*+'\"#@]+", "_", rn)
    fn = unlist(paste('relbars_', rn, '.png', sep = ''))
    fn2 = unlist(paste('absbars_', rn, '.png', sep = ''))
    png(fn)
    #component.peak.pattern = "^pyruvate.* peak1+-(.*)"
    component.peak.pattern = rnp
    assay.plot(p,
               repeat.inj.pattern,
               component.peak.pattern,
               condition.pattern)
    title(fn)
    dev.off()
    png(fn2)
    component.peak.pattern = rnp
    assay.plot(p,
               repeat.inj.pattern,
               component.peak.pattern,
               condition.pattern,
               FALSE)
    title(fn2)
    dev.off()
  }
  return(p)
}
#

assay.plot = function(q,
                      repeat.inj.pattern = ".*n[[:digit:]]-*(.*)",
                      component.peak.pattern = "^Galactose.* peak4+-(.*)",
                      condition.pattern = "(^.*)-[[:digit:]]+$",
                      relative = TRUE) {
  q = reshape.assay.result2(q, xpatt = repeat.inj.pattern, ypatt = component.peak.pattern, mean)
  if (relative) {
    sums = apply(q, 2, sum)
    s = matrix(sums, nrow(q), ncol(q), byrow = TRUE)
    q = q / s
  }
  #r = reshape.assay.result2(q, xpatt= "([[:digit:]]+min)utepulse.*", ypatt= "^(.*)$", median)
  #s = reshape.assay.result2(q, xpatt= condition.pattern, ypatt= "^(.*)$", sd)
  
  par(las = 2)
  par(mar = c(10, 10, 3, 4))
  
  q <- as.matrix(q)
  q[is.infinite(q)] <- NA
  
  barplot(
    as.matrix(q),
    xlim = c(0, ncol(q) + 5),
    col = brewer.pal(nrow(q), "Paired"),
    legend.text = TRUE,
    args.legend = list(
      x = ncol(q) + 5,
      y = max(colSums(q)),
      bty = "n"
    )
  )
}


reshape.assay.result2 = function(r, xpatt, ypatt, fun) {
  x = grep(xpatt, colnames(r))
  y = grep(ypatt, rownames(r))
  if (length(x) == 0) {
    cn = paste(colnames(r), collapse = " ")
    cat(
      paste(
        "reshape.assay.result2: xpatt {",
        xpatt,
        "} did not match anything in colnames(r):",
        cn
        ,
        "\n"
      )
    )
    return()
  }
  if (length(y) == 0) {
    rn = paste(rownames(r), collapse = " ")
    cat(
      paste(
        "reshape.assay.result2: ypatt {",
        ypatt,
        "} did not match anything in rownames(r):",
        rn,
        "\n"
      )
    )
    return()
  }
  g = r[y, x]
  g$id = gsub(ypatt, '\\1', rownames(g))
  m = melt(g, id.vars = 'id')
  # m$variable = gsub(xpatt,'\\1',m$variable)
  m = m[!is.na(m$value), ]
  d = dcast(m, formula = id ~ variable, fun.aggregate = fun)
  rownames(d) = d$id
  d$id = NULL
  return(d)
}



labels.with.component = function(q) {
  q$id = rownames(q)
  m = melt(q, id.vars = 'id')
  m$id2 = gsub('-(13C|15N|2H)-[[:digit:]]+$|0$', '', m$variable)
  m$variable = gsub('^.*?(-(13C|15N|2H)-[[:digit:]]+$|0$)', '\\1', m$variable)
  m$variable = gsub('^0$', '-0', m$variable)
  r = dcast(m, id + variable ~ id2)
  rownames(r) = paste(r$id, r$variable, sep = '')
  r$variable = NULL
  r$id = NULL
  return(r)
}

## note: calls to the following should be optional!
shorten.names = function(n, sep = "[.]") {
  # remove redundancy from a bunch of names (which must have the same format!)
  splt = as.data.frame(strsplit(n, sep))
  u = apply(splt, 1, unique)
  i = which(lapply(u, length) > 1)
  newn = as.vector(unlist(lapply(splt[i, ], paste, collapse = "_")))
  return(paste("n", newn, sep = ""))
}



standardize.RTs = function(list.of.TICs,
                           pattern = "std|standard",
                           match.mode = "omit",
                           shorten.col.names = TRUE) {
  if (length(list.of.TICs) == 0) {
    print("ERROR: list.of.TICs is empty")
    return()
  }
  # first, let's get some basic info on the data...
  median.diffs = c()
  min.RTs = c()
  max.RTs = c()
  for (tic in list.of.TICs) {
    median.diffs = c(median.diffs, median(diff(tic$rt)))
    min.RTs = c(min.RTs, min(tic$rt))
    max.RTs = c(max.RTs, max(tic$rt))
  }
  print(median.diffs)
  print(names(list.of.TICs))
  med.diff = median(median.diffs)
  min.RT = max(min.RTs)
  max.RT = min(max.RTs)
  
  
  # now generate an x series (retention time with regular intervals)
  # There is an issue 4/7/21 with min.RT going to inf. will set a safe min.RT here and reset it to this if it is set to inf in error
  
  if (is.finite(min.RT)) {
    min.RT.fin <<- min.RT
  }
  else{
    min.RT <- min.RT.fin
    print("min.RT was not finite")
  }
  
  if (is.finite(max.RT)) {
    max.RT.fin <<- max.RT
  }
  else{
    max.RT <- max.RT.fin
    print("max.RT was not finite")
  }
  
  if (is.na(med.diff)) {
    med.diff <- med.diff.fin
    print("med.diff was NA")
  }
  else{
    med.diff.fin <<- med.diff
    print("med.diff was NA")
  }
  
  RTs = seq(min.RT, max.RT, med.diff)
  # now we can interpolate...
  ## note: calls to the following should be optional!
  n <- names(list.of.TICs)
  if (shorten.col.names) {
    n = shorten.names(n)
  }
  chromatogram = data.frame(rt = RTs)
  for (i in 1:length(list.of.TICs)) {
    # this could be apply instead of for, if we first filter out which
    # chromatograms match the omit pattern
    name = n[i]
    ######################### THIS IS WHERE...
    ######################## we can omit the standards from peak picking.
    if ((length(grep(pattern, name)) == 0 && match.mode == "omit")
        ||
        length(grep(pattern, name)) != 0) {
      tic = list.of.TICs[[i]]
      interp = approx(tic$rt, tic$sumIntensity, RTs)
      chromatogram[name] = interp$y
    }
  }
  # now we'll calculate the maximum chromatogram, which will contain
  # all the peaks we can to detect...
  
  if (ncol(chromatogram) < 2) {
    print("ERROR: chromatogram table contains no data")
    return()
  }
  
  chromatogram$max = apply(chromatogram[2:ncol(chromatogram)], 1, max)
  ### (though maybe this shouldn't be part of this function as the function name suggests nothing about it)
  return(chromatogram)
}

baseline = function(x, d) {
  # d is the peak-detection return value
  m = c()
  if (length(names(d)) > 1) {
    m = d[, 2]
  }
  if (length(d) > 2) {
    m = apply(d[, 2:length(d)], 1, max)
  }
  l = m > 0
  dl = diff(l)
  starts = which(dl == 1)
  ends = which(dl == -1)
  # probably we should be checking here
  # for peaks that start or end outside our window
  # actually, we should in fact fail if this is the
  # case so the user knows to adjust the window...
  if (length(starts) == 0 || length(ends) == 0) {
    cat("no peaks found!")
    #fail()
  } else if (length(starts) != length(ends) ||
             starts[1] > ends[1] ||
             starts[length(starts)] > ends[length(ends)]) {
    cat("peaks overlap window!")
    #fail()
  }
  
  lastsample = length(x) - 1
  n = names(x)
  for (j in 2:lastsample) {
    np = x[, j]
    np[l] = NA
    if (length(ends) > 0) {
      for (i in 1:length(ends)) {
        #cat(paste("i =", i, ", length(ends)=", length(ends), "\n"))
        #cat(paste("starts[i] =", starts[i], ", ends[i]=", ends[i], "\n"))
        start = starts[i]
        end = ends[i] + 1
        L = end - start
        prestart = start - L
        postend = end + L
        if (prestart < 1) {
          prestart = 1
        }
        if (postend > length(np)) {
          postend = length(np)
        }
        
        # fix NA problem...
        if (is.na(sum(np[end:postend]))) {
          # there are NAs
          nais = which(is.na(np[end:postend]))
          postend = nais[1] + end
          #cat("postend fixed: ", postend, "\n")
        }
        
        
        #cat(prestart, start, end, postend,length(np),"\n")
        
        med1i = round(mean(c(prestart, start), na.rm = TRUE))
        med2i = round(mean(c(end, postend), na.rm = TRUE))
        # for the max ...
        ########################### THIS IS WHERE...
        ######################### we could use min in stead of median to get neighbour-peak-resistant baseline
        min1 = min(np[prestart:start], na.rm = TRUE)
        min2 = min(np[end:postend], na.rm = TRUE)
        #med1 = median(np[prestart:start], na.rm = TRUE)
        #med2 = median(np[end:postend], na.rm = TRUE)
        if (length(min1) < 1) {
          min1 = 0
        }
        if (length(min2) < 1) {
          min2 = 0
        }
        #cat(prestart, start, end, postend, med1i, med2i, length(np),"\n")
        #cat(np[med1i],np[med2i],"\n")
        a = approx(c(med1i, med2i), c(min1, min2), start:end)
        # then e.g.
        np[start:end] = a$y
      }
    }
    newname = unlist(paste(n[j], ".baseline", sep = ""))
    x[newname] = np
  }
  return(x)
}

limit.tic = function(tic, rt1, rt2) {
  i = tic$rt > rt1 & tic$rt < rt2
  return (tic[i, ])
}

read.tic = function(file) {
  #tic = read.delim(file,skip=1)
  tic = read.delim(file)
  tic$rt = tic$rt / 60
  return(tic)
}

mexican.hat = function(peaksamples, normlim = 5) {
  # 1,1 give the original data.
  # 5 is a good limit (5 standard deviations), as the hat levels off
  # at zero nicely.  The sign-change happens at about 1 sd.
  
  #totalsamples = peaksamples * normlim
  #normstep = normlim / totalsamples
  #which can be simplified as
  peaksamples = peaksamples / 2
  normstep = 1 / peaksamples
  
  # grab a normal distribution...
  mynorm = dnorm(seq(
    from = -normlim,
    to = normlim,
    by = normstep
  ))
  # calculate the 1st derivative
  d1mynorm = mynorm[2:length(mynorm)] - mynorm[1:length(mynorm) - 1]
  # calculate the 2nd derivative
  d2mynorm = d1mynorm[2:length(d1mynorm)] - d1mynorm[1:length(d1mynorm) -
                                                       1]
  # return the normalized result
  midstart = peaksamples * normlim * 0.8
  midend = peaksamples * normlim * 1.2
  mysum = sum(d2mynorm[midstart:midend])
  return(d2mynorm / mysum)
}

# THEN DO SOMETHING LIKE:
# plot(tic$sumIntensity,type="l",col="red")
# lines(filter(tic$sumIntensity, mexican.hat(20)))

logical.peaks = function(x, wavethresh = 1e5) {
  #cat("==================================================\nx:\n------------\n")
  #print(x)
  
  ######################## THIS IS WHERE
  #### we'll start to fix some of the one-ended peak errors...
  
  maxima = which.peaks(x, decreasing = FALSE)
  minima = which.peaks(x, decreasing = TRUE)
  minima = c(50, minima, length(x) - 50)   #### these are not the correct default end values!
  
  
  #cat("==================================================\nmaxima:\n------------\n")
  #print(maxima)
  #cat("==================================================\nminima:\n------------\n")
  #print(minima)
  
  threshmaxima = maxima[x[maxima] > wavethresh]
  
  # this is failing now because there are no peaks that make the threshold,
  # so threshmaxima is empty...
  # so starts and ends are empty... so whoever calls this funciton should be aware
  
  #cat("==================================================\nthreshmaxima:\n------------\n")
  #print(threshmaxima)
  
  if (length(threshmaxima) == 0) {
    return(data.frame(starts = c(), ends = c()))
  }
  
  starts = c()
  ends = c()
  
  for (i in threshmaxima) {
    starti = max(which(minima < i))
    endi = min(which(minima > i))
    starts = c(starts, minima[starti])
    ends = c(ends, minima[endi])
  }
  
  if (is.na(ends[length(ends)])) {
    ends[length(ends)] = length(x) - 1
  }
  if (starts < 0) {
    starts = 0
  }
  if (ends > length(x)) {
    ends = length(x)
  }
  #cat("==================================================\nstarts:\n------------\n")
  #print(starts)
  
  #cat("==================================================\nends:\n------------\n")
  #print(ends)
  
  return(data.frame(start = starts, end = ends))
  
}

which.peaks <- function(x,
                        partial = TRUE,
                        decreasing = FALSE) {
  if (decreasing) {
    if (partial) {
      which(diff(c(FALSE, diff(x) > 0, TRUE)) > 0)
    } else {
      which(diff(diff(x) > 0) > 0) + 1
    }
  } else {
    if (partial) {
      which(diff(c(TRUE, diff(x) >= 0, FALSE)) < 0)
    } else {
      which(diff(diff(x) >= 0) < 0) + 1
    }
  }
}

detect.peaks = function(tic,
                        samples = 20,
                        wavethresh = 2e4,
                        normlim = 5) {
  ########################## THIS IS WHERE
  ##### we can pad tic so hat always fits...
  
  firstvalue <- tic[1]
  lastvalue <- tic[length(tic)]
  buffer <- 2 * samples
  tic <- c(rep(firstvalue, buffer), tic, rep(lastvalue, buffer))
  f = filter(tic, mexican.hat(samples, normlim))
  f[is.na(f)] <- 0
  # I just don't understand how this can still crash with filter longer than sample.
  
  A <- buffer + 1
  B <- length(tic) - buffer
  #lpc = lp
  tic <- tic[A:B]
  f <- f[A:B]
  
  lp = logical.peaks(f, wavethresh)
  r = data.frame("tic" = tic)
  ### now lp's indices are all out by 'samples'
  #lp <- lp - buffer # should correct it
  
  # lp is now a data.frame with start and end
  
  # can use the minima that define peak edges to say
  # whether two peaks are actually a split peak.
  # this means the intensity considered one that
  # part of the peak as background will be the average
  # of the non-shared minima... if you know what I mean...
  #
  
  
  #peakn = 1
  
  lastpeak = ""
  lastminimum = 0
  
  
  # lp can be empty...
  if (length(lp) > 0) {
    for (peakn in 1:length(lp[, 1])) {
      peakstart = lp[peakn, 1]
      peakend = lp[peakn, 2]
      peaktic = tic
      L = length(peaktic)
      peaktic[1:peakstart] = 0
      peaktic[peakend:L] = 0
      name = paste("peak", peakn, sep = "")
      r[name] = peaktic
      if (lp[peakn, 1] == lastminimum) {
        cat(lastpeak, " is joined to ", name, "\n")
      }
      lastminimum = lp[peakn, 2]
      lastpeak = name
      peakn = peakn + 1
    }
  }
  # r might only have tic in it!
  return (r)
}

get.peaks = function(path = ".",
                     pattern = ".tsv",
                     rt.min = 0,
                     rt.max = Inf,
                     seconds = 10,
                     threshold = 1e7,
                     chrom.pattern = "std|standard",
                     chrom.match.mode = "omit",
                     shorten.col.names = TRUE) {
  initialpath = getwd()
  setwd(path)
  files = list.files(
    path = '.',
    pattern = pattern,
    full.names = FALSE,
    ignore.case = TRUE,
    recursive = TRUE
  )
  tics = list()
  max = 0
  filemax = ""
  for (file in files) {
    tic = read.tic(file = file)
    tic = limit.tic(tic, rt.min, rt.max)
    thismax = max(tic$sumIntensity)
    if (thismax > max) {
      max = thismax
      filemax = file
    }
    # smoothing?
    tic$sumIntensity = stats::smooth(tic$sumIntensity)
    tics[[file]] = tic
  }
  setwd(initialpath)
  
  # shortening names messes things up 8/10/2020
  x = standardize.RTs(
    tics,
    pattern = chrom.pattern,
    match.mode = chrom.match.mode,
    shorten.col.names = shorten.col.names
  )
  
  
  ################## HERE IS WHERE
  ###### we can convert seconds into number of samples for mex hat
  # x$rt has fixed period.  Can use this to convert.
  
  
  td <-
    (x$rt[2] - x$rt[1]) * 60 ## time difference between two samples.
  
  samples <- 1 +     # the point at the start
    (seconds / td)   # the points after to make up the correct number of intervals that add up to seconds
  
  samples <- round(samples)
  
  if (samples < 6) {
    seconds <- 5 * td
    print (
      paste(
        "WARNING:",
        seconds,
        "seconds for",
        samples,
        "data points, minimum is 6s.  Resetting to 6 data points for",
        seconds,
        "s"
      )
    )
    samples <- 1 +     # the point at the start
      (seconds / td)   # the points after to make up the correct number of intervals that add up to seconds
  }
  
  d = detect.peaks(x$max, samples, threshold)
  
  
  # d might have only $tics and no peaks defined...
  if (length(d) == 1) {
    plot(x$rt,
         x$max,
         type = 'l',
         xlab = "RT (min)",
         ylab = "intensity")
    return (x)
  }
  
  
  b = baseline(x, d)
  
  xmax = max(x$max)
  
  plot(b$rt,
       b$max,
       type = 'n',
       xlab = "RT (min)",
       ylab = "intensity")
  
  #cols = heat.colors(length(d))
  cols = topo.colors(length(d))
  
  
  for (i in 2:length(d)) {
    ii = which(d[, i] > 0)
    start = x$rt[min(ii)]
    end = x$rt[max(ii)]
    thisx = c(start, start, end, end)
    thisy = c(0, xmax, xmax, 0)
    #polygon(x$rt,d[,i], col=cols[i])
    polygon(thisx, thisy, col = cols[i], border = NA)
  }
  
  cols = rainbow(length(x))
  
  lastsample = length(x) - 1
  for (i in 2:lastsample) {
    col = cols[i]
    lines(x$rt, x[, i], col = col)
    k = i + lastsample
    lines(x$rt, b[, k], col = col)
    
  }
  lines(x$rt, x$max, col = "black")
  
  
  
  r = data.frame(rt = x$rt)
  sumsummary = data.frame(row.names = names(d)[2:length(d)])
  for (j in 2:lastsample) {
    k = j + lastsample
    x_b = x[, j] - b[, k]
    x_b[x_b < 0] = 0
    peaksi = c()
    for (i in 2:length(d)) {
      x_b_copy = x_b
      newname = unlist(paste(names(x)[j], names(d)[i]))
      x_b_copy[d[, i] == 0] = 0
      r[newname] = x_b_copy
      si = sum(x_b_copy)
      peaksi = c(peaksi, si)
    }
    sumsummary[names(x)[j]] = peaksi
  }
  
  ## ADD SOME KIND OF BAR PLOT HERE?
  ## WOULD NEED A PLOT FOR EACH PEAK...
  
  result = c()
  result$data = x
  result$baseline = b
  result$peaks = d
  result$summary = sumsummary
  return(result)
}


ms.access.tics = function(dir = "D:/Processing/Jimi/Marcus1/mzXML",
                          MZS = c(145.0142, 115.0037, 147.0299, 117.0193),
                          PPMs = c(2.5, 2.5, 2.5, 2.5),
                          C13 = c(0, 0, 0, 0),
                          N15 = c(0, 0, 0, 0),
                          H2 = c(0, 0, 0, 0),
                          EVENTS) {
  olddir = getwd()
  setwd(dir)
  
  # path to msaccess (if it's not already in your path)
  msaccess <- c("msaccess")
  if (file.exists("/home/jicawi/bin/msaccess")) {
    msaccess <- c("/home/jicawi/bin/msaccess")
  }
  # it's quicker if these are centroided
  FILES <-
    list.files(recursive = TRUE,
               full.names = TRUE,
               pattern = "\\.mzML")
  cat("FILES:")
  cat (FILES)
  show(FILES)
  
  stamp = strftime(Sys.time(), format = "%Y%m%d-%H%M%S")
  outputbase = unlist(paste(olddir, stamp, sep = "/"))
  dir.create(outputbase)
  
  # TODO: Here we need to do something about isotope labels
  C13delta = 13.0033548378 - 12
  N15delta = 15.0001088982 - 14.003074004
  H2delta = 2.0141017778 - 1.00782503207
  
  for (j in 1:length(MZS)) {
    PPM = PPMs[j]
    mlo = 1 - PPM / 1000000
    mhi = 1 + PPM / 1000000
    mzlo = MZS[j] * mlo
    mzhi = MZS[j] * mhi
    outputdir = unlist(paste(outputbase, MZS[j], sep = "/"))
    dir.create(outputdir)
    
    C13max = C13[j]
    N15max = N15[j]
    H2max  = H2[j]
    EVENT = EVENTS[j]
    
    # C13 LOOP
    for (i_C13 in 0:C13max) {
      for (i_N15 in 0:N15max) {
        for (i_H2 in 0:H2max) {
          cat("\n------------------------------------------------------\n")
          cat("mz =", mzlo, "..", mzhi, "\n")
          cat("13C:", i_C13, "15N:", i_N15, "2H:", i_H2, "\n")
          
          delta = i_C13 * C13delta + i_N15 * N15delta + i_H2 *
            H2delta
          
          cat("delta =", delta, "\n")
          
          imzlo = mzlo + delta
          imzhi = mzhi + delta
          
          ########################## HERE IS WHERE...
          ##################### we can track which XICs are already extracted and skip over them
          
          cat("isotope mz =", imzlo, "..", imzhi, "\n")
          
          
          # warning: shorten.names will fail if you try to use . or _ as sep...
          isoname = ""
          if (i_C13 > 0) {
            isoname = paste(isoname, "13C", i_C13, sep = '-')
          }
          if (i_N15 > 0) {
            isoname = paste(isoname, "15N", i_N15, sep = '-')
          }
          if (i_H2 > 0) {
            isoname = paste(isoname, "2H", i_H2, sep = '-')
          }
          if (isoname == "") {
            isoname = "0"
          }
          cat("name", isoname, "\n")
          
          for (i in 1:length(FILES)) {
            cat(FILES[i], "\n")
            l1 = list.files(outputdir) # list files
            system (
              paste(
                msaccess,
                "-x \"",
                # use the actual m/z measured,
                #paste("tic mz=",imzlo,",",imzhi,isoname,sep=""),
                # or that dervied from?
                paste("tic mz=", imzlo, ",", imzhi, sep =
                        ""),
                "delimiter=tab",
                "\"",
                "--filter \"",
                paste("scanEvent ", EVENT),
                "\"",
                "-o",
                outputdir,
                "-v",
                FILES[i]
              )
            )
            l2 = list.files(outputdir) # list files again
            
            # now grab the name of the new file and set about renaming it!
            newlygeneratedfile = setdiff(l2, l1) # filename
            isofilename = sub("\\.tsv$",
                              paste(".", isoname, ".tsv", sep = ""),
                              newlygeneratedfile)
            src = paste(outputdir, newlygeneratedfile, sep =
                          "/")
            dest = paste(outputdir, isofilename, sep = "/")
            file.rename(src, dest)
            cat(
              "file",
              newlygeneratedfile,
              "renamed to",
              isofilename,
              "in",
              outputdir,
              "\n"
            )
          }
        }
      }
    }
  }
  setwd(olddir)
  return(outputbase)
}

mzR.tics = function(fileName, config, inputDir, outputDir) {
  MZS = config$mz
  PPMs = config$ppm
  C13 = config$C13
  N15 = config$N15
  H2 = config$H2
  
  print(basename(fileName))
  
  # TODO: Here we need to do something about isotope labels
  C13delta = 13.0033548378 - 12
  N15delta = 15.0001088982 - 14.003074004
  H2delta = 2.0141017778 - 1.00782503207
  
  aa = openMSfile(fileName)
  hdr = getTicHeader(aa)
  
  for (j in 1:length(MZS)) {
    PPM = PPMs[j]
    mlo = 1 - PPM / 1000000
    mhi = 1 + PPM / 1000000
    mzlo = MZS[j] * mlo
    mzhi = MZS[j] * mhi
    outputdir = unlist(paste(outputDir, MZS[j], sep = "/"))
    
    C13max = C13[j]
    N15max = N15[j]
    H2max  = H2[j]
    
    # C13 LOOP
    for (i_C13 in 0:C13max) {
      for (i_N15 in 0:N15max) {
        for (i_H2 in 0:H2max) {
          delta = i_C13 * C13delta + i_N15 * N15delta + i_H2 * H2delta
          imzlo = mzlo + delta
          imzhi = mzhi + delta
          
          ########################## HERE IS WHERE...
          ##################### we can track which XICs are already extracted and skip over them
          
          # warning: shorten.names will fail if you try to use . or _ as sep...
          isoname = ""
          if (i_C13 > 0) {
            isoname = paste(isoname, "13C", i_C13, sep = '-')
          }
          if (i_N15 > 0) {
            isoname = paste(isoname, "15N", i_N15, sep = '-')
          }
          if (i_H2 > 0) {
            isoname = paste(isoname, "2H", i_H2, sep = '-')
          }
          if (isoname == "") {
            isoname = "0"
          }
          
          p = getTickPeak(aa, hdr, imzlo, imzhi)
          outputFile = writePeak(basename(fileName),
                                 isoname,
                                 imzlo,
                                 imzhi,
                                 p,
                                 outputDir = outputdir)
        }
      }
    }
  }
}

writePeak <-
  function(mzMLfile,
           isoname,
           mzlo,
           mzhi,
           res,
           outputDir = '.') {
    outFile <- paste(
      outputDir,
      '/',
      mzMLfile,
      '.tic.',
      format(mzlo, nsmall = 2, digits = 2),
      '-',
      format(mzhi, nsmall = 2, digits = 2),
      '.',
      isoname,
      '.tsv',
      sep = ''
    )
    write.table(
      res,
      row.names = FALSE,
      col.names = TRUE,
      file = outFile,
      sep = "\t"
    )
    return(outFile)
  }

# Can this be vectorized?
getTicHeader <- function(aa) {
  v <- 1:runInfo(aa)$scanCount
  dim(v) <- length(v)
  
  h <- t(apply(v, 1, function(i) {
    h <- header(aa, i)
    return(c(h$lowMZ, h$highMZ, h$retentionTime))
  }))
  
  return(h)
}

getTickPeak = function(aa, h, mzlo, mzhi) {
  # Vectorization would happen here
  s <- which(h[, 1] <= mzlo & mzhi <= h[, 2])
  if (length(s) == 0) {
    print(paste(
      "Error: ( mzlo mzhi ) = (",
      mzlo,
      mzhi,
      ") did no match any scan ranges"
    ))
  }
  dim(s) <- length(s)
  
  x <- t(apply(s, 1, function(i) {
    p <- as.data.frame(peaks(aa, i))
    t <- sum(p[mzlo <= p[, 1] & p[, 1] <= mzhi, 2])
    return(c(h[i, 3], t))
  }))
  colnames(x) <- c("rt", "sumIntensity")
  return (x)
}

run.config.tics = function(sampleSheet,
                           inputDir,
                           outputDir,
                           BPPARAM = SerialParam()) {
  #
  config = read.delim(sampleSheet)
  config$C13[is.na(config$C13) | config$C13 == ""] <- 0
  config$N15[is.na(config$N15) | config$N15 == ""] <- 0
  config$H2[is.na(config$H2) | config$H2 == ""] <- 0
  
  # Get a list of files in inputDir
  flist = list.files(
    recursive = TRUE,
    full.names = TRUE,
    path = inputDir,
    pattern = "\\.mzML"
  )
  print(flist)
  if (!dir.exists(outputDir))
    dir.create(outputDir, recursive = TRUE)
  lapply(config$mz, function(x) {
    odir = sprintf("%s/%s", outputDir, x)
    if (!dir.exists(odir))
      dir.create(odir, recursive = TRUE)
  })
  bplapply(flist, function(x)
    mzR.tics(x, config, inputDir, outputDir), BPPARAM = BPPARAM)
  outputDir
}


run.config.peaks = function(path.to.config.tsv = 'config.tsv',
                            path.to.tics,
                            Interactive = TRUE,
                            chrom.pattern = "std|standard",
                            chrom.match.mode = "omit",
                            shorten.col.names = TRUE) {
  config = read.delim (path.to.config.tsv)
  
  summ = data.frame()
  for (i in 1:length(config[, 1])) {
    s <- run.config.peak(
      path.to.config.tsv = path.to.config.tsv,
      path.to.tics = path.to.tics,
      index = i,
      Interactive = Interactive,
      chrom.pattern = chrom.pattern,
      chrom.match.mode = chrom.match.mode,
      shorten.col.names = shorten.col.names
    )
    ################### THIS IS WHERE...
    ############### we can try to fix that crash-on-empty-first-XIC error.
    if (length(s) == 0 ||
        nrow(s) == 0) {
      # if s is empty, do nothin
      # do nothing!
    } else if (nrow(summ) == 0) {
      #if summ is yet empty, set it to s
      summ = s
    } else {
      # otherwise, add missing columns and bind rows
      s = addabsentcols(s, summ)
      summ = addabsentcols(summ, s)
      summ = rbind(summ, s)
    }
  }
  write.csv(summ, "result1.csv")
  return(summ)
}


run.config.peak = function(path.to.config.tsv = 'config.tsv',
                           path.to.tics,
                           index = 1,
                           Interactive = TRUE,
                           chrom.pattern = "std|standard",
                           chrom.match.mode = "omit",
                           shorten.col.names = TRUE) {
  config = read.delim (path.to.config.tsv)
  config$C13[is.na(config$C13) | config$C13 == ""] <- 0
  config$N15[is.na(config$N15) | config$N15 == ""] <- 0
  config$H2[is.na(config$H2) | config$H2 == ""] <- 0
  config$interactive[grep("Yes|Inter|Okay|1|True|T",
                          config$interactive,
                          ignore.case = TRUE)] <- "yes"
  # check if all the files necessary exist in this place?
  notfound = 0
  for (i in 1:length(config[, 1])) {
    path = unlist(paste(path.to.tics, config$mz[i], sep = "/"))
    if (!file.exists(path)) {
      notfound = notfound + 1
      cat(path, " not found!\n")
    }
  }
  if (notfound > 0) {
    cat ("SOME FILES NOT FOUND, TRY CREATING THEM!\n")
    return()
  }
  # OK SO FAR? THE LET'S PROCEED!
  i <- index
  if (i <= length(config[, 1])) {
    path = unlist(paste(path.to.tics, config$mz[i], sep = "/"))
    #pngs in the current directory!
    
    # TODO: Here we can optionally enter a loop
    # in which we run get peaks and show the result,
    # and wait for the user to accept the result, or
    # update the parameters.
    # we should print the parameters and then ask which
    # needs to be updated, eg:
    #
    #     1. rt.min (12)
    #     2. rt.max (18)
    #     3. samples (20)
    #     4. threshold (1e6)
    #     0. Yes, it looks OK. Continue.
    #
    #     Pick one... > _
    #
    if (Interactive == TRUE && config$interactive[i] == "yes") {
      while (1) {
        gp = get.peaks(
          path = path,
          pattern = '.tsv',
          rt.min = config$rt.min[i],
          rt.max = config$rt.max[i],
          seconds = config$seconds[i],
          threshold = config$threshold[i],
          chrom.pattern = chrom.pattern,
          chrom.match.mode = chrom.match.mode,
          shorten.col.names = shorten.col.names
        )
        title(main = config$name[i], sub = "shading represents detected peak")
        
        cat(
          paste(
            "\n---------------------------------------\nPeak picking for ",
            config$target[i],
            "\n\n",
            sep = ""
          )
        )
        
        if (length(grep("peaks", names(gp))) == 0) {
          cat("\nThese setting did not give any peaks! Please try again... \n\n")
        }
        
        cat(paste("  1. rt.min(", config$rt.min[i], ")\n", sep = ""))
        cat(paste("  2. rt.max (", config$rt.max[i], ")\n", sep = ""))
        cat(paste("  3. seconds (", config$seconds[i], ")\n", sep = ""))
        cat(paste("  4. threshold (", config$threshold[i], ")\n", sep =
                    ""))
        cat("  a. Set all.\n")
        cat("  y. Yes, it looks OK. Continue.\n")
        cat("  q. quit.\n")
        
        usersays = readline("Pick one: ")
        if (usersays == "y") {
          if (length(grep("peaks", names(gp))) == 0) {
            cat("\nNo peaks.  Skipping\n\n")
          }
          break
        } else if (usersays == "q") {
          return()
        } else if (usersays == "a") {
          # set all at once
          line = readline(paste("enter new value for rt.min [", config$rt.min[i], "]: "))
          if (nchar(line) > 0) {
            config$rt.min[i] = as.double(line)
          }
          line = readline(paste("enter new value for rt.max [", config$rt.max[i], "]: "))
          if (nchar(line) > 0) {
            config$rt.max[i] = as.double(line)
            ############################## THIS IS WHERE
            ### we can immunize against rt min and max being inverted
            if (config$rt.max[i] <= config$rt.min[i]) {
              print("WARNING: rt.max < rt.min, setting rt.max = rt.min + 2")
              config$rt.max[i] <- config$rt.min[i] + 2
            }
          }
          line = readline(paste("enter new value for seconds [", config$seconds[i], "]: "))
          if (nchar(line) > 0) {
            config$seconds[i] = as.double(line)
          }
          line = readline(paste(
            "enter new value for threshold [",
            config$threshold[i],
            "]: "
          ))
          if (nchar(line) > 0) {
            config$threshold[i] = as.double(line)
          }
          ### set individually
        } else if (usersays == "1") {
          line = readline(paste("enter new value for rt.min [", config$rt.min[i], "]: "))
          if (nchar(line) > 0) {
            config$rt.min[i] = as.double(line)
            ############################## THIS IS WHERE
            ### we can immunize against rt min and max being inverted
            if (config$rt.max[i] <= config$rt.min[i]) {
              print("WARNING: rt.max < rt.min, setting rt.min = rt.max - 2")
              config$rt.min[i] <- max((config$rt.max[i] - 2), 0)
            }
          }
        } else if (usersays == "2") {
          line = readline(paste("enter new value for rt.max [", config$rt.max[i], "]: "))
          if (nchar(line) > 0) {
            config$rt.max[i] = as.double(line)
            ############################## THIS IS WHERE
            ### we can immunize against rt min and max being inverted
            if (config$rt.max[i] <= config$rt.min[i]) {
              print("WARNING: rt.max < rt.min, setting rt.max = rt.min + 2")
              config$rt.max[i] <- config$rt.min[i] + 2
            }
          }
        } else if (usersays == "3") {
          line = readline(paste("enter new value for seconds [", config$seconds[i], "]: "))
          if (nchar(line) > 0) {
            config$seconds[i] = as.double(line)
          }
        } else if (usersays == "4") {
          line = readline(paste(
            "enter new value for threshold [",
            config$threshold[i],
            "]: "
          ))
          if (nchar(line) > 0) {
            config$threshold[i] = as.double(line)
          }
        }
      }
      ############################### HERE IS WHERE
      ############## we can update config with parameters used for current row
    }
    else {
      gp = get.peaks(
        path = path,
        pattern = '.tsv',
        rt.min = config$rt.min[i],
        rt.max = config$rt.max[i],
        seconds = config$seconds[i],
        threshold = config$threshold[i],
        chrom.pattern = chrom.pattern,
        chrom.match.mode = chrom.match.mode,
        shorten.col.names = shorten.col.names
      )
    }
    
    if (length(grep("peaks", names(gp))) > 0) {
      ##################### HERE IS WHERE
      ########### we can immunize against illegal filename characters
      filename <- config$name[i]
      filename <- gsub("\\?", "(maybe)", filename)
      filename <- gsub("\\\\", "(bslash)", filename)
      filename <- gsub("/", "(fslash)", filename)
      filename <- gsub("\\*", "(aster)", filename)
      filename <- gsub("\\+", "(pos)", filename)
      png(unlist(paste("xic_", filename, ".png")))
      gp = get.peaks(
        path = path,
        pattern = '.tsv',
        rt.min = config$rt.min[i],
        rt.max = config$rt.max[i],
        seconds = config$seconds[i],
        threshold = config$threshold[i],
        chrom.pattern = chrom.pattern,
        chrom.match.mode = chrom.match.mode,
        shorten.col.names = shorten.col.names
      )
      
      
      title(main = config$name[i], sub = "shading represents detected peak")
      dev.off()
      s = gp$summary
      row.names(s) = paste(config$name[i], row.names(s))
      
      names(s) = sub("_[[:digit:]]+_[[:digit:]]+-[[:digit:]]+_[[:digit:]]+_",
                     "",
                     names(s))
      names(s) = sub("_[[:digit:]]+_[[:digit:]]+-[[:digit:]]+_",
                     "",
                     names(s))
      
      #print(config)
      #write.table (config, file=gsub(".tsv","-updated.tsv",path.to.config.tsv), sep="\t", row.names=FALSE)
      write.table (config,
                   file = path.to.config,
                   sep = "\t",
                   row.names = FALSE)
      return(s)
      
    }
  }
}

addabsentcols = function(addto, addfrom) {
  n1 = names(addto)
  L = length(addto[, 1])
  n2 = names(addfrom)
  for (name in n2) {
    if (length(grep(name, n1, fixed = TRUE)) == 0) {
      addto[name] = rep(NA, L)
    }
  }
  return(addto)
}


run.config = function(path.to.config.tsv = 'config.tsv',
                      path.to.mzMLs = 'mzMLneg/',
                      chrom.pattern = "std|standard",
                      chrom.match.mode = "omit",
                      shorten.col.names = TRUE) {
  starttime = Sys.time()
  ticpath = run.config.tics(path.to.config.tsv, path.to.mzMLs)
  print (Sys.time() - starttime)
  mysummary = run.config.peaks(
    path.to.config.tsv,
    ticpath,
    chrom.pattern = chrom.pattern,
    chrom.match.mode = chrom.match.mode,
    shorten.col.names = shorten.col.names
  )
  print (Sys.time() - starttime)
  return(mysummary)
}



# # biocLite("KEGGREST") to install
#
#
# mz = function(name=NULL,mass=NULL,z=-1,C13=0,N15=0){
#   if(is.null(mass)){
#     if(is.null(name)){
#       cat("You must specify either mass or name")
#       return(-2)
#     }
#     mass = unlist(lapply(name, kegg_exact_mass))
#   }
#   iC13 = 1.003355
#   iN15 = 0.9970348934
#   H = 1.0073
#   mass = mass + z*H
#   mass = mass + C13*iC13
#   mass = mass + N15*N15
#   if(z==0) z=1
#   mz = mass / abs(z)
#   mz[mz < 0] = round(mz[mz<0])
#   return(mz)
# }
#
# kegg_exact_mass = function(name){
#   r = keggFind("compound",name)
#   patt = unlist(paste("(^|; )",name,"(;|$)",sep=""))
#   g = grep(patt, r, ignore.case=TRUE)
#   if(length(g) < 1){
#     cat("No results :(\n")
#     return(-1)
#   }
#   k = unlist(lapply(names(r)[g], function(n){
#     k=keggGet(n)
#     return(as.numeric(k[[1]]$EXACT_MASS))
#   }))
#   if(min(k) != max(k)){
#     cat("More than one result, with different masses :(\n")
#     return()
#   }
#   return(k[1])
# }



### remember to add this to the main repo
msconvert_all = function() {
  # current director, extract all from raw... divide into pos and neg
  system(
    'msconvert --mzML --simAsSpectra --filter "peakPicking true 1-" -o mzMLpos --filter "polarity +" *.raw'
  )
  system(
    'msconvert --mzML --simAsSpectra --filter "peakPicking true 1-" -o mzMLneg --filter "polarity -" *.raw'
  )
}

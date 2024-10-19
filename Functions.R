#' @title Forest plot of MR results
#' @param mr_res mr_res or harmonized data of exposure(s) and outcome(s)
#' @param combine logical, combine multiple exposures and outcomes
#' @param x_trans Change axis scale, Allowed values are one of c("none", "log", "log2", "log10"). Default is "none",
#' @param wrap_str length of wrap labels
#' @param nSNP minimum number of SNPs.
#' @param col_pal a vector of colors
#' @param xlim xlims used in forest plot.
#' @param rmID logical, remove ID of exposure and outcome.
#' @param addCol add additional columns in mr_res, e.g., Q for Cochrane's Q, R2, F, Power, etc.
#' @param main title of plot.
#' @param arrow_lab Labels for the arrows, string vector of length two (left and right).
#' @param genType combine mr_res by `method` or `exposure`.
#' @param append logical, append results of `metagen` to `mr_res`.
#' @param append.color color of appended rows.
#' @param exponentiate Convert effects to OR? Default is TRUE.
#' @param ... other parameters used in \link[forestploter]{forest_theme} or \link[meta]{forest}.
#'
#' @family MR analysis
#' @export
#'
#' @return a data.frame of MR result, and associated data sets of xxx_mr_res.rda and xxx_mr_dat.rda will be saved to current working directory.
#'
forest_mr <- function(mr_res, combine = F, x_trans = 'none', wrap_str = 50, nSNP = 2,
                      col_pal = NULL, xlim = NULL, rmID = T, addCol = '', main = NULL, arrow_lab = NULL,
                      metagen = F, genType = 'exposure', append = F, append.color = 'blue', exponentiate = T,
                      footnotes = NULL, vert_line = NULL, ...) {
  suppressPackageStartupMessages({library(forestploter); library(ggplot2)})
  if (is.null(col_pal))
    col_pal <- c('#2BCE48FF', '#FF5005FF', "#F0A0FFFF", "#0075DCFF", '#993F00FF', '#4C005CFF')

  #
  {
    mr_res = mr_res[mr_res$nsnp > nSNP,]
    mr_res = generate_odds_ratios(mr_res)
    if (!exponentiate) {
      mr_res$or = mr_res$b
      mr_res$or_lci95 = mr_res$lo_ci
      mr_res$or_uci95 = mr_res$up_ci
    }
    dups = duplicated(paste0(mr_res$exposure, mr_res$outcome))
    addCol = addCol[addCol %in% names(mr_res)]
    if (length(addCol) > 0) mr_res[dups, addCol] = ''
    if (rmID) {
      mr_res$outcome = gsub(' \\|.*$', '', mr_res$outcome)
      mr_res$exposure = gsub(' \\|.*$', '', mr_res$exposure)
    }
    # combine with meta
    suppressPackageStartupMessages({library(meta)})
    mr_res$expo_out = paste0(mr_res$exposure, ' -> ', mr_res$outcome)
    sm = ifelse(exponentiate, 'OR', 'SMD')
    if (genType == 'method') {
      if (length(unique(mr_res$outcome)) == 1) {
        metaB = metagen(b, se, data = mr_res, studlab = exposure, sm = sm, byvar = method)
      } else {
        metaB = metagen(b, se, data = mr_res, studlab = expo_out, sm = sm, byvar = method)
        }
    } else {
      if (length(unique(mr_res$outcome)) == 1) {
        metaB = metagen(b, se, data = mr_res, studlab = method, sm = sm, byvar = exposure)
      } else {
        metaB = metagen(b, se, data = mr_res, studlab = method, sm = sm, byvar = expo_out)
      }
    }

    # combine with forestploter
    # mr_res = subset(mr_res, mr_res$exposure %in% unique(mr_res$exposure)[1:3])
    {
      mr_res$` ` <- paste0(rep(' ', 20), collapse = ' ')
      if (exponentiate) {
        mr_res$or = ifelse(mr_res$or > 1e3 | mr_res$or < 1e-3,
                           format(mr_res$or, scientific = T, digits = 3),
                           round(mr_res$or, 3))
        mr_res$or[which(as.numeric(mr_res$or) == 0)] = 0
        mr_res$or_lci95 = ifelse(mr_res$or_lci95 > 1e3 | mr_res$or_lci95 < 1e-3,
                                 format(mr_res$or_lci95, scientific = T, digits = 3),
                                 round(mr_res$or_lci95, 3))
        mr_res$or_lci95[which(as.numeric(mr_res$or_lci95) == 0)] = 0
        mr_res$or_uci95 = ifelse(mr_res$or_uci95 > 1e3 | mr_res$or_uci95 < 1e-3,
                                 format(mr_res$or_uci95, scientific = T, digits = 3),
                                 round(mr_res$or_uci95, 3))
        mr_res$or_uci95[which(as.numeric(mr_res$or_uci95) == 0)] = 0
      } else {
        mr_res$or = round(mr_res$or, 3)
        mr_res$or_lci95 = round(mr_res$or_lci95, 3)
        mr_res$or_uci95 = round(mr_res$or_uci95, 3)
      }
      mr_res$`OR (95% CI)` = paste0(mr_res$or, " (", mr_res$or_lci95, ", ", mr_res$or_uci95, ")")
      mr_res$or = as.numeric(mr_res$or)
      mr_res$or_lci95 = as.numeric(mr_res$or_lci95)
      mr_res$or_uci95 = as.numeric(mr_res$or_uci95)
      # mr_res$`OR (95% CI)` <- ifelse(is.na(mr_res$se), "",
      #                            sprintf("%.2f (%.2f to %.2f)",
      #                                    mr_res$or, mr_res$or_lci95, mr_res$or_uci95))
      mr_res$pSig = cut(as.numeric(mr_res$pval),
                        c(0, 0.001, 0.01, 0.05, 1), c('***', '**', '*', ''), right = F)
      mr_res$pval = format(mr_res$pval, scientific = T, digits = 3)
      mr_res$pSig[which(is.na(mr_res$pSig))] = ''
      mr_res$pval = paste0(mr_res$pval, mr_res$pSig)
      mr_res$`exposure(n)` = paste0(mr_res$exposure, ' (', mr_res$nsnp, ')')
      # mr_res$`exposure(n)`[which(mr_res$method != unique(mr_res$method)[1])] = ' '
      mr_res$`exposure(n)`[which(dups)] = ''
      mr_res$`outcome(n)` = paste0(mr_res$outcome, ' (', mr_res$nsnp, ')')
      mr_res$`outcome(n)`[which(dups)] = ''
      #
      if (x_trans %in% c("log", "log2", "log10")) {
        # do.call('log', list(10))
        # eval(call('log', 10))
        mr_res$or = ifelse(mr_res$or > 0, eval(call(x_trans, mr_res$or)), -eval(call(x_trans, mr_res$or)))
        mr_res$or_lci95 = ifelse(mr_res$or_lci95 > 0,
                                 eval(call(x_trans, mr_res$or_lci95)),
                                 -eval(call(x_trans, mr_res$or_lci95)))
        mr_res$or_uci95 = ifelse(mr_res$or_uci95 > 0,
                                 eval(call(x_trans, mr_res$or_uci95)),
                                 -eval(call(x_trans, mr_res$or_uci95)))
      } else x_trans = 'none'
    }
    # append results of meta
    if (append) {
      meteRes = meta.res(metaB, save = F, showPlot = T)
      meteRes = meteRes[meteRes$levels != '',]
      meteRes$or = lapply(1:nrow(meteRes), function(x) strsplit(meteRes$`ES[LCI; UCI]`[x], '\\[')[[1]][1]) %>% unlist() %>% as.numeric()
      meteRes$LUCI = lapply(1:nrow(meteRes), function(x) strsplit(meteRes$`ES[LCI; UCI]`[x], '\\[|\\]')[[1]][2]) %>% unlist()
      meteRes$or_lci95 = lapply(1:nrow(meteRes), function(x) strsplit(meteRes$LUCI[x], ';')[[1]][1]) %>% unlist() %>% as.numeric()
      meteRes$or_uci95 = lapply(1:nrow(meteRes), function(x) strsplit(meteRes$LUCI[x], ';')[[1]][2]) %>% unlist() %>% as.numeric()
      meteRes$`ES[LCI; UCI]` = paste0(meteRes$or, ' (', meteRes$or_lci95, ', ', meteRes$or_uci95, ')')
      #
      app.dat = data.frame(matrix(nrow = nrow(meteRes), ncol = ncol(mr_res)))
      names(app.dat) = names(mr_res)
      app.dat$or = meteRes$or
      app.dat$or_lci95 = meteRes$or_lci95
      app.dat$or_uci95 = meteRes$or_uci95
      app.dat$pval = as.character(meteRes$P.val)
      if (length(unique(mr_res$outcome)) == 1) {
        app.dat$exposure = meteRes$levels
        app.dat$outcome = mr_res$outcome[match(app.dat$exposure, mr_res$exposure)]
      } else {
        app.dat$exposure = unlist(lapply(meteRes$levels, function(x) strsplit(x, ' -> ')[[1]][1]))
        app.dat$outcome = unlist(lapply(meteRes$levels, function(x) strsplit(x, ' -> ')[[1]][2]))
      }
      app.dat$`OR (95% CI)` = meteRes$`ES[LCI; UCI]`
      # if ('P.het' %in% names(app.dat)) app.dat$P.het = meteRes$P.val.Q
      # app.dat$`exposure(n)` = mr_res$`exposure(n)`[match(app.dat$exposure, mr_res$exposure)]
      app.dat$`exposure(n)` = paste0("(I-square: ", meteRes$`I2(%)`, '%, pval.Q: ', meteRes$P.val.Q, ')')
      app.dat$`outcome(n)` = app.dat$`exposure(n)`
      app.dat$method = 'Pooled (I-V)'
      app.dat$se = meteRes$seTE
      #
      app.dat[is.na(app.dat)] = ''
      # joinVars = intersect(names(mr_res), names(app.dat))
      # mr_res = rbind(mr_res[,joinVars], app.dat[,joinVars])
      # or
      mr_res = merge(mr_res, app.dat, all = T, sort = F)
    }

    # combine forest ----
    mr_list = split(mr_res, ~outcome)
    if (length(mr_list) > 1 & combine) {
      if (append) {
        mr_res = mr_res[order(factor(mr_res$exposure, levels = unique(mr_res$exposure)),
                              factor(mr_res$method, levels = unique(mr_res$method))),]
        summaryRows = ifelse(mr_res$method == 'Pooled (I-V)', TRUE, FALSE)
        # mr_res$summary = summaryRows
        summaryRows = split(summaryRows, mr_res$outcome)
      } else summaryRows = NULL
      # combine with ggplot
      mr_res$ypos = 1:nrow(mr_res)
      labels = unique(mr_res$exposure)
      # breaks
      breaks <- 1:length(unique(mr_res$exposure)) * length(unique(mr_res$outcome)) * length(unique(mr_res$method))
      if (max(breaks) != nrow(mr_res)) {
        breaks <- which(!duplicated(mr_res$exposure))
        breaks <- c(breaks[-1] - 1, nrow(mr_res))
        # breaks[length(breaks)] <- nrow(mr_res)
      }
      # diff(which(!duplicated(mr_res$exposure)))
      # table(mr_res$exposure)
      # unique(mr_res$outcome)
      # showColor()
      mr_res$sig = ifelse(grepl('\\*', mr_res$pSig), '"*"', '')
      mr_res$method = factor(mr_res$method, unique(mr_res$method))
      mr_res$outcome = factor(mr_res$outcome, names(mr_list))
      p = ggplot(mr_res, aes(x = or, y = ypos, color = outcome)) +
        geom_vline(xintercept = ifelse(x_trans == 'none', 1, 0), size = 2, alpha = .5, color = "grey50") +
        geom_segment(aes(x = or_lci95,
                         xend = or_uci95,
                         yend = ypos), size = 1, alpha = .5) +
        geom_point(aes(shape = method), size = 4) +
        scale_color_manual(values = col_pal) +
        scale_y_continuous(breaks = breaks,
                           labels = stringr::str_wrap(labels, width = wrap_str), expand = c(.05, .05)) +
        theme_minimal(base_family = "Times New Roman", base_size = 12) +
        theme(legend.position = "right",
              panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_line(size = 4, color = "grey90"),
              axis.text.y = element_text(vjust = .3, size = 12)) +
        ggtitle('Estimated causal effects of exposure(s) on outcome(s)') +
        ylab(' ') + xlab(ifelse(exponentiate, 'OR (95% CI)', 'ES (95% CI)'))

      p <- p + annotate('text', x = mr_res$or, y = mr_res$ypos, col = 'black',
                        label = mr_res$sig, vjust = .7, size = 5, parse = T)
      if (x_trans != 'none') p <- p + xlab(paste0(x_trans, '_OR (95% CI)'))
      if (any(c(mr_res$or > 100, mr_res$or_uci95 > 100))) p <- p + scale_x_log10() + xlab('log10_OR (95% CI)')
      print(p)

      # check paired list
      NumLen = sapply(names(mr_list), function(x) nrow(mr_list[[x]]))
      if (length(unique(NumLen)) != 1) {
        message('Length of list were not equal. \nYou can set "combine = F" to see all results.')
        # Nc = names(which.max(NumLen))
        # for (i in Nc) mr_list[[i]] = mr_list[[i]][mr_list[[2]]$method %in% mr_list[[1]]$method,]
        Allm = intersect(mr_list[[which.max(NumLen)]]$method, mr_list[[which.min(NumLen)]]$method)
        mr_res = mr_res[mr_res$method %in% Allm,]
      }
      # combine multiple plots
      tm <- forest_theme(base_size = 10,
                         refline_lty = "solid",
                         ci_col = col_pal[1:length(mr_list)],
                         footnote_col = "blue",
                         legend_name = "Group",
                         legend_value = names(mr_list),
                         ...)
      mr_dat = unique(mr_res[,c('exposure', 'method', ' ', 'nsnp')])
      mr_dat$`exposure(n)` = paste0(mr_dat$`exposure`, ' (', mr_dat$nsnp, ')')
      mr_dat = mr_dat[!duplicated(mr_dat[,c('exposure', 'method')]),]
      mr_dat$`exposure(n)`[which(duplicated(mr_dat$`exposure`))] = ' '
      mr_dat$`exposure(n)` = stringr::str_wrap(mr_dat$`exposure(n)`, width = wrap_str)
      if (any(duplicated(mr_dat$`exposure(n)`[mr_dat$`exposure(n)` != ' ' & mr_dat$method != 'Pooled (I-V)']))) Colecho('Duplicates found, set `rmID = F` and run it again.')
      if (append) {
        mr_dat = mr_dat[order(factor(mr_dat$exposure, levels = unique(mr_dat$exposure)),
                              factor(mr_dat$method, levels = unique(mr_dat$method))),]
      }
      # mr_res$Sig = cut(as.numeric(mr_res$pval),
      #                   c(0, 0.001, 0.01, 0.05, 1), c('***', '**', '*', '-'), right = F)
      #
      # if (any(c(mr_res$or_lci95, mr_res$or_uci95, mr_res$or) == 0)) {
      #   message('"zore" should not be included in transformated outcome. Used default parameter "none".')
      #   x_trans = 'none'
      # } else x_trans = trans
      if (is.null(xlim)) {
        xlim = c(floor(range(mr_res$or)[1]), ceiling(range(mr_res$or)[2]))
        diffOR = max(mr_res$or_uci95 - mr_res$or_lci95)
        rangeOr = range(c(mr_res$or_lci95, mr_res$or_uci95))
        if (diffOR > 10) {
          ORs = c(mr_res$or_lci95, mr_res$or_uci95)
          rangeOr = range(ORs[!ORs %in% rangeOr])
          xlim = c(floor(rangeOr[1]), ceiling(rangeOr[2]))
        } else if (diffOR < 1) {
          xlim = c(rangeOr[1] * .95, rangeOr[2] * 1.05)
        }
      }
      if (is.null(footnotes)) {
        footnotes = ifelse(x_trans == 'none', ' ', paste0(x_trans, ' transformed'))
      } else {
        transX = ifelse(x_trans == 'none', ' ', paste0(x_trans, ' transformed'))
        footnotes = paste0(footnotes, '\n', transX)
      }
      #
      # sizes <- sqrt(1/mr_res$se)
      # sizes[is.infinite(sizes)] = NA
      # mr_res$sizes = sizes / max(sizes, na.rm = TRUE)
      mr_res$sizes = scale(abs(mr_res$or), F)
      if (is.null(main)) main = 'Estimated causal effects of exposure(s) on outcome(s)'
      if (is.null(arrow_lab)) arrow_lab = c("Negative", "Positive")
      forestploter::forest(mr_dat[,c('exposure(n)', 'method', ' ')],
                           est = split(mr_res$or, mr_res$outcome),
                           lower = split(mr_res$or_lci95, mr_res$outcome),
                           upper = split(mr_res$or_uci95, mr_res$outcome),
                           sizes = split(mr_res$sizes, mr_res$outcome),
                           # is_summary = summaryRows,
                           arrow_lab = arrow_lab, xlim = xlim,
                           title = main, ci_column = 3,
                           ref_line = ifelse(x_trans == 'none' & exponentiate, 1, 0),
                           vert_line = vert_line,
                           footnote = footnotes,
                           #x_trans = x_trans, # only suit to number > 0
                           theme = tm) %>%
        add_border(part = "header", row = 1, where = "bottom") %>%
        # Edit background of multiple exposures
        edit_plot(row = which(mr_dat$`exposure(n)` != ''), which = "background",
                  gp = grid::gpar(fill = "grey80")) %>% print()
    } else {
      # single forest ----
      tm <- forest_theme(...)
      # check length of exposure
      if (length(unique(mr_res$exposure)) == 1 & length(unique(mr_res$outcome)) > 1) {
        Colecho('One exposure compared with more than one outcomes.')
        mr_dat = mr_res
        mr_dat$`outcome(n)` = stringr::str_wrap(mr_dat$`outcome(n)`, width = wrap_str)
        if (any(duplicated(mr_dat$`outcome(n)`[mr_dat$`outcome(n)` != '' & mr_dat$method != 'Pooled (I-V)']))) Colecho('Duplicates found, set `rmID = F` and run it again.')
        # xlim
        if (is.null(xlim)) {
          xlim = c(floor(range(mr_dat$or)[1]), ceiling(range(mr_dat$or)[2]))
          diffOR = max(mr_dat$or_uci95 - mr_dat$or_lci95)
          rangeOr = range(c(mr_dat$or_lci95, mr_dat$or_uci95))
          if (diffOR > 10) {
            ORs = c(mr_res$or_lci95, mr_res$or_uci95)
            rangeOr = range(ORs[!ORs %in% rangeOr])
            xlim = c(floor(rangeOr[1]), ceiling(rangeOr[2]))
          } else if (diffOR < 1) {
            xlim = c(rangeOr[1] * .9, rangeOr[2] * 1.1)
          }
        }
        # highlight
        if (append) {
          mr_dat = mr_dat[order(factor(mr_dat$outcome, levels = unique(mr_dat$outcome)),
                                factor(mr_dat$method, levels = unique(mr_dat$method))),]
          summaryRows = ifelse(mr_dat$method == 'Pooled (I-V)', TRUE, FALSE)
        } else summaryRows = NULL
        mrk = grep('\\*', mr_dat$pval)
        if (length(mrk) != 0) mrk = mrk else mrk = NULL
        # plot
        # sizes <- sqrt(1/mr_dat$se)
        # sizes[is.infinite(sizes)] = NA
        # sizes <- sizes/max(sizes, na.rm = TRUE)
        if (is.null(main)) title = paste0('Estimated causality of outcomes for exposure\n', unique(mr_res$exposure))
        if (is.null(arrow_lab)) arrow_lab = c("Negative", "Positive")
        sel_cols = c('outcome(n)', 'method', ' ', 'OR (95% CI)', 'pval', addCol)
        if (!exponentiate) {
          mr_dat$`ES (95% CI)` = mr_dat$`OR (95% CI)`
          sel_cols = c('outcome(n)', 'method', ' ', 'ES (95% CI)', 'pval', addCol)
        }
        if (is.null(footnotes)) {
          footnotes = ifelse(x_trans == 'none', ' ', paste0(x_trans, ' transformed'))
        } else {
          transX = ifelse(x_trans == 'none', ' ', paste0(x_trans, ' transformed'))
          footnotes = paste0(footnotes, '\n', transX)
        }
        #
        forestploter::forest(mr_dat[,sel_cols],
                             est = mr_dat$or,
                             lower = mr_dat$or_lci95,
                             upper = mr_dat$or_uci95,
                             sizes = scale(abs(mr_dat$or), F),
                             ci_column = 3,
                             is_summary = summaryRows,
                             # x_trans = x_trans,
                             arrow_lab = arrow_lab,
                             title = ifelse(is.null(main), title, main),
                             ref_line = ifelse(x_trans == 'none' & exponentiate, 1, 0),
                             vert_line = vert_line, xlim = xlim,
                             footnote = footnotes,
                             ticks_at = NULL, theme = tm) %>%
          add_border(part = "header", row = 1, where = "bottom") %>%
          # Edit fontface of results
          edit_plot(row = which(mr_dat$method == 'Pooled (I-V)'), col = c(1, 2, 4, 5), which = "text",
                    gp = grid::gpar(fontface = "italic", col = append.color)) %>%
          edit_plot(row = mrk, col = ifelse(is.null(mrk), 6, 5), which = "text",
                    gp = grid::gpar(fontface = "bold")) %>%
          # Edit background of multiple exposures
          edit_plot(row = which(mr_dat$`outcome(n)` != '' & mr_dat$method != 'Pooled (I-V)'),
                    which = "background", gp = grid::gpar(fill = "grey80")) %>% print()
        #
        return('Reverse them automatically.')
      }
      # seperate exposures by more outcomes
      for (i in seq_along(mr_list)) {
        mr_dat = mr_list[[i]]
        mr_dat$`exposure(n)` = stringr::str_wrap(mr_dat$`exposure(n)`, width = wrap_str)
        if (any(duplicated(mr_dat$`exposure(n)`[mr_dat$`exposure(n)` != '' & mr_dat$method != 'Pooled (I-V)']))) Colecho('Duplicates found, set `rmID = F` and run it again.')
        # if (any(c(mr_dat$or_lci95, mr_dat$or_uci95, mr_dat$or) == 0)) {
        #   message('"zore" should not be included in transformated outcome. Used default parameter "none".')
        #   x_trans = 'none'
        # } else x_trans = trans
        # if (x_trans %in% c("log", "log2", "log10")) {
        #   breakN = base::pretty(c(mr_dat$or_lci95, mr_dat$or_uci95), 3) #[2:4]
        #   xlim = NULL
        # } else breakN = NULL
        # {
        #   breakN = base::pretty(c(mr_dat$or_lci95, mr_dat$or_uci95), 3) #[2:4]
        # } else breakN = NULL
        if (is.null(xlim)) {
          xlim = c(floor(range(mr_dat$or)[1]), ceiling(range(mr_dat$or)[2]))
          diffOR = max(mr_dat$or_uci95 - mr_dat$or_lci95)
          rangeOr = range(c(mr_dat$or_lci95, mr_dat$or_uci95))
          if (diffOR > 10) {
            ORs = c(mr_res$or_lci95, mr_res$or_uci95)
            rangeOr = range(ORs[!ORs %in% rangeOr])
            xlim = c(floor(rangeOr[1]), ceiling(rangeOr[2]))
          } else if (diffOR < 1) {
            xlim = c(rangeOr[1] * .9, rangeOr[2] * 1.1)
          }
        }
        # rangeOr = range(c(mr_dat$or_lci95, mr_dat$or_uci95))
        # if (rangeOr[1] > xlim[1] & rangeOr[2] < xlim[2]) xlim = rangeOr
        # highlight
        if (append) {
          mr_dat = mr_dat[order(factor(mr_dat$exposure, levels = unique(mr_dat$exposure)),
                                factor(mr_dat$method, levels = unique(mr_dat$method))),]
          # mrk = grep('\\*', mr_dat$pval)
          # mrk = c(mrk, which(mr_dat$pval != 'ns' & mr_dat$method == 'Pooled (I-V)'))
          # not work on ci_*
          # ci_pchs = rep(15, length(unique(mr_dat$method)))
          # ci_pchs = rep(15, nrow(mr_dat))
          # ci_pchs = ci_pchs[which(unique(mr_dat$method) == 'Pooled (I-V)')] = 16
          # ci_cols = rep('black', length(unique(mr_dat$method)))
          # ci_cols = rep('black', nrow(mr_dat))
          # ci_cols = ci_cols[which(unique(mr_dat$method) == 'Pooled (I-V)')] = 'red'
          # tm <- forest_theme(ci_pch = ci_pchs, ci_col = ci_cols, ...)
          summaryRows = ifelse(mr_dat$method == 'Pooled (I-V)', TRUE, FALSE)
        } else summaryRows = NULL
        mrk = grep('\\*', mr_dat$pval)
        if (length(mrk) != 0) mrk = mrk else mrk = NULL
        # plot
        # sizes <- sqrt(1/mr_dat$se)
        # sizes[is.infinite(sizes)] = NA
        # sizes <- sizes/max(sizes, na.rm = TRUE)
        if (is.null(main)) title = paste0('Estimated causal effects of exposure(s) on\n', names(mr_list)[i])
        if (is.null(arrow_lab)) arrow_lab = c("Negative", "Positive")
        sel_cols = c('exposure(n)', 'method', ' ', 'OR (95% CI)', 'pval', addCol)
        if (!exponentiate) {
          mr_dat$`ES (95% CI)` = mr_dat$`OR (95% CI)`
          sel_cols = c('exposure(n)', 'method', ' ', 'ES (95% CI)', 'pval', addCol)
        }
        if (is.null(footnotes)) {
          footnotes = ifelse(x_trans == 'none', ' ', paste0(x_trans, ' transformed'))
        } else {
          transX = ifelse(x_trans == 'none', ' ', paste0(x_trans, ' transformed'))
          footnotes = paste0(footnotes, '\n', transX)
        }
        #
        forestploter::forest(mr_dat[,sel_cols],
                             est = mr_dat$or,
                             lower = mr_dat$or_lci95,
                             upper = mr_dat$or_uci95,
                             sizes = scale(abs(mr_dat$or), F),
                             ci_column = 3,
                             is_summary = summaryRows,
                             # x_trans = x_trans,
                             # arrow_lab = c("Outcome", "Exposure"),
                             arrow_lab = arrow_lab,
                             title = ifelse(is.null(main), title, main),
                             ref_line = ifelse(x_trans == 'none' & exponentiate, 1, 0),
                             vert_line = vert_line, xlim = xlim,
                             footnote = footnotes,
                             ticks_at = NULL, theme = tm) %>%
          add_border(part = "header", row = 1, where = "bottom") %>%
          # Edit fontface of results
          edit_plot(row = which(mr_dat$method == 'Pooled (I-V)'), col = c(1, 2, 4, 5), which = "text",
                    gp = grid::gpar(fontface = "italic", col = append.color)) %>%
          edit_plot(row = mrk, col = ifelse(is.null(mrk), 6, 5), which = "text",
                    gp = grid::gpar(fontface = "bold")) %>%
          # Edit background of multiple exposures
          edit_plot(row = which(mr_dat$`exposure(n)` != '' & mr_dat$method != 'Pooled (I-V)'),
                    which = "background", gp = grid::gpar(fill = "grey80")) %>% print()
          # insert_text('Beta', col = 4, part = 'header', just = 'left',
          #             gp = grid::gpar(fontface = 'bold')) %>% print()
        # gpar(col = "red", fontface = "bold", fill = "darkolivegreen1")
        # median(as.numeric(c(mr_dat$or_lci95, mr_dat$or_uci95)))
        # fivenum(as.numeric(c(mr_dat$or_lci95, mr_dat$or_uci95)))
        # breakN = base::pretty(c(mr_dat$or_lci95, mr_dat$or_uci95), 5)[2:4]
      }
    }
  }
}


#' @title MR plot
#' @description Integrate plots of MR: scatter_plot, forest_plot, funnel_plot, leaveoneout_plot, density_plot, rucker_jackknife, plot_radial, etc.
#' @param dat,mr_res harmonized `dat` and result of `mr()`.
#' @param type `scatter_plot`, `forest_plot`, `funnel_plot`, `leaveoneout_plot`, `density_plot`, `radial_plot`, `jackknife`.
#' @param all Combine all plots to one file if study was one-arm. Or all plots will be merged by group of plot types.
#' @param rmID logical, remove ID of exposure and outcome.
#' @param add.genes logical, add nearest genes of SNP.
#' @param outliers logical, show outliers name on plot.
#' @param scatter,interactive Plot of scatter. Default plot by \link[TwoSampleMR]{mr_scatter_plot}, or set to `MR` to use \link[MendelianRandomization]{mr_plot}, `MRall` to use all methods. `interactive`: logical, interactive plot.
#' @param radial type of plot: `egger` or `ivw`. Default if `egger`.
#' @param ... other parameters used in \link[cowplot]{plotgrid}.
#'
#' @family MR analysis
#' @export
#'
MR_plot <- function(dat, mr_res, type = "scatter_plot", all = F, rmID = T,
                    add.genes = F, outliers = F,
                    scatter = "", interactive = F, radial = "egger", ...) {
    if (missing(dat) | missing(mr_res)) stop("Please provide dat and mr_res.")
    p1 <- p2 <- p3 <- p4 <- p5 <- p6 <- NULL
    if (rmID) {
        dat$outcome <- gsub(" \\|.*$", "", dat$outcome)
        dat$exposure <- gsub(" \\|.*$", "", dat$exposure)
        mr_res$outcome <- gsub(" \\|.*$", "", mr_res$outcome)
        mr_res$exposure <- gsub(" \\|.*$", "", mr_res$exposure)
    } else {
        dat$outcome <- gsub(" \\|\\|", "\n", dat$outcome)
        dat$exposure <- gsub(" \\|\\|", "\n", dat$exposure)
        mr_res$outcome <- gsub(" \\|\\|", "\n", mr_res$outcome)
        mr_res$exposure <- gsub(" \\|\\|", "\n", mr_res$exposure)
    }
    # add.genes
    if (add.genes) {
        geneRes <- chromSNP(dat, dat$SNP, Genes = T)
        dat$genes <- sapply(dat$SNP, function(s) geneRes$transcripts[[s]][[1]][1]) %>% unlist()
        dat$SNP <- paste0(dat$SNP, "(", dat$genes, ")")
    }
    # scatter
    if (type == "scatter_plot" | all) {
        Colecho("scatter_plot:", cols = "green")
        if (scatter %in% c("MR", "MRall")) {
            onLoadpkg("MendelianRandomization")
            mr_dat <- dat_to_MRInput(dat, get_correlations = F)
            #
            if (scatter == "MR") {
                p1 <- lapply(seq_along(mr_dat), function(x) {
                    mr_plot(mr_input(
                        bx = mr_dat[[x]]@betaX, bxse = mr_dat[[x]]@betaXse,
                        by = mr_dat[[x]]@betaY, byse = mr_dat[[x]]@betaYse,
                        # correlation = mr_dat[[1]]@correlation,
                        exposure = mr_dat[[x]]@exposure,
                        outcome = mr_dat[[x]]@outcome,
                        snps = mr_dat[[x]]@snps
                    ),
                    interactive = interactive, labels = T, line = "egger"
                    )
                })
            }
            #
            if (scatter == "MRall") {
                p1 <- lapply(seq_along(mr_dat), function(x) {
                    mr_plot(mr_allmethods(mr_input(
                        bx = mr_dat[[x]]@betaX, bxse = mr_dat[[x]]@betaXse,
                        by = mr_dat[[x]]@betaY, byse = mr_dat[[x]]@betaYse,
                        exposure = mr_dat[[x]]@exposure,
                        outcome = mr_dat[[x]]@outcome
                    ), method = "all"))
                })
            }
        } else {
            p1 <- mr_scatter_plot(mr_res, dat)
        }
        #
        # if (!all) print(p1[[1]])
        if (length(p1) > 1) {
            message("The result including ", length(p1), " plots.")
            resAsk <- readline("Do you want show all of them? (Y/N)")
            if (resAsk %in% c("y", "Y")) showAll <- T else showAll <- F
            if (all | showAll) cowplot::plot_grid(plotlist = p1, ...) %>% print() else print(p1[[1]])
        } else {
            print(p1[[1]])
        }
    }
    #
    if (type %in% c("forest_plot", "density_plot", "funnel_plot") | all) {
        res_single <- mr_singlesnp(dat)
        # forest_plot
        if (type == "forest_plot" | all) {
            Colecho("scatter_plot:", cols = "green")
            p2 <- mr_forest_plot(res_single)
            # if (!all) print(p2[[1]])
            if (length(p2) > 1) {
                message("The result including ", length(p2), " plots.")
                resAsk <- readline("Do you want show all of them? (Y/N)")
                if (resAsk %in% c("y", "Y")) showAll <- T else showAll <- F
                if (all | showAll) cowplot::plot_grid(plotlist = p2, ...) %>% print() else print(p2[[1]])
            } else {
                print(p2[[1]])
            }
        }
        # funnel_plot
        if (type == "funnel_plot" | all) {
            Colecho("funnel_plot:", cols = "green")
            p4 <- mr_funnel_plot(res_single)
            # if (!all) print(p4[[1]])
            if (length(p4) > 1) {
                message("The result including ", length(p4), " plots.")
                resAsk <- readline("Do you want show all of them? (Y/N)")
                if (resAsk %in% c("y", "Y")) showAll <- T else showAll <- F
                if (all | showAll) cowplot::plot_grid(plotlist = p4, ...) %>% print() else print(p4[[1]])
            } else {
                print(p4[[1]])
            }
        }
        # density_plot
        if (type == "density_plot" | all) {
            Colecho("density_plot:", cols = "green")
            p5 <- mr_density_plot(res_single, mr_res)
            # if (!all) print(p5[[1]])
            if (length(p5) > 1) {
                message("The result including ", length(p5), " plots.")
                resAsk <- readline("Do you want show all of them? (Y/N)")
                if (resAsk %in% c("y", "Y")) showAll <- T else showAll <- F
                if (all | showAll) cowplot::plot_grid(plotlist = p5, ...) %>% print() else print(p5[[1]])
            } else {
                print(p5[[1]])
            }
        }
    }
    # leave-one-out analysis
    if (type == "leaveoneout_plot" | all) {
        Colecho("leaveoneout_plot:", cols = "green")
        res_loo <- mr_leaveoneout(dat)
        p3 <- mr_leaveoneout_plot(res_loo)
        # if (!all) print(p3[[1]])
        if (length(p3) > 1) {
            message("The result including ", length(p3), " plots.")
            resAsk <- readline("Do you want show all of them? (Y/N)")
            if (resAsk %in% c("y", "Y")) showAll <- T else showAll <- F
            if (all | showAll) cowplot::plot_grid(plotlist = p3, ...) %>% print() else print(p3[[1]])
        } else {
            print(p3[[1]])
        }
    }
    # plot_radial
    if (type == "radial_plot" | all) {
        Colecho("radial_plot:", cols = "green")
        if (length(unique(dat$exposure)) == 1) dat_list <- split(dat, ~outcome) # outcome as list name
        if (length(unique(dat$outcome)) == 1) dat_list <- split(dat, ~exposure) # exposure as list name
        if (length(unique(dat$outcome)) > 1 & length(unique(dat$exposure)) > 1) {
            dat$expo_out <- paste0(dat$exposure, " on\n", dat$outcome)
            dat_list <- split(dat, ~expo_out)
            # dat_list = split(dat, ~exposure + outcome)
        }
        #
        if (radial == "ivw") {
            p6 <- lapply(seq_along(dat_list), function(x) {
                Colecho(names(dat_list)[x])
                RadialP <- RadialMR::plot_radial(RadialMR::ivw_radial(dat_list[[x]])) +
                    labs(subtitle = names(dat_list)[x])
                # tag outliers
                if (outliers) {
                    RadialP$data$labSNP <- RadialP$data$SNP
                    RadialP$data$labSNP[which(RadialP$data$Outliers != "Outlier")] <- ""
                    # RadialP + geom_label(aes(label = RadialP$data$labSNP))
                    # RadialP + geom_label_repel(aes(label = RadialP$data$labSNP))
                    RadialP <- RadialP + annotate("text",
                        x = RadialP$data$Wj, y = RadialP$data$BetaWj,
                        label = RadialP$data$labSNP, vjust = 1.5, size = 3
                    )
                }
                return(RadialP)
            })
        } else {
            p6 <- lapply(seq_along(dat_list), function(x) {
                Colecho(names(dat_list)[x])
                RadialP <- RadialMR::plot_radial(RadialMR::egger_radial(dat_list[[x]])) +
                    labs(subtitle = names(dat_list)[x])
                # tag outliers
                if (outliers) {
                    RadialP$data$labSNP <- RadialP$data$SNP
                    RadialP$data$labSNP[which(RadialP$data$Outliers != "Outlier")] <- ""
                    # RadialP + geom_label(aes(label = RadialP$data$labSNP))
                    # RadialP + geom_label_repel(aes(label = RadialP$data$labSNP))
                    RadialP <- RadialP + annotate("text",
                        x = RadialP$data$Wj, y = RadialP$data$BetaWj,
                        label = RadialP$data$labSNP, vjust = 1.5, size = 3
                    )
                }
                return(RadialP)
            })
        }
        # if (!all) print(p6[[1]])
        if (length(p6) > 1) {
            message("The result including ", length(p6), " plots.")
            resAsk <- readline("Do you want show all of them? (Y/N)")
            if (resAsk %in% c("y", "Y")) showAll <- T else showAll <- F
            if (all | showAll) cowplot::plot_grid(plotlist = p6, ...) %>% print() else print(p6[[1]])
        } else {
            print(p6[[1]])
        }
    }
    # mr_rucker_jackknife
    if (type == "jackknife" | all) {
        Colecho("rucker_jackknife_plot:", cols = "green")
        jk <- mr_rucker_jackknife(dat)
        # jk[[1]]$q_plot
        if (!all) print(jk[[1]]$e_plot)
        if (length(jk) > 1) message("The result including ", length(jk), " plots.")
    }
    #
    if (all & length(p1) == 1 & length(p4) == 1) {
        cowplot::plot_grid(p1[[1]], p4[[1]], p5[[1]], p6[[1]], rel_heights = c(3, 2)) %>% print()
        cowplot::plot_grid(p2[[1]], p3[[1]]) %>% print()
    }
}



#' @title Pan-survival analysis in TCGA database
#' @param tumorS user specified tumors abbreviation, default is 33 tumor types.
#' @param genes user specified gene names.
#' @param Surv type of survival analysis: OS, PFS, or DFS.
#' @param KMplot plot style of 'gg' or 'base', show Kaplan-Meier survival plot.
#' @param ncol number columns to show in `KMplot`.
#' @param forestplot whether show forestplot.
#' @param time survival time limited days. Default is 10 years.
#' @param matchVars Matching for Causal Inference \link[MatchIt]{MatchIt}: 'Age' and/or 'Gender'.
#' @param ExprSurv logical, only return a list of `gExprSurv` data.
#' @param subDat select subset group of *Var1* in definition table of sample codes. Default is c('01', '03', '06').
#' * 01	Primary Solid Tumor
#' * 03	Primary Blood Derived,Cancer-Peripheral Blood
#' * 06	Metastatic
#' @param subPheno add Phenotypes in data.
#' @param ... other parameters used in \link[forestplot]{forestplot}
#' @export
#'
panSurv <- function(tumorS, genes, Surv = 'OS',
                    KMplot = 'base', forestplot = F,
                    ncol = 3, time = 365 * 50, matchVars = NULL,
                    ExprSurv = F, subDat = c('01', '03', '06'),
                    subPheno = NULL, ...) {

  onLoadpkg(c("survival", "forestplot", "survminer"))

  genes = unique(genes)
  tumorN = list.dirs("E:/Scientist career/TCGA_database", recursive = F, full.names = F)
  tumorN = tumorN[nchar(tumorN) <= 4]
  if (missing(tumorS)) tumorS = tumorN
  tumorS = toupper(tumorS)
  if (!any(tumorS %in% tumorN)) {
    cat("\n", tumorS[!tumorS %in% tumorN], "were not in database.")
    tumorS = tumorS[tumorS %in% tumorN]
  }
  #
  if (length(tumorS) == 0) stop("None tumor available!!")
  surDat = list()
  for (tm in tumorS) {
    Datlist = LoadTCGA(tm, matchPat = T)
    Phenotype = Datlist[["Phenotype"]]
    genes.f = genes[genes %in% names(Datlist[["Expr"]])]
    if (length(genes.f) > 0) {
      tmDat <- subset(Datlist[["Expr"]], select = genes.f)
      # subset Expr and Phenotype data
      if (!is.null(subPheno)) {
        #Phenotype = subset(Phenotype, subset = subPheno)
        #tmDat = tmDat[match(Phenotype$sampleID, rownames(tmDat)), ]
        ssubPheno = subPheno[subPheno %in% names(Phenotype)]
        tmDat[, ssubPheno] = Phenotype[, ssubPheno]
      }
      #
      if (Surv == "OS") {
        gExprSurv = cbind(tmDat, Phenotype[, c("OS", "OS.time")])
        gExprSurv$Stime = gExprSurv$`OS.time`
        gExprSurv$Sx = gExprSurv$OS
        main = "Overall Survival"
      }
      if (Surv == "PFS" | Surv == 'PFI') {
        gExprSurv = cbind(tmDat, Phenotype[, c("PFI", "PFI.time")])
        gExprSurv$Stime = gExprSurv$`PFI.time`
        gExprSurv$Sx = gExprSurv$PFI
        main = 'Progression Free Interval'
      }
      # match age/gender
      if (!is.null(matchVars)) {
        if (!require(MatchIt)) install.packages("MatchIt")
        matchVars = tolower(matchVars)
        if ("age" %in% matchVars) {
          library(dplyr)
          gExprSurv$age = Phenotype[, grep("age", names(Phenotype))[1]] %>% as.numeric()
        }
        if ("gender" %in% matchVars) {
          gExprSurv$gender = Phenotype[, grep("gender", names(Phenotype))[1]] %>% as.factor()
        }
        if ("age" %in% matchVars) gExprSurv = gExprSurv[!is.na(gExprSurv$age), ]
        if ("gender" %in% matchVars) gExprSurv = gExprSurv[!is.na(gExprSurv$gender), ]
        gExprSurv = gExprSurv[!is.na(gExprSurv$Sx), ]
        model = as.formula(paste0("Sx ~", paste0(matchVars, collapse = " + ")))
        f <- matchit(model, data = gExprSurv, method = "nearest", caliper = 0.1, ratio = 1)
        gExprSurv <- match.data(f)
      }
      #
      gExprSurv = gExprSurv[gExprSurv$`Stime` <= time,]
      if (!is.null(subDat)) gExprSurv = gExprSurv[substring(rownames(gExprSurv), 14, 15) %in% subDat, ]
      # tmDat$OS = Phenotype$OS[match(rownames(tmDat), Phenotype$sampleID)]
      # tmDat$`OS.time` = Phenotype$`OS.time`[match(rownames(tmDat), Phenotype$sampleID)]
      surDat[[tm]] = gExprSurv
    }
  }

  # return gExprSurv
  if (ExprSurv) return(surDat)

  # single tumor
  if (length(tumorS) == 1) {
    if (forestplot) {
      univ_res = NULL
      gExprSurv = surDat[[tumorS]]
      for (g in genes) {
        if (g %in% names(gExprSurv)) {
          univ_f = as.formula(paste("Surv(Stime, Sx)~", g))
          univ_m = coxph(univ_f, data = gExprSurv)
          x <- summary(univ_m)
          HR <- signif(x$coef[2], 2)
          lower <- signif(x$conf.int[, "lower .95"], 2)
          upper <- signif(x$conf.int[, "upper .95"], 2)
          univ_res = rbind.data.frame(univ_res, data.frame(
            Symbol = g,
            p.value = signif(x$wald["pvalue"], 3),
            wald.test = signif(x$wald["test"], 2),
            beta = signif(x$coef[1], 2),
            lower = lower, upper = upper,
            HR = paste0(HR, " (", lower, "-", upper, ")")
          ))
        }
      }
      #
      test_data <- data.frame(
        coef = exp(as.numeric(univ_res$beta)),
        lower = as.numeric(univ_res$lower),
        upper = as.numeric(univ_res$upper)
      )
      test_data <- rbind(rep(NA, 3), test_data)
      row_names <- cbind(
        c("Symbol", univ_res$Symbol),
        c("Hazard Ratio (95% CI)", univ_res$HR),
        c("P value", univ_res$p.value)
      )

      forestplot::forestplot(
        labeltext = row_names,
        test_data[, c("coef", "lower", "upper")],
        graph.pos = 4, hrzl_lines = T, # zero = 1,
        title = paste0(tumorS, " (Cox proportional hazards regression model)"),
        is.summary = c(T, rep(F, nrow(test_data))), ...
      ) %>% print()
    } else {
      # K-M survival plot
      if (KMplot == "gg") {
        gExprSurv = surDat[[tumorS]]
        genes.f = genes[genes %in% names(gExprSurv)]
        ps = lapply(genes.f, function(g) {
          gExprSurv$group = ifelse(gExprSurv[, g] > median(gExprSurv[, g], na.rm = T), "High", "Low")
          data.survdiff = survdiff(Surv(Stime, Sx) ~ group, data = gExprSurv)
          sfit <- survfit(Surv(Stime, Sx) ~ group, data = gExprSurv)
          p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
          ggp = ggsurvplot(sfit,
            data = gExprSurv, palette = "npg", risk.table = F, conf.int = F,
            surv.median.line = c("v"),
            title = paste0(Surv, " of gene ", g),
            font.main = c(12, "bold", "black"),
            subtitle = "Kaplan-Meier estimates",
            font.submain = c(10, "italic", "gray"),
            xlab = "Survival time(days)",
            font.x = c(10, "bold", "black"),
            font.y = c(10, "bold", "black"),
            font.tickslab = c(8, "plain", "black"),
            legend = c(0.7, 0.8), # "top",#
            legend.title = "",
            legend.labs = c("High", "Low"), size = 0.5,
            pval = round(p.val, 3),
            pval.method = T
          )
          return(ggp$plot)
        })
        #
        nrow = length(ps) %/% ncol
        nrow = ifelse(length(ps) %% ncol == 0, nrow, nrow + 1)
        layoutM = matrix(1:(nrow * ncol), byrow = T, nrow = nrow, ncol = ncol)
        mv = gridExtra::marrangeGrob(ps, ncol = ncol, nrow = nrow, top = tumorS, layout_matrix = layoutM)
        print(mv)
      }
      #
      genes.f <- genes[genes %in% names(gExprSurv)]
      nrow <- length(genes.f) %/% ncol
      nrow <- ifelse(length(genes.f) %% ncol == 0, nrow, nrow + 1)
      if (KMplot == "base") GepiaSurv(gExprSurv, genes.f, main = main, ncol = ncol, nrow = nrow)
    }
    #
    return(NULL)
  }

  # multiple tumors
  for (g in genes) {
    if (forestplot) {
      univ_res = NULL
      for (tm in names(surDat)) {
        gExprSurv = surDat[[tm]]
        if (g %in% names(gExprSurv)) {
          univ_f = as.formula(paste("Surv(Stime, Sx)~", g))
          univ_m = coxph(univ_f, data = gExprSurv)
          x <- summary(univ_m)
          HR <- signif(x$coef[2], 2)
          lower <- signif(x$conf.int[, "lower .95"], 2)
          upper <- signif(x$conf.int[, "upper .95"], 2)
          univ_res = rbind.data.frame(univ_res, data.frame(
            tumor = tm,
            p.value = signif(x$wald["pvalue"], 3),
            wald.test = signif(x$wald["test"], 2),
            beta = signif(x$coef[1], 2),
            lower = lower, upper = upper,
            HR = paste0(HR, " (", lower, "-", upper, ")")
          ))
        }
      }
      #
      test_data <- data.frame(
        coef = exp(as.numeric(univ_res$beta)),
        lower = as.numeric(univ_res$lower),
        upper = as.numeric(univ_res$upper)
      )
      test_data <- rbind(rep(NA, 3), test_data)
      row_names <- cbind(
        c("Cancer", univ_res$tumor),
        c("Hazard Ratio (95% CI)", univ_res$HR),
        c("P value", univ_res$p.value)
      )
      #
      forestplot::forestplot(
        labeltext = row_names,
        test_data[, c("coef", "lower", "upper")],
        graph.pos = 4, hrzl_lines = T,
        title = paste0(g, " (Cox proportional hazards regression model)"),
        is.summary = c(T, rep(F, nrow(test_data))), ...
      ) %>% print()
      return(NULL)
    }
  }

  # K-M survival plot
  if (KMplot == 'gg') {
    #  multiple tumors
    for (g in genes) {
      ps = lapply(tumorS, function(tm) {
        gExprSurv = surDat[[tm]]
        if (g %in% names(gExprSurv)) {
          gExprSurv$group = ifelse(gExprSurv[, g] > median(gExprSurv[, g], na.rm = T), "High", "Low")
          data.survdiff = survdiff(Surv(Stime, Sx) ~ group, data = gExprSurv)
          sfit <- survfit(Surv(Stime, Sx) ~ group, data = gExprSurv)
          p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
          ggp = ggsurvplot(sfit,
            data = gExprSurv, palette = "npg", risk.table = F, conf.int = F,
            surv.median.line = c("v"),
            title = paste0(Surv, " of gene ", g),
            font.main = c(12, "bold", "black"),
            subtitle = "Kaplan-Meier estimates",
            font.submain = c(10, "italic", "gray"),
            xlab = "Survival time(days)",
            font.x = c(10, "bold", "black"),
            font.y = c(10, "bold", "black"),
            font.tickslab = c(8, "plain", "black"),
            legend = c(0.7, 0.8), # "top",#
            legend.title = "",
            legend.labs = c("High", "Low"), size = 0.5,
            pval = round(p.val, 3),
            pval.method = T
          )
          return(ggp$plot)
        }
      })
      #
      nrow = length(ps) %/% ncol
      nrow = ifelse(length(ps) %% ncol == 0, nrow, nrow + 1)
      layoutM = matrix(1:(nrow * ncol), byrow = T, nrow = nrow, ncol = ncol)
      mv = gridExtra::marrangeGrob(ps, ncol = ncol, nrow = nrow, top = tm, layout_matrix = layoutM)
      print(mv)
      #ggsave(filename, mv, width = nc * 3, height = nr * 3)
    }
  }
  #
  if (KMplot == 'base') {
    for (tm in tumorS) {
      gExprSurv = surDat[[tm]]
      genes.f = genes[genes %in% names(gExprSurv)]
      nrow = length(genes.f) %/% ncol
      nrow = ifelse(length(genes.f) %% ncol == 0, nrow, nrow + 1)
      GepiaSurv(gExprSurv, genes.f, main = main, ncol = ncol, nrow = nrow)
    }
  }
}


#' @title Add percent
#' @param var variable name
#' @param data data.frame
#' @param top label top groups
#'
Addper <- function(var, data, top = 1, digit = 0, sepL = "\n") {
    pos <- which(names(data) == var)
    topTar <- na.omit(names(sort(table(data[, pos]), T))[1:top])
    data[, pos][!data[, pos] %in% topTar] <- "Other"
    for (Tar in c(topTar, "Other")) {
        ADper <- round(length(which(data[, pos] == Tar)) / nrow(data) * 100, digit)
        data[, pos][which(data[, pos] == Tar)] <- paste0(Tar, sepL, "(", ADper, "%)")
    }
    return(data)
}

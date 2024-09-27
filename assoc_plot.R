#!/usr/bin/env Rscript
###
#  The MIT License
#
#  Copyright (C) 2021-2024 Giulio Genovese
#
#  Author: Giulio Genovese <giulio.genovese@gmail.com>
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
###

options(error = function() {traceback(3); q()})

assoc_plot_version <- '2024-09-27'

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
if (capabilities()[['cairo']]) options(bitmapType = 'cairo')

parser <- OptionParser('usage: assoc_plot.R [options] --genome <GRCh37|GRCh38>|--cytoband <cytoband.txt.gz> --vcf|--tbx <file>')
parser <- add_option(parser, c('--genome'), type = 'character', help = 'genome assembly (e.g. GRCh38)', metavar = '<assembly>')
parser <- add_option(parser, c('--cytoband'), type = 'character', help = 'cytoband file', metavar = '<cytoband.txt.gz>')
parser <- add_option(parser, c('--nauto'), type = 'integer', default = 22, help = 'number of autosomes', metavar = '<integer>')
parser <- add_option(parser, c('--vcf'), type = 'character', help = 'input VCF file', metavar = '<file.vcf>')
parser <- add_option(parser, c('--pheno'), type = 'character', help = 'phenotype to select from GWAS-VCF file', metavar = '<string>')
parser <- add_option(parser, c('--as'), action = 'store_true', default = FALSE, help = 'input VCF file has allelic shift information')
parser <- add_option(parser, c('--csq'), action = 'store_true', default = FALSE, help = 'whether coding variant should be flagged as red')
parser <- add_option(parser, c('--tbx'), type = 'character', help = 'input REGENIE/PLINK summary statistics', metavar = '<file.gz>')
parser <- add_option(parser, c('--region'), type = 'character', help = 'region to plot', metavar = '<region>')
parser <- add_option(parser, c('--min-af'), type = 'double', default = 0.0, help = 'minimum minor allele frequency [0]', metavar = '<float>')
parser <- add_option(parser, c('--min-lp'), type = 'integer', help = 'minimum -log10 p-val [2]', metavar = '<integer>')
parser <- add_option(parser, c('--loglog-pval'), type = 'integer', default = 10, help = '-log10 p-val threshold for using log-log scale in manhattan plot', metavar = '<integer>')
parser <- add_option(parser, c('--cyto-ratio'), type = 'integer', default = 25, help = 'plot to cytoband ratio [25]', metavar = '<integer>')
parser <- add_option(parser, c('--max-height'), type = 'integer', help = 'phred score ceiling', metavar = '<integer>')
parser <- add_option(parser, c('--spacing'), type = 'integer', default = 20, help = 'spacing between chromosomes [10]', metavar = '<integer>')
parser <- add_option(parser, c('--pdf'), type = 'character', help = 'output PDF file', metavar = '<file.pdf>')
parser <- add_option(parser, c('--png'), type = 'character', help = 'output PNG file', metavar = '<file.png>')
parser <- add_option(parser, c('--width'), type = 'double', default = 7.0, help = 'inches width of the output file [7.0]', metavar = '<float>')
parser <- add_option(parser, c('--height'), type = 'double', help = 'inches height of the output file [3.5/7.0]', metavar = '<float>')
parser <- add_option(parser, c('--fontsize'), type = 'integer', default = 12, help = 'font size [12]', metavar = '<integer>')

args <- parse_args(parser, commandArgs(trailingOnly = TRUE), convert_hyphens_to_underscores = TRUE)

write(paste('assoc_plot.R', assoc_plot_version, 'http://github.com/freeseek/score'), stderr())

if (is.null(args$genome) && is.null(args$cytoband)) {print_help(parser); stop('either --genome or --cytoband is required\nTo download the cytoband file run:\nwget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz')}
if (!is.null(args$genome) && !is.null(args$cytoband)) {print_help(parser); stop('cannot use --genome and --cytoband at the same time')}
if (!is.null(args$genome) && args$genome != 'GRCh37' && args$genome != 'GRCh38') {print_help(parser); stop('--genome accepts only GRCh37 or GRCh38')}
if (is.null(args$vcf) && is.null(args$tbx)) {print_help(parser); stop('either --vcf or --tbx is required')}
if (!is.null(args$vcf) && !is.null(args$tbx)) {print_help(parser); stop('cannot use --vcf and --tbx at the same time')}
if (is.null(args$vcf) && !is.null(args$pheno))  {print_help(parser); stop('--pheno requires --vcf')}
if (is.null(args$vcf) && args$as) {print_help(parser); stop('--as requires --vcf')}
if (args$as && args$min_af > 0) {print_help(parser); stop('--as cannot be used with --min-af')}
if (is.null(args$vcf) && is.null(args$as_vcf) && args$csq)  {print_help(parser); stop('--csq requires --vcf')}
if (is.null(args$pdf) && is.null(args$png)) {print_help(parser); stop('either --pdf or --png is required')}
if (!is.null(args$pdf) && !is.null(args$png)) {print_help(parser); stop('cannot use --pdf and --png at the same time')}
if (!is.null(args$png) && !capabilities('png')) {print_help(parser); stop('unable to start device PNG: no png support in this version of R\nyou need to reinstall R with support for PNG to use the --png option')}

if (is.null(args$min_lp)) {
  if (is.null(args$region)) {
    args$min_lp <- 2
  } else {
    args$min_lp <- 0
  }
}

if (!is.null(args$cytoband)) {
  df_cyto <- setNames(read.table(args$cytoband, sep = '\t', header = FALSE), c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))
  df_cyto$chrom <- gsub('chr', '', df_cyto$chrom)
  chrlen <- tapply(df_cyto$chromEnd, df_cyto$chrom, max)
  chrs <- unique(df_cyto$chrom)
  modified_chrs <- gsub('^M[T]?$', args$nauto + 4, gsub('^Y$', args$nauto + 2, gsub('^X$', args$nauto + 1, chrs)))
  ord <- order(suppressWarnings(as.numeric(modified_chrs)))
  chrs <- chrs[ord]

  idx_p <- df_cyto$gieStain == 'acen' & substr(df_cyto$name, 1, 3) == 'p11'
  idx_q <- df_cyto$gieStain == 'acen' & substr(df_cyto$name, 1, 3) == 'q11'
  if (sum(idx_p) > 0 && sum(idx_q)) {
    df_cen <- rbind(cbind(setNames(df_cyto[idx_p, c('chrom', 'name', 'chromStart')], c('chrom', 'name', 'x')), y = -1),
                    cbind(setNames(df_cyto[idx_p, c('chrom', 'name', 'chromEnd')], c('chrom', 'name', 'x')), y = -1/2),
                    cbind(setNames(df_cyto[idx_p, c('chrom', 'name', 'chromStart')], c('chrom', 'name', 'x')), y = 0),
                    cbind(setNames(df_cyto[idx_q, c('chrom', 'name', 'chromEnd')], c('chrom', 'name', 'x')), y = -1),
                    cbind(setNames(df_cyto[idx_q, c('chrom', 'name', 'chromStart')], c('chrom', 'name', 'x')), y = -1/2),
                    cbind(setNames(df_cyto[idx_q, c('chrom', 'name', 'chromEnd')], c('chrom', 'name', 'x')), y = 0))
  }

  # idx <- df_cyto$gieStain %in% c('acen', 'gvar', 'stalk')
  # if (sum(idx) > 0) {
  #   cen_beg <- tapply(df_cyto$chromEnd[idx], df_cyto$chrom[idx], min)
  #   cen_end <- tapply(df_cyto$chromEnd[idx], df_cyto$chrom[idx], max)
  #   df_chrs <- data.frame(chrlen = chrlen[chrs], cen_beg = cen_beg[chrs], cen_end = cen_end[chrs], CHROM = chrs)
  # }

  chrlen <- chrlen[c(1:args$nauto, 'X')]
} else if ( args$genome == 'GRCh37' ) {
  chrlen <- setNames(c(249251621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63026520, 48129895, 51305566, 155270560), c(1:22,'X'))
} else if ( args$genome == 'GRCh38' ) {
  chrlen <- setNames(c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895), c(1:22,'X'))
}

if ( !is.null(args$vcf) ) {
  fmt <- '%CHROM\\t%POS'
  names <- c('chrom', 'pos')
  if (args$as) {
    fmt <- paste0(fmt, '\\t%AS{0}\\t%AS{1}')
    names <- c(names, 'as0', 'as1')
  } else {
    fmt <- paste0(fmt, '[\\t%LP{0}]')
    names <- c(names, 'lp')
  }
  if (!args$as && (args$min_af>0 || args$min_lp>0)) {
    opt_filter <- ' --include \''
    if (args$min_af>0) opt_filter <- paste0(opt_filter, 'AF>', args$min_af, ' & AF<1-', args$min_af)
    if (args$min_af>0 && args$min_lp>0) opt_filter <- paste0(opt_filter, ' & ')
    if (args$min_lp>0) opt_filter <- paste0(opt_filter, 'LP>=', args$min_lp)
    opt_filter <- paste0(opt_filter, '\'')
  } else if (args$as) {
    opt_filter <- ' --include \'sum(AS)>0\''
  } else opt_filter <- ''
  if (!is.null(args$region)) opt_regions <- paste0(' --regions ', args$region) else opt_regions <- ''
  if (!is.null(args$pheno)) opt_samples <- paste0(' --samples ', args$pheno) else opt_samples <- ''
  if (args$csq) {
    fmt <- paste0(fmt, '\\t%Consequence')
    names <- c(names, 'consequence')
    cmd <- paste0('bcftools +split-vep --output-type u --columns Consequence', opt_filter, opt_regions, ' "', args$vcf, '" | bcftools query --format "', fmt, '\\n"', opt_samples)
  } else {
    cmd <- paste0('bcftools query --format "', fmt, '\\n"', opt_filter, opt_regions, opt_samples, ' "', args$vcf, '"')
  }
} else {
  if ( !is.null(args$region) ) {
    cmd <- past0e('tabix --print-header "', args$tbx, '" ', strsplit(args$region,','))
  } else {
    cmd <- paste0('zcat "', args$tbx, '"')
  }
  if (!is.null(args$min_af)) {
    filter <- paste0('$af>', args$min_af, ' && $af<', 1-args$min_af, ' && $lp!="NA" && $lp>=', args$min_lp)
  } else {
    filter <- paste0('$lp!="NA" && $lp>=', args$min_lp)
  }
  cmd <- paste(cmd, '| awk \'NR==1 {for (i=1; i<=NF; i++) f[$i] = i; if ("CHROM" in f) chrom=f["CHROM"]; else chrom=f["#CHROM"]; if ("GENPOS" in f) pos=f["GENPOS"]; else pos=f["POS"]; if ("A1FREQ" in f) af=f["A1FREQ"]; else af=f["A1_FREQ"]; if ("NEG_LOG10_P" in f) lp=f["NEG_LOG10_P"]; else if ("LOG10P" in f) lp=f["LOG10P"]; else lp=f["LOG10_P"]} NR==1 || NR>1 &&', filter , '{print $chrom"\\t"$pos"\\t"$lp}\'')
  names <- c('chrom', 'pos', 'lp')
}

write(paste('Command:', cmd), stderr())
if (packageVersion('data.table') < '1.11.6') {
  df <- setNames(fread(cmd, sep = '\t', header = TRUE, na.strings = '.', colClasses = list(character = c(1)), data.table = FALSE), names)
} else {
  df <- setNames(fread(cmd = cmd, sep = '\t', header = TRUE, na.strings = '.', colClasses = list(character = c(1)), data.table = FALSE), names)
}

if (args$as) {
  df$lp = -(log(2) + pbinom(pmin(df$as0, df$as1), df$as0 + df$as1, .5, log.p = TRUE)) / log(10)
  df <- df[df$lp >= args$min_lp & !is.na(df$lp),]
} else {
  df <- df[!is.na(df$lp),]
}

if (nrow(df) == 0) stop('nothing to be plotted')

df$chrom <- as.factor(gsub('^chr', '', gsub('^chrM', 'MT', df$chrom)))
ord <- order(as.numeric(gsub('MT', '26', gsub('Y', '24', gsub('X', '23', levels(df$chrom))))))
df$chrom <- factor(df$chrom, levels(df$chrom)[ord])

if (is.null(args$height)) {
  if (length(unique(df$chrom)) == 1) {
    args$height <- 7.0
  } else {
    args$height <- 3.5
  }
}

if (!is.null(args$max_height)) {
  max_height <- max(args$max_height, -log10(5e-8))
} else {
  max_height <- max(df$lp[!is.infinite(df$lp)], -log10(5e-8), na.rm = TRUE)
}

# see http://github.com/FINNGEN/saige-pipelines/blob/master/scripts/qqplot.R
max_loglog <- args$loglog_pval * log10(max_height) / log10(args$loglog_pval)
df$lp[df$lp > args$loglog & !is.na(df$lp)] = args$loglog_pval * log10(df$lp[df$lp > args$loglog_pval & !is.na(df$lp)]) / log10(args$loglog_pval)
tick_pos <- round(seq(1, max_loglog, length.out = 10))
tick_lab <- sapply(tick_pos, function(x) { round(ifelse(x < args$loglog_pval, x, args$loglog_pval^(x/args$loglog_pval))) })

if (length(unique(df$chrom)) == 1) {
  p <- ggplot(df, aes(x = pos/1e6, y = lp)) +
    geom_hline(yintercept = -log10(5e-8), color = 'gray', lty = 'longdash') +
    geom_point(size = 1/2) +
    scale_x_continuous(paste('Chromosome', unique(df$chrom), '(Mbp position)'), expand = c(.01,.01)) +
    scale_y_continuous('-log10(p-value)', breaks = tick_pos, labels = tick_lab, expand = c(.01,.01)) +
    theme_bw(base_size = args$fontsize)
} else {
  df$chrompos <- cumsum(c(0,args$spacing * 1e6 + chrlen))[as.numeric(df$chrom)] + df$pos
  p <- ggplot(df, aes(x = chrompos/1e6, y = lp, color = as.numeric(chrom)%%2 == 1)) +
    geom_hline(yintercept = -log10(5e-8), color = 'gray', lty = 'longdash') +
    geom_point(size = 1/2) +
    scale_x_continuous(NULL, breaks = (cumsum(args$spacing * 1e6 + chrlen) - chrlen/2 - args$spacing * 1e6) / 1e6, labels = names(chrlen), expand = c(.01,.01)) +
    scale_y_continuous('-log10(p-value)', breaks = tick_pos, labels = tick_lab, expand = c(.01,.01)) +
    scale_color_manual(guide = 'none', values = c('FALSE' = 'dodgerblue', 'TRUE' = 'gray')) +
    theme_bw(base_size = args$fontsize)
}

# this errors out when df is empty
if (args$csq) {
  df$coding <- FALSE
  for (annotation in c('missense', 'inframe', 'protein_altering', 'transcript_amplification', 'exon_loss', 'disruptive', 'start_lost', 'stop_lost', 'stop_gained', 'frameshift', 'splice_acceptor', 'splice_donor', 'transcript_ablation')) {
    df$coding <- df$coding | grepl(annotation, df$consequence)
  }
  p <- p + geom_point(data = df[df$coding, ], size = 1/3, color = 'red')
}

if (!is.null(args$cytoband) && length(unique(df$chrom)) == 1) {
  cyto_height <- (max_loglog - args$min_lp) / args$cyto_ratio
  p <- p  +
    geom_rect(data = df_cyto[df_cyto$chrom == unique(df$chrom) & df_cyto$gieStain != 'acen',], aes(x = NULL, y = NULL, xmin = chromStart/1e6, xmax = chromEnd/1e6, fill = gieStain, shape = NULL), ymin = args$min_lp - cyto_height, ymax = args$min_lp, color = 'black', linewidth = 1/4, show.legend = FALSE) +
    scale_fill_manual(values = c('gneg' = 'white', 'gpos25' = 'lightgray', 'gpos50' = 'gray50', 'gpos75' = 'darkgray', 'gpos100' = 'black', 'gvar' = 'lightblue', 'stalk' = 'slategrey'))
  if (exists('df_cen')) {
    p <- p + geom_polygon(data = df_cen[df_cen$chrom == unique(df$chrom),], aes(x = x/1e6, y = args$min_lp + cyto_height * y, shape = NULL, group = name), color = 'black', fill = 'red', linewidth = 1/8)
  }
} else {
  cyto_height <- 0
}

xlim <- FALSE
if ( !is.null(args$region) ) {
  if (grepl(':', args$region) && grepl('-', args$region)) {
    left <- as.numeric(gsub('^.*:', '', gsub('-.*$', '', args$region)))
    right <- as.numeric(gsub('^.*-', '', args$region))
    xlim <- TRUE
  }
}

if (!is.null(args$max_height) && xlim) {
  p <- p + coord_cartesian(xlim = c(left, right)/1e6, ylim = c(args$min_lp - cyto_height, max_loglog)/10)
} else if (xlim) {
  p <- p + coord_cartesian(xlim = c(left, right)/1e6)
} else if (!is.null(args$max_height)) {
  p <- p + coord_cartesian(ylim = c(args$min_lp - cyto_height, max_loglog)/10)
}

if (!is.null(args$pdf)) {
  pdf(args$pdf, width = args$width, height = args$height)
} else {
  png(args$png, width = args$width, height = args$height, units = 'in', res = 150)
}
print(p)
invisible(dev.off())

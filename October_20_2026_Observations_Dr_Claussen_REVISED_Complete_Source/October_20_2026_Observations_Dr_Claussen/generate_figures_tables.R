# October 20, 2026 Observations from Dr. Claussen
# FIGURE/TABLE REPRODUCTION SCRIPT (Base R graphics and calculations only)
# R packages: none beyond R's installed base distributions.
# Works from a self-contained source_data directory (recommended) or
# the original Windows GCV_Full_Field_Project stage folders.
# This script READS SAVED STAGE 1-3C RESULTS; it does not rerun GAM fits.

options(stringsAsFactors = FALSE)
local_data <- 'source_data'
windows_root <- 'C:/Users/pbbs/Claussen/GCV_Full_Field_Project'
output_dir <- 'generated_figures_tables'
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

get_input <- function(name, stage) {
  p <- file.path(local_data, name)
  if (!file.exists(p)) p <- file.path(windows_root, stage, name)
  if (!file.exists(p)) stop('Missing data file: ', name,
      '\nPut the saved CSV in source_data/ or the original stage directory.')
  p
}
read_stage <- function(name, stage) {
  read.csv(get_input(name, stage), check.names=FALSE)
}
manual_grid <- read_stage('Manual_Full_Field_Grid.csv','02_Manual_Full_Field')
manual_summary <- read_stage('Manual_Full_Field_Summary.csv','02_Manual_Full_Field')
manual_minima <- read_stage('Manual_Candidate_Minima.csv','02_Manual_Full_Field')
mgcv_grid <- read_stage('Stage3C_Full_Grid.csv','03C_MGCV_Fixed_Grid_Full_Field')
mgcv_summary <- read_stage('Stage3C_Auto_vs_Grid.csv','03C_MGCV_Fixed_Grid_Full_Field')
mgcv_minima <- read_stage('Stage3C_Interior_Grid_Dips.csv','03C_MGCV_Fixed_Grid_Full_Field')
pilot <- read_stage('Stage3A_Auto_vs_Grid.csv','03A_MGCV_Pilot')

stopifnot(nrow(manual_grid)==79600L, nrow(manual_summary)==199L,
          nrow(mgcv_grid)==119400L, nrow(mgcv_summary)==199L,
          all(sort(unique(manual_grid$Easting))==seq(2,398,2)),
          all(sort(unique(mgcv_grid$Easting))==seq(2,398,2)))

# Choose any available original taxonomy file. Otherwise reconstruct the
# directly supported flags without altering the original stage results.
taxname <- 'GCV_Stage3C_Anomaly_Taxonomy_All_Eastings.csv'
taxpath <- file.path(local_data,taxname)
if(file.exists(taxpath)) {
  taxonomy <- read.csv(taxpath,stringsAsFactors=FALSE)
} else {
  taxonomy <- data.frame(
    Easting=mgcv_summary$Easting,
    Multiple_Minima=mgcv_summary$Interior_Grid_Dip_Count>=2,
    Auto_Worse_Grid=mgcv_summary$Auto_Minus_Best_Grid_GCV>1e-4,
    Auto_Worse_Linear=mgcv_summary$Auto_Worse_Than_Linear,
    Auto_Minus_Best_Grid_Beta=mgcv_summary$Auto_Minus_Best_Grid_Beta)
}

# Run each base-graphics plot twice for print and preview formats.
render_figure <- function(stem, draw, width=11, height=6) {
  grDevices::pdf(file.path(output_dir,paste0(stem,'.pdf')),
                 width=width,height=height,onefile=FALSE)
  draw(); grDevices::dev.off()
  grDevices::png(file.path(output_dir,paste0(stem,'.png')),
                 width=round(width*135),height=round(height*135),res=135)
  draw(); grDevices::dev.off()
}

# Build rectangular grids with rows=Easting and columns=log(lambda).
as_grid <- function(d, xcol, ycol='GCV') {
  d <- d[order(d$Easting,d[[xcol]]),,drop=FALSE]
  xs <- sort(unique(d[[xcol]])); ys <- sort(unique(d$Easting))
  if(nrow(d) != length(xs)*length(ys)) stop('Grid not rectangular')
  list(xs=xs,ys=ys,Z=matrix(d[[ycol]],nrow=length(ys),
                              ncol=length(xs),byrow=TRUE))
}
plot_heatmap <- function(d, xcol, summary_data, linear_col, stem, method) {
  z <- as_grid(d,xcol)
  idx <- match(z$ys,summary_data$Easting)
  lin <- summary_data[[linear_col]][idx]
  delta <- sweep(z$Z,1,lin,'-')
  draw <- function() {
    par(mfrow=c(1,2),mar=c(4.5,4.2,4,3))
    image(z$xs,z$ys,t(z$Z),xlab=paste0('log(lambda): ',method),
          ylab='Easting',main='Absolute GCV',
          col=hcl.colors(96,'Viridis'),useRaster=TRUE)
    image(z$xs,z$ys,t(delta),xlab=paste0('log(lambda): ',method),
          ylab='Easting',main='GCV minus exact linear limit',
          col=hcl.colors(96,'Blue-Red 3'),useRaster=TRUE)
    par(mfrow=c(1,1))
  }
  render_figure(stem,draw,width=12.4,height=5.1)
}
plot_heatmap(manual_grid,'log_lambda',manual_summary,
             'GCV_Linear_Limit','F01_manual_absolute_relative_heatmaps','manual')
plot_heatmap(mgcv_grid,'LogLambda',mgcv_summary,
             'Exact_Linear_GCV','F02_mgcv_absolute_relative_heatmaps','mgcv')

# Four complete sampled GCV curves, with automatic, grid, and linear comparators.
render_figure('F03_four_mgcv_GCV_profiles',function() {
  par(mfrow=c(2,2),mar=c(4,4,3,1))
  for(e in c(214,216,218,374)) {
    x <- mgcv_grid[mgcv_grid$Easting==e,]
    s <- mgcv_summary[mgcv_summary$Easting==e,]
    plot(x$LogLambda,x$GCV,type='l',lwd=1.8,
         xlab='log(lambda_mgcv)',ylab='GCV',main=paste('Easting',e))
    points(s$Automatic_LogLambda,s$Automatic_GCV,pch=4,cex=1.5)
    points(s$Best_Grid_LogLambda,s$Best_Grid_GCV,pch=19,cex=.9)
    abline(h=s$Exact_Linear_GCV,lty=2)
    legend('topright',c('Automatic (x)','Lowest grid (dot)','Linear limit (dash)'),
           pch=c(4,19,NA),lty=c(NA,NA,2),bty='n',cex=.63)
  }
  par(mfrow=c(1,1))
},width=11.5,height=7.3)

# Individual Easting 374 plot with the two independently detected minima.
z374 <- mgcv_grid[mgcv_grid$Easting==374,]
d374 <- mgcv_minima[mgcv_minima$Easting==374,]
s374 <- mgcv_summary[mgcv_summary$Easting==374,]
render_figure('F04_Easting_374_two_minima',function() {
  plot(z374$LogLambda,z374$GCV,type='l',lwd=2,
       xlab='log(lambda_mgcv) = log(sp)',ylab='Absolute GCV',
       main='Easting 374: two competing GCV minima')
  points(d374$LogLambda,d374$GCV,pch=19,cex=1.3)
  text(d374$LogLambda,d374$GCV,
       labels=sprintf('GCV %.3f, log(lambda) %.3f',d374$GCV,d374$LogLambda),
       pos=c(1,3)[seq_len(nrow(d374))],cex=.75)
  points(s374$Automatic_LogLambda,s374$Automatic_GCV,pch=4,cex=1.4,lwd=2)
  abline(h=s374$Exact_Linear_GCV,lty=2)
  legend('top',c('Interior minima','Automatic x','Linear limit dashed'),
         pch=c(19,4,NA),lty=c(NA,NA,2),bty='n',cex=.72)
},width=10,height=4.8)

render_figure('F05_Easting_374_GCV_RSS_EDF',function() {
  par(mfrow=c(2,1),mar=c(4,4,3,4))
  plot(z374$LogLambda,z374$GCV,type='l',lwd=2,
       xlab='',ylab='GCV',main='GCV rises toward the linear-limit plateau')
  abline(h=s374$Exact_Linear_GCV,lty=2)
  plot(z374$LogLambda,z374$RSS,type='l',lwd=2,
       xlab='log(lambda_mgcv)',ylab='RSS',main='Fit deterioration versus complexity reduction')
  par(new=TRUE)
  plot(z374$LogLambda,z374$Smooth_EDF,type='l',lty=2,axes=FALSE,
       xlab='',ylab='')
  axis(4);mtext('Smooth EDF',side=4,line=2.5)
  par(mfrow=c(1,1))
},width=10,height=7)

render_figure('F06_product_coefficient_trajectories',function() {
  par(mfrow=c(2,1),mar=c(4,4,3,1))
  plot(manual_summary$Easting,manual_summary$Lowest_Candidate_Product_B,
       type='l',xlab='',ylab='Product B',main='Stage 2: manual spline')
  lines(manual_summary$Easting,manual_summary$Product_B_Linear_Limit,lty=2)
  legend('topright',c('Lowest listed candidate','Exact linear'),lty=1:2,cex=.8,bty='n')
  plot(mgcv_summary$Easting,mgcv_summary$Automatic_Product_B,
       type='l',xlab='Easting',ylab='Product B',main='Stage 3C: mgcv')
  lines(mgcv_summary$Easting,mgcv_summary$Best_Grid_Product_B,lty=2)
  lines(mgcv_summary$Easting,mgcv_summary$Exact_Linear_Product_B,lty=3)
  legend('topright',c('Automatic','Lowest grid','Exact linear'),lty=1:3,cex=.8,bty='n')
  par(mfrow=c(1,1))
},width=12,height=8)

render_figure('F07_optimizer_gaps_vs_coefficient_sensitivity',function() {
  par(mfrow=c(2,1),mar=c(4,5,3,1))
  plot(mgcv_summary$Easting,mgcv_summary$Auto_Minus_Best_Grid_GCV,
       type='b',pch='.',xlab='',ylab='Automatic - grid GCV',main='GCV optimization gap')
  abline(h=1e-4,lty=2)
  plot(mgcv_summary$Easting,mgcv_summary$Auto_Minus_Best_Grid_Beta,
       type='b',pch='.',xlab='Easting',ylab='Automatic - grid Product B',
       main='Coefficient consequence')
  abline(h=0,lty=2);par(mfrow=c(1,1))
},width=11.6,height=7.5)

render_figure('F08_manual_structural_gap',function() {
  plot(manual_summary$Easting,
       manual_summary$Structural_Gap_Interior_Minus_Linear,
       pch=19,cex=.45,xlab='Easting',ylab='Best interior GCV - exact linear GCV',
       main='Manual Stage 2: structural criterion preference')
  abline(h=0,lty=2)
},width=10.8,height=4.2)

render_figure('F09_detected_minima_histogram',function() {
  am=table(factor(manual_summary$Interior_Dip_Count,levels=0:3))
  ag=table(factor(mgcv_summary$Interior_Grid_Dip_Count,levels=0:3))
  barplot(rbind(am,ag),beside=TRUE,names.arg=0:3,legend.text=c('Manual','mgcv'),
          xlab='Detected interior minima per Easting',ylab='Number of Eastings',
          main='Manual versus mgcv minimum counts')
},width=8.2,height=4.5)

render_figure('F10_regime_preference_map',function() {
  mmanual <- grepl('linear',manual_summary$Lowest_Candidate,ignore.case=TRUE)
  mmgcv <- mgcv_summary$Structural_Gap_Interior_Minus_Linear > 0
  plot(manual_summary$Easting,rep(1,199),pch=15,
       col=ifelse(mmanual,'firebrick','steelblue'),ylim=c(.7,2.3),
       yaxt='n',xlab='Easting',ylab='',main='Observed lowest candidate regime')
  points(mgcv_summary$Easting,rep(2,199),pch=15,
         col=ifelse(mmgcv,'firebrick','steelblue'))
  axis(2,at=c(1,2),labels=c('Manual','mgcv'),las=1)
  legend('topright',c('Linear','Interior'),fill=c('firebrick','steelblue'),bty='n',cex=.85)
},width=11,height=3.5)

bad <- mgcv_summary[mgcv_summary$Auto_Minus_Best_Grid_GCV>1e-4,]
render_figure('F11_sixteen_inferior_automatic_selections',function() {
  par(mfrow=c(2,1),mar=c(5,5,3,1))
  barplot(bad$Auto_Minus_Best_Grid_GCV,names.arg=bad$Easting,
          las=2,ylab='Automatic - grid GCV',main='Inferior automatic selections')
  barplot(bad$Auto_Minus_Best_Grid_Beta,names.arg=bad$Easting,
          las=2,ylab='Automatic - grid Product B',main='Coefficient effects')
  par(mfrow=c(1,1))
},width=11,height=7.5)

m176 <- manual_summary[manual_summary$Easting==176,]
g176 <- mgcv_summary[mgcv_summary$Easting==176,]
render_figure('F12_easting_176_basis_sensitivity',function() {
  par(mfrow=c(1,2),mar=c(7,4,3,1))
  barplot(c(m176$Best_Interior_GCV,g176$Best_Interior_Grid_GCV,m176$GCV_Linear_Limit),
          names.arg=c('Manual interior','mgcv interior','Linear limit'),las=2,
          ylim=c(145.9,146.6),ylab='GCV',main='Easting 176: near-tie')
  barplot(c(m176$Lowest_Candidate_Product_B,g176$Best_Grid_Product_B),
          names.arg=c('Manual selected','mgcv grid selected'),las=2,
          ylim=c(12,15.4),ylab='Product B',main='Different preferred coefficients')
  par(mfrow=c(1,1))
},width=10.5,height=4.5)

render_figure('F13_taxonomy_counts',function() {
  counts <- c('Multiple minima'=sum(taxonomy$Multiple_Minima),
              'Automatic worse than grid'=sum(taxonomy$Auto_Worse_Grid),
              'Automatic worse than linear'=sum(taxonomy$Auto_Worse_Linear),
              'Abs(Product B difference)>=0.25'=
                 sum(abs(taxonomy$Auto_Minus_Best_Grid_Beta)>=0.25))
  par(mar=c(4,11,3,1))
  barplot(rev(counts),horiz=TRUE,las=1,
          xlab='Number of Eastings; categories may overlap',
          main='Nonexclusive taxonomy flags')
},width=10,height=4.4)

# Every new table is saved as a plain CSV, then displayed in the accompanying
# Quarto report. No original Stage 1--3C CSV is modified.
tab <- function(name, d) write.csv(d,file.path(output_dir,name),row.names=FALSE)

stage_summary <- data.frame(
 Stage=c('1A','1B','2','3A','3B','3C'),
 Result=c('4,776 observations; 199 Eastings; 24 passes','72 historical slice records match',
          '79,600 GCV evaluations; 27 multi-minimum Eastings',
          '3,600 fixed-grid pilot GAM fits',
          '199 automatic fits; 6 worse than linear limit',
          '119,400 fixed fits; 16 inferior automatic; 28 multi-minimum Eastings'))
tab('T01_stage_summary.csv',stage_summary)
pilot_cols <- c('Easting','Automatic_GCV','Best_Grid_GCV','Optimizer_GCV_Gap','Product_B_Difference')
tab('T02_six_Easting_pilot.csv',pilot[,pilot_cols])
orig <- data.frame(Easting=c(214,216,218))
orig$Manual_Interior_GCV <- manual_summary$Best_Interior_GCV[match(orig$Easting,manual_summary$Easting)]
orig$Linear_GCV <- manual_summary$GCV_Linear_Limit[match(orig$Easting,manual_summary$Easting)]
orig$Automatic_mgcv_GCV <- mgcv_summary$Automatic_GCV[match(orig$Easting,mgcv_summary$Easting)]
orig$Best_mgcv_Grid_GCV <- mgcv_summary$Best_Grid_GCV[match(orig$Easting,mgcv_summary$Easting)]
tab('T03_original_slice_comparison.csv',orig)
tab('T04_sixteen_inferior_automatic.csv',bad[,c('Easting','Auto_Minus_Best_Grid_GCV','Auto_Minus_Best_Grid_Beta','Interior_Grid_Dip_Count')])
multi <- mgcv_summary[mgcv_summary$Interior_Grid_Dip_Count>=2,]
tab('T05_twenty_eight_multi_minima.csv',multi[,c('Easting','Interior_Grid_Dip_Count','Auto_Worse_Than_Grid','Auto_Minus_Best_Grid_GCV')])
tab('T06_sensitivity_contrasts.csv',mgcv_summary[mgcv_summary$Easting %in% c(174,374),
     c('Easting','Auto_Minus_Best_Grid_GCV','Auto_Minus_Best_Grid_Beta')])
tab('T07_Easting374_minima.csv',rbind(d374[,c('Easting','LogLambda','GCV','Smooth_EDF','Product_B')],
   data.frame(Easting=374,LogLambda=s374$Automatic_LogLambda,
              GCV=s374$Automatic_GCV,Smooth_EDF=s374$Automatic_Smooth_EDF,
              Product_B=s374$Automatic_Product_B)))
tab('T08_minimum_counts.csv',data.frame(Dips=0:3,
      Manual=as.integer(table(factor(manual_summary$Interior_Dip_Count,levels=0:3))),
      mgcv=as.integer(table(factor(mgcv_summary$Interior_Grid_Dip_Count,levels=0:3)))))
tab('T09_Easting176_contrast.csv',data.frame(Model=c('Manual','mgcv'),
 Interior_GCV=c(m176$Best_Interior_GCV,g176$Best_Interior_Grid_GCV),
 Linear_GCV=c(m176$GCV_Linear_Limit,g176$Exact_Linear_GCV),
 Selected_Product_B=c(m176$Lowest_Candidate_Product_B,g176$Best_Grid_Product_B)))
tab('T10_Easting374_grid_landmarks.csv',z374[c(1,
  mgcv_minima$Grid_Index[mgcv_minima$Easting==374],nrow(z374)),
  c('Easting','LogLambda','GCV','Smooth_EDF')])
codebook <- data.frame(Code=c('M0-M3','A','B','C','D','E','F','G','H'),
 Rule=c('Count 0--3 interior minima','At least 2 interior minima',
 'Automatic GCV - best grid GCV > 0.0001',
 'Automatic GCV - exact linear GCV > 0.0001',
 'Best interior GCV - linear GCV > 0.0001',
 'Best interior GCV - linear GCV < -0.0001',
 'Absolute coefficient difference >= 0.25',
 'Multiple minima score spread >= 5',
 'Best grid GCV at right boundary'))
tab('T11_taxonomy_codebook.csv',codebook)
cat('Generated 13 figures in PDF/PNG and 11 source-linked CSV tables in:',
    normalizePath(output_dir), '\n')

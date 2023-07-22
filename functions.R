############fun for fitting panel data expectile model with random effects for diff tau


fit_eeffects<-function(formul,#formula as lm formula
                       dat,#dataframe
                       cluster="no",#how to cluster the panel data
                       taus=seq(0.05,0.95,by=0.05),#expectile levels
                       iters=100,#iterations,
                       ses=T
){
  coefs0<-ses0<-NULL
  mod0<-list(0)
  r2<-NULL
  library(sandwich)
  environment(formul) <- environment()#for the lm function: env of formula
  
  for(i in 1:length(taus)){
    tau<-taus[i]
    w<-rep(0.5,nrow(dat))
    conv<-F
    iter<-0
    while(!conv){
      mod<-lm(formula=formul,data=dat,weights=w)
      wn<-ifelse(as.numeric(resid(mod))>0,tau,1-tau)
      if(all(wn==w)|iter>iters){
        conv<-T
      }
      iter<-iter+1
      w<-wn
      
    }
    
    mod<-lm(formula=formul,data=dat,weights=w)
    mod0[[i]]<-mod
    coefs0<-cbind(coefs0,coef(mod))
    if(ses){
      if(cluster=="no"){
        cov0<-vcov(mod)
      } else if (cluster=="all"){
        cov0<-vcovPL(mod)
      } else{
        cov0<-vcovPL(mod,cluster=cluster)
      }
      ses0<-cbind(ses0, sqrt(diag(cov0)))
    }
    
    r2<-c(r2,summary(mod)$r.squared)
  }
  return(list(coefs=coefs0,ses=ses0, mods=mod0,r2=r2))
}

##########function for plotting the coefficients of the ERRE models with asymptotic $95\%$-CI for varying $\tau$.-->

plot_eefects<-function(mod,#fitted model
                       al=0.05,#alpha (CI error)
                       which_ones,#which coefs to plot
                       taus,#which taus to plot
                       lbls=NA, #labels fro the coefs on the plot
                       special=NULL, #supply user defined limits for the graph
                       sc=rep(1, length(which_ones)),
                       cex.axis=0.5
){
  tau_axis<-taus
  coefs<-mod$coefs[which_ones,which_tau]
  ses<-mod$ses[which_ones,which_tau]
  
  if(is.na(lbls[1])){
    lbls<-rownames(coefs)
  }
  par(mfrow=c(ceiling(length(which_ones)/2),2))
  
  for(j in 1:nrow(coefs)){
    # if(!is.null(special)){
    #   ylims<-c(special[j,1],special[j,2])
    #   ats <- c(round(ylims[1],3),0,round(ylims[2]/2,3),round(ylims[2],3))
    #   labelss<- c(ats[1],0,"",ats[4])
    # } else{
    #   ylims<-c(-0.01,0.01)*sc[j]
    #   ats<-c(-0.01,-0.005,0,0.005,0.01)*sc[j]
    #   labelss<-c(-0.01*sc[j],"",0,"",0.01*sc[j])
    #   
    # }
    # plot(coefs[j,],type="l",main=lbls[j],xlab=bquote(tau), ylab=bquote(beta[tau]),col=4,lwd=2,axes=F)
    # axis(2,at=ats,labels=labelss,...)
    # axis(1,at=seq(2,ncol(coefs),2),labels=seq(0.1,0.9,0.1),...)
    # lines(coefs[j,] - qnorm(1-al/2)*ses[j,],lty=2,col=4)
    # lines(coefs[j,] + qnorm(1-al/2)*ses[j,],lty=2,col=4)
    # abline(h=0,lwd=1,col=2)
    plot(taus,coefs[j,],type="l",main=lbls[j],xlab=bquote(tau), ylab=bquote(beta[tau]),col=4,lwd=2,axes=F,
         ylim=special[c(1,3),j])
    axis(2,at=special[,j])
    axis(1,at=taus, labels=round(taus,2),cex.axis=cex.axis)
    lines(taus,coefs[j,] - qnorm(1-al/2)*ses[j,],lty=2,col=4)
    lines(taus,coefs[j,] + qnorm(1-al/2)*ses[j,],lty=2,col=4)
    abline(h=0,lwd=1,col=2)
    abline(h=coefs[j,which(taus==0.5)],lwd=1,col=3)
  }
}

##################make the coefficient table function-->

make_knitr_table_all<-function(mod, #fitted model
                               which_tau=c(1,2,5,10,15,18,19), #taus to put in the table
                               caption="Coefficients and $R^2$ of the ERRE models with investor opinion proxies and different $\\tau$s.",
                               taus=seq(0.05,0.95,0.05),
                               which_ones=(1:14),
                               row.names=c("Constant","$turn_{i,t-1}$","$pre\\_ean_{i,t}$","$ean_{i,t}$","$post\\_ean_{i,t}$","$news_{i,t}$",
                                           "$\\log(size_{i,t})$","$btm_{i,t-1}$", "$turn_{i,t-1}*days\\_ean_{i,t}$", "$(btm_{i,t-1}\\leq 1)*post\\_ean_{i,t}$",
                                           "$supr_{i,t}$",  "$lm\\_tone_{i,t}$", "$hv\\_tone_{i,t}$", "$ml\\_tone_{i,t}$","$R^2$")
                               
){
  library(knitr)
  #library(kableExtra)
  #https://stackoverflow.com/questions/53341155/coloring-rows-with-kableextra-based-on-cell-values
  #https://haozhu233.github.io/kableExtra/awesome_table_in_pdf.pdf
  
  #coefficient table
  options(scipen=999)
  tab<-round(mod$coefs[,which_tau],5)
  
  
  #pvals and stars
  pvals<-2*(1-pt(abs(mod$coefs[rownames(mod$coefs)%in%rownames(mod$ses),]/mod$ses)[,which_tau], mod$mods[[1]]$df.residual))
  stars<-matrix(rep("",nrow(mod$coefs)),dim(pvals)[1],dim(pvals)[2])
  stars[pvals<0.01]<-"***"; stars[pvals>0.01&pvals<0.05]<-"**";  stars[pvals>0.01&pvals<0.1]<-"*";
  #add stars
  tab[tab>0]<-paste0(" ",tab[tab>0])
  tab<-matrix(paste0(tab,stars,sep=""),dim(pvals)[1],dim(pvals)[2])
  tab<-rbind(tab,round(mod$r2[which_tau],5))
  if(is.null(row.names)){
    row.names =c(names(coef(mod$mods[[1]])),"R^2")
  }
  rownames(tab)<-row.names
  coln<-paste0("$\\hat b_{",round(taus[which_tau],2),"}$")
  
  # knitr::kable(tab, align = "lllllll",
  #              col.names = coln,
  #              row.names = TRUE,
  #              #table.attr = "id=\"table1\"",
  #              digits = 5,
  #              caption = caption
  # )%>%kable_styling() %>%
  #   row_spec(2:6, bold = T, color = "black", background = "#D2F0FA")%>%
  #   row_spec(7:10, bold = T, color = "black", background = "#D2FAE4")%>%
  #   row_spec(11:14, bold = T, color = "black", background = "#F9E8D2")%>%
  #   row_spec(15, bold = T, color = "black", background = "#D9D6DA")
  knitr::kable(tab, align = "lllllll",
               col.names = coln,
               row.names = TRUE,
               #table.attr = "id=\"table1\"",
               digits = 5,
               caption = caption)
}

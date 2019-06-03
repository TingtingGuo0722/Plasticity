JGRA=function(LOC,pheno,envir,enp,tt.line,tt.e,mets)
{
  
  colnames(pheno)=as.character(unlist(pheno[1,]));
  envir.name=colnames(pheno)[-1];
  
  pheno=pheno[-1,];
  m=pheno[,-1];
  m=as.numeric(as.character(unlist(m)));m <- matrix(data=m, ncol=dim(pheno)[2]-1, nrow=dim(pheno)[1]);
  colnames(m)=colnames(pheno)[-1];
  pheno_=data.frame(line_code=pheno$line_code,m);
  colnames(pheno_)=c("line_code",envir.name);
  pheno=pheno_;
  
  envir=envir[envir$env_code%in%envir.name,];
  LOC=LOC[LOC$env_code%in%envir.name,];
  
  order=match(envir$env_code,LOC$env_code);
  LOC=LOC[order,];
  #pheno=pheno[,as.character(envir$env_code)];
  
  rm.env=colnames(pheno)[colSums(is.na(pheno))>=tt.line];
  LOC=LOC[which(!(LOC$env_code%in%rm.env)),];
  envir=envir[which(!(envir$env_code%in%rm.env)),];
  pheno=pheno[,which(!(colnames(pheno)%in%rm.env))];
  
  pheno.1=pheno
  keep=dim(pheno.1)[2]-rowSums(is.na(pheno.1))
  pheno=pheno.1[which(keep>=tt.e),];
  n.line=dim(pheno)[1];
  n.envir=dim(pheno)[2]-1;
  
  library("yarrr")
  coloo=piratepal(palette="basel",trans=0)[1:n.envir];
  
  if(mets=="RM.E") 
  {
    pheno.hat=matrix(999,n.line,dim(envir)[1]);
    cor.whole=numeric();
    for(k in 1:n.envir)
    {
      for(j in 1:n.line)
      {
        x1=envir[,enp][-k];
        y1=as.vector(t(pheno[j,-c(1,1+k)]));
        
        coe=lm(y~x,data=data.frame(x=x1,y=y1));
        y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*envir[,enp][k];
        pheno.hat[j,k]=y.hat;
      }
      cor.each=cor(pheno.hat[,k],pheno[,k+1], use = "complete.obs");
      cor.whole=c(cor.whole,cor.each);
    }
    observe=as.vector(as.matrix(pheno[,-1]));
    predict=as.vector(pheno.hat);
    outforfigure=data.frame(obs=observe,pre=predict,col=rep(coloo,times=rep(n.line,n.envir)))
    r_within=cor.whole;
    names(r_within)=colnames(pheno)[-1];
  }
  
  
  
  
  
  if(mets=="RM.G") 
  {
    intercept=numeric();
    slope=numeric();
    for(j in 1:n.line)
    {
      x1=envir[,enp];
      y1=as.vector(t(pheno[j,-c(1)]));
      
      coe=lm(y~x,data=data.frame(x=x1,y=y1));
      inter=summary(coe)$coefficients[1,1]
      slop=summary(coe)$coefficients[2,1];
      intercept=c(intercept,inter);
      slope=c(slope,slop);
    }
    
    genotype.match=match(pheno[,1],geno[,1])
    genotype=geno[genotype.match,];
    Marker=genotype[,-1];
    
    intercept.hat=numeric();slope.hat=numeric();r_within.1=numeric();
    for(i in 1:5)
    {
      r1 = rbinom(intercept, 1, p=0.5);
      id.T=c(1:n.line)[r1==1]; id.V=c(1:n.line)[r1==0];
      ##Intercept###
      y0=intercept; 
      ans<-mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
      e=as.matrix(ans$u)
      G.pred=Marker[id.V,]
      y_pred=as.matrix(G.pred) %*% e
      GEBV.inter=c(y_pred[,1])+ans$beta;
      
      ##Slope###
      y0=slope; 
      ans<-mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
      e=as.matrix(ans$u)
      G.pred=Marker[id.V,]
      y_pred=as.matrix(G.pred) %*% e
      GEBV.slope=c(y_pred[,1])+ans$beta;
      ###All the predicted slope and intercept
      cor.whole=numeric();yhat.whole=numeric();yobs.whole=numeric();
      for(j in 1:n.envir)
      {
        yhat=GEBV.inter+GEBV.slope*envir[j,enp];
        yobs=pheno[id.V,j+1];
        cor.each=cor(yhat,yobs, use = "complete.obs");
        cor.whole=c(cor.whole,cor.each);
        yhat.whole=c(yhat.whole,yhat);
        yobs.whole=c(yobs.whole,yobs);
      }
      r_within.1=rbind(r_within.1,cor.whole);
    }
    outforfigure=data.frame(obs=yobs.whole,pre=yhat.whole,col=rep(coloo,times=rep(sum(r1==0,na.rm=FALSE),n.envir)))
    r_within=apply(r_within.1,2,mean);
    names(r_within)=colnames(pheno)[-1];
  }   
  
  
  
  if(mets=="RM.GE") 
  { 
    r_within=numeric();yhat.whole=numeric();yobs.whole=numeric();record=numeric();
    for(k in 1:n.envir)
    {
      intercept=numeric();
      slope=numeric();
      for(j in 1:n.line)
      {
        x1=envir[-k,enp];
        y1=as.vector(t(pheno[j,-c(1,1+k)]));
        
        coe=lm(y~x,data=data.frame(x=x1,y=y1));
        inter=summary(coe)$coefficients[1,1]
        slop=summary(coe)$coefficients[2,1];
        intercept=c(intercept,inter);
        slope=c(slope,slop);
      }
      
      genotype.match=match(pheno[,1],geno[,1])
      genotype=geno[genotype.match,];
      Marker=genotype[,-1];
      
      intercept.hat=numeric();slope.hat=numeric();cor_50=numeric();
      for(i in 1:5)
      {
        r1 = rbinom(intercept, 1, p=0.5);
        id.T=c(1:n.line)[r1==1]; id.V=c(1:n.line)[r1==0];
        ##Intercept###
        y0=intercept; 
        ans<-mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
        e=as.matrix(ans$u)
        G.pred=Marker[id.V,]
        y_pred=as.matrix(G.pred) %*% e
        GEBV.inter=c(y_pred[,1])+ans$beta;
        
        ##Slope###
        y0=slope; 
        ans<-mixed.solve(y=y0[id.T],Z=genotype[id.T,-1])
        e=as.matrix(ans$u)
        G.pred=Marker[id.V,]
        y_pred=as.matrix(G.pred) %*% e
        GEBV.slope=c(y_pred[,1])+ans$beta;
        ###All the predicted slope and intercept
        yhat=GEBV.inter+GEBV.slope*envir[k,enp];
        yobs=pheno[id.V,k+1];
        cor.each=cor(yhat,yobs, use = "complete.obs");
        cor_50=c(cor_50,cor.each);
        
      }
      r_within=c(r_within,mean(cor_50));       
      yhat.whole=c(yhat.whole,yhat);
      yobs.whole=c(yobs.whole,yobs);
      record=c(record,length(yhat));
    }    
    
    outforfigure=data.frame(obs=yobs.whole,pre=yhat.whole,
                            col=rep(coloo,times=record))
    names(r_within)=colnames(pheno)[-1];
  }
  
  return(list(outforfigure,r_within));
}

JGRA.marker=function(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets)
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

genotype.match=match(pheno[,1],geno[,1])
genotype=geno[genotype.match,];
Marker=genotype[,-1];
n.marker=dim(Marker)[2];

if(mets=="RM.E") 
{
 
  effect=matrix(999,n.marker,n.envir);
  intercept=matrix(999,1,n.envir)
  
  for(i in 1:n.envir)
  {
    fit=mixed.solve(pheno[,1+i],Z=Marker)
    effect[,i]=fit$u
    intercept[,i]=fit$beta
  }
  
  
  
  pheno.hat=matrix(999,n.line,dim(envir)[1]);
  cor.whole=numeric();
  for(k in 1:n.envir)
  {
    
    effect.hat=numeric();
    for(j in 1:n.marker)
    {
      x1=envir[,enp][-k];
      y1=effect[j,-k];
      
      coe=lm(y~x,data=data.frame(x=x1,y=y1));
      y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*envir[,enp][k];
      effect.hat=c(effect.hat,y.hat);
    }
    
    ##Environment intercept####
    reg.intercept=data.frame(x=as.numeric(envir[,enp][-k]),y=as.vector(intercept[-k]));
    coe.intercept=lm(y~x,data=reg.intercept)
    y.intercept=summary(coe.intercept)$coefficients[1,1]+summary(coe.intercept)$coefficients[2,1]*as.numeric(envir[,enp][k]);
    
    pheno.hat[,k]=y.intercept+as.matrix(Marker)%*%as.matrix(effect.hat);
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
  r_within.1=numeric();
  for(i in 1:50)
  {
    r1 = rbinom(n.line, 1, p=0.5);
    id.T=c(1:n.line)[r1==1]; id.V=c(1:n.line)[r1==0];

    ##marker effect###
    effect=matrix(999,n.marker,n.envir);
    intercept=matrix(999,1,n.envir)
    
    for(k in 1:n.envir)
    {
      fit=mixed.solve(pheno[id.T,1+k],Z=Marker[id.T,])
      effect[,k]=fit$u
      intercept[,k]=fit$beta
    }
    
    ##Slope###
    effect.hat=numeric();
      for(j in 1:n.marker)
      {
        x1=envir[,enp];
        y1=effect[j,];
        coe=lm(y~x,data=data.frame(x=x1,y=y1));
        y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*envir[,enp];
        effect.hat=rbind(effect.hat,y.hat);
      }
      
    ##Environment intercept####
    reg.intercept=data.frame(x=as.numeric(envir[,enp]),y=as.vector(intercept));
    coe.intercept=lm(y~x,data=reg.intercept)
    y.intercept=summary(coe.intercept)$coefficients[1,1]+summary(coe.intercept)$coefficients[2,1]*as.numeric(envir[,enp]);
    
    ###All the predicted slope and intercept
    cor.whole=numeric();yhat.whole=numeric();yobs.whole=numeric();
    for(j in 1:n.envir)
    {
      yhat=y.intercept[j]+(as.matrix(Marker[id.V,])%*%as.matrix(effect.hat))[,j];
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
  r_within.1=numeric();

    for(i in 1:50)
    {
      r1 = rbinom(n.line, 1, p=0.5);
      id.T=c(1:n.line)[r1==1]; id.V=c(1:n.line)[r1==0];
      
      ##marker effect###
      effect=matrix(999,n.marker,n.envir);
      intercept=matrix(999,1,n.envir)
      
      for(k in 1:n.envir)
      {
        fit=mixed.solve(pheno[id.T,1+k],Z=Marker[id.T,])
        effect[,k]=fit$u
        intercept[,k]=fit$beta
      }
      
      ##predict marker effect###
      cor.whole=numeric();yhat.whole=numeric();yobs.whole=numeric();
      for(kk in 1:n.envir)
      {
      effect.hat=numeric();
      for(j in 1:n.marker)
      {
        x1=envir[-kk,enp];
        y1=effect[j,-kk];
        coe=lm(y~x,data=data.frame(x=x1,y=y1));
        y.hat=summary(coe)$coefficients[1,1]+summary(coe)$coefficients[2,1]*envir[kk,enp];
        effect.hat=rbind(effect.hat,y.hat);
      }
      ##Environment intercept####
      reg.intercept=data.frame(x=as.numeric(envir[-kk,enp]),y=as.vector(intercept[-kk]));
      coe.intercept=lm(y~x,data=reg.intercept)
      y.intercept=summary(coe.intercept)$coefficients[1,1]+summary(coe.intercept)$coefficients[2,1]*as.numeric(envir[kk,enp]);
      
      ###All the predicted slope and intercept

        yhat=y.intercept+(as.matrix(Marker[id.V,])%*%as.matrix(effect.hat));
        yobs=pheno[id.V,kk+1];
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

return(list(outforfigure,r_within));
}

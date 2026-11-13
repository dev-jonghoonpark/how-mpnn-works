/* dft.html 핵심 코드: 1차원 장난감 Kohn–Sham 계산기.
   원자 단위(ħ = mₑ = e = 1). 격자 위 유한차분 운동에너지 + 국소 퍼텐셜이라 해밀토니안이 삼중대각 행렬이 된다.
   전자 간 반발은 연화 쿨롱 1/√((x−x')²+1), 교환-상관은 3차원 LDA 교환식을 그대로 빌려 쓴 장난감 근사다. */
(function(root){
'use strict';

/* 삼중대각 대칭 행렬의 고윳값·고유벡터 (암묵적 QL, Numerical Recipes tqli) */
function tqli(d,e,n,z){
  var m,l,iter,i,k,s,r,p,g,f,dd,c,b;
  for(i=1;i<n;i++)e[i-1]=e[i];e[n-1]=0;
  for(l=0;l<n;l++){iter=0;
    do{
      for(m=l;m<n-1;m++){dd=Math.abs(d[m])+Math.abs(d[m+1]);if(Math.abs(e[m])<=1e-15*dd)break;}
      if(m!==l){
        if(iter++>60)break;
        g=(d[l+1]-d[l])/(2*e[l]);r=Math.hypot(g,1);g=d[m]-d[l]+e[l]/(g+(g>=0?Math.abs(r):-Math.abs(r)));
        s=c=1;p=0;
        for(i=m-1;i>=l;i--){
          f=s*e[i];b=c*e[i];e[i+1]=(r=Math.hypot(f,g));
          if(r===0){d[i+1]-=p;e[m]=0;break;}
          s=f/r;c=g/r;g=d[i+1]-p;r=(d[i]-g)*s+2*c*b;d[i+1]=g+(p=s*r);g=c*r-b;
          for(k=0;k<n;k++){f=z[k][i+1];z[k][i+1]=s*z[k][i]+c*f;z[k][i]=c*z[k][i]-s*f;}
        }
        if(r===0&&i>=l)continue;
        d[l]-=p;e[l]=g;e[m]=0;
      }
    }while(m!==l);
  }
}
/* 대각 V와 격자 간격 h로 H = −½∇² + V 의 낮은 고유상태 nOrb개 */
function solve(V,h,nOrb){
  var n=V.length,d=new Float64Array(n),e=new Float64Array(n),z=[],i,j;
  for(i=0;i<n;i++){d[i]=1/(h*h)+V[i];e[i]=i?-0.5/(h*h):0;var row=new Float64Array(n);row[i]=1;z.push(row);}
  tqli(d,e,n,z);
  var idx=[];for(i=0;i<n;i++)idx.push(i);idx.sort(function(a,b){return d[a]-d[b];});
  var out=[];
  for(var k=0;k<nOrb;k++){var c=idx[k],psi=new Float64Array(n),nn=0;
    for(i=0;i<n;i++){psi[i]=z[i][c];nn+=psi[i]*psi[i]*h;}
    nn=Math.sqrt(nn);var sg=0;for(i=0;i<n;i++){psi[i]/=nn;sg+=psi[i];}if(sg<0)for(i=0;i<n;i++)psi[i]=-psi[i];
    out.push({eps:d[c],psi:psi});}
  return out;
}

var CX=Math.pow(3/Math.PI,1/3);
function System(opts){
  this.L=opts.L||10;this.n=opts.n||161;this.h=2*this.L/(this.n-1);
  this.x=new Float64Array(this.n);for(var i=0;i<this.n;i++)this.x[i]=-this.L+i*this.h;
  this.Ne=opts.Ne||2;this.kind=opts.kind||'harmonic';this.inter=opts.inter!==false;this.xc=opts.xc!==false;
  this.vext=new Float64Array(this.n);
  for(i=0;i<this.n;i++){var x=this.x[i];
    if(this.kind==='harmonic')this.vext[i]=0.5*0.25*x*x;
    else if(this.kind==='double')this.vext[i]=0.012*(x*x-16)*(x*x-16)/4-2;
    else /* 분자: 핵 두 개의 연화 쿨롱 인력 */ this.vext[i]=-2/Math.sqrt((x-1.6)*(x-1.6)+1)-2/Math.sqrt((x+1.6)*(x+1.6)+1);}
  /* 연화 쿨롱 핵 */
  this.K=[];for(i=0;i<this.n;i++){var row=new Float64Array(this.n);for(var j=0;j<this.n;j++){var dx=this.x[i]-this.x[j];row[j]=1/Math.sqrt(dx*dx+1);}this.K.push(row);}
  this.dens=new Float64Array(this.n);this.iter=0;this.hist=[];
}
System.prototype.hartree=function(n){var o=new Float64Array(this.n);for(var i=0;i<this.n;i++){var s=0,K=this.K[i];for(var j=0;j<this.n;j++)s+=K[j]*n[j];o[i]=this.inter?s*this.h:0;}return o;};
System.prototype.vxc=function(n){var o=new Float64Array(this.n);if(!this.xc||!this.inter)return o;for(var i=0;i<this.n;i++)o[i]=-CX*Math.cbrt(Math.max(n[i],0));return o;};
System.prototype.exc=function(n){if(!this.xc||!this.inter)return 0;var s=0;for(var i=0;i<this.n;i++)s+=-0.75*CX*Math.pow(Math.max(n[i],0),4/3);return s*this.h;};
/* 입력 밀도로 유효 퍼텐셜을 만들고 → 오비탈을 풀고 → 새 밀도를 만든다 */
System.prototype.step=function(alpha){
  var n=this.dens,vH=this.hartree(n),vX=this.vxc(n),veff=new Float64Array(this.n),i;
  for(i=0;i<this.n;i++)veff[i]=this.vext[i]+vH[i]+vX[i];
  var nOcc=Math.ceil(this.Ne/2),orbs=solve(veff,this.h,nOcc+2),nout=new Float64Array(this.n);
  for(var k=0;k<nOcc;k++){var occ=Math.min(2,this.Ne-2*k);for(i=0;i<this.n;i++)nout[i]+=occ*orbs[k].psi[i]*orbs[k].psi[i];}
  var diff=0;for(i=0;i<this.n;i++)diff+=Math.abs(nout[i]-n[i])*this.h;
  /* 에너지: E = Σ fᵢεᵢ − E_H[n] − ∫v_xc n + E_xc[n] (입력 밀도 기준) */
  var sumEps=0;for(k=0;k<nOcc;k++)sumEps+=Math.min(2,this.Ne-2*k)*orbs[k].eps;
  var EH=0,Vxc=0;for(i=0;i<this.n;i++){EH+=0.5*vH[i]*n[i]*this.h;Vxc+=vX[i]*n[i]*this.h;}
  var E=sumEps-EH-Vxc+this.exc(n);
  var mixed=new Float64Array(this.n);for(i=0;i<this.n;i++)mixed[i]=(1-alpha)*n[i]+alpha*nout[i];
  this.last={veff:veff,vH:vH,vX:vX,orbs:orbs,nOcc:nOcc,nin:n,nout:nout,diff:diff,E:E};
  this.dens=mixed;this.iter++;this.hist.push({diff:diff,E:E});
  return this.last;
};
System.prototype.reset=function(){this.dens=new Float64Array(this.n);this.iter=0;this.hist=[];this.last=null;};

var API={tqli:tqli,solve:solve,System:System};
if(typeof module!=='undefined'&&module.exports)module.exports=API;else root.DFT=API;
})(this);

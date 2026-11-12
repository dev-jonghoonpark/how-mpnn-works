/* how-mpnn-works 핵심 코드: 분자 그래프, 장난감 MPNN, 브라우저 학습기.
   index.html이 <script src="mpnn.js">로 불러 쓰고, node에서도 require로 검증할 수 있다. */
(function(root){
'use strict';

function rng(seed){return function(){seed|=0;seed=seed+0x6D2B79F5|0;var t=Math.imul(seed^seed>>>15,1|seed);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
function gauss(r){var u=1-r(),v=r();return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);}

/* ── 분자: 무거운 원자만 노드로 두고 수소는 개수 특징으로 (논문 기본 설정, 암시적 수소) ──
   x, y는 그림용 2차원 좌표(Å 단위 근사). b: [i, j, 결합 종류] — 1 단일, 2 이중, 3 삼중, 4 방향족 */
var MOLS={
  ethanol:{name:'에탄올',formula:'C₂H₅OH',
    atoms:[{el:'C',h:3,hyb:'sp3',x:0,y:0},{el:'C',h:2,hyb:'sp3',x:1.26,y:0.73},{el:'O',h:1,hyb:'sp3',x:2.5,y:0}],
    b:[[0,1,1],[1,2,1]]},
  acetic:{name:'아세트산',formula:'CH₃COOH',
    atoms:[{el:'C',h:3,hyb:'sp3',x:0,y:0},{el:'C',h:0,hyb:'sp2',x:1.3,y:0.75},{el:'O',h:0,hyb:'sp2',x:1.3,y:2.0},{el:'O',h:1,hyb:'sp2',x:2.55,y:0}],
    b:[[0,1,1],[1,2,2],[1,3,1]]},
  acetonitrile:{name:'아세토나이트릴',formula:'CH₃CN',
    atoms:[{el:'C',h:3,hyb:'sp3',x:0,y:0},{el:'C',h:0,hyb:'sp',x:1.46,y:0},{el:'N',h:0,hyb:'sp',x:2.62,y:0}],
    b:[[0,1,1],[1,2,3]]},
  glycine:{name:'글리신',formula:'NH₂CH₂COOH',
    atoms:[{el:'N',h:2,hyb:'sp3',x:0,y:0},{el:'C',h:2,hyb:'sp3',x:1.27,y:0.73},{el:'C',h:0,hyb:'sp2',x:2.55,y:0},{el:'O',h:0,hyb:'sp2',x:2.55,y:-1.25},{el:'O',h:1,hyb:'sp2',x:3.8,y:0.73}],
    b:[[0,1,1],[1,2,1],[2,3,2],[2,4,1]]},
  pyridine:{name:'피리딘',formula:'C₅H₅N',atoms:[],b:[]},
  benzene:{name:'벤젠',formula:'C₆H₆',atoms:[],b:[]},
  octanol:{name:'1-옥탄올',formula:'C₈H₁₇OH',atoms:[],b:[]}
};
(function(){
  /* 고리 두 개는 정육각형, 옥탄올은 지그재그로 좌표를 만든다 */
  function ring(m,els){for(var i=0;i<6;i++){var a=Math.PI/2+i*Math.PI/3;
    m.atoms.push({el:els[i],h:els[i]==='C'?1:0,hyb:'sp2',arom:true,x:1.4*Math.cos(a),y:1.4*Math.sin(a)});
    m.b.push([i,(i+1)%6,4]);}}
  ring(MOLS.benzene,['C','C','C','C','C','C']);
  ring(MOLS.pyridine,['N','C','C','C','C','C']);
  var o=MOLS.octanol;
  for(var i=0;i<9;i++){o.atoms.push({el:i===8?'O':'C',h:i===0?3:(i===8?1:2),hyb:'sp3',x:i*1.26,y:i%2?0.73:0});if(i)o.b.push([i-1,i,1]);}
})();

var EL={H:1,C:6,N:7,O:8,F:9};
var BOND_NAME={1:'단일',2:'이중',3:'삼중',4:'방향족'};
var BOND_LEN={'C-C1':1.54,'C-C4':1.40,'C-N4':1.34,'C-O1':1.43,'C-O2':1.21,'C-N1':1.47,'C-N3':1.16};

/* 논문 표 1의 원자 특징. 받개·주개는 RDKit 정의를 단순화: N·O는 받개, H가 붙은 N·O는 주개 */
function atomFeatures(a){
  var types=['H','C','N','O','F'];
  var hy=['sp','sp2','sp3'];
  return {
    type:types.map(function(t){return a.el===t?1:0;}),
    Z:EL[a.el],
    acceptor:(a.el==='N'||a.el==='O')?1:0,
    donor:((a.el==='N'||a.el==='O')&&a.h>0)?1:0,
    aromatic:a.arom?1:0,
    hyb:hy.map(function(t){return a.hyb===t?1:0;}),
    nH:a.h
  };
}
function featVec(a){var f=atomFeatures(a);return f.type.concat([f.Z,f.acceptor,f.donor,f.aromatic],f.hyb,[f.nH]);}

function adjacency(m){
  var n=m.atoms.length,adj=[];for(var i=0;i<n;i++)adj.push([]);
  m.b.forEach(function(e){adj[e[0]].push({w:e[1],t:e[2]});adj[e[1]].push({w:e[0],t:e[2]});});
  return adj;
}
/* 모든 쌍의 그래프 거리 (BFS) */
function graphDist(n,adj){
  var D=[];for(var s=0;s<n;s++){var d=new Array(n).fill(Infinity);d[s]=0;var q=[s];
    while(q.length){var u=q.shift();adj[u].forEach(function(e){if(d[e.w]===Infinity){d[e.w]=d[u]+1;q.push(e.w);}});}
    D.push(d);}
  return D;
}
function eucl(m,i,j){var a=m.atoms[i],b=m.atoms[j],dx=a.x-b.x,dy=a.y-b.y;return Math.sqrt(dx*dx+dy*dy);}

/* 논문 6절 거리 구간: [0,2], [2,6]을 8등분, [6,∞) — 결합 종류 4개와 합쳐 알파벳 14개 */
function distBin(d){if(d<2)return 0;if(d>=6)return 9;return 1+Math.min(7,Math.floor((d-2)/0.5));}
function binLabel(k){if(k===0)return '[0, 2)';if(k===9)return '[6, ∞)';var lo=2+(k-1)*0.5;return '['+lo.toFixed(1)+', '+(lo+0.5).toFixed(1)+')';}

/* ── 장난감 MPNN (데모용, 가중치 고정) ──
   h⁰ = 특징을 d차원으로 사영, M = A_{결합종류} h_w (행렬 곱 메시지), U = tanh(W h + m + b) (논문은 GRU) */
function ToyMPNN(d,seed){
  var r=rng(seed||7),i;
  this.d=d;
  this.P=[];for(i=0;i<d*13;i++)this.P.push(gauss(r)*0.35);
  this.A={};[1,2,3,4].forEach(function(t){var M=[];for(var k=0;k<d*d;k++)M.push(gauss(r)*0.45/Math.sqrt(d));this.A[t]=M;},this);
  this.W=[];for(i=0;i<d*d;i++)this.W.push(gauss(r)*0.5/Math.sqrt(d));
  this.bias=[];for(i=0;i<d;i++)this.bias.push(gauss(r)*0.1);
}
function matvec(M,v,d){var o=new Array(d).fill(0);for(var i=0;i<d;i++){var s=0;for(var j=0;j<v.length;j++)s+=M[i*v.length+j]*v[j];o[i]=s;}return o;}
ToyMPNN.prototype.h0=function(m){var d=this.d,P=this.P;return m.atoms.map(function(a){
  var x=featVec(a).slice();x[5]/=9;x[12]/=3;return matvec(P,x,d).map(Math.tanh);});};
/* T 단계를 돌리며 단계마다 h와 메시지를 기록한다 */
ToyMPNN.prototype.run=function(m,T){
  var d=this.d,adj=adjacency(m),H=[this.h0(m)],MS=[],self=this;
  for(var t=0;t<T;t++){
    var h=H[t],msgs=[],nh=[];
    for(var v=0;v<h.length;v++){
      var parts=adj[v].map(function(e){return {w:e.w,t:e.t,msg:matvec(self.A[e.t],h[e.w],d)};});
      var mv=new Array(d).fill(0);parts.forEach(function(p){for(var k=0;k<d;k++)mv[k]+=p.msg[k];});
      var wh=matvec(self.W,h[v],d),out=[];
      for(var k=0;k<d;k++)out.push(Math.tanh(wh[k]+mv[k]+self.bias[k]));
      msgs.push({parts:parts,m:mv});nh.push(out);
    }
    MS.push(msgs);H.push(nh);
  }
  return {H:H,M:MS,adj:adj};
};

/* ── 브라우저 학습기 ──
   과제: 무작위 분자 모양 나무(원자 6~9개, C/N/O)에서 "그래프 거리가 정확히 3인 헤테로원자(N·O) 쌍의 수"를 맞히기.
   답을 알려면 원자마다 3칸 떨어진 곳까지 봐야 하므로, 메시지 전달 횟수 T가 3보다 작으면 원리적으로 못 맞힌다.
   모델: h⁰ = E[원자 종류], h^{t+1}_v = tanh(U h_v + A Σ_{w∈N(v)} h_w + b) (가중치는 t에 공유),
         ŷ = Σ_v wᵀ tanh(R h^T_v + c) + b₀ */
function randTree(r){
  var n=6+Math.floor(r()*4),types=[],edges=[],deg=[];
  for(var i=0;i<n;i++){var u=r();types.push(u<0.6?0:(u<0.8?1:2));deg.push(0);}
  for(i=1;i<n;i++){var p,tries=0;do{p=Math.floor(r()*i);tries++;}while(deg[p]>=3&&tries<20);edges.push([p,i]);deg[p]++;deg[i]++;}
  var adj=[];for(i=0;i<n;i++)adj.push([]);edges.forEach(function(e){adj[e[0]].push(e[1]);adj[e[1]].push(e[0]);});
  return {n:n,types:types,edges:edges,adj:adj};
}
function treeTarget(g,K){
  var adjE=g.adj.map(function(l){return l.map(function(w){return {w:w};});});
  var D=graphDist(g.n,adjE),c=0;
  for(var i=0;i<g.n;i++)for(var j=i+1;j<g.n;j++)if(g.types[i]&&g.types[j]&&D[i][j]===K)c++;
  return c;
}
function makeData(seed,N,K){var r=rng(seed),out=[];for(var i=0;i<N;i++){var g=randTree(r);g.y=treeTarget(g,K);out.push(g);}return out;}

function Param(n,s,r){this.v=new Float64Array(n);this.g=new Float64Array(n);this.m=new Float64Array(n);this.s=new Float64Array(n);
  for(var i=0;i<n;i++)this.v[i]=r?gauss(r)*s:0;}
function Net(T,d,seed){
  var r=rng(seed);this.T=T;this.d=d;
  this.E=new Param(3*d,0.5,r);
  this.U=new Param(d*d,0.6/Math.sqrt(d),r);
  this.A=new Param(d*d,0.6/Math.sqrt(d),r);
  this.b=new Param(d,0,null);
  this.R=new Param(d*d,1/Math.sqrt(d),r);
  this.c=new Param(d,0,null);
  this.w=new Param(d,0.3,r);
  this.b0=new Param(1,0,null);
  this.ps=[this.E,this.U,this.A,this.b,this.R,this.c,this.w,this.b0];
  this.step=0;
}
Net.prototype.forward=function(g){
  var d=this.d,n=g.n,T=this.T,H=[],v,k,j;
  var h0=[];for(v=0;v<n;v++){var e=new Float64Array(d);for(k=0;k<d;k++)e[k]=this.E.v[g.types[v]*d+k];h0.push(e);}
  H.push(h0);
  var S=[];
  for(var t=0;t<T;t++){
    var h=H[t],nh=[],ss=[];
    for(v=0;v<n;v++){
      var s=new Float64Array(d);g.adj[v].forEach(function(w){for(var q=0;q<d;q++)s[q]+=h[w][q];});ss.push(s);
      var o=new Float64Array(d);
      for(k=0;k<d;k++){var a=this.b.v[k];for(j=0;j<d;j++)a+=this.U.v[k*d+j]*h[v][j]+this.A.v[k*d+j]*s[j];o[k]=Math.tanh(a);}
      nh.push(o);
    }
    S.push(ss);H.push(nh);
  }
  var Z=[],y=this.b0.v[0],hT=H[T];
  for(v=0;v<n;v++){var z=new Float64Array(d);for(k=0;k<d;k++){var a2=this.c.v[k];for(j=0;j<d;j++)a2+=this.R.v[k*d+j]*hT[v][j];z[k]=Math.tanh(a2);y+=this.w.v[k]*z[k];}Z.push(z);}
  return {H:H,S:S,Z:Z,y:y};
};
/* 손실 ½(ŷ−y)²의 기울기를 누적한다 */
Net.prototype.backward=function(g,f){
  var d=this.d,n=g.n,T=this.T,v,k,j,dy=f.y-g.y;
  this.b0.g[0]+=dy;
  var dH=[];for(v=0;v<n;v++)dH.push(new Float64Array(d));
  var hT=f.H[T];
  for(v=0;v<n;v++)for(k=0;k<d;k++){
    this.w.g[k]+=dy*f.Z[v][k];
    var da=dy*this.w.v[k]*(1-f.Z[v][k]*f.Z[v][k]);
    this.c.g[k]+=da;
    for(j=0;j<d;j++){this.R.g[k*d+j]+=da*hT[v][j];dH[v][j]+=da*this.R.v[k*d+j];}
  }
  for(var t=T-1;t>=0;t--){
    var h=f.H[t],o=f.H[t+1],ndH=[],dS=[];
    for(v=0;v<n;v++){ndH.push(new Float64Array(d));dS.push(new Float64Array(d));}
    for(v=0;v<n;v++)for(k=0;k<d;k++){
      var dl=dH[v][k]*(1-o[v][k]*o[v][k]);if(dl===0)continue;
      this.b.g[k]+=dl;
      for(j=0;j<d;j++){
        this.U.g[k*d+j]+=dl*h[v][j];this.A.g[k*d+j]+=dl*f.S[t][v][j];
        ndH[v][j]+=dl*this.U.v[k*d+j];dS[v][j]+=dl*this.A.v[k*d+j];
      }
    }
    /* s_v = Σ_{w∈N(v)} h_w 이므로 dS_v를 이웃들에게 돌려준다 */
    for(v=0;v<n;v++)g.adj[v].forEach(function(w){for(var q=0;q<d;q++)ndH[w][q]+=dS[v][q];});
    dH=ndH;
  }
  for(v=0;v<n;v++)for(k=0;k<d;k++)this.E.g[g.types[v]*d+k]+=dH[v][k];
  return 0.5*dy*dy;
};
Net.prototype.adam=function(lr,bs){
  this.step++;var b1=0.9,b2=0.999,c1=1-Math.pow(b1,this.step),c2=1-Math.pow(b2,this.step);
  this.ps.forEach(function(p){for(var i=0;i<p.v.length;i++){var gr=p.g[i]/bs;
    p.m[i]=b1*p.m[i]+(1-b1)*gr;p.s[i]=b2*p.s[i]+(1-b2)*gr*gr;
    p.v[i]-=lr*(p.m[i]/c1)/(Math.sqrt(p.s[i]/c2)+1e-8);p.g[i]=0;}});
};
Net.prototype.trainBatch=function(data,r,bs,lr){
  var L=0;for(var i=0;i<bs;i++){var g=data[Math.floor(r()*data.length)];L+=this.backward(g,this.forward(g));}
  this.adam(lr,bs);return L/bs;
};
Net.prototype.evalMAE=function(data){var s=0;for(var i=0;i<data.length;i++)s+=Math.abs(this.forward(data[i]).y-data[i].y);return s/data.length;};

var API={rng:rng,gauss:gauss,MOLS:MOLS,EL:EL,BOND_NAME:BOND_NAME,atomFeatures:atomFeatures,featVec:featVec,adjacency:adjacency,
  graphDist:graphDist,eucl:eucl,distBin:distBin,binLabel:binLabel,ToyMPNN:ToyMPNN,matvec:matvec,
  randTree:randTree,treeTarget:treeTarget,makeData:makeData,Net:Net};
if(typeof module!=='undefined'&&module.exports)module.exports=API;else root.MP=API;
})(this);

#!/usr/bin/env python3
"""Eight bilingual figures for the SR manuscript, drawn as editable SVG.

Theoretical panels are deterministic illustrations. Figure 7 reads archived
empirical summary tables. Inkscape exports vector PDFs (outlined text for robust
CJK rendering) and PNG previews. SVG text is editable and no raster is embedded.
Run from any directory: python code/redraw_figures.py
"""
from pathlib import Path
from xml.sax.saxutils import escape
import argparse, csv, json, math, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'source_data'
for folder in ('svg','pdf','png'):(ROOT/'figures'/folder).mkdir(parents=True,exist_ok=True)
W=680
C={'ink':'#20272D','muted':'#65717B','line':'#B8C1C7','grid':'#E6EBEE',
   'F':'#0072B2','B':'#D55E00','teal':'#009E73','purple':'#8B6BA8',
   'grey':'#87929A','lightF':'#E5F2F8','lightB':'#FBECE3',
   'lightP':'#F0EAF5','pale':'#F4F6F7','white':'#FFFFFF'}
ANALYTIC={}
REGISTRY=[]

def val(col):return C.get(col,col)
def blend(col,t):
    h=val(col); a=np.array([int(h[i:i+2],16) for i in (1,3,5)])
    b=np.rint(255*(1-t)+a*t).astype(int)
    return '#'+''.join(f'{v:02X}' for v in b)

class SVG:
    def __init__(self,number,stem,lang,height,title):
        self.number,self.stem,self.lang,self.h=number,stem,lang,height
        self.font='Noto Sans CJK SC' if lang=='zh' else 'DejaVu Sans'
        self.p=[f'<svg xmlns="http://www.w3.org/2000/svg" width="170mm" height="{height/4}mm" viewBox="0 0 680 {height}" role="img" aria-labelledby="title desc">',
            f'<title id="title">{escape(title)}</title>',
            f'<desc id="desc">Figure {number}. Original editable vector figure; panel descriptions and assumptions are in the accompanying manuscript caption.</desc>',
            '<defs>']
        for k in ['ink','muted','F','B','teal','purple','grey']:
            self.p.append(f'<marker id="arr-{k}" markerWidth="7" markerHeight="7" refX="6.3" refY="3.5" orient="auto" markerUnits="userSpaceOnUse"><path d="M0 0 L7 3.5 L0 7 Z" fill="{C[k]}"/></marker>')
        self.p+=['</defs>','<rect width="680" height="100%" fill="#FFFFFF"/>']
    def t(self,x,y,txt,size=12,color='ink',weight=400,anchor='start',rotate=None):
        trans=f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ''
        self.p.append(f'<text id="txt-{len(self.p)}" x="{x:.3f}" y="{y:.3f}" font-family="{self.font}, DejaVu Sans, sans-serif" font-size="{size}" font-weight="{weight}" fill="{val(color)}" text-anchor="{anchor}"{trans}>{escape(str(txt))}</text>')
    def lines(self,x,y,lines,size=12,color='ink',weight=400,anchor='start',lead=17):
        for i,z in enumerate(lines):self.t(x,y+i*lead,z,size,color,weight,anchor)
    def label(self,label,x,y,title):
        self.t(x,y,label,17,weight=700)
        self.t(x+22,y,title,14,weight=550)
    def line(self,x1,y1,x2,y2,color='muted',sw=1,dash=None,arrow=False,opacity=1):
        d=f' stroke-dasharray="{dash}"' if dash else ''
        a=f' marker-end="url(#arr-{color})"' if arrow else ''
        self.p.append(f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" stroke="{val(color)}" stroke-width="{sw}" opacity="{opacity}"{d}{a}/>')
    def rect(self,x,y,w,h,fill='none',stroke=None,sw=1,rx=0,opacity=1,dash=None):
        st=f' stroke="{val(stroke)}" stroke-width="{sw}"' if stroke else ''
        ds=f' stroke-dasharray="{dash}"' if dash else ''
        self.p.append(f'<rect x="{x:.3f}" y="{y:.3f}" width="{w:.3f}" height="{h:.3f}" rx="{rx}" fill="{val(fill)}" opacity="{opacity}"{st}{ds}/>')
    def circle(self,x,y,r,fill='white',stroke=None,sw=1,opacity=1):
        st=f' stroke="{val(stroke)}" stroke-width="{sw}"' if stroke else ''
        self.p.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{r}" fill="{val(fill)}" opacity="{opacity}"{st}/>')
    def path(self,pts,color='F',sw=1.8,fill='none',opacity=1,dash=None,close=False):
        d='M '+' L '.join(f'{x:.3f} {y:.3f}' for x,y in pts)+(' Z' if close else '')
        ds=f' stroke-dasharray="{dash}"' if dash else ''
        self.p.append(f'<path d="{d}" stroke="{val(color)}" stroke-width="{sw}" fill="{val(fill)}" opacity="{opacity}" stroke-linejoin="round" stroke-linecap="round"{ds}/>')
    def diamond(self,x,y,r,color='ink'):
        self.path([(x,y-r),(x+r,y),(x,y+r),(x-r,y)],color,1,color,close=True)
    def heat(self,M,x,y,w,h,high='F',signed=False,vmax=None,annot=False,annotation_format='.2g'):
        M=np.asarray(M);n,k=M.shape
        mx=float(np.max(np.abs(M))) if vmax is None else vmax
        mx=mx or 1
        for i in range(n):
            for j in range(k):
                z=float(M[i,j]);color='B' if signed and z<0 else high
                self.rect(x+j*w/k,y+i*h/n,w/k+.03,h/n+.03,blend(color,min(abs(z)/mx,1)))
                if annot:self.t(x+(j+.5)*w/k,y+(i+.5)*h/n+4,format(z,annotation_format),11.5,'white' if abs(z)/mx>.65 else 'ink',anchor='middle')
        self.rect(x,y,w,h,'none','line',.7)
    def colorbar(self,x,y,w,high='F',lo='0',hi='1'):
        for j in range(60):self.rect(x+j*w/60,y,w/60+.1,5,blend(high,j/59))
        self.t(x,y+19,lo,11.5,'muted');self.t(x+w,y+19,hi,11.5,'muted',anchor='end')
    def save(self):
        name=f'{self.stem}_{self.lang}'
        p=ROOT/'figures/svg'/f'{name}.svg'
        p.write_text('\n'.join(self.p+['</svg>']),encoding='utf-8')
        REGISTRY.append({'number':self.number,'language':self.lang,'stem':name,'width_mm':170,'height_mm':self.h/4,'min_font_pt':11.5*170/680*72/25.4,'svg':f'figures/svg/{name}.svg','pdf':f'figures/pdf/{name}.pdf','png':f'figures/png/{name}.png'})
        return p

def text(s,en,zh):return zh if s.lang=='zh' else en
def ring(s,cx,cy,r,n=8,step=1,color='F',labels=True):
    pts=[(cx+r*math.cos(-math.pi/2+2*math.pi*j/n),cy+r*math.sin(-math.pi/2+2*math.pi*j/n)) for j in range(n)]
    for j in range(n):
        a=np.array(pts[j]);b=np.array(pts[(j+step)%n]);u=(b-a)/np.linalg.norm(b-a)
        s.line(*(a+11*u),*(b-12*u),color,.9,arrow=True)
    for j,(x,y) in enumerate(pts):
        s.circle(x,y,9.4,'white',color,.9)
        if labels:s.t(x,y+4,str(j),11.5,color,anchor='middle')

def axes(s,x,y,w,h,xlim,ylim,xticks,yticks,xlabel=None,ylabel=None):
    X=lambda a:x+(a-xlim[0])/(xlim[1]-xlim[0])*w
    Y=lambda a:y+h-(a-ylim[0])/(ylim[1]-ylim[0])*h
    for a in yticks:
        s.line(x,Y(a),x+w,Y(a),'grid',.6)
        s.t(x-8,Y(a)+4,f'{a:g}',11.5,'muted',anchor='end')
    for a in xticks:
        s.line(X(a),y+h,X(a),y+h+4,'muted',.8)
        s.t(X(a),y+h+19,f'{a:g}',11.5,'muted',anchor='middle')
    s.line(x,y,x,y+h,'muted',.9);s.line(x,y+h,x+w,y+h,'muted',.9)
    if xlabel:s.t(x+w/2,y+h+38,xlabel,12,'muted',anchor='middle')
    if ylabel:s.t(x-40,y+h/2,ylabel,12,'muted',anchor='middle',rotate=-90)
    return X,Y

def rotation(s,cx,cy,r,ang,col='F',label=None):
    s.circle(cx,cy,r,'none','line',.8)
    s.line(cx,cy,cx+r,cy,'grey',1,arrow=True)
    a=math.radians(ang)
    s.line(cx,cy,cx+r*math.cos(a),cy-r*math.sin(a),col,1.7,arrow=True)
    aa=np.linspace(0,a,30)
    s.path([(cx+r*.43*math.cos(t),cy-r*.43*math.sin(t)) for t in aa],col,1.2)
    if label:s.t(cx,cy+r+22,label,12,col,anchor='middle')

def fig01(lang):
    s=SVG(1,'fig01_foundations',lang,585,'Predictive occupancy and two exact observational equivalences')
    n=8;P=np.roll(np.eye(n),1,axis=1)
    MF=np.linalg.inv(np.eye(n)-.8*P);MB=MF.T
    DF=np.linalg.norm(MF[:,None,:]-MF[None,:,:],axis=-1)
    DB=np.linalg.norm(MB[:,None,:]-MB[None,:,:],axis=-1)
    assert np.allclose(DF,DB)
    s.label('a',18,24,text(s,'From transitions to occupancy','从转移到占据'))
    ring(s,87,101,47)
    s.t(87,171,'P',13,anchor='middle')
    s.line(155,103,188,103,'muted',1,arrow=True)
    terms=np.array([(.7**k)*np.linalg.matrix_power(P,k)[0] for k in range(4)])
    s.heat(terms,233,52,192,96,'F',vmax=1)
    for k,l in enumerate(['I','γP','γ²P²','γ³P³']):s.t(222,68+24*k,l,12,anchor='end')
    for j in range(8):s.t(245+j*24,166,str(j),11.5,'muted',anchor='middle')
    s.t(329,185,text(s,'First four terms','前四项'),12,'muted',anchor='middle')
    s.t(444,107,'…',19,'muted',anchor='middle');s.line(460,103,483,103,'muted',1,arrow=True)
    M=np.linalg.inv(np.eye(n)-.7*P)
    for j,z in enumerate(M[0]):
        s.rect(510+j*18,145-72*z/M[0].max(),11,72*z/M[0].max(),'F')
        s.t(515.5+j*18,165,str(j),11.5,'muted',anchor='middle')
    s.line(502,145,651,145,'muted',.8)
    s.t(580,53,'Mγ(P)',13,anchor='middle')
    s.t(580,185,'γ = 0.7',12,'muted',anchor='middle')

    s.label('b',18,225,text(s,'Different patterns, identical distances','模式不同，距离相同'))
    s.t(106,251,text(s,'Forward','前瞻'),12,'F',anchor='middle')
    s.t(251,251,text(s,'Backward','后瞻'),12,'B',anchor='middle')
    s.heat(MF,58,263,96,96,'F');s.heat(MB,203,263,96,96,'B')
    s.t(178,317,'≠',23,'ink',anchor='middle')
    s.line(325,289,420,308,'F',1,arrow=True)
    s.line(325,336,420,317,'B',1,arrow=True)
    s.heat(DF,451,263,96,96,'grey')
    s.t(499,251,text(s,'One RDM','同一 RDM'),12,anchor='middle')
    s.t(106,383,'M₀.₈(P)',12,'F',anchor='middle')
    s.t(251,383,'M₀.₈(Pᵀ)',12,'B',anchor='middle')
    s.t(499,383,'DF = DB',13,anchor='middle')
    s.lines(577,298,['PPᵀ','=', 'PᵀP'],13,'muted',anchor='middle',lead=21)

    s.label('c',18,425,text(s,'A learned reference changes the operator','已学习参考改变相对算子'))
    P7=np.zeros((7,7))
    for i in range(7):
        for d in [1,2,4]:P7[i,(i+d)%7]=1/3
    Q=np.zeros((7,6))
    for j in range(6):
        Q[:j+1,j]=1/math.sqrt((j+1)*(j+2));Q[j+1,j]=-(j+1)/math.sqrt((j+1)*(j+2))
    A=Q.T@P7@Q
    T=np.linalg.inv(np.eye(6)-.3*A)@(np.eye(6)-.8*A)
    a=634/519;beta=75/173;Tb=a*np.linalg.inv(np.eye(6)-beta*A.T)
    assert np.allclose(T,Tb)
    s.heat(P7,51,455,77,77,'grey')
    s.t(89.5,554,text(s,'7 states','7 个状态'),12,'muted',anchor='middle')
    s.line(149,493,214,493,'muted',1,arrow=True)
    s.heat(T,248,455,90,90,'F',True);s.heat(Tb,491,455,90,90,'F',True)
    s.t(414,508,'=',25,anchor='middle')
    s.t(293,566,'M₀.₃(A) M₀.₈(A)⁻¹',12,anchor='middle')
    s.t(536,566,'aMβ(Aᵀ)',12,anchor='middle')
    ANALYTIC['fig01']={'gamma_basics':.7,'weighted_first_four_rows':terms.tolist(),'complete_state_zero_row':M[0].tolist(),'gamma_geometry':.8,'max_distance_difference':float(np.max(np.abs(DF-DB))),'paley_gain':a,'paley_backward_discount':beta,'paley_error':float(np.max(np.abs(T-Tb)))}
    return s.save()

def fig02(lang):
    s=SVG(2,'fig02_families',lang,440,'Continuous prediction families retain properties while parameters vary')
    s.label('a',18,24,text(s,'Continuous prediction families','连续预测族'))
    x,y,w,h=58,60,262,303
    X=lambda a:x+(a+1.5)/3*w;Y=lambda a:y+(1.5-a)/3*h
    s.rect(x,y,w,h/2,'lightF');s.rect(x,y+h/2,w,h/2,'lightB')
    s.line(x,Y(0),x+w+9,Y(0),'muted',1,arrow=True)
    s.line(X(0),y+h,X(0),y-8,'muted',1,arrow=True)
    s.t(x+w+13,Y(0)+4,'u',12);s.t(X(0)-8,y-10,'v',12)
    for k in [-1,1]:
        s.line(X(k),Y(0)-3,X(k),Y(0)+3,'muted',.8);s.t(X(k),Y(0)+18,str(k),11.5,'muted',anchor='middle')
        s.line(X(0)-3,Y(k),X(0)+3,Y(k),'muted',.8);s.t(X(0)-8,Y(k)+4,str(k),11.5,'muted',anchor='end')
    s.t(x+12,y+22,text(s,'Forward family','前向族'),13,'F',550)
    s.t(x+12,y+h-14,text(s,'Backward family','后向族'),13,'B',550)
    g=np.linspace(0,.999,200)
    uu=1/(1+g*g);vv=g/(1+g*g)
    s.path([(X(a),Y(b)) for a,b in zip(uu,vv)],'F',2.2)
    s.path([(X(a),Y(-b)) for a,b in zip(uu,vv)],'B',2.2)
    s.circle(X(.5),Y(.5),3,'white','F',1)
    s.circle(X(.5),Y(-.5),3,'white','B',1)
    s.circle(X(1),Y(1),4.4,'F');s.t(X(1)-6,Y(1)-12,'I + J',12,'F',anchor='end')
    s.circle(X(1),Y(0),3.3,'white','ink',1)
    s.t(189,396,'G = uI + vJ',14,anchor='middle')
    s.t(189,420,text(s,'Direct reference','直接参考'),12,'muted',anchor='middle')

    s.label('b',359,24,text(s,'Exact compensation','精确补偿'))
    Xb,Yb=axes(s,400,62,241,122,(0,1),(1,2),[0,.5,1],[1,1.5,2],xlabel='γ')
    gs=np.linspace(0,1,201)
    s.path([(Xb(g),Yb(1+g)) for g in gs],'F',1.9)
    s.path([(Xb(g),Yb(1+g*g)) for g in gs],'purple',1.9)
    s.circle(Xb(1),Yb(2),3.1,'white','F',1.2)
    s.t(419,80,'c = 1 + γ',11.5,'F')
    s.t(542,177,'α = 1 + γ²',11.5,'purple')
    s.t(520,248,'Gγ = I + J',14,'F',550,anchor='middle')

    s.label('c',359,273,text(s,'One observable map','同一可见地图'))
    J=np.array([[0.,-1],[1,0]])
    mix=.5*np.linalg.inv(np.eye(2)-.2*J)+.5*np.linalg.inv(np.eye(2)-.85*J)
    ge=mix[1,0]/mix[0,0];ae=mix[0,0]*(1+ge*ge)
    eff=ae*np.linalg.inv(np.eye(2)-ge*J)
    assert np.allclose(mix,eff)
    s.t(429,302,text(s,'Two scales','两个时域'),12,'purple',anchor='middle')
    s.t(590,302,text(s,'Effective SR','有效 SR'),12,'F',anchor='middle')
    s.heat(mix,394,315,70,70,'F',True,vmax=mix.max())
    s.heat(eff,555,315,70,70,'F',True,vmax=mix.max())
    s.t(510,359,'=',24,anchor='middle')
    s.t(429,409,'γ = 0.20, 0.85',11.5,'muted',anchor='middle')
    s.t(590,409,f'γeff = {ge:.3f}',11.5,'muted',anchor='middle')
    ANALYTIC['fig02']={'mixture_discounts':[.2,.85],'weights':[.5,.5],'mixture_matrix':mix.tolist(),'effective_matrix':eff.tolist(),'effective_discount':float(ge),'effective_gain':float(ae),'max_effective_error':float(np.max(np.abs(mix-eff)))}
    return s.save()

def fig03(lang):
    s=SVG(3,'fig03_reciprocity',lang,450,'Reciprocal task structure and the complete four-family angle ranges')
    s.label('a',18,24,text(s,'A shared reciprocal plane','共同互反平面'))
    ring(s,85,104,49,step=1,color='F');ring(s,252,104,49,step=3,color='B')
    s.t(85,183,'P',14,'F',anchor='middle');s.t(252,183,'P³',14,'B',anchor='middle')
    s.line(324,103,375,103,'muted',1,arrow=True)
    rotation(s,444,102,40,90,'F','J')
    rotation(s,592,102,40,-90,'B','−J')
    s.t(444,196,'CP = JC',12.5,anchor='middle')
    s.t(592,196,'CP³ = −JC',12.5,anchor='middle')
    s.label('b',18,238,text(s,'The full competing families','全部竞争族'))
    x0,x1=140,648;A=lambda phi:x0+(phi+90)/180*(x1-x0)
    s.rect(A(-90),260,A(-45)-A(-90),145,'lightB')
    s.rect(A(45),260,A(90)-A(45),145,'lightF')
    s.t((A(-90)+A(-45))/2,255,text(s,'BB only','仅 BB'),11.5,'B',anchor='middle')
    s.t((A(45)+A(90))/2,255,text(s,'FF only','仅 FF'),11.5,'F',anchor='middle')
    rows=[('FF',0,90,'F'),('BB',-90,0,'B'),('FB',-45,45,'purple'),('BF',-45,45,'purple')]
    for j,(name,a,b,col) in enumerate(rows):
        yy=278+j*34
        s.t(112,yy+4,name,12,col,550,anchor='end')
        s.line(A(a),yy,A(b),yy,col,5)
        for phi in [a,b]:s.circle(A(phi),yy,3.3,'white' if abs(phi) in [45,90] else col,col,1)
    for phi in [-90,-45,0,45,90]:
        s.line(A(phi),404,A(phi),409,'muted',.8)
        s.t(A(phi),426,f'{phi}°',11.5,'muted',anchor='middle')
    s.line(x0,404,x1,404,'muted',.9)
    for phi,col in [(30,'grey'),(70,'F')]:
        s.line(A(phi),262,A(phi),393,col,.8,'3 3')
    s.t(A(30),445,'30°',11.5,'grey',anchor='middle');s.t(A(70),445,'70°',11.5,'F',anchor='middle')
    ANALYTIC['fig03']={'discount_domain':'0 <= gamma, eta < 1','closures':{'FF':[0,90],'BB':[-90,0],'FB':[-45,45],'BF':[-45,45]},'open_endpoints_deg':[-90,-45,45,90],'state_task':'eight-cycle steps 1 and 3; common k=2 plane'}
    return s.save()

def fig04(lang):
    s=SVG(4,'fig04_measurement',lang,505,'Compatible readout and joint uncertainty preserve a scientific property')
    s.label('a',18,24,text(s,'Preserve the task action','保留任务作用'))
    R=np.array([[.8,.2,-.4,.1],[-.1,.5,.3,-.6]])
    X=np.array([[.1,.7,-.2,.4],[.4,-.2,.6,.8]])
    J=np.array([[0.,-1],[1,0]]);T=.5*np.eye(2)+.8*J
    z=lambda x:np.array([np.sum(R*x),np.sum((J@R)*x)])/np.linalg.norm(R)
    assert np.allclose(z(T@X),T@z(X))
    s.heat(X,56,59,96,48,'F',True);s.heat(T@X,269,59,96,48,'F',True)
    s.t(104,48,'X',12,anchor='middle');s.t(317,48,'GX',12,anchor='middle')
    s.line(172,82,245,82,'muted',1,arrow=True);s.t(209,68,'G',13,anchor='middle')
    s.line(104,120,104,161,'teal',1,arrow=True);s.line(317,120,317,161,'teal',1,arrow=True)
    s.t(119,147,'z',13,'teal');s.t(332,147,'z',13,'teal')
    # Two actual coordinate vectors, scaled jointly.
    bz=z(X);yz=z(T@X);norm=max(np.linalg.norm(bz),np.linalg.norm(yz));scale=31/norm
    for cx,v in [(104,bz),(317,yz)]:
        s.line(cx-35,198,cx+40,198,'line',.7);s.line(cx,224,cx,170,'line',.7)
        s.line(cx,198,cx+v[0]*scale,198-v[1]*scale,'F',1.8,arrow=True)
    s.line(172,195,245,195,'muted',1,arrow=True);s.t(209,181,'G',13,anchor='middle')
    s.t(104,246,'z(X)',12,anchor='middle');s.t(317,246,'z(GX)',12,anchor='middle')
    s.t(210,272,'z(GX) = Gz(X)',13,'teal',anchor='middle')

    s.label('b',424,24,text(s,'Joint error','联合误差'))
    S=np.loadtxt(DATA/'fig06_joint_covariance.csv',delimiter=',',skiprows=1)
    s.heat(S,472,66,128,128,'F',True)
    labs=['b₁','b₂','y₁','y₂']
    for j,l in enumerate(labs):
        s.t(488+j*32,55,l,11.5,anchor='middle');s.t(460,86+j*32,l,11.5,anchor='end')
    s.line(536,66,536,194,'ink',.9);s.line(472,130,600,130,'ink',.9)
    s.t(536,222,'rᵢ(T) = yᵢ − Tbᵢ',13,anchor='middle')
    s.t(536,247,text(s,'Reference + target','参考 + 目标'),12,'muted',anchor='middle')

    s.label('c',18,315,text(s,'Retain all compatible explanations','保留全部相容解释'))
    meta=json.loads((DATA/'fig06_metadata.json').read_text())
    center=np.array(meta['joint_acceptance_circle_center']);radius=meta['joint_acceptance_radius'];mean=np.array(meta['target_mean'])
    origin=(118,469);scale=85;U=lambda v:origin[0]+v*scale;V=lambda v:origin[1]-v*scale
    angles=np.linspace(0,math.pi/2,91)
    s.path([origin]+[(U(1.28*math.cos(a)),V(1.28*math.sin(a))) for a in angles[45:]]+[origin],'none',0,'lightF',close=True)
    s.line(U(0),V(0),U(1.65),V(0),'muted',.8,arrow=True);s.line(U(0),V(0),U(0),V(1.35),'muted',.8,arrow=True)
    s.line(U(0),V(0),U(1.2),V(1.2),'grey',.8,'3 3')
    s.t(U(1.65)+4,V(0)+4,'u',11.5);s.t(U(0)-7,V(1.35)-4,'v',11.5)
    s.circle(U(center[0]),V(center[1]),scale*radius,'lightF','F',1.6)
    s.circle(U(mean[0]),V(mean[1]),2.8,'ink')
    s.t(235,368,text(s,'Compatible set','相容集合'),11.5,'F');s.line(230,373,U(center[0])+scale*radius,V(center[1]),'F',.8)
    s.t(224,416,'45°',11.5,'muted')
    s.t(170,489,text(s,'Prescribed moments','指定矩示例'),11.5,'muted',anchor='middle')
    labels=[('F','Shared property','共同性质'),('purple','Ambiguous','仍有歧义'),('B','No compatible model','无相容模型'),('grey','Ineligible measurement','测量条件未满足')]
    for j,(col,en,zh) in enumerate(labels):
        yy=350+j*40
        s.circle(430,yy-4,4.3,col)
        s.t(447,yy,text(s,en,zh),12,col)
    ANALYTIC['fig04']={'R':R.tolist(),'X':X.tolist(),'G':T.tolist(),'z_X':z(X).tolist(),'z_GX':z(T@X).tolist(),'readout_equivariance_error':float(np.max(np.abs(z(T@X)-T@z(X)))),'acceptance_set_source':'source_data/fig06_metadata.json','statistical_scope':'candidate-specific two-dimensional residual inversion; not four-dimensional mean confidence region','n':24,'alpha':.05}
    return s.save()

def fig05(lang):
    s=SVG(5,'fig05_quantitative',lang,500,'From a parameter set to point identification through independent propagation')
    s.label('a',18,24,text(s,'A sharp parameter set','尖锐参数集合'))
    X,Y=axes(s,67,57,220,153,(0,1),(0,1),[0,.5,1],[0,.5,1],xlabel='γ',ylabel='η')
    low=math.tan(math.radians(25));gs=np.linspace(low+1e-5,1-1e-5,301);etas=np.tan(math.radians(70)-np.arctan(gs))
    s.rect(X(0),Y(1),X(low)-X(0),153,'pale')
    s.rect(X(low),Y(low),X(1)-X(low),Y(0)-Y(low),'pale')
    s.path([(X(g),Y(e)) for g,e in zip(gs,etas)],'F',2)
    s.circle(X(low),Y(1),3,'white','F',1);s.circle(X(1),Y(low),3,'white','F',1)
    s.t(83,81,'φ = 70°',12,'F')
    s.t(181,273,'γ, η > 0.466',12,'F',anchor='middle')

    s.label('b',359,24,text(s,'Independent propagation','独立传播形状'))
    M=np.eye(3)
    for j,label in enumerate(['B','PB','P²B']):
        xx=398+j*82
        s.t(xx+22,63,label,13,anchor='middle')
        s.heat(M[:,j:j+1],xx,80,44,114,'F',vmax=1)
    s.t(520,230,text(s,'rank [B, PB, P²B] = 3','秩 [B, PB, P²B] = 3'),12,'teal',anchor='middle')

    s.label('c',18,296,text(s,'Two visible task planes','两个可见任务平面'))
    rotation(s,105,371,41,90,'F','90°')
    rotation(s,253,371,41,45,'teal','45°')
    s.t(179,453,text(s,'Shared c, α, γ','共同 c、α、γ'),12,anchor='middle')
    s.t(179,479,'c₁ = c₂',13,'teal',anchor='middle')

    s.label('d',359,296,text(s,'Recover the discount','恢复折扣'))
    Xd,Yd=axes(s,405,333,237,117,(0,1),(.25,1),[0,.5,1],[.3,.6,.9],xlabel='γ',ylabel='q')
    gs=np.linspace(0,1,251);q=1-np.sqrt(2)*gs/(1+gs*gs)
    s.path([(Xd(g),Yd(v)) for g,v in zip(gs,q)],'F',1.9)
    q0=1-math.sqrt(2)*.5/1.25
    s.line(Xd(.5),Yd(.25),Xd(.5),Yd(q0),'grey',.8,'3 3');s.line(Xd(0),Yd(q0),Xd(.5),Yd(q0),'grey',.8,'3 3')
    s.circle(Xd(.5),Yd(q0),3.8,'F');s.circle(Xd(1),Yd(q[-1]),3,'white','F',1)
    s.t(550,355,'γ = 0.50',11.5,'F')
    ANALYTIC['fig05']={'partial_phi_deg':70,'discount_lower_infimum':low,'ratio_gamma':.5,'ratio_q':q0,'three_shape_matrix':M.tolist()}
    return s.save()

def fig06(lang):
    s=SVG(6,'fig06_kernels',lang,535,'Finite memory makes temporal weights separately visible')
    s.label('a',18,24,text(s,'Four-bit memory','四位记忆'))
    rows=[['0','1','0','1'],['1','0','1','b₁'],['0','1','b₁','b₂'],['1','b₁','b₂','b₃'],['b₁','b₂','b₃','b₄']]
    for i,row in enumerate(rows):
        yy=55+i*32
        s.t(53,yy+20,str(i),11.5,'muted',anchor='end')
        for j,v in enumerate(row):
            col='teal' if v.startswith('b') else 'grey'
            s.rect(75+j*40,yy,32,26,blend(col,.13),col,.7)
            s.t(91+j*40,yy+18,v,12,col,anchor='middle')
    s.t(53,40,'t',12,'muted',anchor='end')
    s.t(155,234,text(s,'16 states · 2 successors','16 个状态 · 2 个后继'),12,'muted',anchor='middle')
    s.label('b',359,24,text(s,'Each lag has its own position','每个滞后分别可见'))
    cols=['grey','F','teal','purple'];mx,my,cell=433,58,38
    for i in range(4):
        for j in range(4):
            k=i-j
            s.rect(mx+j*cell,my+i*cell,cell,cell,blend(cols[k],.17) if k>=0 else 'white','line',.6)
            if k>=0:s.t(mx+(j+.5)*cell,my+(i+.5)*cell+4,['w₀','w₁','w₂','w₃'][k],12,cols[k],anchor='middle')
    s.t(509,235,'f(N)   ·   N⁴ = 0',13,anchor='middle')

    s.label('c',18,287,text(s,'Test the weight sequence','检验权重序列'))
    X,Y=axes(s,67,325,225,137,(1,3),(0,.65),[1,2,3],[0,.2,.4,.6],xlabel=text(s,'Lag','滞后'),ylabel=text(s,'Weight','权重'))
    lag=np.arange(1,4);single=.6**lag;mixture=.5*(.2**lag+.85**lag)
    for vs,col in [(single,'F'),(mixture,'purple')]:
        s.path([(X(k),Y(v)) for k,v in zip(lag,vs)],col,1.6)
        for k,v in zip(lag,vs):s.circle(X(k),Y(v),3.1,col)
    s.line(83,312,102,312,'F',1.6);s.t(109,316,text(s,'Single','单尺度'),11.5,'F')
    s.line(190,312,209,312,'purple',1.6);s.t(216,316,text(s,'Mixture','混合'),11.5,'purple')
    s.t(179,520,'Δ = w₁w₃ − w₂²',12,anchor='middle')

    s.label('d',359,287,text(s,'Separate a learned reference','分离已学习参考'))
    N=np.diag(np.ones(3),-1);gamma=.6;eta=.4
    T=np.linalg.inv(np.eye(4)-gamma*N)@(np.eye(4)-eta*N.T)
    assert np.allclose(T[:,0],gamma**np.arange(4))
    s.heat(T,413,324,168,168,'F',True,vmax=1,annot=True,annotation_format='.3g')
    s.rect(412,323,43,170,'none','teal',2)
    s.line(590,407,629,407,'teal',1,arrow=True)
    s.lines(642,354,['f₀','f₁','f₂','f₃'],12,'teal',anchor='middle',lead=42)
    s.t(498,520,'T = af(N)q(Nᵀ)⁻¹',12,anchor='middle')
    ANALYTIC['fig06']={'single_discount':.6,'mixture_discounts':[.2,.85],'mixture_weights':[.5,.5],'single_weights':single.tolist(),'mixture_weights_lags':mixture.tolist(),'moment_gap_single':float(single[0]*single[2]-single[1]**2),'moment_gap_mixture':float(mixture[0]*mixture[2]-mixture[1]**2),'reverse_reference_gamma':gamma,'reverse_reference_eta':eta,'relative_matrix':T.tolist()}
    return s.save()

def loadcsv(name):
    with (DATA/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def fig07(lang):
    s=SVG(7,'fig07_calibration',lang,535,'What repeated neural measurements support in calibration')
    s.label('a',18,24,text(s,'Reference prediction','参考预测'))
    rows=loadcsv('fig07_reference_blocks.csv');agg=loadcsv('fig07_reference_aggregate.csv')
    X,Y=axes(s,67,59,220,150,(0,1),(-12,64),[],[0,20,40,60],ylabel=text(s,'Gain (%)','改善（%）'))
    xx=[115,239]
    byblock={}
    for r in rows:byblock.setdefault(r['block_label'],{})[r['window']]=float(r['gain_percent'])
    for k,(name,r) in enumerate(byblock.items()):
        jit=(k-1.5)*5
        s.line(xx[0]+jit,Y(r['early']),xx[1]+jit,Y(r['late']),'line',.9)
        for j,win in enumerate(['early','late']):s.circle(xx[j]+jit,Y(r[win]),3.2,'white','grey',1.3)
    for j,win in enumerate(['early','late']):
        a=next(r for r in agg if r['window']==win);gain=float(a['gain_percent'])
        s.diamond(xx[j],Y(gain),4.2,'F');s.t(xx[j]+12,Y(gain)+4,f'{gain:.2f}',11.5,'F')
        s.t(xx[j],232,text(s,win.capitalize(),'早窗' if win=='early' else '晚窗'),12,'muted',anchor='middle')
    s.t(177,254,text(s,'4 acquisition blocks','4 个实际采集块'),11.5,'muted',anchor='middle')

    s.label('b',359,24,text(s,'Use the relevant error space','使用相关误差空间'))
    Xb,Yb=axes(s,404,59,238,150,(0,3),(0,85),[],[0,20,40,60,80],ylabel=text(s,'Conditions','条件数'))
    dfb=loadcsv('dFB_error_space_counts.csv')
    order=['full','common_prediction','common_prediction_current_removed']
    vals=[int(next(r for r in dfb if r['space']==k)['multiclass']) for k in order]
    assert vals==[0,6,77] and all(int(r['total'])==864 for r in dfb)
    labels_en=[['Full'],['Predictive'],['Current','profiled']];labels_zh=[['全空间'],['预测空间'],['消去','当前项']]
    for j,v in enumerate(vals):
        cx=443+78*j
        if v:s.rect(cx-18,Yb(v),36,Yb(0)-Yb(v),'F' if j==2 else 'grey')
        else:s.circle(cx,Yb(0),2.6,'grey')
        s.t(cx,Yb(v)-9,str(v),12,'F' if j==2 else 'muted',550,anchor='middle')
        s.lines(cx,232,(labels_zh if lang=='zh' else labels_en)[j],11.5,'muted',anchor='middle',lead=15)
    s.t(641,49,text(s,'864 conditions','864 个条件'),11.5,'muted',anchor='end')

    s.label('c',18,300,text(s,'Known-transform calibration','已知变换校准'))
    cond=loadcsv('fig07_conditional_sr.csv')
    Xc,Yc=axes(s,67,337,220,125,(0,3),(0.8,1.42),[],[.8,1,1.2,1.4],ylabel=text(s,'Boundary / threshold','边界量 / 阈值'))
    s.line(67,Yc(1),287,Yc(1),'grey',1,'4 3')
    for j,r in enumerate(cond):
        xx0=104+j*73;rr=float(r['ratio'])
        s.line(xx0,Yc(1),xx0,Yc(rr),'B',2.2);s.circle(xx0,Yc(rr),4,'B')
        s.t(xx0,Yc(rr)-10,f'{rr:.2f}',11.5,'B',anchor='middle')
        labs=(['Native','γ = .95'] if j==0 else ['Designed',f'γ = {float(r["gamma"]):.2f}']) if lang=='en' else (['原生任务','γ = .95'] if j==0 else ['设计任务',f'γ = {float(r["gamma"]):.2f}'])
        s.lines(xx0,485,labs,11.5,'muted',anchor='middle',lead=16)
    s.t(288,326,'n = 17',11.5,'muted',anchor='end')

    s.label('d',359,300,text(s,'Acquisition and transfer','采集与跨记录'))
    blocks=loadcsv('fig07_within_physical_block.csv');transport=loadcsv('fig07_transport_summary.csv')
    assert [int(r['n_pairs']) for r in blocks]==[6,7]
    assert all(r['decision']=='abstain' for r in blocks)
    assert sum(int(r['computational_conditions']) for r in transport)==1200
    assert all(int(r['declared'])==0 for r in transport)
    s.t(411,339,text(s,'Evaluation','评价单位'),11.5,'muted')
    s.t(641,339,text(s,'Declarations','声明数'),11.5,'muted',anchor='end')
    labs=[('Block 1','实际块 1','6 pairs','6 对'),('Block 2','实际块 2','7 pairs','7 对'),('Across records','跨记录','1,200 conditions','1,200 个条件')]
    for j,(en,zh,ne,nz) in enumerate(labs):
        yy=370+j*55
        s.line(402,yy+23,647,yy+23,'grid',.8)
        s.t(411,yy,text(s,en,zh),12)
        s.t(411,yy+17,text(s,ne,nz),11.5,'muted')
        s.circle(603,yy+2,4.1,'white','grey',1.2)
        s.t(641,yy+6,'0',12,'grey',anchor='end')
    ANALYTIC['fig07']={'sources':['fig07_reference_blocks.csv','fig07_reference_aggregate.csv','fig07_conditional_sr.csv','fig07_within_physical_block.csv','fig07_transport_summary.csv'],'dfb_count_source':'SR_dFB_report_zh.pdf p. 9 and SR_dFB_key_results.csv; 2026-11-24 post hoc deterministic calibration','dfb_multiclass_counts':vals,'dfb_total':864,'conditional_ratios':[float(r['ratio']) for r in cond],'physical_block_counts':[6,7],'transport_conditions':1200}
    return s.save()

def fig08(lang):
    s=SVG(8,'fig08_programme',lang,530,'General and particular models and the identification of common properties')
    s.label('a',18,24,text(s,'General models and particular fits','一般模型与特殊拟合'))
    # The author's two-level diagram: general claim above, actual fit below.
    boxes=[(52,58,211,57,'pale','line','General model','一般模型'),
           (415,58,211,57,'pale','line','Psychological mechanism','心理机制'),
           (52,218,211,57,'lightF','F','Particular model','特殊模型'),
           (415,218,211,57,'lightF','F','Observed data','实验数据')]
    for x,y,w,h,fill,stroke,en,zh in boxes:
        s.rect(x,y,w,h,fill,stroke,.9,rx=3)
        s.t(x+w/2,y+33,text(s,en,zh),13.2 if en!='Psychological mechanism' or lang=='zh' else 12.1,anchor='middle')
    s.line(278,86,398,86,'B',1.4,'4 3',True)
    s.t(339,69,text(s,'General claim','一般机制断言'),11.5,'B',anchor='middle')
    s.line(278,246,398,246,'F',1.6,arrow=True)
    s.t(339,267,text(s,'Fit','拟合'),12,'F',anchor='middle')
    s.line(157.5,127,157.5,205,'muted',1,arrow=True)
    s.line(520.5,127,520.5,205,'muted',1,arrow=True)
    s.t(170,171,text(s,'Specify','具体化'),12,'muted')
    s.t(533,171,text(s,'Observe','观测'),12,'muted')
    s.line(339,225,339,105,'B',1.5,'4 3',True)
    s.rect(287,146,104,43,'white')
    s.lines(339,163,[text(s,'Identify','辨识'),text(s,'the property','机制性质')],12,'B',550,anchor='middle',lead=17)

    s.label('b',18,327,text(s,'Parameter compensation','参数补偿'))
    X,Y=axes(s,64,364,197,105,(0,1),(1,2),[0,.5,1],[1,1.5,2],xlabel='γ')
    gs=np.linspace(0,1,201)
    s.path([(X(g),Y(1+g)) for g in gs],'F',1.8)
    s.path([(X(g),Y(1+g*g)) for g in gs],'purple',1.8)
    s.circle(X(1),Y(2),3,'white','F',1)
    s.t(85,380,'c',12,'F');s.t(194,448,'α',12,'purple')
    s.line(273,413,307,413,'muted',1,arrow=True)
    s.t(313,418,'I + J',12,'F')

    s.label('c',359,327,text(s,'Properties to quantities','从性质到定量'))
    s.t(518,366,text(s,'Direction · structure · weights','方向 · 结构 · 权重'),12,anchor='middle')
    s.line(518,377,518,392,'teal',1,arrow=True)
    s.rect(401,402,234,39,'lightF','F',.8,rx=3)
    s.t(518,426,text(s,'q(θ) = constant on C(y)','q(θ) 为常数，θ ∈ C(y)'),12,'F',anchor='middle')
    s.line(518,450,518,465,'teal',1,arrow=True)
    s.t(518,486,text(s,'Parameter sets · kernels','参数域 · 时间核'),12,anchor='middle')
    s.line(413,508,438,508,'purple',1,'3 3',True)
    s.t(447,512,text(s,'Lie symmetries: future work','李群对称：后续研究'),11.5,'purple')
    ANALYTIC['fig08']={'panel_a':'author original two-level model/psychological-mechanism and fitted-instance/observed-data argument preserved','panel_b':'c=1+gamma, alpha=1+gamma^2 yields I+J for 0<=gamma<1, direct-reference quarter-turn block','panel_c':'property constancy on compatible set, additional independent constraints support quantitative targets; Lie-symmetry inverse programme is future work'}
    return s.save()

def export(p):
    stem=p.stem
    subprocess.run(['inkscape',str(p),'--export-type=pdf','--export-text-to-path',f'--export-filename={ROOT/"figures/pdf"/(stem+".pdf")}'],check=True,capture_output=True)
    subprocess.run(['inkscape',str(p),'--export-type=png','--export-width=2040',f'--export-filename={ROOT/"figures/png"/(stem+".png")}'],check=True,capture_output=True)
    return stem

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--figures',nargs='*',type=int);ap.add_argument('--no-export',action='store_true');args=ap.parse_args()
    funcs=[fig01,fig02,fig03,fig04,fig05,fig06,fig07,fig08]
    requested=args.figures or list(range(1,9));paths=[]
    for n in requested:
        for lang in ['en','zh']:paths.append(funcs[n-1](lang))
    if not args.no_export:
        with ThreadPoolExecutor(max_workers=3) as pool:
            for stem in pool.map(export,paths):print(stem)
    (DATA/'redrawn_figure_registry.json').write_text(json.dumps(REGISTRY,ensure_ascii=False,indent=2),encoding='utf-8')
    (DATA/'redrawn_analytic_values.json').write_text(json.dumps(ANALYTIC,ensure_ascii=False,indent=2),encoding='utf-8')

if __name__=='__main__':main()

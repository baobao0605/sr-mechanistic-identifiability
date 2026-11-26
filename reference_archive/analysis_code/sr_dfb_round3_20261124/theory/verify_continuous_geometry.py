"""Independent algebraic/numerical checks; these are not synthetic subjects."""
from pathlib import Path
import json,math,time
import numpy as np
from scipy.linalg import block_diag
from continuous_geometry import SpectralSRGeometry

ROOT=Path(__file__).resolve().parent
rng=np.random.default_rng(20261124)
records=[]

def check(name,value,tolerance,**extra):
    value=float(value);passed=bool(value<=tolerance)
    records.append(dict(name=name,error=value,tolerance=tolerance,passed=passed,**extra))
    if not passed:raise AssertionError((name,value,tolerance))

def direct(A,B,g,kind):
    K=A if kind=="forward" else A.T if kind=="backward" else (A+A.T)/2
    X=K@np.linalg.solve(np.eye(len(A))-g*K,B)
    return X-B*np.vdot(B,X).real/np.vdot(B,B).real

def dense_point(geo,kt,gt,kr,lo=0,hi=.95,count=12001):
    u=geo.curve(kt,gt)["unit"]
    lam=geo._mode_values[kr]
    gs=np.linspace(lo,hi,count)
    z=lam[None,:]/(1-gs[:,None]*lam)
    z-=np.real(z@geo.p)[:,None]
    norm=np.sqrt(np.maximum(0,abs(z)**2@geo.p))
    co=np.real(z.conj()@(geo.p*u))/norm
    co=np.clip(co,-1,1)
    c=float(co.max());return math.sqrt(max(0,1-max(0,c)**2)),float(gs[np.argmax(co)])

for n in (2,4,7,13):
    A=np.array([[0.,-1.],[1.,0.]]) if n==2 else np.roll(np.eye(n),-1,axis=1)
    B=rng.normal(size=(n,6))
    if n>2:B-=B.mean(0)
    geo=SpectralSRGeometry.from_reference(A,B)
    for kind in ("forward","backward","symmetric_transition"):
        for gamma in (0.,.25,.5,.8,.95):
            q=direct(A,B,gamma,kind)
            check(f"direct_residual_norm_n{n}_{kind}_{gamma}",
                  abs(np.linalg.norm(q)-geo.residual_norm(kind,gamma)),1e-10)
    for gt in (.25,.5,.8):
        for kr in ("backward","symmetric_transition"):
            bound=geo.point_margin("forward",gt,kr,atol=1e-5,max_evals=10001)
            if geo.is_current_equivalent(kr):
                check(f"null_family_n{n}_{gt}_{kr}",max(abs(bound["lower"]-1),abs(bound["upper"]-1)),1e-12)
            else:
                dense,g=dense_point(geo,"forward",gt,kr)
                violation=max(0,bound["lower"]-dense,dense-bound["upper"]-3e-6)
                check(f"continuous_dense_enclosure_n{n}_{gt}_{kr}",violation,2e-10,
                      lower=bound["lower"],upper=bound["upper"],dense=dense,
                      evaluations=bound["evaluations"],status=bound["status"])
    # Direct real Krylov basis, independent from the spectral projection.
    cols=[];x=B.copy()
    for _ in range(n):cols.append(x.ravel());x=A@x
    U,s,_=np.linalg.svd(np.stack(cols,axis=1),full_matrices=False)
    U=U[:,s>s[0]*1e-11]
    E=rng.normal(size=B.shape)
    if n>2:E-=E.mean(0)
    direct_projection=(U@(U.T@E.ravel())).reshape(B.shape)
    projected=geo.project_krylov(E)
    check(f"krylov_svd_projection_n{n}",np.linalg.norm(projected-direct_projection),2e-10)
    f=.8*B+.7*(A@np.linalg.solve(np.eye(n)-.3*A,B))
    b=1.1*B+.2*(A.T@np.linalg.solve(np.eye(n)-.8*A.T,B))
    raw=(np.linalg.norm(E-f)**2-np.linalg.norm(E-b)**2)
    pr=(np.linalg.norm(projected-f)**2-np.linalg.norm(projected-b)**2)
    check(f"sse_difference_preserved_n{n}",abs(raw-pr),2e-10)

# Repeated eigenvalues must be pooled, otherwise the prediction subspace
# would be incorrectly enlarged to arbitrary within-eigenspace changes.
th=.37;R=np.array([[np.cos(th),-np.sin(th)],[np.sin(th),np.cos(th)]])
A=block_diag(R,R);B=rng.normal(size=(4,5));E=rng.normal(size=B.shape)
geo=SpectralSRGeometry.from_reference(A,B)
D=np.stack([np.linalg.matrix_power(A,k).dot(B).ravel() for k in range(4)],axis=1)
U,s,_=np.linalg.svd(D,full_matrices=False);U=U[:,s>s[0]*1e-11]
check("repeated_eigenspace_projection",np.linalg.norm(geo.project_krylov(E).ravel()-U@(U.T@E.ravel())),2e-10,
      krylov_rank=int(U.shape[1]))

# Analytic 90-degree block: all positive F/B rays oppose each other.
A=np.array([[0.,-1.],[1.,0.]]);B=np.array([[1.,2.],[.3,-1.]])
geo=SpectralSRGeometry.from_reference(A,B)
z=geo.family_margin("forward","backward",atol=1e-4,max_evals=50001)
check("rotation_family_ray_exact",max(abs(z["lower"]-1),abs(z["upper"]-1)),1e-12)
check("rotation_family_chord_encloses_two",max(0,z["chord_lower"]-2,2-z["chord_upper"]),1e-12,
      lower=z["chord_lower"],upper=z["chord_upper"],evaluations=z["evaluations"])

# A two-sided strength floor gives a chord, whereas allowing the rival to
# shrink to zero gives a ray residual. At opposite directions: 2 versus 1.
check("opposite_ray_and_chord_distinction",abs(z["lower"]-1)+abs(z["chord_upper"]-2),1e-12)

# API boundary states are explicit; empty domains are empty sets, not a
# positive certificate concerning a nonempty scientific alternative.
z=geo.family_margin("current","current")
check("boundary_current_vs_current",abs(z["lower"])+abs(z["upper"]),0.)
z=geo.family_margin("forward","current")
check("boundary_forward_vs_current",abs(z["lower"]-1)+abs(z["upper"]-1),0.)
z=geo.family_margin("forward","symmetric_transition")
check("boundary_rotation_symmetric_is_current",abs(z["lower"]-1)+abs(z["upper"]-1),0.)
scalar=SpectralSRGeometry.from_reference(.4*np.eye(3),rng.normal(size=(3,2)))
z=scalar.family_margin("forward","backward")
check("boundary_two_degenerate_families",abs(z["lower"])+abs(z["upper"]),0.)
z=geo.point_margin("forward",.5,"backward",intervals=[])
check("boundary_empty_point_domain",0. if math.isinf(z["lower"]) and math.isinf(z["upper"]) and z["status"]=="empty_rival_domain" else 1.,0.)
z=geo.family_margin("forward","backward",intervals_b=[])
check("boundary_empty_family_domain",0. if math.isinf(z["lower"]) and math.isinf(z["upper"]) and z["status"]=="empty_family_domain" else 1.,0.)

# Independent derivative-bound sanity checks inside adaptive intervals.
A=np.roll(np.eye(11),-1,axis=1);B=rng.normal(size=(11,4));B-=B.mean(0)
geo=SpectralSRGeometry.from_reference(A,B)
for kind in ("forward","backward","symmetric_transition"):
    worst=0.
    for a,b in ((0,.1),(.2,.3),(.7,.75),(.9,.95)):
        bound=geo._interval_bounds(kind,a,b)
        if not bound["valid"]:continue
        for g in np.linspace(a,b,17):
            c=geo.curve(kind,g);up=geo.norm(c["unit_prime"])
            worst=max(worst,up-bound["up"])
    check("unit_derivative_bound_"+kind,max(0,worst),1e-10)

# Unknown-reference ambiguity with one calibration center.
Bhat=1.;e0=1/3;TF=1.;TB=2.;BF=4/3;BB=2/3
check("two_latent_references_same_observation",abs(TF*BF-TB*BB),1e-14,
      reference_error_F=abs(BF-Bhat),reference_error_B=abs(BB-Bhat),
      common_B_min_distance=2/3)

# Unbounded current gain can collapse unknown-reference robust family
# separation even when both predictive residuals remain nonzero.
n=5;A=np.roll(np.eye(n),-1,axis=1);Bhat=np.eye(n)-np.ones((n,n))/n
NF=A@np.linalg.inv(np.eye(n)-.5*A);NB=NF.T
collapse=[]
for c in (10.,100.,1000.):
    target=c*Bhat
    tf=c*np.eye(n)+NF;tb=c*np.eye(n)+NB
    bf=np.linalg.solve(tf,target);bb=np.linalg.solve(tb,target)
    qf=direct(A,bf,.5,"forward");qb=direct(A,bb,.5,"backward")
    err=max(np.linalg.norm(tf@bf-target),np.linalg.norm(tb@bb-target))
    check("unbounded_current_gain_"+str(c),err,1e-10)
    collapse.append(dict(current_gain=c,reference_F=float(np.linalg.norm(bf-Bhat)),
                         reference_B=float(np.linalg.norm(bb-Bhat)),
                         predictive_F=float(np.linalg.norm(qf)),predictive_B=float(np.linalg.norm(qb)),
                         target_residual=err))

out={"status":"PASS","scope":"Algebraic and numerical verification, not simulated participants or biological evidence",
     "checks":records,"n_checks":len(records),
     "max_error":max(x["error"] for x in records),
     "unknown_reference_unbounded_current_example":collapse,
     "numeric_scope":"Continuous bounds are analytic inequalities evaluated in float64, not outward-rounded arithmetic."}
(ROOT/"geometry_verification.json").write_text(json.dumps(out,indent=2))
print(json.dumps({k:out[k] for k in ("status","n_checks","max_error")},indent=2))

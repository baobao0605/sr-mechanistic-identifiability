"""SR-specific observable geometry. No synthetic participants.

The main family is Y = c B + alpha N_gamma(P) B, where c is an
unrestricted current-only nuisance gain, alpha >= 0, and
N_gamma(P) = P (I - gamma P)^(-1).  gamma=0 is a one-step kernel,
not the gamma=0 member of the unaugmented M_gamma SR family.
All arrays must already use a common, training-only feature metric.
The optional state map C must be fixed by the design, not fitted to Y.
Grid minima are descriptive numerical searches, not global certificates.
"""
from __future__ import annotations

import numpy as np

KINDS = ("forward", "backward", "symmetric_transition", "balanced_fb")

def templates(P, B, gammas, kind="forward", C=None, backward=None):
    """Return pure-future templates with shape (n_gamma, states, features).

    Default backward=P.T is valid for doubly-stochastic transitions or an
    explicitly defined finite-chain predecessor operator. Nonuniform
    stationary chains require backward=diag(pi)^-1 P.T diag(pi).
    """
    P=np.asarray(P,dtype=float); B=np.asarray(B,dtype=float)
    R=P.T if backward is None else np.asarray(backward,dtype=float)
    I=np.eye(len(P)); out=[]
    for g in np.asarray(gammas,dtype=float):
        if kind=="forward": X=P@np.linalg.solve(I-g*P,B)
        elif kind=="backward": X=R@np.linalg.solve(I-g*R,B)
        elif kind=="symmetric_transition":
            S=(P+R)/2; X=S@np.linalg.solve(I-g*S,B)
        elif kind=="balanced_fb":
            X=(P@np.linalg.solve(I-g*P,B)+R@np.linalg.solve(I-g*R,B))/2
        else: raise ValueError(kind)
        out.append(X if C is None else np.asarray(C)@X)
    return np.stack(out)

def residualize_current(X, B, nuisance=None):
    """Project feature vectors off current B and optional common nuisances.

    nuisance: observation-space columns (B.size, n_nuisance). Whitening
    must precede this operation. Orthogonal projection is numerically
    computed by SVD; the same projection is used for all hypotheses.
    """
    v=np.asarray(B,dtype=float).ravel()
    Z=v[:,None] if nuisance is None else np.c_[v,nuisance]
    U,s,_=np.linalg.svd(Z,full_matrices=False)
    tol=(s[0] if s.size else 0)*max(Z.shape)*np.finfo(float).eps
    U=U[:,s>tol]
    X=np.asarray(X,dtype=float)
    flat=X.reshape((-1,v.size))
    residual=flat-(flat@U)@U.T
    # Preserve exact model degeneracies against machine-precision leftovers.
    zero=np.linalg.norm(residual,axis=1)<=1e-12*np.maximum(
        np.linalg.norm(flat,axis=1),np.finfo(float).tiny)
    residual[zero]=0
    return residual.reshape(X.shape)

def normalized_rays(X, tolerance=1e-12):
    X=np.asarray(X,float).reshape((len(X),-1))
    norms=np.linalg.norm(X,axis=1)
    threshold=tolerance*max(float(norms.max(initial=0)),np.finfo(float).tiny)
    valid=norms>threshold
    out=np.zeros_like(X); out[valid]=X[valid]/norms[valid,None]
    return out,norms,valid

def ray_geometry(X,Y,signed_gain=False):
    """Pairwise angular geometry after the same nuisance projection.

    residual_per_unit = distance from unit X to the positive cone of Y.
    Invalid zero templates are flagged, never converted into evidence.
    """
    ux,nx,vx=normalized_rays(X); uy,ny,vy=normalized_rays(Y)
    cos=np.clip(ux@uy.T,-1,1)
    if signed_gain: cos=np.abs(cos)
    invalid=~(vx[:,None]&vy[None,:])
    dist=np.sqrt(np.maximum(0,1-np.maximum(cos,0)**2))
    chord=np.sqrt(np.maximum(0,2-2*cos))
    cos[invalid]=np.nan;dist[invalid]=np.nan;chord[invalid]=np.nan
    return {"cosine":cos,"residual_per_unit":dist,"unit_chord":chord,
            "x_norm":nx,"y_norm":ny,"x_valid":vx,"y_valid":vy}

def fit_affine(B, X, Y, alpha_nonnegative=True):
    """Fit current gain and predictive strength for every candidate template.

    Unrestricted c is a conservative nuisance. X is (n_candidate,n,p).
    Return all candidate SSEs and coefficients, without choosing a gamma.
    """
    b=np.asarray(B,float).ravel(); x=np.asarray(X,float).reshape((len(X),-1))
    y=np.asarray(Y,float).ravel(); bb=float(b@b)
    if bb<=0: raise ValueError("zero current template")
    by=float(b@y)/bb; bx=(x@b)/bb
    r=x-bx[:,None]*b; yr=y-by*b
    rr=np.sum(r*r,axis=1); r_y=r@yr
    alpha=np.divide(r_y,rr,out=np.zeros_like(r_y),where=rr>np.finfo(float).eps*max(float(np.sum(x*x,axis=1).max()),np.finfo(float).tiny))
    if alpha_nonnegative: alpha=np.maximum(alpha,0)
    c=by-alpha*bx
    err=y[None,:]-c[:,None]*b-alpha[:,None]*x
    return {"sse":np.sum(err*err,axis=1),"current_gain":c,"predictive_gain":alpha}

def kernel_design(P,B,n_lags=3,C=None,drop_current=True):
    """Columns vec(C P^k B), k=1..n_lags after current projection by default."""
    P=np.asarray(P,float);B=np.asarray(B,float)
    X=B.copy(); out=[]
    for k in range(n_lags+1):
        out.append(X if C is None else C@X);X=P@X
    out=np.stack(out)
    if drop_current: out=residualize_current(out[1:],out[0])
    return out.reshape((len(out),-1)).T

def reference_ball_residual(T,Bhat,Yhat,e0):
    """Exact fixed-T Frobenius-ball compatibility, up to numeric root tolerance.

    min_{||Delta||_F <= e0} ||Yhat - T(Bhat+Delta)||_F.
    Output is a numerical solution, not an interval-arithmetic certificate.
    This function does not attach statistical coverage to e0 or e1.
    """
    from scipy.optimize import brentq
    from scipy.linalg import norm as stable_norm
    T=np.asarray(T,float);Bhat=np.asarray(Bhat,float);Yhat=np.asarray(Yhat,float)
    if e0<0 or not np.isfinite(e0): raise ValueError("error radius must be finite and nonnegative")
    fnorm=lambda x:float(stable_norm(np.asarray(x).ravel(),check_finite=False))
    r=Yhat-T@Bhat
    # This endpoint has no finite positive multiplier and needs no SVD.
    if e0==0:
        return {"residual":fnorm(r),"delta":np.zeros_like(Bhat),
                "reference_error_used":0.,"lambda":float("inf"),
                "lambda_scaled":float("inf"),"feasibility_roundoff_adjusted":False}
    U,s,Vh=np.linalg.svd(T,full_matrices=False)
    z=U.T@r
    # Evaluate only strictly positive SVD modes. Zero modes otherwise produce
    # 0/0 at the lower endpoint of the constrained multiplier equation.
    # Tiny positive modes are not dropped merely because of a rank tolerance.
    use=s>0
    delta_coords=np.zeros_like(z)
    delta_coords[use]=z[use]/s[use,None]
    unconstrained_norm=fnorm(delta_coords)
    if unconstrained_norm<=e0:
        delta=Vh.T@delta_coords
        return {"residual":fnorm(r-T@delta),"delta":delta,
                "reference_error_used":fnorm(delta),"lambda":0.,
                "lambda_scaled":0.,"feasibility_roundoff_adjusted":False}
    # Scale lambda by the largest singular value squared. The root equation
    # uses ||Delta|| itself and ratios, not squared denominators that can
    # underflow. Its lambda=0 value is the explicit minimum-norm solution.
    smax=float(s[use].max());sn=s[use]/smax;zn=z[use]/smax
    def coords(eta):
        return (sn/(sn*sn+eta))[:,None]*zn
    def radius(eta):
        return unconstrained_norm if eta==0 else fnorm(coords(eta))
    hi=1.
    while radius(hi)>e0:
        hi*=2
        if not np.isfinite(hi): raise FloatingPointError("could not bracket reference-ball multiplier")
    eta=brentq(lambda v:radius(v)-e0,0.,hi,
               xtol=np.nextafter(0.,1.),rtol=4*np.finfo(float).eps,maxiter=500)
    delta=Vh[use].T@coords(eta)
    used=fnorm(delta);adjusted=used>e0
    if adjusted:
        delta*=e0/used
        used=fnorm(delta)
    # Direct substitution retains arbitrarily small components perpendicular
    # to the range of T; subtracting two total energies can erase them.
    return {"residual":fnorm(r-T@delta),"delta":delta,
            "reference_error_used":used,"lambda":float(eta*smax*smax),
            "lambda_scaled":float(eta),"feasibility_roundoff_adjusted":adjusted}

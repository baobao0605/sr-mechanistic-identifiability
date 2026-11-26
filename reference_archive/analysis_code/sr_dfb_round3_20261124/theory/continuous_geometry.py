"""Continuous SR ray separation for a known, normal state operator.

The analytic interval bounds are evaluated in ordinary double precision;
they are not outward-rounded interval-arithmetic certificates.  Returned
lower/upper brackets cover the continuous gamma domains in exact arithmetic.
The implementation exposes its numerical tolerances and never silently
promotes a finite grid minimum to a continuous lower bound.

Observation metric must already be fixed (e.g. training-only whitening).
Only the common current-reference nuisance is removed.  More nuisance
columns require a different geometry.  This module does not attach
statistical coverage to empirical repeat differences.
"""
from __future__ import annotations

import heapq
import itertools
import math
import numpy as np
from scipy.linalg import schur

ALIASES={"F":"forward","B":"backward","S":"symmetric_transition",
         "symmetric":"symmetric_transition","current":"current"}


def _kind(k):
    k=ALIASES.get(k,k)
    if k not in ("forward","backward","symmetric_transition","current"):
        raise ValueError("Use forward/backward/symmetric_transition/current; balanced_fb is the aligned symmetric class.")
    return k


def _ray_distance(c):
    return math.sqrt(max(0.,1.-max(0.,min(1.,float(c)))**2))


class SpectralSRGeometry:
    """Normal-operator SR geometry using weighted spectral coefficients.

    Use from_reference(A,B), where A is a real normal contraction and B is
    an already whitened reference. F uses A, B uses A.T, S uses (A+A.T)/2.
    Functions of A share a spectral basis, so only modal reference energies
    are needed for geometry. The stored basis also provides the exact
    common Krylov-prediction-space projection (within numerical tolerance).
    """
    def __init__(self,lambdas,weights,*,vectors=None,reference_modes=None,
                 normality_error=0.,spectral_error=0.,zero_tol=1e-12):
        self.lambdas=np.asarray(lambdas,dtype=complex)
        self.weights=np.asarray(weights,dtype=float)
        if len(self.lambdas)!=len(self.weights) or np.any(self.weights<0):
            raise ValueError("Mismatched modes or negative weights")
        self.total_weight=float(self.weights.sum())
        if not self.total_weight>0:raise ValueError("Reference is zero")
        self.p=self.weights/self.total_weight
        if np.max(abs(self.lambdas))>1+1e-9:
            raise ValueError("Operator is not a contraction in this normal basis")
        self.vectors=vectors;self.reference_modes=reference_modes
        self.normality_error=float(normality_error)
        self.spectral_error=float(spectral_error)
        self.zero_tol=float(zero_tol)
        self._curve_cache={};self._interval_cache={}
        self._mode_values={
            "forward":self.lambdas,
            "backward":self.lambdas.conj(),
            "symmetric_transition":self.lambdas.real.astype(complex)}
        # Group repeated eigenvalues, since a polynomial must share its
        # coefficient across every vector in the same eigenspace.
        self.groups=[]
        for i,lam in enumerate(self.lambdas):
            for group in self.groups:
                if abs(lam-self.lambdas[group[0]])<=1e-10:
                    group.append(i);break
            else:self.groups.append([i])

    @classmethod
    def from_reference(cls,A,B,**kwargs):
        A=np.asarray(A,dtype=float);B=np.asarray(B,dtype=float)
        if A.ndim!=2 or A.shape[0]!=A.shape[1] or B.shape[0]!=len(A):
            raise ValueError("A must be square and match B rows")
        normality=np.linalg.norm(A@A.T-A.T@A)
        if normality>1e-9*max(1.,np.linalg.norm(A)**2):
            raise ValueError("Spectral shortcut requires a normal operator")
        S,V=schur(A,output="complex")
        off=S-np.diag(np.diag(S));err=np.linalg.norm(off)
        if err>1e-9*max(1.,np.linalg.norm(A)):
            raise ValueError("Schur form is not diagonal")
        modes=V.conj().T@B
        weights=np.sum(abs(modes)**2,axis=1)
        obj=cls(np.diag(S),weights,vectors=V,reference_modes=modes,
                normality_error=normality,spectral_error=err,**kwargs)
        # Computing the symmetric operator first preserves exact zero for
        # skew-symmetric rotation blocks. Re(lambda) from a floating Schur
        # factorization alone could invent a tiny nonzero symmetric signal.
        sym=(A+A.T)/2
        obj._mode_values["symmetric_transition"]=np.diag(V.conj().T@sym@V)
        return obj

    def inner(self,x,y):
        return float(np.real(np.vdot(x*self.p,y)))

    def norm(self,x):
        return math.sqrt(max(0.,self.inner(x,x)))

    def _center(self,z):
        # Real scalar nuisance: the imaginary weighted mean cancels for
        # a real normal A and real B, up to decomposition roundoff.
        return z-float(np.real(np.dot(self.p,z)))

    def curve(self,kind,gamma):
        kind=_kind(kind);g=float(gamma)
        if not (0<=g<1):raise ValueError("gamma must be in [0,1)")
        key=(kind,g)
        if key in self._curve_cache:return self._curve_cache[key]
        if kind=="current":
            out={"norm":0.,"unit":np.zeros_like(self.lambdas),
                 "unit_prime":np.zeros_like(self.lambdas),"q_prime_norm":0.,"valid":False}
        else:
            lam=self._mode_values[kind];den=1-g*lam
            z=lam/den;qp=self._center(lam**2/den**2)
            q=self._center(z);s=self.norm(q)
            valid=s>self.zero_tol*max(self.norm(z),np.finfo(float).tiny)
            if valid:
                u=q/s;up=(qp-u*self.inner(u,qp))/s
            else:u=np.zeros_like(q);up=np.zeros_like(q)
            out={"norm":s,"unit":u,"unit_prime":up,
                 "q_prime_norm":self.norm(qp),"valid":valid}
        self._curve_cache[key]=out
        return out

    def residual_norm(self,kind,gamma):
        return self.curve(kind,gamma)["norm"]*math.sqrt(self.total_weight)

    def is_current_equivalent(self,kind):
        kind=_kind(kind)
        if kind=="current":return True
        lam=self._mode_values[kind]
        return self.norm(self._center(lam))<=self.zero_tol*max(self.norm(lam),np.finfo(float).tiny)

    def _interval_bounds(self,kind,a,b):
        """Bound ||u'|| and ||u''|| throughout a gamma interval."""
        kind=_kind(kind);a=float(a);b=float(b);key=(kind,a,b)
        if key in self._interval_cache:return self._interval_cache[key]
        m=(a+b)/2;h=(b-a)/2;c=self.curve(kind,m)
        if kind=="current":return {"up":0.,"upp":0.,"valid":False}
        lam=self._mode_values[kind];ab=abs(lam)
        # Exact minimum of |1-gamma*lambda| on [a,b].
        t=np.divide(lam.real,ab**2,out=np.zeros_like(ab),where=ab>0)
        t=np.clip(t,a,b);den_min=abs(1-t*lam)
        raw1=math.sqrt(float(np.dot(self.p,ab**4/den_min**4)))
        raw2=2*math.sqrt(float(np.dot(self.p,ab**6/den_min**6)))
        L1=min(raw1,c["q_prime_norm"]+h*raw2)
        qlow=c["norm"]-h*L1
        if qlow<=0 or not c["valid"]:
            out={"up":float("inf"),"upp":float("inf"),"valid":False}
        else:
            up=L1/qlow
            out={"up":up,"upp":raw2/qlow+3*up**2,"valid":True}
        self._interval_cache[key]=out
        return out

    @staticmethod
    def _domains(intervals):
        intervals=[(0.,.95)] if intervals is None else [(float(a),float(b)) for a,b in intervals]
        if any(a<0 or a>b or b>=1 for a,b in intervals):raise ValueError("Domains must satisfy 0<=a<=b<1")
        return intervals

    def point_margin(self,kind_true,gamma_true,kind_rival,*,intervals=None,
                     atol=1e-4,max_evals=2049):
        """Unit true-residual distance to an entire continuous positive-ray family.

        Current nuisance is free in both models; rival gain is nonnegative
        and includes zero.  Multiply bounds by the actual true predictive
        residual norm in this reference to get a measurement-space margin.
        """
        kt=_kind(kind_true);kr=_kind(kind_rival);dom=self._domains(intervals)
        true=self.curve(kt,gamma_true)
        if not true["valid"]:
            return {"lower":0.,"upper":0.,"status":"no_true_residual",
                    "evaluations":0,"gamma_witness":None,"cosine_lower":None,"cosine_upper":None}
        if not dom:
            return {"lower":float("inf"),"upper":float("inf"),
                    "status":"empty_rival_domain","evaluations":0,"gamma_witness":None,
                    "cosine_lower":None,"cosine_upper":None}
        if kr=="current":
            return {"lower":1.,"upper":1.,"status":"exact_current",
                    "evaluations":0,"gamma_witness":None,"cosine_lower":0.,"cosine_upper":0.}
        u=true["unit"];best=-float("inf");bestg=None;calls=0;serial=itertools.count();heap=[]
        # All spectral values equal implies a current-only rival for every
        # gamma, proven by its scalar action on the observed eigenspaces.
        if self.is_current_equivalent(kr):
            return {"lower":1.,"upper":1.,"status":"numerical_current_equivalence",
                    "evaluations":0,"gamma_witness":None,"cosine_lower":0.,"cosine_upper":0.}
        def ev(g):
            nonlocal best,bestg,calls
            v=self.curve(kr,g);calls+=1
            if not v["valid"]:return None
            c=float(np.clip(self.inner(u,v["unit"]),-1.,1.))
            if c>best:best=c;bestg=float(g)
            cp=self.inner(u,v["unit_prime"])
            return c,cp
        def add(a,b):
            m=(a+b)/2;h=(b-a)/2;value=ev(m)
            if value is None:upper=1.
            else:
                c,cp=value;bound=self._interval_bounds(kr,a,b)
                if bound["valid"]:
                    upper=min(1.,c+h*bound["up"],c+abs(cp)*h+.5*bound["upp"]*h*h)
                    upper=min(1.,upper+128*np.finfo(float).eps)
                else:upper=1.
            heapq.heappush(heap,(-upper,next(serial),a,b))
        for a,b in dom:
            ev(a)
            if b!=a:ev(b);add(a,b)
            else:heapq.heappush(heap,(-best,next(serial),a,b))
        status="max_evaluations"
        while heap:
            cub=max(best,-heap[0][0]);lo=_ray_distance(cub);hi=_ray_distance(best)
            if hi-lo<=atol:
                status="analytic_bracket_tolerance";break
            if calls+2>max_evals:break
            neg,_,a,b=heapq.heappop(heap)
            if -neg<=best:continue
            m=(a+b)/2
            if m==a or m==b:
                heapq.heappush(heap,(neg,next(serial),a,b));status="machine_resolution";break
            add(a,m);add(m,b)
        cub=max(best,-heap[0][0]) if heap else best
        return {"lower":_ray_distance(cub),"upper":_ray_distance(best),
                "cosine_lower":best,"cosine_upper":cub,
                "gamma_witness":bestg,"evaluations":calls,"status":status,
                "numeric_scope":"analytic interval bounds evaluated in float64; no outward rounding"}

    def family_margin(self,kind_a,kind_b,*,intervals_a=None,intervals_b=None,
                      atol=1e-4,max_evals=20001):
        """Continuous two-gamma maximum correlation and separation bounds.

        lower/upper: true residual >=1 versus rival amplitude >=0.
        chord_lower/chord_upper: BOTH residual amplitudes >=1.
        Scale by a common positive amplitude floor for the chosen question.
        """
        ka=_kind(kind_a);kb=_kind(kind_b)
        da=self._domains(intervals_a);db=self._domains(intervals_b)
        if not da or not db:
            return {"lower":float("inf"),"upper":float("inf"),
                    "chord_lower":float("inf"),"chord_upper":float("inf"),
                    "evaluations":0,"status":"empty_family_domain"}
        za=self.is_current_equivalent(ka);zb=self.is_current_equivalent(kb)
        if za or zb:
            dist=0. if za and zb else 1.
            return {"lower":dist,"upper":dist,"chord_lower":dist,"chord_upper":dist,
                    "evaluations":0,"status":"current_equivalence_merged",
                    "scope":"Zero-residual family merged with current; a positive residual floor for that family is infeasible."}
        best=-float("inf");witness=None;calls=0;serial=itertools.count();heap=[]
        def ev(g,h):
            nonlocal best,witness,calls
            x=self.curve(ka,g);y=self.curve(kb,h);calls+=1
            if not(x["valid"] and y["valid"]):return None
            c=float(np.clip(self.inner(x["unit"],y["unit"]),-1,1))
            if c>best:best=c;witness=[float(g),float(h)]
            return c,self.inner(x["unit_prime"],y["unit"]),self.inner(x["unit"],y["unit_prime"])
        def add(a,b,c,d):
            g=(a+b)/2;h=(c+d)/2;hx=(b-a)/2;hy=(d-c)/2;v=ev(g,h)
            bx=self._interval_bounds(ka,a,b);by=self._interval_bounds(kb,c,d)
            if v is None or not(bx["valid"] and by["valid"]):upper=1.
            else:
                f,fx,fy=v
                upper=min(1.,f+bx["up"]*hx+by["up"]*hy,
                          f+abs(fx)*hx+abs(fy)*hy+.5*bx["upp"]*hx*hx+
                          bx["up"]*by["up"]*hx*hy+.5*by["upp"]*hy*hy)
                upper=min(1.,upper+256*np.finfo(float).eps)
            # Relative gamma widths avoid always splitting one dimension
            # while the other remains too wide to bound its normalization.
            heapq.heappush(heap,(-upper,next(serial),a,b,c,d))
        for a,b in da:
            for c,d in db:
                for g,h in ((a,c),(a,d),(b,c),(b,d)):ev(g,h)
                add(a,b,c,d)
        status="max_evaluations"
        while heap:
            cub=max(best,-heap[0][0])
            gap=max(_ray_distance(best)-_ray_distance(cub),
                    math.sqrt(max(0.,2-2*best))-math.sqrt(max(0.,2-2*cub)))
            if gap<=atol:status="analytic_bracket_tolerance";break
            if calls+2>max_evals:break
            neg,_,a,b,c,d=heapq.heappop(heap)
            if -neg<=best:continue
            if (b-a)>=(d-c):
                m=(a+b)/2
                if m==a or m==b:
                    heapq.heappush(heap,(neg,next(serial),a,b,c,d));status="machine_resolution";break
                add(a,m,c,d);add(m,b,c,d)
            else:
                m=(c+d)/2
                if m==c or m==d:
                    heapq.heappush(heap,(neg,next(serial),a,b,c,d));status="machine_resolution";break
                add(a,b,c,m);add(a,b,m,d)
        cub=max(best,-heap[0][0]) if heap else best
        return {"lower":_ray_distance(cub),"upper":_ray_distance(best),
                "chord_lower":math.sqrt(max(0.,2-2*cub)),
                "chord_upper":math.sqrt(max(0.,2-2*best)),
                "cosine_lower":best,"cosine_upper":cub,"gamma_witness":witness,
                "evaluations":calls,"status":status,
                "numeric_scope":"analytic interval bounds evaluated in float64; no outward rounding"}

    def project_krylov(self,E):
        """Orthogonal projection onto span vec(A^k B), including current.

        This preserves every pairwise squared-error difference among the
        current and all continuous SR predictions on this SAME reference.
        E may have shape (states,channels) or (...,states,channels).
        """
        if self.vectors is None or self.reference_modes is None:
            raise ValueError("Construct with from_reference to project measurements")
        E=np.asarray(E,float);single=E.ndim==2
        shape=E.shape;batch=E.reshape((-1,)+self.reference_modes.shape)
        modes=np.einsum('ij,bjk->bik',self.vectors.conj().T,batch)
        out=np.zeros_like(modes)
        for group in self.groups:
            ref=self.reference_modes[group]
            weight=float(np.sum(abs(ref)**2))
            # Zero-energy eigenmodes are genuinely absent from this
            # reference. Tiny positive modes are not thresholded away.
            if weight==0:continue
            coeff=np.einsum('ik,bik->b',ref.conj(),modes[:,group])/weight
            out[:,group]=coeff[:,None,None]*ref[None]
        result=np.einsum('ij,bjk->bik',self.vectors,out)
        imaginary=float(np.linalg.norm(result.imag))
        if imaginary>1e-8*max(1.,float(np.linalg.norm(result.real))):
            raise FloatingPointError("Projection did not preserve real data")
        return result.real.reshape(shape)

    def krylov_norm(self,E):
        p=self.project_krylov(E)
        return np.linalg.norm(p,axis=(-2,-1))


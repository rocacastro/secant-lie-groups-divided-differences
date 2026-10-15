"""Algebraic, domain, counting, and local consistency tests."""
from __future__ import annotations
import unittest
from unittest.mock import patch
import mpmath as mp
import numpy as np
from scipy.linalg import expm
from threadpoolctl import threadpool_limits
from codes.lie_divdiff import norm, adapted_frame, reduced_divdiff, rank_one_correction, central_divdiff
from codes.so3_utils import hat, vee, exp, log
from codes.verify_semilocal import verify as verify_semilocal
from experiments.bessel import benchmark
from verify_reference import verify_metadata, verify_so3


class HighPrecisionTests(unittest.TestCase):
    def setUp(self):
        self.old_dps=mp.mp.dps
        mp.mp.dps=100
        self.eps=mp.mpf('1e-85')

    def tearDown(self):
        mp.mp.dps=self.old_dps

    def test_householder_all_dimensions_and_signs(self):
        for m in [1,2,3,4,8,45]:
            for sign in [1,-1]:
                v=mp.matrix([sign*mp.mpf(i+1)/(m+1) for i in range(m)])
                frame=adapted_frame(v)
                self.assertLess(norm(frame.T*frame-mp.eye(m)),self.eps)
                self.assertLess(norm(frame[:,0]-v/norm(v)),self.eps)

    def test_householder_equator_and_axis(self):
        for values in [[0,1,0],[1,0,0],[-1,0,0],[0,0,1]]:
            v=mp.matrix(values);Q=adapted_frame(v)
            self.assertLess(norm(Q.T*Q-mp.eye(3)),self.eps)
            self.assertLess(norm(Q[:,0]-v),self.eps)

    def test_householder_zero_rejected(self):
        with self.assertRaises(ValueError):
            adapted_frame(mp.zeros(3,1))

    def test_hat_vee(self):
        v=mp.matrix([mp.mpf('.1'),mp.mpf('-.2'),mp.mpf('.3')])
        self.assertEqual(list(vee(hat(v))),list(v))
        self.assertLess(abs(norm(hat(v))-mp.sqrt(2)*norm(v)),self.eps)

    def test_exp_log_round_trip(self):
        for values in [['0','0','0'],['.1','-.2','.3'],['1e-70','-2e-70','3e-70']]:
            v=mp.matrix(list(map(mp.mpf,values)));R=exp(v)
            self.assertLess(norm(log(R)-v),self.eps)
            self.assertLess(norm(R.T*R-mp.eye(3)),self.eps)

    def test_cut_locus_rejected(self):
        with self.assertRaises(ValueError):
            log(mp.diag([1,-1,-1]))

    def test_rank_one_correction_and_transverse_action(self):
        A=mp.matrix([[2,3,1],[0,2,4],[3,-1,1]])
        v=mp.matrix([1,2,-1]);df=mp.matrix([3,1,2])
        L=rank_one_correction(A,v,df)
        self.assertLess(norm(L*v-df),self.eps)
        Q=adapted_frame(v)
        self.assertLess(norm((L-A)*Q[:,1]),self.eps)
        self.assertLess(norm((L-A)*Q[:,2]),self.eps)

    def test_rank_one_invalid_arguments(self):
        with self.assertRaises(ValueError):
            rank_one_correction(mp.eye(3),mp.zeros(3,1),mp.ones(3,1))
        with self.assertRaises(ValueError):
            rank_one_correction(mp.eye(2),mp.ones(3,1),mp.ones(3,1))

    def dense_field(self):
        C=mp.matrix([[2,mp.mpf('.2'),mp.mpf('-.1')],
                     [mp.mpf('-.3'),mp.mpf('1.5'),mp.mpf('.25')],
                     [mp.mpf('.1'),mp.mpf('-.2'),mp.mpf('1.8')]])
        calls=[0]
        def field(X):
            calls[0]+=1
            g=vee(X-X.T)
            return C*g+mp.matrix([mp.mpf('.1')*mp.sin(g[i]) for i in range(3)])
        return field,calls,2*C+mp.mpf('.2')*mp.eye(3)

    def test_reduced_and_central_counts_and_identity(self):
        for builder,expected in [(reduced_divdiff,4),(central_divdiff,6)]:
            field,calls,_=self.dense_field()
            X=exp(mp.matrix(list(map(mp.mpf,['.02','-.03','.01']))))
            Y=X*exp(mp.matrix(list(map(mp.mpf,['.015','.006','-.004']))))
            fx,fy=field(X),field(Y);before=calls[0]
            L,v,Q,h=builder(field,X,Y,fx,fy,exp,log)
            self.assertEqual(calls[0]-before,expected)
            self.assertLess(norm(L*v-(fy-fx)),self.eps)
            self.assertGreater(h,0)

    def test_second_order_midpoint_consistency(self):
        for builder in [reduced_divdiff,central_divdiff]:
            field,_,J=self.dense_field()
            direction=mp.matrix([1,2,-1]);direction/=norm(direction)
            errors=[]
            for r in map(mp.mpf,['.01','.005','.0025']):
                X=exp(-r*direction/2);Y=exp(r*direction/2)
                L,*_=builder(field,X,Y,field(X),field(Y),exp,log)
                errors.append(norm(L-J))
            for i in range(2):
                self.assertGreater(errors[i]/errors[i+1],mp.mpf('3.999'))
                self.assertLess(errors[i]/errors[i+1],mp.mpf('4.001'))

    def test_reduced_additive_group_one_dimension(self):
        field=lambda x:mp.matrix([x[0]**2+1])
        X=mp.matrix([1]);Y=mp.matrix([2])
        L,v,Q,h=reduced_divdiff(field,X,Y,field(X),field(Y),lambda x:x,lambda x:x,
                               multiply=lambda x,y:x+y,inverse=lambda x:-x)
        self.assertEqual(L[0,0],3)
        self.assertEqual(v[0],1)

    def test_reduced_invalid_tau(self):
        field=lambda X:mp.ones(3,1)
        X=mp.eye(3);Y=exp(mp.matrix([1,0,0]))
        for tau in [0,-1]:
            with self.assertRaises(ValueError):
                reduced_divdiff(field,X,Y,field(X),field(Y),exp,log,tau=tau)

    def test_explicit_transverse_frame_dependence(self):
        e1,e2,e3=mp.eye(3)[:,0],mp.eye(3)[:,1],mp.eye(3)[:,2]
        h=mp.mpf('.05')
        field=lambda x:x+e1*x[0]**3
        def transverse(a,b):
            return (e3*e3.T+(field(h*a)-field(-h*a))*a.T/(2*h)
                    +(field(h*b)-field(-h*b))*b.T/(2*h))
        L1=transverse(e1,e2)
        L2=transverse((e1+e2)/mp.sqrt(2),(-e1+e2)/mp.sqrt(2))
        self.assertLess(norm(L1-L2-h*h/2*e1*e1.T),self.eps)

    def test_exact_semilocal_bounds(self):
        self.assertTrue(all(verify_semilocal()['checks'].values()))

    def test_reference_so3_precision_agreement(self):
        self.assertEqual(len(verify_so3()),2)

    def test_metadata(self):
        self.assertTrue(verify_metadata()['metadata_consistent'])


class BesselTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pool=threadpool_limits(limits=1)
        with np.load(benchmark.BASE,allow_pickle=False) as data:
            cls.field=benchmark.Field(data)

    @classmethod
    def tearDownClass(cls):
        cls.pool.restore_original_limits()

    def test_dataset_shape_and_prescribed_root(self):
        p=self.field
        self.assertEqual((p.n,p.m),(10,45))
        np.testing.assert_array_equal(p.g(p.Rstar),np.zeros(45))
        np.testing.assert_array_equal(p.A,-(p.Rstar-p.Rstar.T))

    def test_geometric_linear_map(self):
        p=self.field
        v=np.linspace(-.1,.2,p.m)
        V=p.Xstart@p.hat(v)
        gamma=V[p.i,p.j]-V[p.j,p.i]
        np.testing.assert_allclose(p.G(p.Xstart)@v,gamma,atol=1e-15,rtol=1e-14)
        np.testing.assert_allclose(p.G(np.eye(p.n)),2*np.eye(p.m),atol=1e-15)

    def test_cached_secant_no_evaluations(self):
        p=self.field;v=.01*np.ones(p.m)/np.sqrt(p.m)
        X=p.Xstart@expm(-p.hat(v));Y=p.Xstart
        fx,wx,px,_=p.eval(X);fy,wy,py,_=p.eval(Y)
        before=p.F_calls
        with patch.object(benchmark,'jv',wraps=benchmark.jv) as calls:
            L=p.alg_secant(X,Y,fx,wx,px,fy,wy,py,v)
        self.assertEqual(calls.call_count,0)
        self.assertEqual(p.F_calls,before)
        self.assertLess(np.linalg.norm(L@v-(fy-fx)),1e-13)

    def test_scalar_fallback_is_counted_and_finite(self):
        p=self.field;v=1e-14*np.ones(p.m)/np.sqrt(p.m)
        X=p.Xstart;Y=X@expm(p.hat(v))
        fx,wx,px,_=p.eval(X);fy,wy,py,_=p.eval(Y)
        before=p.fallbacks
        with patch.object(benchmark,'jv',wraps=benchmark.jv) as calls:
            L=p.alg_secant(X,Y,fx,wx,px,fy,wy,py,v)
        self.assertEqual(p.fallbacks-before,45)
        self.assertEqual(calls.call_count,2)
        self.assertTrue(np.all(np.isfinite(L)))
        self.assertLess(np.linalg.norm(L@v-(fy-fx)),1e-25)

    def test_analytic_jacobian_directional_check(self):
        self.assertLess(benchmark.validate_jacobian(self.field),1e-7)

    def test_increment_reuse(self):
        p=self.field;s=.01*np.linspace(-1,1,p.m)
        Y=p.Xstart@expm(p.hat(s))
        from scipy.linalg import logm
        V=logm(p.Xstart.T@Y)
        self.assertLess(np.max(np.abs(np.imag(V))),1e-12)
        self.assertLess(np.linalg.norm(V.real[p.i,p.j]-s),1e-12)

    def test_reference_timing_statistics(self):
        from experiments.bessel.verify import verify
        result=verify(False)
        self.assertEqual(result['secant_faster_blocks'],50)
        self.assertAlmostEqual(result['median_paired_ratio'],0.8807853631014997,14)


if __name__=='__main__':
    unittest.main()

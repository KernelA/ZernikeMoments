// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
/*

                          3D Zernike Moments
    Copyright (C) 2003 by Computer Graphics Group, University of Bonn
           http://www.cg.cs.uni-bonn.de/project-pages/3dsearch/

Code by Marcin Novotni:     marcin@cs.uni-bonn.de

for more information, see the paper:

@inproceedings{novotni-2003-3d,
    author = {M. Novotni and R. Klein},
    title = {3{D} {Z}ernike Descriptors for Content Based Shape Retrieval},
    booktitle = {The 8th ACM Symposium on Solid Modeling and Applications},
    pages = {216--225},
    year = {2003},
    month = {June},
    institution = {Universit\"{a}t Bonn},
    conference = {The 8th ACM Symposium on Solid Modeling and Applications, June 16-20, Seattle, WA}
}
 *---------------------------------------------------------------------------*
 *                                                                           *
 *                                License                                    *
 *                                                                           *
 *  This library is free software; you can redistribute it and/or modify it  *
 *  under the terms of the GNU Library General Public License as published   *
 *  by the Free Software Foundation, version 2.                              *
 *                                                                           *
 *  This library is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        *
 *  Library General Public License for more details.                         *
 *                                                                           *
 *  You should have received a copy of the GNU Library General Public        *
 *  License along with this library; if not, write to the Free Software      *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                *
 *                                                                           *
\*===========================================================================*/

#pragma once

// ----- local program includes -----
#include "ScaledGeometricMoments.hpp"

/**
 * Struct representing a complex coefficient of a moment
 * of order (p_,q_,r_)
 */
template<class T>
struct ComplexCoeff
{
    typedef     std::complex<T>     ValueT;

    /**
 * Constructor with scalar args
 */
    ComplexCoeff(int _p, int _q, int _r, const ValueT & _value) :
        p_(_p), q_(_q), r_(_r), value_(_value)
    {
    }

    /**
 * Copy constructor
 */
    ComplexCoeff(const ComplexCoeff<T> & _cc) : p_(_cc.p_), q_(_cc.q_), r_(_cc.r_), value_(_cc.value_)
    {
    }

    ComplexCoeff() : p_(0), q_(0), r_(0)
    {
    }

    int                 p_, q_, r_;
    ValueT     value_;
};

/**
 * Class representing the Zernike moments
 */
template<class InputVoxelIterator, class MomentT>
class ZernikeMoments
{
public:
    // ---- public typedefs ----
    typedef MomentT             T;
    typedef vector<T>           T1D;        // vector of scalar type
    typedef vector<T1D>         T2D;        // 2D array of scalar type
    typedef vector<T2D>         T3D;        // 3D array of scalar type
    typedef vector<T3D>         T4D;        // 3D array of scalar type

    typedef std::complex<T>                      ComplexT;       // complex type
    typedef vector<vector<vector<ComplexT> > >   ComplexT3D;     // 3D array of complex type

    typedef ComplexCoeff<T>                      ComplexCoeffT;
    typedef vector<vector<vector<vector<ComplexCoeffT> > > >    ComplexCoeffT4D;

public:
    // ---- public member functions ----
    ZernikeMoments(int _order, ScaledGeometricalMoments<InputVoxelIterator, MomentT> & _gm)
    {
        Init(_order, _gm);
    }

    ZernikeMoments() :
        order_(0)
    {
    }

    /* Todo: Init() doesn't depend on the geomentrical moments.
     * Removing the  parameter and moving it to Compute() will make the class more efficient to reuse to
     * compute ZM for several objects at once.

     * Computes all coefficients that are input data independent
     */
    void Init(int _order, ScaledGeometricalMoments<InputVoxelIterator, MomentT> & _gm)
    {
        static_assert(std::is_floating_point<T>::value, "MomentT must be float, double or long double");
        /* Todo: Init() doesn't depend on the geomentrical moments.
         * Removing the  parameter and moving it to Compute() will make the class more efficient to reuse to
         * compute ZM for several objects at once.
         */
        gm_ = _gm;
        order_ = _order;

        ComputeCs();
        ComputeQs();
        ComputeGCoefficients();
    }

    /**
 * Computes the Zernike moments. This computation is data dependent
 * and has to be performed for each new object and/or transformation.
 */
    void Compute()
    {
        // geometrical moments have to be computed first
        if (!order_)
        {
            throw std::runtime_error("ZernikeMoments<InputVoxelIterator,MomentT>::ComputeZernikeMoments (): attempting to \
                     compute Zernike moments without setting valid geometrical \
                     moments first.");
        }

        /*
         indexing:
           n goes 0..order_
           l goes 0..n so that n-l is even
           m goes -l..l
        */

        T nullMoment;

        constexpr T three_quarters_div_pi = boost::math::constants::three_quarters<T>() * 1 / boost::math::constants::pi<T>();

        zernikeMoments_.resize(order_ + 1);
        for (int n = 0; n <= order_; ++n)
        {
            zernikeMoments_[n].resize(n / 2 + 1);

            int l0 = n % 2, li = 0;
            for (int l = l0; l <= n; ++li, l += 2)
            {
                zernikeMoments_[n][li].resize(l + 1);
                for (int m = 0; m <= l; ++m)
                {
                    // Zernike moment of according indices [nlm]
                    ComplexT zm(static_cast<T>(0), static_cast<T>(0));

                    int nCoeffs = (int)gCoeffs_[n][li][m].size();
                    for (int i = 0; i < nCoeffs; ++i)
                    {
                        ComplexCoeffT cc = gCoeffs_[n][li][m][i];
                        //T scale = gm_.GetScale ();
                        //T fact = std::pow (scale, cc.p_+cc.q_+cc.r_+3);

                        //zm +=  std::conj (cc.value_) * gm_.GetMoment(cc.p_, cc.q_, cc.r_) * fact;
                        zm += std::conj(cc.value_) * gm_.GetMoment(cc.p_, cc.q_, cc.r_);
                    }

                    zm *= three_quarters_div_pi;

                    if (n == 0 && l == 0 && m == 0)
                    {
                        // FixMe Unused variable! What is it for?
                        nullMoment = zm.real();
                    }
                    zernikeMoments_[n][li][m] = zm;

                    std::cout << "zernike moment[nlm]: " << n << "\t" << l << "\t" << m << "\t" << zernikeMoments_[n][li][m] << "\n";
                }
            }
        }
    }

    inline ComplexT GetMoment(int _n, int _l, int _m)
    {
        if (_m >= 0)
        {
            return zernikeMoments_[_n][_l / 2][_m];
        }
        else
        {
            T sign{};

            if (_m % 2)
            {
                sign = static_cast<T>(-1);
            }
            else
            {
                sign = static_cast<T>(1);
            }
            return sign * std::conj(zernikeMoments_[_n][_l / 2][abs(_m)]);
        }
    }

    // ---- debug functions/arguments ----
    /**
 * The function previously encoded as complex valued Zernike
 * moments, is reconstructed. _grid is the output grid containing
 * the reconstructed function.
 */
    void Reconstruct(ComplexT3D & _grid,                // grid containing the reconstructed function
        T             _xCOG,                // center of gravity
        T             _yCOG,
        T             _zCOG,
        T             _scale,               // scaling factor to map into unit ball
        int           _minN = 0,            // min value for n freq index
        int           _maxN = 100,          // min value for n freq index
        int           _minL = 0,            // min value for l freq index
        int           _maxL = 100)          // max value for l freq index
    {
        int dimX = _grid.size();
        int dimY = _grid[0].size();
        int dimZ = _grid[0][0].size();

        //scaling
        T scale = _scale;

        //translation
        T vx = _xCOG;
        T vy = _yCOG;
        T vz = _zCOG;

        T point[3];

        if (_maxN == 100)
        {
            _maxN = order_;
        }

        for (int x = 0; x < dimX; ++x)
        {
            std::cout << x << ". layer being processed\n";
            for (int y = 0; y < dimY; ++y)
            {
                for (int z = 0; z < dimZ; ++z)
                {
                    // the origin is in the middle of the grid, all voxels are
                    // projected into the unit ball
                    point[0] = (static_cast<T>(x) - vx) * scale;
                    point[1] = (static_cast<T>(y) - vy) * scale;
                    point[2] = (static_cast<T>(z) - vz) * scale;

                    if (point[0] * point[0] + point[1] * point[1] + point[2] * point[2] > 1.0)
                    {
                        continue;
                    }

                    // function value at point
                    ComplexT fVal{ 0, 0 };

                    for (int n = _minN; n <= _maxN; ++n)
                        //for (int n=14; n<=14; ++n)
                    {
                        int maxK = n / 2;
                        for (int k = 0; k <= maxK; ++k)
                        {
                            for (int nu = 0; nu <= k; ++nu)
                            {
                                int l = n - 2 * k;
                                // check whether l is within bounds
                                if (l < _minL || l > _maxL)
                                {
                                    continue;
                                }

                                for (int m = -l; m <= l; ++m)
                                {
                                    // zernike polynomial evaluated at point
                                    ComplexT zp(0, 0);

                                    int absM = std::abs(m);

                                    int nCoeffs = gCoeffs_[n][l / 2][absM].size();
                                    for (int i = 0; i < nCoeffs; ++i)
                                    {
                                        ComplexCoeffT cc = gCoeffs_[n][l / 2][absM][i];
                                        ComplexT cvalue = cc.value_;

                                        // conjugate if m negative
                                        if (m < 0)
                                        {
                                            cvalue = std::conj(cvalue);

                                            // take care of the sign
                                            if (m % 2)
                                            {
                                                cvalue *= (T)(-1);
                                            }
                                        }

                                        zp += cvalue *
                                            std::pow(point[0], (T)cc.p_) *
                                            std::pow(point[1], (T)cc.q_) *
                                            std::pow(point[2], (T)cc.r_);
                                    }

                                    fVal += zp * GetMoment(n, l, m);
                                }
                            }
                        }
                    }
                    _grid[x][y][z] = fVal;
                }
            }
        }

        //NormalizeGridValues (_grid);
    }

    void NormalizeGridValues(ComplexT3D & _grid)
    {
        int xD = _grid.size();
        int yD = _grid[0].size();
        int zD = _grid[0][0].size();

        T max{ 0 };
        for (int k = 0; k < zD; ++k)
        {
            for (int j = 0; j < yD; ++j)
            {
                for (int i = 0; i < xD; ++i)
                {
                    if (_grid[i][j][k].real() > max)
                    {
                        max = _grid[i][j][k].real();
                    }
                }
            }
        }

        std::cout << "\nMaximal value in grid: " << max << "\n";

        T invMax = static_cast<T>(1) / max;

        for (int k = 0; k < zD; ++k)
        {
            for (int j = 0; j < yD; ++j)
            {
                for (int i = 0; i < xD; ++i)
                {
                    _grid[i][j][k] *= invMax;
                }
            }
        }
    }

    void CheckOrthonormality(int _n1, int _l1, int _m1, int _n2, int _l2, int _m2)
    {
        int li1 = _l1 / 2;
        int li2 = _l2 / 2;
        int dim = 64;

        // the total sum of the scalar product
        ComplexT sum(static_cast<T>(0), static_cast<T>(0));

        int nCoeffs1 = (int)gCoeffs_[_n1][li1][_m1].size();
        int nCoeffs2 = (int)gCoeffs_[_n2][li2][_m2].size();

        for (int i = 0; i < nCoeffs1; ++i)
        {
            ComplexCoeffT cc1 = gCoeffs_[_n1][li1][_m1][i];
            for (int j = 0; j < nCoeffs2; ++j)
            {
                ComplexCoeffT cc2 = gCoeffs_[_n2][li2][_m2][j];

                T temp{ 0 };

                int p = cc1.p_ + cc2.p_;
                int q = cc1.q_ + cc2.q_;
                int r = cc1.r_ + cc2.r_;

                sum += cc1.value_ *
                    std::conj(cc2.value_) *
                    EvalMonomialIntegral(p, q, r, dim);
            }
        }

        std::cout << "\nInner product of [" << _n1 << "," << _l1 << "," << _m1 << "]";
        std::cout << " and [" << _n2 << "," << _l2 << "," << _m2 << "]: ";
        std::cout << sum << "\n\n";
    }

private:
    // ---- private member functions ----

/**
 * Computes all the normalizing factors $c_l^m$ for harmonic polynomials e
 */
    void ComputeCs()
    {
        using namespace boost::math;

        /*
         indexing:
           l goes from 0 to n
           m goes from -l to l, in fact from 0 to l, since c(l,-m) = c (l,m)
        */

        cs_.resize(order_ + 1);

        for (size_t l = 0; l <= order_; ++l)
        {
            cs_[l].resize(l + 1);
            for (size_t m = 0; m <= l; ++m)
            {
                /*         T n_sqrt = ((T)2 * l + (T)1) *
                             Factorial<T>::Get(l + 1, l + m);
                         T d_sqrt = Factorial<T>::Get(l - m + 1, l);*/

                T n_sqrt = static_cast<T>((2 * l + 1) * rising_factorial(l + 1, m));
                T d_sqrt = static_cast<T>(rising_factorial(l - m + 1, m));

                cs_[l][m] = std::sqrt(n_sqrt / d_sqrt);
            }
        }
    }

    /**
 * Computes all coefficients q for orthonormalization of radial polynomials
 * in Zernike polynomials.
 */
    void ComputeQs()
    {
        using namespace boost::math;

        /*
         indexing:
           n goes 0..order_
           l goes 0..n, so that n-l is even
           mu goes 0..(n-l)/2
        */

        qs_.resize(order_ + 1);            // there is order_ + 1 n's

        for (size_t n = 0; n <= order_; ++n)
        {
            qs_[n].resize(n / 2 + 1);      // there is floor(n/2) + 1 l's

            size_t l0 = n % 2;
            for (size_t l = l0; l <= n; l += 2)
            {
                size_t k = (n - l) / 2;

                qs_[n][l / 2].resize(k + 1);   // there is k+1 mu's

                for (size_t mu = 0; mu <= k; ++mu)
                {
                    T nom = binomial_coefficient<T>(2 * k, k) * // nominator of straight part
                        binomial_coefficient<T>(k, mu) * binomial_coefficient<T>(2 * (k + l + mu) + 1, 2 * k);

                    if ((k + mu) % 2)
                    {
                        nom *= static_cast<T>(-1);
                    }

                    T den = std::pow(static_cast<T>(2), static_cast<T>(2 * k)) *     // denominator of straight part
                        binomial_coefficient<T>(k + l + mu, k);

                    T n_sqrt = static_cast<T>(2 * l + 4 * k + 3);      // nominator of sqrt part
                    T d_sqrt = static_cast<T>(3);                        // denominator of sqrt part

                    qs_[n][l / 2][mu] = nom / den * sqrt(n_sqrt / d_sqrt);
                }
            }
        }
    }

    /**
 * Computes the coefficients of geometrical moments in linear combinations
 * yielding the Zernike moments for each applicable [n,l,m] for n<=order_.
 * For each such combination the coefficients are stored with according
 * geom. moment order (see ComplexCoeff).
 */
    void ComputeGCoefficients()
    {
        using namespace boost::math;

        //DD
        size_t countCoeffs = 0;
        //DD
        gCoeffs_.resize(order_ + 1);

        for (size_t n = 0; n <= order_; ++n)
        {
            gCoeffs_[n].resize(n / 2 + 1);
            size_t li = 0, l0 = n % 2;
            for (size_t l = l0; l <= n; ++li, l += 2)
            {
                gCoeffs_[n][li].resize(l + 1);
                for (size_t m = 0; m <= l; ++m)
                {
                    T w = cs_[l][m] / std::pow(static_cast<T>(2), static_cast<T>(m));

                    size_t k = (n - l) / 2;
                    for (size_t nu = 0; nu <= k; ++nu)
                    {
                        T w_Nu = w * qs_[n][li][nu];
                        for (size_t alpha = 0; alpha <= nu; ++alpha)
                        {
                            T w_NuA = w_Nu * binomial_coefficient<T>(nu, alpha);
                            for (size_t beta = 0; beta <= nu - alpha; ++beta)
                            {
                                T w_NuAB = w_NuA * binomial_coefficient<T>(nu - alpha, beta);
                                for (size_t p = 0; p <= m; ++p)
                                {
                                    T w_NuABP = w_NuAB * binomial_coefficient<T>(m, p);
                                    for (size_t mu = 0; mu <= (l - m) / 2; ++mu)
                                    {
                                        T w_NuABPMu = w_NuABP *
                                            binomial_coefficient<T>(l, mu) *
                                            binomial_coefficient<T>(l - mu, m + mu) /
                                            static_cast<T>(std::pow(2.0, (double)(2 * mu)));
                                        for (size_t q = 0; q <= mu; ++q)
                                        {
                                            // the absolute value of the coefficient
                                            T w_NuABPMuQ = w_NuABPMu * binomial_coefficient<T>(mu, q);

                                            // the sign
                                            if ((m - p + mu) % 2)
                                            {
                                                w_NuABPMuQ *= static_cast<T>(-1);
                                            }

                                            // * i^p
                                            size_t rest = p % 4;
                                            ComplexT c;
                                            switch (rest)
                                            {
                                                case 0: c = ComplexT(w_NuABPMuQ, static_cast <T>(0)); break;
                                                case 1: c = ComplexT(static_cast<T>(0), w_NuABPMuQ); break;
                                                case 2: c = ComplexT(static_cast<T>(-1) * w_NuABPMuQ, static_cast<T>(0)); break;
                                                case 3: c = ComplexT(static_cast<T>(0), static_cast<T>(-1) * w_NuABPMuQ); break;
                                            }

                                            // determination of the order of according moment
                                            int z_i = l - m + 2 * (nu - alpha - beta - mu);
                                            int y_i = 2 * (mu - q + beta) + m - p;
                                            int x_i = 2 * q + p + 2 * alpha;
                                            //DD
                                                                                    //std::cout << x_i << " " << y_i << " " << z_i;
                                                                                    //std::cout << "\t" << n << " " << l << " " << m;
                                                                                    //std::cout << "\t" << c.real () << " " << c.imag () << std::endl;
                                            //DD
                                            ComplexCoeffT cc(x_i, y_i, z_i, c);
                                            gCoeffs_[n][li][m].push_back(cc);
                                            //DD
                                            countCoeffs++;
                                            //DD
                                        } // q
                                    } // mu
                                } // p
                            } // beta
                        } // alpha
                    } // nu
                } // m
            } // l
        } // n
    //DD
        //std::cout << countCoeffs << std::endl;
        //DD
    }

    // ---- private attributes -----
    ComplexCoeffT4D     gCoeffs_;           // coefficients of the geometric moments
    ComplexT3D          zernikeMoments_;    // nomen est omen
    T3D                 qs_;                // q coefficients (radial polynomial normalization)
    T2D                 cs_;                // c coefficients (harmonic polynomial normalization)

    ScaledGeometricalMoments<InputVoxelIterator, MomentT> gm_;
    int                 order_;             // := max{n} according to indexing of Zernike polynomials

    // ---- debug functions/arguments ----
    void PrintGrid(ComplexT3D & _grid)
    {
        int xD = _grid.size();
        int yD = _grid[0].size();
        int zD = _grid[0][0].size();

        std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);

        T max = static_cast<T>(0);

        for (int k = 0; k < zD; ++k)
        {
            for (int j = 0; j < yD; ++j)
            {
                for (int i = 0; i < xD; ++i)
                {
                    if (fabs(_grid[i][j][k].real()) > max)
                    {
                        max = fabs(_grid[i][j][k].real());
                    }
                }
            }
        }

        for (int k = 0; k < zD; ++k)
        {
            std::cout << k << ". layer:\n";
            for (int j = 0; j < yD; ++j)
            {
                for (int i = 0; i < xD; ++i)
                {
                    //std::cout << _grid[i][j][k].real () / max << "\t";
                    std::cout << _grid[i][j][k].real() << "\t";
                }
                std::cout << "\n";
            }
            std::cout << "\n";
        }
        std::cout.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
    }

    /**
 * Evaluates the integral of a monomial x^p*y^q*z^r within the unit sphere
 * Attention : a very stupid implementation, thus it's accordingly very slow
 */
    T EvalMonomialIntegral(int _p, int _q, int _r, int _dim)
    {
        T radius = static_cast<T>(_dim - 1) / static_cast <T>(2);
        T scale = std::pow(static_cast <T>(1) / radius, 3);
        T center = static_cast<T>(_dim - 1) / static_cast <T>(2);

        T result = static_cast<T>(0);
        T point[3];

        constexpr T three_quarters_div_pi = boost::math::constants::three_quarters<T>() * 1 / boost::math::constants::pi<T>();

        for (int x = 0; x < _dim; ++x)
        {
            point[0] = (static_cast<T>(x) - center) / radius;
            for (int y = 0; y < _dim; ++y)
            {
                point[1] = (static_cast<T>(y) - center) / radius;
                for (int z = 0; z < _dim; ++z)
                {
                    point[2] = (static_cast<T>(z) - center) / radius;

                    if (point[0] * point[0] + point[1] * point[1] + point[2] * point[2] > (T)1)
                    {
                        continue;
                    }

                    result += std::pow(point[0], (T)_p) *
                        std::pow(point[1], (T)_q) *
                        std::pow(point[2], (T)_r);
                }
            }
        }

        result *= three_quarters_div_pi * scale;
        return result;
    }
};

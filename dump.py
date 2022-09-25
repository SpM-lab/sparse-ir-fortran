import argparse
import sparse_ir
import xprec # Make sure xprec is install for good accuracy
import numpy as np

def run():
    parser = argparse.ArgumentParser(
        prog='irdump',
        description='Dump sampling points of IR',
    )
    parser.add_argument('lambda_', action='store', default=None,
                        type=float, help='beta * wmax')
    parser.add_argument('eps', action='store', default=None,
                        type=float, help='eps')
    parser.add_argument('output', action='store', default=None,
                        type=str, help='output filename')
    args = parser.parse_args()

    print(args.lambda_)
    print(args.eps)

    beta = 1.0
    lambda_ = args.lambda_
    wmax_ = lambda_/beta
    eps = args.eps

    basis = sparse_ir.FiniteTempBasisSet(1.0, lambda_, eps)
    basis_f = basis.basis_f
    basis_b = basis.basis_b
    s = basis_f.s / np.sqrt(0.5 * lambda_)

    # For the logistic kernel, the fermionic and bosonic basis functions are equivalent.
    with open(args.output, "w") as f:
        print(f"#version 1", file=f)
        print(f"#lambda {lambda_}", file=f)
        print(f"#eps {eps}", file=f)
        print(f"#singular_values", file=f)
        print(basis_f.size, file=f)
        for s in s:
            print('{:.16e}'.format(s), file=f)

        print(f"# sampling times", file=f)
        #   Add tau = 0, beta to sampling points (for Matsubara summation)
        xs = np.unique(np.hstack([-1, 2*basis_f.default_tau_sampling_points()/beta - 1, 1]))
        times = 0.5 * (xs + 1) * beta
        print(xs.size, file=f)
        for x in xs:
            print('{:.16e}'.format(x), file=f)

        print(f"# u", file=f)
        uval = np.sqrt(0.5 * beta) * basis_f.u(times)
        for l in range(basis_f.size):
            for t in range(times.size):
                print('{:.16e}'.format(uval[l, t]), file=f)

        for stat in ["F", "B"]:
            print(f"# sampling frequencies {stat}", file=f)
            basis = {"F": basis_f, "B": basis_b}[stat]
            freqs = basis.default_matsubara_sampling_points()
            print(freqs.size, file=f)
            for x in freqs:
                print(x, file=f)

            print(f"# uhat {stat}", file=f)
            uhatval = (1/np.sqrt(beta)) * basis.uhat(freqs)
            for l in range(basis.size):
                for i in range(freqs.size):
                    print('{:.16e} {:.16e}'.
                        format(uhatval[l, i].real, uhatval[l, i].imag), file=f)
        
        print(f"# sampling poles", file=f)
        #   Add w_p = -wmax, +wmax to sampling poles
        omega = np.unique(np.hstack([-1, basis_f.default_omega_sampling_points()/wmax_, 1]))
        print(omega.size, file=f)
        for x in omega:
            print('{:.16e}'.format(x), file=f)

        print(f"# v", file=f)
        vval = np.sqrt(wmax_) * basis_f.v(omega)
        for l in range(basis_f.size):
            for i in range(omega.size):
                print('{:.16e}'.format(vval[l, i]), file=f)

if __name__ == '__main__':
    run()

import argparse
import sparse_ir

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

    lambda_ = args.lambda_
    eps = args.eps

    kernel = sparse_ir.LogisticKernel(lambda_)
    sve_result = sparse_ir.compute_sve(kernel, eps)
    basis_f = sparse_ir.IRBasis("F", lambda_, eps, kernel=kernel, sve_result=sve_result)
    basis_b = sparse_ir.IRBasis("B", lambda_, eps, kernel=kernel, sve_result=sve_result)

    # For the logistic kernel, the fermionic and bosonic basis functions are equivalent.
    with open(args.output, "w") as f:
        print(f"#version 1", file=f)
        print(f"#lambda {lambda_}", file=f)
        print(f"#eps {eps}", file=f)
        print(f"#singular_values", file=f)
        print(basis_f.size, file=f)
        for s in basis_f.s:
            print('{:.16e}'.format(s), file=f)

        print(f"# sampling times", file=f)
        times = basis_f.default_tau_sampling_points()
        print(times.size, file=f)
        for x in times:
            print('{:.16e}'.format(x), file=f)

        print(f"# u", file=f)
        uval = basis_f.u(times)
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
            uhatval = basis.uhat(freqs)
            for l in range(basis.size):
                for i in range(freqs.size):
                    print('{:.16e} {:.16e}'.
                        format(uhatval[l, i].real, uhatval[l, i].imag), file=f)

if __name__ == '__main__':
    run()

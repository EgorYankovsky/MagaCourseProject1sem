public class LU_LOS(int maxIter = 100_000, double eps = 1E-15) : ISolver
{
    private int _maxIter = maxIter;
    private double _eps = eps;

    private static void PartitialLU(GlobalMatrix A)
    {
        for (int i = 0; i < A.Size; i++)
        {
            for (int j = A._ig[i]; j < A._ig[i + 1]; j++)
            {
                int jCol = A._jg[j];
                int jk = A._ig[jCol];
                int k = A._ig[i];
                int sdvig = A._jg[A._ig[i]] - A._jg[A._ig[jCol]];
                if (sdvig > 0) jk += sdvig;
                else k -= sdvig;
                double sumL = 0.0;
                double sumU = 0.0;
                for (; k < j && jk < A._ig[jCol + 1]; k++, jk++)
                {
                    sumL += A._al[k] * A._au[jk];
                    sumU += A._au[k] * A._al[jk];
                }
                A._al[j] -= sumL;
                A._au[j] -= sumU;
                A._au[j] /= A._diag[jCol];
            }
            double sumD = 0.0;
            for (int j = A._ig[i]; j < A._ig[i + 1]; j++)
                sumD += A._al[j] * A._au[j];
            A._diag[i] -= sumD;
        }
    }

    private static GlobalVector Forward(GlobalMatrix Matrix, GlobalVector b)
    {
        var result = new GlobalVector(b);
        for (int i = 0; i < Matrix.Size; i++)
        {
            for (int j = Matrix._ig[i]; j < Matrix._ig[i + 1]; j++)
                result[i] -= Matrix._al[j] * result[Matrix._jg[j]];
            result[i] /= Matrix._diag[i];
        }
        return result;
    }

    private static GlobalVector Backward(GlobalMatrix A, GlobalVector b)
    {
        var result = new GlobalVector(b);
        for (int i = A.Size - 1; i >= 0; i--)
            for (int j = A._ig[i + 1] - 1; j >= A._ig[i]; j--)
                result[A._jg[j]] -= A._au[j] * result[i];
        return result;
    }

    public (GlobalVector, GlobalVector) Solve(GlobalMatrix A, GlobalVector b)
    {
        GlobalVector x = new(b.Size);
        GlobalVector x_ = new(b.Size);
        GlobalMatrix LU = new(A);
        PartitialLU(LU);
        GlobalVector r = Forward(LU, b - A * x);
        var r0 = new GlobalVector(r);
        GlobalVector z = Backward(LU, r);
        GlobalVector p = Forward(LU, A * z);
        GlobalVector tmp = new(b.Size);
        GlobalVector r_ = new(b.Size);
        GlobalVector z_ = new(b.Size);
        GlobalVector p_ = new(b.Size);
        double alph = 0.0D;
        double beta = 0.0D;
        int iter = 0;
        do
        {
            x_ = x;
            z_ = z;
            r_ = r;
            p_ = p;
            alph = (p_ * r_) / (p_ * p_);
            x = x_ + alph * z_;
            r = r_ - alph * p_;
            tmp = Forward(LU, A * Backward(LU, r));
            beta = -1.0D * (p_ * tmp) / (p_ * p_);
            z = Backward(LU, r) + beta * z_;
            p = tmp + beta * p_;
            iter++;
            if (iter % 10 == 0)
                Console.WriteLine($"{(r.Norma() * r.Norma()) / (r0.Norma() * r0.Norma()):E15}");
        } while (iter < _maxIter && (r.Norma() * r.Norma()) / (r0.Norma() * r0.Norma()) >= _eps * _eps);
        Console.WriteLine(
        $@"Computing finished!
Total iterations: {iter}
Relative residuality: {r.Norma() / b.Norma():E15}");
        return (x, x_);
    }
}
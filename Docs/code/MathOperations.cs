using System.Threading.Tasks;

namespace MathObjects;

public static class MathOperations
{
    public static GlobalVector Diff(GlobalVector gv1, GlobalVector gv2)
    {
        if (gv1.Size != gv2.Size)
            throw new Exception("Найти разность векторов не возможно, т.к. они имеют разные размеры.");
        
        GlobalVector gv = new(gv1.Size);
        for (int i = 0; i < gv.Size; i++)
            gv[i] = gv1[i] - gv2[i];

        return gv;
    }

    public static GlobalVector DiffPar(GlobalVector gv1, GlobalVector gv2)
    {
        if (gv1.Size != gv2.Size)
            throw new Exception("Найти разность векторов не возможно, т.к. они имеют разные размеры.");
        GlobalVector gv = new(gv1.Size);
        Parallel.For(0, gv1.Size, i => {
            gv[i] = gv1[i] - gv2[i];
        });
        return gv;
    }

    public static GlobalVector Sum(GlobalVector gv1, GlobalVector gv2)
    {
        if (gv1.Size != gv2.Size)
            throw new Exception("Найти разность векторов не возможно, т.к. они имеют разные размеры.");
        
        GlobalVector gv = new(gv1.Size);
        for (int i = 0; i < gv.Size; i++)
            gv[i] = gv1[i] + gv2[i];

        return gv;
    }

    public static GlobalVector SumPar(GlobalVector gv1, GlobalVector gv2)
    {
        if (gv1.Size != gv2.Size)
            throw new Exception("Найти разность векторов не возможно, т.к. они имеют разные размеры.");
        
        GlobalVector gv = new(gv1.Size);
        Parallel.For(0, gv1.Size, i => {
            gv[i] = gv1[i] + gv2[i];
        });
        return gv;
    }

    public static GlobalVector Multiply(double a, GlobalVector gv)
    {
        GlobalVector _gv = new(gv.Size);
        for (int i = 0; i < gv.Size; i++)
            _gv[i] = a * gv[i];
        return _gv;
    }

    public static GlobalVector MultiplyPar(double a, GlobalVector gv)
    {
        GlobalVector _gv = new(gv.Size);
        Parallel.For(0, gv.Size, i => {
            _gv[i] = a * gv[i];
        });
        return _gv;
    }

    public static double Multiply(GlobalVector gv1, GlobalVector gv2)
    {
        if (gv1.Size != gv2.Size) throw new Exception("Невозможно найти скалярное умножение векторов, из-за разности в размерах.");
        double ans = 0.0D;
        for (int i = 0; i < gv1.Size; i++)
            ans += gv1[i] * gv2[i];
        return ans;
    }

    public static double MultiplyPar(GlobalVector gv1, GlobalVector gv2)
    {
        if (gv1.Size != gv2.Size) throw new Exception("Невозможно найти скалярное умножение векторов, из-за разности в размерах.");
        double ans = 0.0D;
        Parallel.For(0, gv1.Size, i => {
            ans = gv1[i] * gv2[i];
        });
        return ans;
    }


    public static GlobalVector Multiply(GlobalMatrix _gm, GlobalVector _gv)
    {
        if (_gm.Size != _gv.Size)
            throw new Exception("Невозможно перемножить матрицу на вектор.");
        GlobalVector ans = new(_gv.Size);

        for (int i = 0; i < _gv.Size; i++)
        {
            for (int j = 0; j < _gm._ig[i + 1] - _gm._ig[i]; j++)
            {
                ans[i] += _gm._al[_gm._ig[i] + j] * _gv[_gm._jg[_gm._ig[i] + j]];
                ans[_gm._jg[_gm._ig[i] + j]] += _gm._au[_gm._ig[i] + j] * _gv[i];
            }
            ans[i] += _gm._diag[i] * _gv[i];
        }
        return ans;
    }

    public static GlobalVector CustomMultiply(GlobalMatrix _gm, GlobalVector _gv)
    {
        if (_gm.Size != _gv.Size)
            throw new Exception("Невозможно перемножить матрицу на вектор.");
        GlobalVector ans = new(_gv.Size);

        for (int i = 0; i < _gv.Size; i++)
        {
            for (int j = _gm._ig[i]; j < _gm._ig[i + 1]; j++)
            {
                ans[i] += _gm._al[j] * _gv[_gm._jg[j]];
                ans[_gm._jg[j]] += _gm._au[j] * _gv[i];
            }
            ans[i] += _gm._diag[i] * _gv[i];
        }
        return ans;
    }


    public static GlobalVector MultiplyPar(GlobalMatrix _gm, GlobalVector _gv)
    {
        if (_gm.Size != _gv.Size)
            throw new Exception("Невозможно перемножить матрицу на вектор.");
        GlobalVector ans = new(_gv.Size);

        for (int i = 0; i < _gv.Size; i++)
        {
            Parallel.For(0, _gm._ig[i + 1] - _gm._ig[i], j => {
                ans[i] += _gm._al[_gm._ig[i] + j] * _gv[_gm._jg[_gm._ig[i] + j]];
            });
            Parallel.For(0, _gm._ig[i + 1] - _gm._ig[i], j => {
                ans[_gm._jg[_gm._ig[i] + j]] += _gm._au[_gm._ig[i] + j] * _gv[i];
            });
            ans[i] += _gm._diag[i] * _gv[i];
        }
        return ans;
    }


    public static GlobalMatrix Multiply(double a, GlobalMatrix gm)
    {
        var ans = new GlobalMatrix(gm);
        for (int i = 0; i < ans.Size; i++)
        {
            for (int j = 0; j < ans._ig[i + 1] - ans._ig[i]; j++)
            {
                ans._al[ans._ig[i] + j] *= a;
                ans._au[ans._ig[i] + j] *= a;
            }
            ans._diag[i] *= a;
        }
        return ans;
    }


    public static GlobalMatrix MultiplyPar(double a, GlobalMatrix gm)
    {
        var ans = new GlobalMatrix(gm);
        Parallel.For(0, ans.Size, i => {
            for (int j = 0; j < ans._ig[i + 1] - ans._ig[i]; j++)
            {
                ans._al[ans._ig[i] + j] *= a;
                ans._au[ans._ig[i] + j] *= a;
            }
            ans._diag[i] *= a;
        });
        return ans;
    }


    public static GlobalMatrix Sum(GlobalMatrix gm1, GlobalMatrix gm2)
    {
        if (!gm1.CheckPortrait(gm2)) throw new ArgumentException("Different matrixes portrait!");
        GlobalMatrix ans = new(gm1);
        for (int i = 0; i < ans.Size; i++)
        {
            for (int j = 0; j < ans._ig[i + 1] - ans._ig[i]; j++)
            {
                ans._al[ans._ig[i] + j] += gm2._al[gm2._ig[i] + j];
                ans._au[ans._ig[i] + j] += gm2._au[gm2._ig[i] + j];
            }
            ans._diag[i] += gm2._diag[i];
        }
        return ans;
    }


    public static GlobalMatrix SumPar(GlobalMatrix gm1, GlobalMatrix gm2)
    {
        if (!gm1.CheckPortrait(gm2)) throw new ArgumentException("Different matrixes portrait!");
        GlobalMatrix ans = new(gm1);

        Parallel.For(0, ans.Size, i => {
            for (int j = 0; j < ans._ig[i + 1] - ans._ig[i]; j++)
            {
                ans._al[ans._ig[i] + j] += gm2._al[gm2._ig[i] + j];
                ans._au[ans._ig[i] + j] += gm2._au[gm2._ig[i] + j];
            }
            ans._diag[i] += gm2._diag[i];
        });
        return ans;
    }
}
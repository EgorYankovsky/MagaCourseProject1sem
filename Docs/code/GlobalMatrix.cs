using DataStructs;
namespace MathObjects;

public enum TypeOfMatrixM
{
    Mr,
    Mrr
}


public class GlobalMatrix : Matrix
{
    public readonly int[]? _ig;
    public readonly List<int> _jg;
    public readonly double[]? _diag;
    public double[]? _al;
    public double[]? _au;
 
    public static GlobalVector operator *(GlobalMatrix _gm, GlobalVector _gv) => MathOperations.Multiply(_gm, _gv);
    public static GlobalMatrix operator *(double a, GlobalMatrix gm) => MathOperations.Multiply(a, gm);
    public static GlobalMatrix operator +(GlobalMatrix gm1, GlobalMatrix gm2) => MathOperations.Sum(gm1, gm2);

    public int Size
    {
        get
        {
            if (_diag is null) throw new Exception("_diag is null");
            return _diag.Length;
        }
    }

    public override double this[int i, int j]
    {
        get 
        {
            if (i > _diag.Length || j > _diag.Length)
                throw new Exception("Index ran out of matrix.");

            switch (i - j)
            {
                case 0: return _diag[i];
                case < 0: return ReturnValueAU(j, i);
                case > 0: return ReturnValueAL(i, j);
            }
        }
        set
        {
            if (i > _diag.Length || j > _diag.Length)
                throw new Exception("Index ran out of matrix.");
            switch (i - j)
            {
                case 0: _diag[i] = value; break;
                case < 0: SetValueAU(j, i, value); break;
                case > 0: SetValueAL(i, j, value); break;
            }
        }
    }

    private void SetValueAL(int i, int j, double val)
    {
        for (int ii = 0; ii < _ig[i + 1] - _ig[i]; ii++)
            if (_jg[_ig[i] + ii] == j)
                _al[_ig[i] + ii] = val;
    }

    private void SetValueAU(int i, int j, double val)
    {
        for (int ii = 0; ii < _ig[i + 1] - _ig[i]; ii++)
            if (_jg[_ig[i] + ii] == j)
                _au[_ig[i] + ii] = val;
    }
    
    private double ReturnValueAL(int i, int j)
    {
        for (int ii = 0; ii < _ig[i + 1] - _ig[i]; ii++)
            if (_jg[_ig[i] + ii] == j)
                return _al[_ig[i] + ii];
        return 0.0D;  
    }

    private double ReturnValueAU(int i, int j)
    {
        for (int ii = 0; ii < _ig[i + 1] - _ig[i]; ii++)
            if (_jg[_ig[i] + ii] == j)
                return _au[_ig[i] + ii];
        return 0.0D;  
    }   

    public GlobalMatrix Transpose() => new(_ig, _jg, _diag, _au, _al);

    public bool CheckPortrait(GlobalMatrix gm)
    {
        if (Size != gm.Size) return false;
        if (_ig.Length != gm._ig.Length) return false;
        if (_jg.Count != gm._jg.Count) return false;
        
        for (int i = 0; i < _ig.Length; i++)
            if (_ig[i] != gm._ig[i])
                return false;
        for (int i = 0; i < _jg.Count; i++)
            if (_jg[i] != gm._jg[i])
                return false;
        return true;
    }

    public GlobalMatrix(int[] ig, List<int> jg, double[] diag, double[] al, double[] au)
    {
        _ig = ig;
        _jg = jg;
        _diag = diag;
        _al = al;
        _au = au; 
    }

    public GlobalMatrix(GlobalMatrix gm)
    {
        _ig = (int[])gm._ig.Clone();
        _jg = gm._jg;
        _diag = (double[])gm._diag.Clone();
        _al = (double[])gm._al.Clone();
        _au = (double[])gm._au.Clone();
    }

    public GlobalMatrix(int arrOfPntLen)
    {        
        _jg = [];
        _ig = new int[arrOfPntLen + 1];
        _diag = new double[arrOfPntLen];
    }

    public GlobalMatrix(ArrayOfElems _arrOfElms, ArrayOfPoints _arrOfPnt, ArrayOfBorders _arrOfBord, List<double> mu0, double koef = 0.0)
    {
        _jg = [];
        _ig = new int[_arrOfPnt.GetLength() + 1];
        _diag = new double[_arrOfPnt.GetLength()];   
    }
}
using System.Drawing;
using DataStructs;
using Functions;

namespace MathObjects;

public class LocalVector : Vector
{
    
    private readonly double _r0;
    private readonly double _r1;
    private readonly double _z0;
    private readonly double _z1;
    private readonly double _hr;
    private readonly double _hz;
    private readonly double _t;
    public override int Size => 4;

    private readonly double[,] _M2R = {{1.0D, 1.0D},
                                       {1.0D, 3.0D}};    
    private readonly double[,] _Mz = {{2.0D, 1.0D},
                                      {1.0D, 2.0D}};
    private readonly double[,] _M1R = {{2.0D, 1.0D}, 
                                       {1.0D, 2.0D}}; 

    public override double this[int i] => i switch
    {
        0 => _Mr[0, 0] * _Mz[0, 0] * Function.F(_r0, _z0, _t) +
             _Mr[0, 1] * _Mz[0, 0] * Function.F(_r1, _z0, _t) +
             _Mr[0, 0] * _Mz[0, 1] * Function.F(_r0, _z1, _t) +
             _Mr[0, 1] * _Mz[0, 1] * Function.F(_r1, _z1, _t),

        1 => _Mr[1, 0] * _Mz[0, 0] * Function.F(_r0, _z0, _t) +
             _Mr[1, 1] * _Mz[0, 0] * Function.F(_r1, _z0, _t) +
             _Mr[1, 0] * _Mz[0, 1] * Function.F(_r0, _z1, _t) +
             _Mr[1, 1] * _Mz[0, 1] * Function.F(_r1, _z1, _t),

        2 => _Mr[0, 0] * _Mz[1, 0] * Function.F(_r0, _z0, _t) +
             _Mr[0, 1] * _Mz[1, 0] * Function.F(_r1, _z0, _t) +
             _Mr[0, 0] * _Mz[1, 1] * Function.F(_r0, _z1, _t) +
             _Mr[0, 1] * _Mz[1, 1] * Function.F(_r1, _z1, _t),

        3 => _Mr[1, 0] * _Mz[1, 0] * Function.F(_r0, _z0, _t) +
             _Mr[1, 1] * _Mz[1, 0] * Function.F(_r1, _z0, _t) +
             _Mr[1, 0] * _Mz[1, 1] * Function.F(_r0, _z1, _t) +
             _Mr[1, 1] * _Mz[1, 1] * Function.F(_r1, _z1, _t),

        _ => throw new IndexOutOfRangeException("Vector out of index"),
    };

    private double[,] _Mr = new double[2, 2];
    
    public LocalVector(List<int> elem, ArrayOfPoints2D arrPt, double t)
    {
        _r0 = arrPt[elem[0]].R;
        _r1 = arrPt[elem[1]].R;
        _t = t;
        
        _z0 = arrPt[elem[0]].Z;
        _z1 = arrPt[elem[2]].Z;
        _hr = _r1 - _r0;
        _hz = _z1 - _z0;
    
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                _Mr[i, j] = (_hr / 6.0) * (_r0 * _M1R[i, j] + (_hr / 2.0) * _M2R[i, j]);
                _Mz[i, j] = (_hz / 6.0) * _Mz[i, j];
            }
        }
    }

    private void WriteVector()
    {
        for (int i = 0; i < 4; i++)
            Console.WriteLine($"{this[i]:E5}");
    }

    public LocalVector(double r0, double r1, double z0, double z1)
    {
        _r0 = r0;
        _r1 = r1;
        _z0 = z0;
        _z1 = z1;
        _hr = _r1 - _r0;
        _hz = _z1 - _z0;
    } 
}
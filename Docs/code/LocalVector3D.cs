using static Functions.Function;

namespace MathObjects;


public class LocalVector3D(Egetter egetter, double x0, double x1, double y0, double y1, double z0, double z1, double t) : Vector
{
    private readonly double x0 = x0;
    private readonly double x1 = x1;
    private readonly double y0 = y0;
    private readonly double y1 = y1;
    private readonly double z0 = z0;
    private readonly double z1 = z1;
    private readonly double xm = 0.5D * (x1 + x0);
    private readonly double ym = 0.5D * (y1 + y0);
    private readonly double zm = 0.5D * (z1 + z0);
    private readonly double t = t;

    private readonly LocalMatrixM3D M = new(1.0, x1 - x0, y1 - y0, z1 - z0);
    private readonly Egetter _egetter = egetter;

    public override double this[int i] 
    { 
        get 
        {
            static double ScalarMult((double, double, double) a, (double, double, double) b) 
            => a.Item1 * b.Item1 + a.Item2 * b.Item2 + a.Item3 * b.Item3;
            List<double> vect = [
                ScalarMult(_egetter(xm, y0, z0, t), (1.0D, 0.0D, 0.0D)),
                ScalarMult(_egetter(xm, y1, z0, t), (1.0D, 0.0D, 0.0D)),
                ScalarMult(_egetter(xm, y0, z1, t), (1.0D, 0.0D, 0.0D)),
                ScalarMult(_egetter(xm, y1, z1, t), (1.0D, 0.0D, 0.0D)),
                ScalarMult(_egetter(x0, ym, z0, t), (0.0D, 1.0D, 0.0D)),
                ScalarMult(_egetter(x1, ym, z0, t), (0.0D, 1.0D, 0.0D)),
                ScalarMult(_egetter(x0, ym, z1, t), (0.0D, 1.0D, 0.0D)),
                ScalarMult(_egetter(x1, ym, z1, t), (0.0D, 1.0D, 0.0D)),
                ScalarMult(_egetter(x0, y0, zm, t), (0.0D, 0.0D, 1.0D)),
                ScalarMult(_egetter(x1, y0, zm, t), (0.0D, 0.0D, 1.0D)),
                ScalarMult(_egetter(x0, y1, zm, t), (0.0D, 0.0D, 1.0D)),
                ScalarMult(_egetter(x1, y1, zm, t), (0.0D, 0.0D, 1.0D)) 
            ];
            double ans = 0.0D;
            for (int j = 0; j < vect.Count; j++)
                ans += M[i, j] * vect[j];
            return ans;
        }
    }
}
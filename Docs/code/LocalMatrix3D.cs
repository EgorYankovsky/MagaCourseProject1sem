namespace MathObjects;

public class LocalMatrixG3D (double mu, double hx, double hy, double hz) : Matrix
{
    private readonly double _mu = mu;
    private readonly double _hx = hx;
    private readonly double _hy = hy;
    private readonly double _hz = hz;
    private readonly double[,] G1 = {{ 2.0D,  1.0D, -2.0D, -1.0D},
                                     { 1.0D,  2.0D, -1.0D, -2.0D},
                                     {-2.0D, -1.0D,  2.0D,  1.0D},
                                     {-1.0D, -2.0D,  1.0D,  2.0D}};
    private readonly double[,] G2 = {{ 2.0D, -2.0D,  1.0D, -1.0D},
                                     {-2.0D,  2.0D, -1.0D,  1.0D},
                                     { 1.0D, -1.0D,  2.0D, -2.0D},
                                     {-1.0D,  1.0D, -2.0D,  2.0D}};
    private readonly double[,] G3 = {{-2.0D,  2.0D, -1.0D,  1.0D},
                                     {-1.0D,  1.0D, -2.0D,  2.0D},
                                     { 2.0D, -2.0D,  1.0D, -1.0D},
                                     { 1.0D, -1.0D,  2.0D, -2.0D}};

    public override double this[int i, int j]
    {
        get
        {
            return (i / 4, j / 4)
            switch
            {
                (0, 0) => (_hx * _hy / (6.0D * _hz) * G1[i % 4, j % 4] + _hx * _hz / (6.0D * _hy) * G2[i % 4, j % 4]) / _mu,
                (0, 1) or (1, 0) => -1.0D * (_hz / 6.0D) * G2[i % 4, j % 4] / _mu,
                (0, 2) => _hy / 6.0D * G3[i % 4, j % 4] / _mu,
                (1, 1) => (_hx * _hy / (6.0D * _hz) * G1[i % 4, j % 4] + _hy * _hz / (6.0D * _hx) * G2[i % 4, j % 4]) / _mu,
                (1, 2) or (2, 1) => -1.0D * _hx / 6.0D * G1[i % 4, j % 4] / _mu,
                (2, 0) => _hy / 6.0D * G3[j % 4, i % 4] / _mu,
                (2, 2) => (_hx * _hz / (6.0D * _hy) * G1[i % 4, j % 4] + _hy * _hz / (6.0D * _hx) * G2[i % 4, j % 4]) / _mu,
                _ => throw new ArgumentOutOfRangeException("Out of local matrix 3d range"),
            };
        }
        set{}
    }
}

public class LocalMatrixM3D(double gamma, double hx, double hy, double hz) : Matrix
{
    private readonly double _gamma = gamma;
    private readonly double _hx = hx;
    private readonly double _hy = hy;
    private readonly double _hz = hz;
    private readonly double[,] D = {{4.0D, 2.0D, 2.0D, 1.0D},
                                    {2.0D, 4.0D, 1.0D, 2.0D},
                                    {2.0D, 1.0D, 4.0D, 2.0D},
                                    {1.0D, 2.0D, 2.0D, 4.0D}};

    public override double this[int i, int j]
    {
        get
        {
            return (i / 4, j / 4)
            switch
            {
                (0, 0) or (1, 1) or (2, 2) => _gamma * _hx * _hy * _hz / 36.0D * D[i % 4, j % 4],
                (0, 1) or (0, 2) or (1, 0) or (1, 2) or (2, 0) or (2, 1) => 0.0D,
                _ => throw new ArgumentOutOfRangeException("Out of local matrix 3d range"),
            };
        }
        set
        {}
    }
}
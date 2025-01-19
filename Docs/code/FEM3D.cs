using MathObjects;
using Grid;
using DataStructs;
using Functions;
using System.Diagnostics;



namespace Project;

public class FEM3D : FEM
{
    private static readonly double mu0 = 4.0D * Math.PI * Math.Pow(10.0D, -7);
    public ArrayOfPoints3D pointsArr;
    public List<GlobalVector> A;
    public List<GlobalVector> E;
    public List<GlobalVector> B;
    public List<GlobalVector> H;
    public ArrayOfRibs? ribsArr;
    private readonly Mesh3Dim mesh3Dim;
    private GlobalMatrix? G;
    List<FEM3D> additionalFields = [];
    private readonly FEM3D? _originalFEM;
    private List<GlobalVector>? _originalE;
    private GlobalMatrix? M;
    internal static List<int> ConvertGlobalToLocalNumeration(List<int> global) => 
            [global[0], global[3], global[8], global[11],
             global[1], global[2], global[9], global[10],
             global[4], global[5], global[6], global[7]];

    public FEM3D(Mesh3Dim mesh, TimeMesh timeMesh) : base(timeMesh)
    {
        string pointsPath = _3dValuesPath + $"AfterConvertation\\Points.poly";
        string elemsPath = _3dValuesPath + $"AfterConvertation\\Elems.poly";
        string bordersPath = _3dValuesPath + $"AfterConvertation\\Borders.poly";
        pointsArr = new ArrayOfPoints3D(pointsPath);
        elemsArr = new(elemsPath, 3);
        bordersArr = new (bordersPath, 3);
        ribsArr = mesh.arrayOfRibs;
        equationType = timeMesh[0] == timeMesh[^1] ? EquationType.Elliptic : EquationType.Parabolic;
        A = [];
        E = [];
        B = [];
        H = [];
        mesh3Dim = mesh;
    }

    public FEM3D(Mesh3Dim mesh, TimeMesh timeMesh, FEM3D originalFEM, int Num) : base(timeMesh)
    {
        string pointsPath = _3dValuesPath + $"Anomaly{Num}\\Points.poly";
        string elemsPath = _3dValuesPath + $"Anomaly{Num}\\Elems.poly";
        string bordersPath = _3dValuesPath + $"Anomaly{Num}\\Borders.poly";
        pointsArr = new(pointsPath);
        elemsArr = new(elemsPath, 3);
        bordersArr = new(bordersPath, 3);
        ribsArr = mesh.arrayOfRibs;
        _originalFEM = originalFEM;
        _originalE = [];
        equationType = timeMesh[0] == timeMesh[^1] ? EquationType.Elliptic : EquationType.Parabolic;
        A = [];
        E = [];
        B = [];
        H = [];
        mesh3Dim = mesh;
        ConstructMatrixes();
    }

    public (double, double, double) GetAAt(double x, double y, double z, double t)
    {
        if (ribsArr is null) throw new ArgumentNullException("Array of ribs not generated");
        for (int tt = 0; tt < Time.Count; tt++)
            if (Time[tt] == t)
            {
                var arr = GetElem(x, y, z);
                if (arr is null) return (0.0D, 0.0D, 0.0D);
                var elem = ConvertGlobalToLocalNumeration(arr);
                double[] q = new double[12];
                for (int i = 0; i < 12; i++)
                    q[i] = A[tt][elem[i]];
                double x0 = ribsArr[elem[0]].a.X;
                double x1 = ribsArr[elem[0]].b.X;
                double y0 = ribsArr[elem[4]].a.Y;
                double y1 = ribsArr[elem[4]].b.Y;
                double z0 = ribsArr[elem[8]].a.Z;
                double z1 = ribsArr[elem[8]].b.Z;
                double eps = (x - x0) / (x1 - x0);
                double nu =  (y - y0) / (y1 - y0);
                double khi = (z - z0) / (z1 - z0);
                var ans = BasisFunctions3DVec.GetValue(eps, nu, khi, q);
                foreach (var solution in additionalFields)
                {
                    var arrCurr = solution.GetElem(x, y, z);
                    if (arrCurr is null) continue;
                    var elemCurr = ConvertGlobalToLocalNumeration(arrCurr);
                    double[] qCurr = new double[12];
                    for (int ii = 0; ii < 12; ii++)
                        qCurr[ii] = solution.A[tt][elemCurr[ii]];
                    double x0Curr = solution.ribsArr[elemCurr[0]].a.X;
                    double x1Curr = solution.ribsArr[elemCurr[0]].b.X;
                    double y0Curr = solution.ribsArr[elemCurr[4]].a.Y;
                    double y1Curr = solution.ribsArr[elemCurr[4]].b.Y;
                    double z0Curr = solution.ribsArr[elemCurr[8]].a.Z;
                    double z1Curr = solution.ribsArr[elemCurr[8]].b.Z;
                    double epsCurr = (x - x0Curr) / (x1Curr - x0Curr);
                    double nuCurr =  (y - y0Curr) / (y1Curr - y0Curr);
                    double khiCurr = (z - z0Curr) / (z1Curr - z0Curr);
                    var ansCurr = BasisFunctions3DVec.GetValue(epsCurr, nuCurr, khiCurr, qCurr);
                    ans.Item1 += ansCurr.Item1;
                    ans.Item2 += ansCurr.Item2;
                    ans.Item3 += ansCurr.Item3;
                }
                return ans;            
            }
        throw new Exception("Out of mesh borders");
    }

    public (double, double, double) AdditioanalGetEAt(double x, double y, double z, double t)
    {
        if (ribsArr is null) throw new ArgumentNullException("Array of ribs not generated");
        for (int tt = 0; tt < Time.Count; tt++)
            if (Time[tt] == t)
            {
                var ans = (0.0D, 0.0D, 0.0D);
                foreach (var solution in additionalFields)
                {
                    var arrCurr = solution.GetElem(x, y, z);
                    if (arrCurr is null) continue;
                    var elemCurr = ConvertGlobalToLocalNumeration(arrCurr);
                    double[] qCurr = new double[12];
                    for (int ii = 0; ii < 12; ii++)
                        qCurr[ii] = solution.E[tt][elemCurr[ii]];
                    double x0Curr = solution.ribsArr[elemCurr[0]].a.X;
                    double x1Curr = solution.ribsArr[elemCurr[0]].b.X;
                    double y0Curr = solution.ribsArr[elemCurr[4]].a.Y;
                    double y1Curr = solution.ribsArr[elemCurr[4]].b.Y;
                    double z0Curr = solution.ribsArr[elemCurr[8]].a.Z;
                    double z1Curr = solution.ribsArr[elemCurr[8]].b.Z;
                    double epsCurr = (x - x0Curr) / (x1Curr - x0Curr);
                    double nuCurr =  (y - y0Curr) / (y1Curr - y0Curr);
                    double khiCurr = (z - z0Curr) / (z1Curr - z0Curr);
                    var ansCurr = BasisFunctions3DVec.GetValue(epsCurr, nuCurr, khiCurr, qCurr);
                    ans.Item1 += ansCurr.Item1;
                    ans.Item2 += ansCurr.Item2;
                    ans.Item3 += ansCurr.Item3;
                }
                return ans;
            }
        throw new Exception("Out of mesh borders");
    }

    public (double, double, double) GetEAt(double x, double y, double z, double t)
    {
        if (ribsArr is null) throw new ArgumentNullException("Array of ribs not generated");
        for (int tt = 0; tt < Time.Count; tt++)
            if (Time[tt] == t)
            {
                var arr = GetElem(x, y, z);
                if (arr is null) return (0.0D, 0.0D, 0.0D);
                var elem = ConvertGlobalToLocalNumeration(arr);
                double[] q = new double[12];
                for (int i = 0; i < 12; i++)
                    q[i] = E[tt][elem[i]];
                double x0 = ribsArr[elem[0]].a.X;
                double x1 = ribsArr[elem[0]].b.X;
                double y0 = ribsArr[elem[4]].a.Y;
                double y1 = ribsArr[elem[4]].b.Y;
                double z0 = ribsArr[elem[8]].a.Z;
                double z1 = ribsArr[elem[8]].b.Z;
                double eps = (x - x0) / (x1 - x0);
                double nu =  (y - y0) / (y1 - y0);
                double khi = (z - z0) / (z1 - z0);
                var ans = BasisFunctions3DVec.GetValue(eps, nu, khi, q);
                foreach (var solution in additionalFields)
                {
                    var arrCurr = solution.GetElem(x, y, z);
                    if (arrCurr is null) continue;
                    var elemCurr = ConvertGlobalToLocalNumeration(arrCurr);
                    double[] qCurr = new double[12];
                    for (int ii = 0; ii < 12; ii++)
                        qCurr[ii] = solution.E[tt][elemCurr[ii]];
                    double x0Curr = solution.ribsArr[elemCurr[0]].a.X;
                    double x1Curr = solution.ribsArr[elemCurr[0]].b.X;
                    double y0Curr = solution.ribsArr[elemCurr[4]].a.Y;
                    double y1Curr = solution.ribsArr[elemCurr[4]].b.Y;
                    double z0Curr = solution.ribsArr[elemCurr[8]].a.Z;
                    double z1Curr = solution.ribsArr[elemCurr[8]].b.Z;
                    double epsCurr = (x - x0Curr) / (x1Curr - x0Curr);
                    double nuCurr =  (y - y0Curr) / (y1Curr - y0Curr);
                    double khiCurr = (z - z0Curr) / (z1Curr - z0Curr);
                    var ansCurr = BasisFunctions3DVec.GetValue(epsCurr, nuCurr, khiCurr, qCurr);
                    ans.Item1 += ansCurr.Item1;
                    ans.Item2 += ansCurr.Item2;
                    ans.Item3 += ansCurr.Item3;
                }
                return ans;
            }
        throw new Exception("Out of mesh borders");
    }

    public (double, double, double) GetBAt(double x, double y, double z, double t)
    {
        if (ribsArr is null) throw new ArgumentNullException("Array of ribs not generated");
        for (int tt = 0; tt < Time.Count; tt++)
            if (Time[tt] == t)
            {
                var arr = GetElem(x, y, z);
                if (arr is null) return (0.0D, 0.0D, 0.0D);
                var elem = ConvertGlobalToLocalNumeration(arr);
                double[] q = new double[12];
                for (int i = 0; i < 12; i++)
                    q[i] = B[tt][elem[i]];
                double x0 = ribsArr[elem[0]].a.X;
                double x1 = ribsArr[elem[0]].b.X;
                double y0 = ribsArr[elem[4]].a.Y;
                double y1 = ribsArr[elem[4]].b.Y;
                double z0 = ribsArr[elem[8]].a.Z;
                double z1 = ribsArr[elem[8]].b.Z;
                double eps = (x - x0) / (x1 - x0);
                double nu =  (y - y0) / (y1 - y0);
                double khi = (z - z0) / (z1 - z0);
                var ans = BasisFunctions3DVec.GetValue(eps, nu, khi, q);
                foreach (var solution in additionalFields)
                {
                    var arrCurr = solution.GetElem(x, y, z);
                    if (arrCurr is null) continue;
                    var elemCurr = ConvertGlobalToLocalNumeration(arrCurr);
                    double[] qCurr = new double[12];
                    for (int ii = 0; ii < 12; ii++)
                        qCurr[ii] = solution.B[tt][elemCurr[ii]];
                    double x0Curr = solution.ribsArr[elemCurr[0]].a.X;
                    double x1Curr = solution.ribsArr[elemCurr[0]].b.X;
                    double y0Curr = solution.ribsArr[elemCurr[4]].a.Y;
                    double y1Curr = solution.ribsArr[elemCurr[4]].b.Y;
                    double z0Curr = solution.ribsArr[elemCurr[8]].a.Z;
                    double z1Curr = solution.ribsArr[elemCurr[8]].b.Z;
                    double epsCurr = (x - x0Curr) / (x1Curr - x0Curr);
                    double nuCurr =  (y - y0Curr) / (y1Curr - y0Curr);
                    double khiCurr = (z - z0Curr) / (z1Curr - z0Curr);
                    var ansCurr = BasisFunctions3DVec.GetValue(epsCurr, nuCurr, khiCurr, qCurr);
                    ans.Item1 += ansCurr.Item1;
                    ans.Item2 += ansCurr.Item2;
                    ans.Item3 += ansCurr.Item3;
                }
                return ans;
            }
        throw new Exception("Out of mesh borders");
    }

    public (double, double, double) GetHAt(double x, double y, double z, double t)
    {
        if (ribsArr is null) throw new ArgumentNullException("Array of ribs not generated");
        for (int tt = 0; tt < Time.Count; tt++)
            if (Time[tt] == t)
            {
                var arr = GetElem(x, y, z);
                if (arr is null) return (0.0D, 0.0D, 0.0D);
                var elem = ConvertGlobalToLocalNumeration(arr);
                double[] q = new double[12];
                for (int i = 0; i < 12; i++)
                    q[i] = H[tt][elem[i]];
                double x0 = ribsArr[elem[0]].a.X;
                double x1 = ribsArr[elem[0]].b.X;
                double y0 = ribsArr[elem[4]].a.Y;
                double y1 = ribsArr[elem[4]].b.Y;
                double z0 = ribsArr[elem[8]].a.Z;
                double z1 = ribsArr[elem[8]].b.Z;
                double eps = (x - x0) / (x1 - x0);
                double nu =  (y - y0) / (y1 - y0);
                double khi = (z - z0) / (z1 - z0);
                var ans = BasisFunctions3DVec.GetValue(eps, nu, khi, q);
                foreach (var solution in additionalFields)
                {
                    var arrCurr = solution.GetElem(x, y, z);
                    if (arrCurr is null) continue;
                    var elemCurr = ConvertGlobalToLocalNumeration(arrCurr);
                    double[] qCurr = new double[12];
                    for (int ii = 0; ii < 12; ii++)
                        qCurr[ii] = solution.H[tt][elemCurr[ii]];
                    double x0Curr = solution.ribsArr[elemCurr[0]].a.X;
                    double x1Curr = solution.ribsArr[elemCurr[0]].b.X;
                    double y0Curr = solution.ribsArr[elemCurr[4]].a.Y;
                    double y1Curr = solution.ribsArr[elemCurr[4]].b.Y;
                    double z0Curr = solution.ribsArr[elemCurr[8]].a.Z;
                    double z1Curr = solution.ribsArr[elemCurr[8]].b.Z;
                    double epsCurr = (x - x0Curr) / (x1Curr - x0Curr);
                    double nuCurr =  (y - y0Curr) / (y1Curr - y0Curr);
                    double khiCurr = (z - z0Curr) / (z1Curr - z0Curr);
                    var ansCurr = BasisFunctions3DVec.GetValue(epsCurr, nuCurr, khiCurr, qCurr);
                    ans.Item1 += ansCurr.Item1;
                    ans.Item2 += ansCurr.Item2;
                    ans.Item3 += ansCurr.Item3;
                }
                return ans;
            }
        throw new Exception("Out of mesh borders");
    }

    // B = rot A
    public void GenerateVectorB()
    {
        if (A is []) throw new Exception("Vector B isn't generated");
        if (ribsArr is null) throw new Exception("Array of ribs isn't generated");
        B = new(A.Count);
        for (int t = 0; t < Time.Count; t++)
        {
            B.Add(new GlobalVector(ribsArr.Count));
            for (int i = 0; i < A[t].Size; i++)
            {
                double ht = Math.Pow(10, -10);
                var pnt = ribsArr[i].GetMiddlePoint();
                var tan = ribsArr[i].GetTangent();
                var vec1 = GetAAt(pnt.X, pnt.Y, pnt.Z + ht, Time[t]);
                var vec2 = GetAAt(pnt.X, pnt.Y, pnt.Z - ht, Time[t]);
                var Bx = -1.0D * (vec1.Item2 - vec2.Item2) / (2.0D * ht);
                vec1 = GetAAt(pnt.X, pnt.Y, pnt.Z + ht, Time[t]);
                vec2 = GetAAt(pnt.X, pnt.Y, pnt.Z - ht, Time[t]);
                var By = (vec1.Item1 - vec2.Item1) / (2.0D * ht);
                vec1 = GetAAt(pnt.X + ht, pnt.Y, pnt.Z, Time[t]);
                vec2 = GetAAt(pnt.X - ht, pnt.Y, pnt.Z, Time[t]);
                var Bz1 = (vec1.Item2 - vec2.Item2) / (2.0D * ht);
                vec1 = GetAAt(pnt.X, pnt.Y + ht, pnt.Z, Time[t]);
                vec2 = GetAAt(pnt.X, pnt.Y - ht, pnt.Z, Time[t]);
                var Bz2 = (vec1.Item1 - vec2.Item1) / (2.0D * ht);
                B[t][i] = Bx * tan.Item1 + By * tan.Item2 + (Bz1 - Bz2) * tan.Item3;
            }
        }
    }

    // H = B / mu_0
    public void GenerateVectorH()
    {
        if (B is []) throw new Exception("Vector B isn't generated");
        H = new(B.Count);
        for (int t = 0; t < B.Count; t++)
        {
            H.Add(new GlobalVector(B[t].Size));
            for (int i = 0; i < B[i].Size; i++)
                H[t][i] = B[t][i] / mu0;
        }
    }

    // E = - dA / dt
    public void GenerateVectorE()
    {
        E = new(A.Count);   
        for (int i = 0; i < A.Count; i++)
        {
            if (i == 0)
                E.Add(new GlobalVector(A[1].Size));
            else
            {
                if (A[i - 1].Size == 0) A[i - 1] = new GlobalVector(ribsArr.Count);
                if (A[i].Size == 0) A[i] = new GlobalVector(ribsArr.Count);
                E.Add(-1.0D / (Time[i] - Time[i - 1]) * (A[i] - A[i - 1]));
            }
        }
    }

    public void AddSolution(FEM3D fem) => additionalFields.Add(fem);

    public List<int>? GetElem(double x, double y, double z)
    {
        if (x < mesh3Dim.nodesX[0] || mesh3Dim.nodesX[^1] < x ||
            y < mesh3Dim.nodesY[0] || mesh3Dim.nodesY[^1] < y ||
            z < mesh3Dim.nodesZ[0] || mesh3Dim.nodesZ[^1] < z)
            return null;
        int i = 0;
        for (; i < mesh3Dim.nodesX.Count - 1; i++)
            if (mesh3Dim.nodesX[i] <= x && x <= mesh3Dim.nodesX[i + 1])
                break;
        int j = 0;
        for (; j < mesh3Dim.nodesY.Count - 1; j++)
            if (mesh3Dim.nodesY[j] <= y && y <= mesh3Dim.nodesY[j + 1])
                break;
        int k = 0;
        for (; k < mesh3Dim.nodesZ.Count - 1; k++)
            if (mesh3Dim.nodesZ[k] <= z && z <= mesh3Dim.nodesZ[k + 1])
                break;
        return elemsArr[k * (mesh3Dim.nodesX.Count - 1) * (mesh3Dim.nodesY.Count - 1) + j * (mesh3Dim.nodesX.Count - 1) + i].Arr;
    }

    public void ConvertResultTo3Dim(FEM2D fem2d)
    {
        if (fem2d.pointsArr is null) throw new ArgumentNullException();
        if (ribsArr is null) throw new ArgumentNullException();
        for (int i = 0; i < Time.Count; i++)
        {
            A.Add(new GlobalVector(ribsArr.Count));
            E.Add(new GlobalVector(ribsArr.Count));
            for (int j = 0; j < ribsArr.Count; j++)
            {
                var X = 0.5D * (ribsArr[j].a.X + ribsArr[j].b.X);
                var Y = 0.5D * (ribsArr[j].a.Y + ribsArr[j].b.Y);
                var Z = 0.5D * (ribsArr[j].a.Z + ribsArr[j].b.Z);
                var antinormal = ((ribsArr[j].b.X - ribsArr[j].a.X) / ribsArr[j].Length,
                                  (ribsArr[j].b.Y - ribsArr[j].a.Y) / ribsArr[j].Length,
                                  (ribsArr[j].b.Z - ribsArr[j].a.Z) / ribsArr[j].Length);
                double fa = 0.0;
                double fe = 0.0;
                var elem = fem2d.GetElem(Math.Sqrt(X * X + Y * Y), Z);
                if (elem is not null)
                {
                    fa = BasisFunctions2D.GetValue(
                                    fem2d.A_phi[i][elem[0]], fem2d.A_phi[i][elem[1]],
                                    fem2d.A_phi[i][elem[2]], fem2d.A_phi[i][elem[3]],
                                    fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
                                    fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
                                    Math.Sqrt(X * X + Y * Y), Z);
                    fe = BasisFunctions2D.GetValue(
                                    fem2d.E_phi[i][elem[0]], fem2d.E_phi[i][elem[1]],
                                    fem2d.E_phi[i][elem[2]], fem2d.E_phi[i][elem[3]],
                                    fem2d.pointsArr[elem[0]].R, fem2d.pointsArr[elem[1]].R,
                                    fem2d.pointsArr[elem[0]].Z, fem2d.pointsArr[elem[3]].Z, 
                                    Math.Sqrt(X * X + Y * Y), Z);
                }
                var Ax = X == 0 && Y == 0 ? 0.0D : -1.0D * (Y / Math.Sqrt(X * X + Y * Y)) * fa;
                var Ay = X == 0 && Y == 0 ? 0.0D : X / Math.Sqrt(X * X + Y * Y) * fa;
                var Az = 0.0D;
                var Ex = X == 0 && Y == 0 ? 0.0D :  -1.0D * (Y / Math.Sqrt(X * X + Y * Y)) * fe;
                var Ey = X == 0 && Y == 0 ? 0.0D :  X / Math.Sqrt(X * X + Y * Y) * fe;
                var Ez = 0.0D;
                A[i][j] = Ax * antinormal.Item1 + Ay * antinormal.Item2 + Az * antinormal.Item3;
                E[i][j] = Ex * antinormal.Item1 + Ey * antinormal.Item2 + Ez * antinormal.Item3;
            }
        }
    }

    public void CheckSolution(List<Point3D> recivers)
    {
        using var sw = new StreamWriter("C:\\Users\\USER\\Desktop\\Test.txt");
        GenerateVectorB();
        GenerateVectorH();
        for (int t = 0; t < Time.Count - 1; t++)
        {
            sw.WriteLine($"Time: {Time[t]}:E15");
            foreach (var reciver in recivers)
            {
                var BAt = GetBAt(reciver.X, reciver.Y, reciver.Z, Time[t]);
                var BAt_1 = GetBAt(reciver.X, reciver.Y, reciver.Z, Time[t + 1]);
                var dBdt = ((BAt_1.Item1 - BAt.Item1) / (Time[t + 1] - Time[t]),
                            (BAt_1.Item2 - BAt.Item2) / (Time[t + 1] - Time[t]),
                            (BAt_1.Item3 - BAt.Item3) / (Time[t + 1] - Time[t]));                
                var Ex_1 = GetEAt(reciver.X - 1E-10, reciver.Y, reciver.Z, Time[t]);
                var Ex_2 = GetEAt(reciver.X + 1E-10, reciver.Y, reciver.Z, Time[t]);
                var Ey_1 = GetEAt(reciver.X, reciver.Y - 1E-10, reciver.Z, Time[t]);
                var Ey_2 = GetEAt(reciver.X, reciver.Y + 1E-10, reciver.Z, Time[t]);
                var Ez_1 = GetEAt(reciver.X, reciver.Y, reciver.Z - 1E-10, Time[t]);
                var Ez_2 = GetEAt(reciver.X, reciver.Y, reciver.Z + 1E-10, Time[t]);
                var rotE = ((Ey_2.Item3 - Ey_1.Item3) / (2.0 * 1E-10) - (Ez_2.Item2 - Ez_1.Item2) / (2.0 * 1E-10),
                   -1.0D * ((Ex_2.Item3 - Ex_1.Item3) / (2.0 * 1E-10) - (Ez_2.Item1 - Ez_1.Item1) / (2.0 * 1E-10)),
                            (Ex_2.Item2 - Ex_1.Item2) / (2.0 * 1E-10) - (Ey_2.Item1 - Ey_1.Item1) / (2.0 * 1E-10));
                sw.WriteLine($"reciver: ({reciver.X:E15}, {reciver.Y:E15}, {reciver.Z:E15})");
                sw.WriteLine($"  rot E: ({rotE.Item1:E15}, {rotE.Item2:E15}, {rotE.Item3:E15})");
                sw.WriteLine($"   dBdt: ({dBdt.Item1:E15}, {dBdt.Item2:E15}, {dBdt.Item3:E15})");
                sw.WriteLine($"discep.: ({Math.Abs(rotE.Item1 - dBdt.Item1):E15}, {Math.Abs(rotE.Item2 - dBdt.Item2):E15}, {Math.Abs(rotE.Item3 - dBdt.Item3):E15})");
            }
            sw.WriteLine();
        } 
        sw.Close();       
    }

    public void ConstructMatrixes()
    {
        if (ribsArr is null) throw new ArgumentNullException("ribsArr is null");
        var sparceMatrix = new GlobalMatrix(ribsArr.Count);
        Generator.BuildPortait(ref sparceMatrix, ribsArr.Count, elemsArr, true);
        G = new GlobalMatrix(sparceMatrix);
        Generator.FillMatrixG(ref G, ribsArr, elemsArr);
        M = new GlobalMatrix(sparceMatrix);
        Generator.FillMatrixM(ref M, ribsArr, elemsArr, mesh3Dim, _originalFEM.mesh3Dim);
    }

    public void Solve()
    {
        if (solver is null) throw new ArgumentNullException("Solver is null");
        if (ribsArr is null) throw new ArgumentNullException("ribs array is null");
        if (_originalFEM is null) throw new ArgumentNullException("original E is null");
        if (G is null) throw new ArgumentNullException();
        if (M is null) throw new ArgumentNullException();
        Stopwatch solutionStopwatch = new();
        solutionStopwatch.Start();
        Solutions = new GlobalVector[Time.Count];
        Discrepancy = new GlobalVector[Time.Count];
        (Solutions[0], Discrepancy[0]) = (new GlobalVector(ribsArr.Count), new GlobalVector(ribsArr.Count));
        if (Time.Count > 1)
        {
            (Solutions[1], Discrepancy[1]) = (Solutions[0], Discrepancy[0]);
            for (int i = 1; i < Time.Count; i++)
            {
                Console.WriteLine($"\n {i} / {Time.Count - 1}. Time layer: {Time[i]}");
                Thread.Sleep(1500);

                double deltT = Time[i] - Time[i - 2];
                double deltT0 = Time[i] - Time[i - 1];
                double deltT1 = Time[i - 1] - Time[i - 2];
                double tau0 = (deltT + deltT0) / (deltT * deltT0);
                double tau1 = deltT / (deltT1 * deltT0);
                double tau2 = deltT0 / (deltT * deltT1);
                double deltT = Time[i] - Time[i - 1];
                double tau = 1.0D / deltT;
                Matrix = G + tau0 * M;
                var b = new GlobalVector(ribsArr.Count);
                Generator.FillVector3D(ref b, _originalFEM.GetEAt, ribsArr, elemsArr, mesh3Dim, _originalFEM.mesh3Dim, Time[i]);
                Vector = b + (tau1 * (M * Solutions[i - 1])) - (tau2 * (M * Solutions[i - 2]));
                Generator.ConsiderBoundaryConditions(ref Matrix, ref Vector, ribsArr, bordersArr, Time[i]);
                var localStopWatch = Stopwatch.StartNew();
                (Solutions[i], Discrepancy[i]) = solver.Solve(Matrix, Vector);
                localStopWatch.Stop();
                var localMS = localStopWatch.ElapsedMilliseconds;
                Console.WriteLine($"Current iteration for {localMS / 60000} min {(localMS % 60000) / 1000} sec");
            }
        }
        A = [.. Solutions];
        solutionStopwatch.Stop();
        var milseconds = solutionStopwatch.ElapsedMilliseconds;
        Console.WriteLine($"Lin eq solved for {milseconds / 60000} min {(milseconds % 60000) / 1000} sec");
    }

    public void WriteData(string path)
    {
        if (Solutions is null) throw new ArgumentNullException("No solutions");
        for (int t = 0; t < Time.Count; t++)
        {
            using var sw = new StreamWriter(path + $"Answer_{Time[t]}.txt");
            for (int i = 0; i < A[t].Size; i++)
                sw.WriteLine($"{i} {A[t][i]:E8}");
            sw.Close();
        }
    }

    public void WriteDrawingForFirst(string path)
    {
        double hx = (50000.0D) / 300.0D;
        double hy = (25000.0D) / 300.0D;
        double hz = (mesh3Dim.nodesZ[^1] - mesh3Dim.nodesZ[0]) / 300.0D;
    
        for (int t = 0; t < Time.Count; t++)
        {
            using var swe = new StreamWriter(path + $"E\\Answer_E_time_layer_{t}.txt");
            for (int k = 0; k < 300; k++)
            {
                for (int i = 0; i < 300; i++)
                {
                    double zCurr = mesh3Dim.nodesZ[0] + k * hz;
                    double xCurr = 0.0D + i * hx;
                    double yCurr = 0.0D + i * hy;
                    var vec = GetEAt(xCurr, yCurr, zCurr, Time[t]);
                    var ans = Math.Sqrt(vec.Item1 * vec.Item1 + vec.Item2 * vec.Item2);
                    var rCurr = Math.Sqrt(xCurr * xCurr + yCurr * yCurr);
                    swe.WriteLine($"{rCurr:E15} {zCurr:E15} {ans:E15}");              
                }
            }
            swe.Close();
        }
    }

    public void WriteDrawingForSecond(string path)
    {
        double hx = (50000.0D) / 300.0D;
        double hy = -1.0D * (12500.0D) / 300.0D;
        double hz = (mesh3Dim.nodesZ[^1] - mesh3Dim.nodesZ[0]) / 300.0D;
        for (int t = 0; t < Time.Count; t++)
        {
            using var swe = new StreamWriter(path + $"E\\Answer_E_time_layer_{t}.txt");
            for (int k = 0; k < 300; k++)
            {
                for (int i = 0; i < 300; i++)
                {
                    double zCurr = mesh3Dim.nodesZ[0] + k * hz;
                    double xCurr = 0.0D + i * hx;
                    double yCurr = 0.0D + i * hy;
                    var vec = GetEAt(xCurr, yCurr, zCurr, Time[t]);
                    var ans = Math.Sqrt(vec.Item1 * vec.Item1 + vec.Item2 * vec.Item2 + vec.Item3 * vec.Item3);
                    var rCurr = Math.Sqrt(xCurr * xCurr + yCurr * yCurr);
                    swe.WriteLine($"{rCurr:E15} {zCurr:E15} {ans:E15}");                }
            }
            swe.Close();
        }
    }

    public void WriteDataToDraw2DimSolution(string path)
    {
        double hx = (25500.0D) / 300.0D;
        double hy = (25500.0D) / 300.0D;
        double hz = (mesh3Dim.nodesZ[^1] - mesh3Dim.nodesZ[0]) / 300.0D;
        for (int t = 0; t < Time.Count; t++)
        {
            using var swe = new StreamWriter(path + $"E\\Answer_E_time_layer_{t}.txt");
            for (int k = 0; k < 300; k++)
            {
                for (int i = 0; i < 300; i++)
                {
                    double zCurr = mesh3Dim.nodesZ[0] + k * hz;
                    double xCurr = 0.0D + i * hx;
                    double yCurr = 0.0D + i * hy;
                    var vec = GetEAt(xCurr, yCurr, zCurr, Time[t]);
                    var ans = Math.Sqrt(vec.Item1 * vec.Item1 + vec.Item2 * vec.Item2 + vec.Item3 * vec.Item3);
                    var rCurr = Math.Sqrt(xCurr * xCurr + yCurr * yCurr);
                    swe.WriteLine($"{rCurr:E15} {zCurr:E15} {ans:E15}");                
                }
            }
            swe.Close();
        }
    }

    public void WriteDataAtLine(string path)
    {
        double x0 = -1177.5;
        double hy = (mesh3Dim.nodesY[^1] - mesh3Dim.nodesY[0]) / 300.0D;
        double z0 = -1050.0D;
        List<double> times = [Time[0], Time[Time.Count / 2], Time[^1]];
        for (int t = 0; t < times.Count; t++)
        {
            using var swa = new StreamWriter(path + $"A_With_anomaly_{times[t]}.txt");

            for (int i = 0; i < 300; i++)
            {
                double xCurr = x0;
                double yCurr = mesh3Dim.nodesY[0] + i * hy;
                double zCurr = z0;
                var vec = GetAAt(xCurr, yCurr, zCurr, times[t]);
                var ans = vec.Item1;
                swa.WriteLine($"{yCurr:E15} {ans:E15}");
            }
            swa.Close();
        }
        for (int t = 0; t < times.Count; t++)
        {
            using var swe = new StreamWriter(path + $"E_With_anomaly_{times[t]}.txt");
            for (int i = 0; i < 300; i++)
            {
                double xCurr = x0;
                double yCurr = mesh3Dim.nodesY[0] + i * hy;
                double zCurr = z0;
                var vec = GetEAt(xCurr, yCurr, zCurr, times[t]);
                var ans = vec.Item1;
                swe.WriteLine($"{yCurr:E15} {ans:E15}");
            }
            swe.Close();
        }
    }

    public void WriteDataToDraw(string path)
    {
        double hx = (mesh3Dim.nodesX[^1] - mesh3Dim.nodesX[0]) / 10.0D;
        double hy = (mesh3Dim.nodesY[^1] - mesh3Dim.nodesY[0]) / 10.0D;
        double hz = (mesh3Dim.nodesZ[^1] - mesh3Dim.nodesZ[0]) / 10.0D;
        double x0 = mesh3Dim.nodesX[0];
        double y0 = mesh3Dim.nodesY[0];
        double z0 = mesh3Dim.nodesZ[0];
        List<Point3D> points = [];
        int i = 0;
        int j = 0;
        int k = 0;
        while (z0 + hz * k <= mesh3Dim.nodesZ[^1])
        {
            j = 0;
            while (y0 + hy * j <= mesh3Dim.nodesY[^1])
            {
                i = 0;
                while (x0 + hx * i <= mesh3Dim.nodesX[^1])
                {
                    points.Add(new Point3D(x0 + i * hx, y0 + j * hy, z0 + k * hz));
                    i++;
                }
                j++;
            }
            k++;
        }
        for (int t = 0; t < Time.Count; t++)
        {
            using var sw = new StreamWriter(path + $"Answer{t}.txt");
            sw.WriteLine($"{mesh3Dim.nodesX[0]:E15} {mesh3Dim.nodesX[^1]:E15} {mesh3Dim.nodesY[0]:E15} {mesh3Dim.nodesY[^1]:E15} {mesh3Dim.nodesZ[0]:E15} {mesh3Dim.nodesZ[^1]:E15}");
            foreach (var point in points)
            {
                var ans = GetEAt(point.X, point.Y, point.Z, Time[t]);
                sw.WriteLine($"{point.X:E15} {point.Y:E15} {point.Z:E15} {ans.Item1:E15} {ans.Item2:E15} {ans.Item3:E15}");
            }
            sw.Close();
        }
    }

    public void ReadData(string AnswerPath)
    {
        A = new(Time.Count);
        string file = string.Empty;
        for (int t = 0; t < Time.Count; t++)
        {
            file = $"Answer_{Time[t]}.txt";
            var fileData = File.ReadAllText(AnswerPath + file).Split("\n");
            A.Add(new GlobalVector(fileData.Length - 1));
            for (int i = 0; i < fileData.Length - 1; i++)
            {
                var val = fileData[i].Split(" ")[1];
                A[t][i] = double.Parse(val);
            }
        }
    }

    public void TestOutput(string path)
    {
        using var sw = new StreamWriter(path + "/A_phi/Answer3D/Answer_Test.txt");
        var absDiscX = 0.0D;
        var absDiscY = 0.0D;
        var absDiscZ = 0.0D;
        var absDivX = 0.0D;
        var absDivY = 0.0D;
        var absDivZ = 0.0D;
        var relDiscX = 0.0D;
        var relDiscY = 0.0D;
        var relDiscZ = 0.0D;
        var relDivX = 0.0D;
        var relDivY = 0.0D;
        var relDivZ = 0.0D;
        int iter = 0;
        foreach (Elem elem in elemsArr)
        {
            int[] elem_local = [elem[0], elem[3], elem[8], elem[11],
                                elem[1], elem[2], elem[9], elem[10],
                                elem[4], elem[5], elem[6], elem[7]];
            var x = 0.5D * (ribsArr[elem_local[0]].a.X + ribsArr[elem_local[0]].b.X);
            var y = 0.5D * (ribsArr[elem_local[4]].a.Y + ribsArr[elem_local[4]].b.Y);
            var z = 0.5D * (ribsArr[elem_local[8]].a.Z + ribsArr[elem_local[8]].b.Z);
            sw.WriteLine($"Points {x:E15} {y:E15} {z:E15}");
            var eps = (x - ribsArr[elem_local[0]].a.X) / (ribsArr[elem_local[0]].b.X - ribsArr[elem_local[0]].a.X);
            var nu =  (y - ribsArr[elem_local[4]].a.Y) / (ribsArr[elem_local[4]].b.Y - ribsArr[elem_local[4]].a.Y);
            var khi = (z - ribsArr[elem_local[8]].a.Z) / (ribsArr[elem_local[8]].b.Z - ribsArr[elem_local[8]].a.Z);
            double[] q = [Solutions[0][elem_local[0]], Solutions[0][elem_local[1]], Solutions[0][elem_local[2]], Solutions[0][elem_local[3]], 
                          Solutions[0][elem_local[4]], Solutions[0][elem_local[5]], Solutions[0][elem_local[6]], Solutions[0][elem_local[7]],
                          Solutions[0][elem_local[8]], Solutions[0][elem_local[9]], Solutions[0][elem_local[10]], Solutions[0][elem_local[11]]];
            var ans = BasisFunctions3DVec.GetValue(eps, nu, khi, q);
            var theorValue = Function.A(x, y, z, 0.0D);
            sw.WriteLine($"FEM A  {ans.Item1:E15} {ans.Item2:E15} {ans.Item3:E15}");
            sw.WriteLine($"Theor  {theorValue.Item1:E15} {theorValue.Item2:E15} {theorValue.Item3:E15}");
            var currAbsDiscX = Math.Abs(ans.Item1 - theorValue.Item1);
            var currAbsDiscY = Math.Abs(ans.Item2 - theorValue.Item2);
            var currAbsDiscZ = Math.Abs(ans.Item3 - theorValue.Item3);
            var currRelDiscX = currAbsDiscX / Math.Abs(theorValue.Item1);
            var currRelDiscY = currAbsDiscY / Math.Abs(theorValue.Item2);
            var currRelDiscZ = currAbsDiscZ / Math.Abs(theorValue.Item3);
            sw.WriteLine($"CurrAD {currAbsDiscX:E15} {currAbsDiscY:E15} {currAbsDiscZ:E15}");
            sw.WriteLine($"CurrRD {currRelDiscX:E15} {currRelDiscY:E15} {currRelDiscZ:E15}\n");
            absDiscX += currAbsDiscX;
            absDiscY += currAbsDiscY;
            absDiscZ += currAbsDiscZ;
            absDivX += theorValue.Item1;
            absDivY += theorValue.Item2;
            absDivZ += theorValue.Item3;
            relDiscX += currRelDiscX * currRelDiscX;
            relDiscY += currRelDiscY * currRelDiscY;
            relDiscZ += currRelDiscZ * currRelDiscZ;
            relDivX += theorValue.Item1 * theorValue.Item1;
            relDivY += theorValue.Item2 * theorValue.Item2;
            relDivZ += theorValue.Item3 * theorValue.Item3;
            iter++;
        }
        sw.WriteLine($"Avg disc: {absDiscX / iter:E15} {absDiscY / iter:E15} {absDiscZ / iter:E15}");
        sw.WriteLine($"Rel disc: {Math.Sqrt(relDiscX / relDivX):E15} {Math.Sqrt(relDiscY / relDivY):E15} {Math.Sqrt(relDiscZ / relDivZ):E15}");
        sw.Close();
    }

    public void MeasureValuesOnReceiversWithoutNormal(List<Point3D> recivers, string OutputPath)
    {
        using var sw_e = new StreamWriter(OutputPath + "E.txt");
        for (int t = 0; t < Time.Count; t++)
        {
            var a_e = AdditioanalGetEAt(recivers[0].X, recivers[0].Y, recivers[0].Z, Time[t]);
            var b_e = AdditioanalGetEAt(recivers[1].X, recivers[1].Y, recivers[1].Z, Time[t]);
            var c_e = AdditioanalGetEAt(recivers[2].X, recivers[2].Y, recivers[2].Z, Time[t]);
            var d_e = AdditioanalGetEAt(recivers[3].X, recivers[3].Y, recivers[3].Z, Time[t]);
            var f_a_e = Math.Sqrt(a_e.Item1 * a_e.Item1 + a_e.Item2 * a_e.Item2 + a_e.Item3 * a_e.Item3);
            var f_b_e = Math.Sqrt(b_e.Item1 * b_e.Item1 + b_e.Item2 * b_e.Item2 + b_e.Item3 * b_e.Item3);
            var f_c_e = Math.Sqrt(c_e.Item1 * c_e.Item1 + c_e.Item2 * c_e.Item2 + c_e.Item3 * c_e.Item3);
            var f_d_e = Math.Sqrt(d_e.Item1 * d_e.Item1 + d_e.Item2 * d_e.Item2 + d_e.Item3 * d_e.Item3);
            sw_e.WriteLine($"{Time[t]} {f_a_e} {f_b_e} {f_c_e} {f_d_e}");
        }
        sw_e.Close();
    }

    public void MeasureValuesOnReceivers(List<Point3D> recivers, string OutputPath)
    {
        using var sw_a = new StreamWriter(OutputPath + "A.txt");
        using var sw_e = new StreamWriter(OutputPath + "E.txt");
        for (int t = 0; t < Time.Count; t++)
        {
            var a_a = GetAAt(recivers[0].X, recivers[0].Y, recivers[0].Z, Time[t]);
            var b_a = GetAAt(recivers[1].X, recivers[1].Y, recivers[1].Z, Time[t]);
            var c_a = GetAAt(recivers[2].X, recivers[2].Y, recivers[2].Z, Time[t]);
            var d_a = GetAAt(recivers[3].X, recivers[3].Y, recivers[3].Z, Time[t]);
            var f_a_a = Math.Sqrt(a_a.Item1 * a_a.Item1 + a_a.Item2 * a_a.Item2 + a_a.Item3 * a_a.Item3);
            var f_b_a = Math.Sqrt(b_a.Item1 * b_a.Item1 + b_a.Item2 * b_a.Item2 + b_a.Item3 * b_a.Item3);
            var f_c_a = Math.Sqrt(c_a.Item1 * c_a.Item1 + c_a.Item2 * c_a.Item2 + c_a.Item3 * c_a.Item3);
            var f_d_a = Math.Sqrt(d_a.Item1 * d_a.Item1 + d_a.Item2 * d_a.Item2 + d_a.Item3 * d_a.Item3);
            var a_e = GetEAt(recivers[0].X, recivers[0].Y, recivers[0].Z, Time[t]);
            var b_e = GetEAt(recivers[1].X, recivers[1].Y, recivers[1].Z, Time[t]);
            var c_e = GetEAt(recivers[2].X, recivers[2].Y, recivers[2].Z, Time[t]);
            var d_e = GetEAt(recivers[3].X, recivers[3].Y, recivers[3].Z, Time[t]);
            var f_a_e = Math.Sqrt(a_e.Item1 * a_e.Item1 + a_e.Item2 * a_e.Item2 + a_e.Item3 * a_e.Item3);
            var f_b_e = Math.Sqrt(b_e.Item1 * b_e.Item1 + b_e.Item2 * b_e.Item2 + b_e.Item3 * b_e.Item3);
            var f_c_e = Math.Sqrt(c_e.Item1 * c_e.Item1 + c_e.Item2 * c_e.Item2 + c_e.Item3 * c_e.Item3);
            var f_d_e = Math.Sqrt(d_e.Item1 * d_e.Item1 + d_e.Item2 * d_e.Item2 + d_e.Item3 * d_e.Item3);
            sw_a.WriteLine($"{Time[t]} {f_a_a} {f_b_a} {f_c_a} {f_d_a}");
            sw_e.WriteLine($"{Time[t]} {f_a_e} {f_b_e} {f_c_e} {f_d_e}");
        }
        sw_a.Close();
        sw_e.Close();
    }
}
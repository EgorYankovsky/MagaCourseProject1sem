using MathObjects;
using DataStructs;
using Grid;
using System.Diagnostics;
using Functions;

namespace Project;

public class FEM2D : FEM
{
    public ArrayOfPoints2D pointsArr = new(_pointspath2D);

    public FEM2D(Mesh2Dim mesh, TimeMesh timeMesh) : base(timeMesh)
    {
        mesh2Dim = mesh;
        if (timeMesh[0] == timeMesh[^1])
            equationType = EquationType.Elliptic;
        else
            equationType = EquationType.Parabolic;
        
        elemsArr = new(_elemspath2D);
        bordersArr = new(_borderspath2D);

        A_phi = new GlobalVector[Time.Count];
        E_phi = new GlobalVector[Time.Count];
        Debug.WriteLine("Generated data submited");
    }

    private readonly Mesh2Dim mesh2Dim;
    public GlobalVector[] A_phi;
    public GlobalVector[] E_phi;

    public void Solve()
    {
        if (solver is null) throw new ArgumentNullException("solver is null !");
        if (Time is null) throw new ArgumentNullException("Time is null!");

        Stopwatch solutionStopwatch = new();
        solutionStopwatch.Start();

        Debug.WriteLine($"\nTime layer: before BC");
        Thread.Sleep(1500);

        Matrix = new GlobalMatrix(pointsArr.GetLength());
        Generator.BuildPortait(ref Matrix, pointsArr.GetLength(), elemsArr);
        Generator.FillMatrix(ref Matrix, pointsArr, elemsArr, TypeOfMatrixM.Mrr);

        Vector = new GlobalVector(pointsArr.GetLength());
        Generator.FillVector(ref Vector, pointsArr, elemsArr, Time[0]);

        Generator.ConsiderBoundaryConditions(ref Matrix, ref Vector, pointsArr, bordersArr, Time[0]);
        (Solutions[0], Discrepancy[0]) = solver.Solve(Matrix, Vector);

        if (Time.Count > 1)
        {
            (Solutions[1], Discrepancy[1]) = (Solutions[0], Discrepancy[0]);
            if (Time.Count > 2)
            for (int i = 2; i < Time.Count; i++)
            {
                Console.WriteLine($"\n {i} / {Time.Count - 1}. Time layer: {Time[i]}");
                //Thread.Sleep(1500);

                double deltT = Time[i] - Time[i - 2];
                double deltT0 = Time[i] - Time[i - 1];
                double deltT1 = Time[i - 1] - Time[i - 2];
                double tau0 = (deltT + deltT0) / (deltT * deltT0);
                double tau1 = deltT / (deltT1 * deltT0);
                double tau2 = deltT0 / (deltT * deltT1);
                double deltT = Time[i] - Time[i - 1];
                double tau = 1.0D / deltT;

                var matrix1 = new GlobalMatrix(pointsArr.GetLength());
                Generator.BuildPortait(ref matrix1, pointsArr.GetLength(), elemsArr);
                Generator.FillMatrix(ref matrix1, pointsArr, elemsArr, TypeOfMatrixM.Mrr);
                var M = new GlobalMatrix(pointsArr.GetLength()); // ???
                Generator.BuildPortait(ref M, pointsArr.GetLength(), elemsArr);
                Generator.FillMatrix(ref M, pointsArr, elemsArr, TypeOfMatrixM.Mr);

                var bi = new GlobalVector(pointsArr.GetLength());
                Matrix = (tau0 * M) + matrix1;
                Vector = bi + (tau1 * (M * Solutions[i - 1])) - (tau2 * (M * Solutions[i - 2]));
                Generator.ConsiderBoundaryConditions(ref Matrix, ref Vector, pointsArr, bordersArr, Time[i]);        
                (Solutions[i], Discrepancy[i]) = solver.Solve(Matrix, Vector);
            }
        }
        A_phi = Solutions;
        solutionStopwatch.Stop();
        var milseconds = solutionStopwatch.ElapsedMilliseconds;
        Console.WriteLine($"Lin eq solved for {milseconds / 60000} min {(milseconds % 60000) / 1000} sec");
    }

    public void WriteData()
    {
        if (Answer is null)
            throw new Exception("Vector _answer is null");
        for (int i = 0; i < Answer.Size; i++)
            Console.WriteLine($"{Answer[i]:E15}");
    }

    public void WriteData(string _path)
    {
        if (A_phi is null) throw new ArgumentNullException();
        if (E_phi is null) throw new ArgumentNullException();

        if (Time.Count != 1)
        {
            for (int i = 0; i < Time.Count; i++)
            {
                using var sw = new StreamWriter($"{_path}\\A_phi\\Answer\\Answer_Aphi_time={Time[i]}.dat");
                for (int j = 0;   j < A_phi[i].Size; j++)
                    sw.WriteLine($"{A_phi[i][j]:E8}");
                sw.Close();
            }
            for (int i = 0; i < Time.Count; i++)
            {
                using var sw = new StreamWriter($"{_path}\\E_phi\\Answer\\Answer_Ephi_time={Time[i]}.dat");
                for (int j = 0;   j < E_phi[i].Size; j++)
                    sw.WriteLine($"{E_phi[i][j]:E8}");
                sw.Close();
            }
        }
        else
        {
            using var sw = new StreamWriter($"{_path}\\A_phi\\Answer\\Answer.dat");
            for (int j = 0; j < A_phi[0].Size; j++)
                sw.WriteLine($"{A_phi[0][j]:E8}");
            sw.Close();
            using var sw1 = new StreamWriter($"{_path}\\E_phi\\Answer\\Answer.dat");
            for (int j = 0; j < E_phi[0].Size; j++)
                sw1.WriteLine($"{E_phi[0][j]:E8}");
            sw1.Close();
        }
    }

    public void WriteDiscrepancy(string _path)
    {
        if (A_phi is null) throw new ArgumentNullException();
        if (E_phi is null) throw new ArgumentNullException();

        if (Time.Count != 1)
        {
            for (int i = 0; i < Time.Count; i++)
            {
                using var sw_d = new StreamWriter($"{_path}\\A_phi\\Discrepancy\\Discrepancy_Aphi_time={Time[i]}.dat");

                int NotNaNamount = 0;
                double maxDisc = 0.0;
                double avgDisc = 0.0;
                double sumU = 0.0D;
                double sumD = 0.0D;

                List<double> TheorAnswer = [];
                foreach (var Z in mesh2Dim.nodesZ)
                    foreach (var R in mesh2Dim.nodesR)
                        TheorAnswer.Add(Function.U(R, Z, Time[i]));

                for (int j = 0; j < A_phi[i].Size; j++)
                {
                    double absDiff = Math.Abs(A_phi[i][j] - TheorAnswer[j]);
                    double currDisc = Math.Abs((A_phi[i][j] - TheorAnswer[j]) / TheorAnswer[j]);
                    
                    if (Math.Abs(maxDisc) < Math.Abs(currDisc))
                        maxDisc = currDisc;

                    if (!double.IsNaN(currDisc) && currDisc > 1E-14)
                    {
                        avgDisc += currDisc;
                        NotNaNamount++;
                        sumU += absDiff * absDiff;
                        sumD += TheorAnswer[j] * TheorAnswer[j];
                    }
                    sw_d.WriteLine($"{absDiff:E8} {currDisc:E8}");
                }
                avgDisc = Math.Sqrt(sumU) / Math.Sqrt(sumD);
                sw_d.WriteLine($"Средняя невязка: {avgDisc:E15}");
                sw_d.WriteLine($"Максимальная невязка: {maxDisc:E15}");
                sw_d.WriteLine($"С: {avgDisc:E7}");
                sw_d.WriteLine($"М: {maxDisc:E7}");
                sw_d.Close();
            }
        }
        else
        {
            using var sw_d = new StreamWriter($"{_path}\\A_phi\\Discrepancy\\Discrepancy_Aphi.dat");

            int NotNaNamount = 0;
                double maxDisc = 0.0;
                double avgDisc = 0.0;
                double sumU = 0.0D;
                double sumD = 0.0D;
                List<double> TheorAnswer = [];
                foreach (var Z in mesh2Dim.nodesZ)
                    foreach (var R in mesh2Dim.nodesR)
                        TheorAnswer.Add(Function.U(R, Z, 0.0D));
                for (int j = 0; j < A_phi[0].Size; j++)
                {
                    double absDiff = Math.Abs(A_phi[0][j] - TheorAnswer[j]);
                    double currDisc = Math.Abs((A_phi[0][j] - TheorAnswer[j]) / TheorAnswer[j]);
                    
                    if (Math.Abs(maxDisc) < Math.Abs(currDisc))
                        maxDisc = currDisc;

                    if (!double.IsNaN(currDisc) && currDisc > 1E-14)
                    {
                        avgDisc += currDisc;
                        NotNaNamount++;
                        sumU += absDiff * absDiff;
                        sumD += TheorAnswer[j] * TheorAnswer[j];
                    }
                    sw_d.WriteLine($"{absDiff:E8} {currDisc:E8}");
                }
                avgDisc = Math.Sqrt(sumU) / Math.Sqrt(sumD);
                sw_d.WriteLine($"Средняя невязка: {avgDisc:E15}");
                sw_d.WriteLine($"Максимальная невязка: {maxDisc:E15}");
                sw_d.WriteLine($"С: {avgDisc:E7}");
                sw_d.WriteLine($"М: {maxDisc:E7}");
                sw_d.Close();
        }
    }

    public void GenerateVectorEphi()
    {
        E_phi = new GlobalVector[A_phi.Length];   
        for (int i = 0; i < E_phi.Length; i++)
        {
            if (i == 0)
                E_phi[i] = new GlobalVector(A_phi[i].Size);
            else
                E_phi[i] = -1.0D / (Time[i] - Time[i - 1]) * (A_phi[i] - A_phi[i - 1]);
        }
    }

    internal List<int>? GetElem(double r, double z)
    {
        if (r < mesh2Dim.nodesR[0] || mesh2Dim.nodesR[^1] < r || z < mesh2Dim.nodesZ[0] || mesh2Dim.nodesZ[^1] < z)
            return null;
        int i = 0;
        for (; i < mesh2Dim.nodesR.Count - 1 && r >= 0.001; i++)
            if (mesh2Dim.nodesR[i] <= r && r <= mesh2Dim.nodesR[i + 1])
                break;
        int j = 0;
        for (; j < mesh2Dim.nodesZ.Count - 1; j++)
            if (mesh2Dim.nodesZ[j] <= z && z <= mesh2Dim.nodesZ[j + 1])
                break;
        return elemsArr[j * (mesh2Dim.nodesR.Count - 1) + i].Arr;
    }
    
    public double GetA_phiAt(double r, double z, double t)
    {
        for (int tt = 0; tt < Time.Count; tt++)
        {
            if (Time[tt] == t)
            {
                var elem = GetElem(r, z);
                if (elem is null) return 0.0D;
                double[] q = new double[4];
                for (int i = 0; i < 4; i++)
                    q[i] = A_phi[tt][elem[i]];
                double r0 = pointsArr[elem[0]].R;
                double r1 = pointsArr[elem[3]].R;
                double z0 = pointsArr[elem[0]].Z;
                double z1 = pointsArr[elem[3]].Z;
                return BasisFunctions2D.GetValue(q[0], q[1], q[2], q[3], r0, r1, z0, z1, r, z);        
            }
        }
        throw new Exception("Out of mesh borders");
    }

    public double GetE_phiAt(double r, double z, double t)
    {
        for (int tt = 0; tt < Time.Count; tt++)
        {
            if (Time[tt] == t)
            {
                var elem = GetElem(r, z);
                if (elem is null) return 0.0D;
                double[] q = new double[4];
                for (int i = 0; i < 4; i++)
                    q[i] = E_phi[tt][elem[i]];
                double r0 = pointsArr[elem[0]].R;
                double r1 = pointsArr[elem[3]].R;
                double z0 = pointsArr[elem[0]].Z;
                double z1 = pointsArr[elem[3]].Z;
                return BasisFunctions2D.GetValue(q[0], q[1], q[2], q[3], r0, r1, z0, z1, r, z);        
            }
        }
        throw new Exception("Out of mesh borders");
    }

    public void ReadAnswer(string AnswerPath)
    {
        string file = string.Empty;
        if (equationType == EquationType.Elliptic)
        {
            file = "Answer.dat";
            var fileData = File.ReadAllText(AnswerPath + "A_phi/Answer/" + file).Split("\n");
            A_phi[0] = new GlobalVector(fileData.Length - 1);
            for (int i = 0; i < fileData.Length - 1; i++)
                A_phi[0][i] = double.Parse(fileData[i]);
            fileData = File.ReadAllText(AnswerPath + "E_phi/Answer/" + file).Split("\n");
            E_phi[0] = new GlobalVector(fileData.Length - 1);
            for (int i = 0; i < fileData.Length - 1; i++)
                E_phi[0][i] = double.Parse(fileData[i]);
        }
        else
        {
            for (int t = 0; t < Time.Count; t++)
            {
                file = $"Answer_Aphi_time={Time[t]}.dat";
                var fileData = File.ReadAllText(AnswerPath + "A_phi/Answer/" + file).Split("\n");
                A_phi[t] = new GlobalVector(fileData.Length - 1);
                for (int i = 0; i < fileData.Length - 1; i++)
                    A_phi[t][i] = double.Parse(fileData[i]);
                file = $"Answer_Ephi_time={Time[t]}.dat";
                fileData = File.ReadAllText(AnswerPath + "E_phi/Answer/" + file).Split("\n");
                E_phi[t] = new GlobalVector(fileData.Length - 1);
                for (int i = 0; i < fileData.Length - 1; i++)
                    E_phi[t][i] = double.Parse(fileData[i]);
            }
        }
    }

    public void MeasureValuesOnReceivers(List<Point3D> recivers, string OutputPath)
    {
        using var sw_a = new StreamWriter(OutputPath + "A.txt");
        using var sw_e = new StreamWriter(OutputPath + "E.txt");
        List<(double, double)> pnt2D = [];
        foreach (var reciver in recivers)
            pnt2D.Add((Math.Sqrt(reciver.X * reciver.X + reciver.Y * reciver.Y), reciver.Z));
        for (int t = 0; t < Time.Count; t++)
        {
            var a_a = GetA_phiAt(pnt2D[0].Item1, pnt2D[0].Item2, Time[t]);
            var b_a = GetA_phiAt(pnt2D[1].Item1, pnt2D[1].Item2, Time[t]);
            var c_a = GetA_phiAt(pnt2D[2].Item1, pnt2D[2].Item2, Time[t]);
            var d_a = GetA_phiAt(pnt2D[3].Item1, pnt2D[3].Item2, Time[t]);
            var a_e = GetE_phiAt(pnt2D[0].Item1, pnt2D[0].Item2, Time[t]);
            var b_e = GetE_phiAt(pnt2D[1].Item1, pnt2D[1].Item2, Time[t]);
            var c_e = GetE_phiAt(pnt2D[2].Item1, pnt2D[2].Item2, Time[t]);
            var d_e = GetE_phiAt(pnt2D[3].Item1, pnt2D[3].Item2, Time[t]);
            sw_a.WriteLine($"{Time[t]} {a_a} {b_a} {c_a} {d_a}");
            sw_e.WriteLine($"{Time[t]} {a_e} {b_e} {c_e} {d_e}");
        }
        sw_a.Close();
        sw_e.Close();
    }

    public void WritePointsToDraw(string pathA, string pathE)
    {
        double hr = (mesh2Dim.nodesR[^1] - mesh2Dim.nodesR[0]) / 150;
        double hz = (mesh2Dim.nodesZ[^1] - mesh2Dim.nodesZ[0]) / 150;
        for (int t = 0; t < Time.Count; t++)
        {
            using var swa = new StreamWriter(pathA + $"Answer_A_time_layer_{t}.txt");
            for (int j = 0; j < 150; j++)
            {
                for (int i = 0; i < 150; i++)
                {
                    double rCurr = mesh2Dim.nodesR[0] + i * hr;
                    double zCurr = mesh2Dim.nodesZ[0] + j * hz;
                    swa.WriteLine($"{rCurr:E15} {zCurr:E15} {GetA_phiAt(rCurr, zCurr, Time[t]):E15}");
                }
            }
            swa.Close();
        }
        for (int t = 0; t < Time.Count; t++)
        {
            using var swe = new StreamWriter(pathE + $"Answer_E_time_layer_{t}.txt");
            for (int j = 0; j < 150; j++)
            {
                for (int i = 0; i < 150; i++)
                {
                    double rCurr = mesh2Dim.nodesR[0] + i * hr;
                    double zCurr = mesh2Dim.nodesZ[0] + j * hz;
                    swe.WriteLine($"{rCurr:E15} {zCurr:E15} {GetE_phiAt(rCurr, zCurr, Time[t]):E15}");
                }
            }
            swe.Close();
        }
    }
}
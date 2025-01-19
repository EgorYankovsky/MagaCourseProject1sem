using System.Collections.Immutable;
using DataStructs;

namespace Grid;

public static class MeshGenerator
{
    private static readonly string _3dValuesPath = Path.GetFullPath("../../../../Data/Subtotals/3_dim/");
    private static readonly string _points2DPath = Path.GetFullPath("../../../../Data/Subtotals/2_dim/Points.poly");
    private static readonly string _elems2DPath = Path.GetFullPath("../../../../Data/Subtotals/2_dim/Elems.poly");
    private static readonly string _borders2DPath = Path.GetFullPath("../../../../Data/Subtotals/2_dim/Borders.poly");

    public static TimeMesh GenerateTimeMesh(double t0, double t1, int tn, double tk)
    {
        double[] arr = new double[tn + 1];
        double h = t1 - t0;
        double denominator = 0.0;

        for (int j = 0; j < tn; j++)
            denominator += Math.Pow(tk, j);

        double x0 = h / denominator;
        arr[0] = t0;
        for(int j = 0; j < tn; j++)
            arr[j + 1] = arr[j] + x0 * Math.Pow(tk, j);
        arr[^1] = t1;
        return new TimeMesh(arr);
    }
    public static void GenerateMesh(Mesh2Dim mesh)
    {
        int currentPosition = 0;
        string[] kek = mesh.infoAboutR.Split();
        for (int i = 0; i < mesh.NodesRWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesR[1 + currentPosition] - mesh.nodesR[currentPosition];
            double denominator = 0.0;
            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);
            double x0 = h / denominator;
            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesR.Insert(currentPosition + 1, mesh.nodesR[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesR_Refs.Add(currentPosition);
        }
        currentPosition = 0;
        kek = mesh.infoAboutZ.Split();
        for (int i = 0; i < mesh.NodesZWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesZ[1 + currentPosition] - mesh.nodesZ[currentPosition];
            double denominator = 0.0;
            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);
            double x0 = h / denominator;
            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesZ.Insert(currentPosition + 1, mesh.nodesZ[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesZRefs.Add(currentPosition);
        }
    }

    public static void ConstructMesh(ref Mesh2Dim mesh)
    {
        GenerateMesh(mesh);
        RemakeBorders(mesh);
        OutputPoints(mesh);
        OutputListOfElems(mesh);
        OutputListOfBorders(mesh);
    }

    public static void RemakeBorders(Mesh2Dim mesh)
    {
        for (int i = 0; i < mesh.borders.Count; i++)
        {
            var newBorder = new Border2D(mesh.borders[i].BorderType, mesh.borders[i].BorderFormula,
                                         mesh.NodesR_Refs[mesh.borders[i].R0], mesh.NodesR_Refs[mesh.borders[i].R1],
                                         mesh.NodesZRefs[mesh.borders[i].Z0], mesh.NodesZRefs[mesh.borders[i].Z1]);
            mesh.borders[i] = newBorder;
        }
    }

    public static void OutputPoints(Mesh2Dim mesh)
    {
        using var sw = new StreamWriter(_points2DPath);
        sw.WriteLine(mesh.NodesAmountTotal);
        int i = 0;
        foreach (var Z in mesh.nodesZ)
            foreach (var R in mesh.nodesR)
            {
                sw.WriteLine($"{SetPointType(mesh, new Point2D(R, Z))}");
                i++;
            }
    }

    private static Point2D SetPointType(Mesh2Dim mesh, Point2D pnt)
    {
        for (int i = 0; i < mesh.Elems.Count; i++)
        {   
            if (mesh.NodesZWithoutFragmentation[mesh.Elems[i].Arr[4]] <= pnt.Z && pnt.Z <= mesh.NodesZWithoutFragmentation[mesh.Elems[i].Arr[5]])
            {
                int minValue = 4;
                foreach (var arr in mesh.borders)
                {
                    if (mesh.nodesR[arr.R0] <= pnt.R && pnt.R <= mesh.nodesR[arr.R1] &&
                        mesh.nodesZ[arr.Z0] <= pnt.Z && pnt.Z <= mesh.nodesZ[arr.Z1] &&
                        arr.BorderType < minValue)
                        {
                            minValue = arr.BorderType;
                            break;
                        }
                }
                switch (minValue)
                {
                    case 1: {pnt.Type = Location.BoundaryI; break; }
                    case 2: {pnt.Type = Location.BoundaryII; break; }
                    case 3: {pnt.Type = Location.BoundaryIII; break; }
                    default: {pnt.Type = Location.Inside; pnt.SubElemNum = i; break;}
                }
            }
            else if (pnt.Type == Location.NotStated)
                pnt.Type = Location.OutSide;
        }
        return pnt;
    }

    public static void OutputListOfElems(Mesh2Dim mesh)
    {
        using var sw = new StreamWriter(_elems2DPath);
        int nr = mesh.NodesAmountR;
        int nz = mesh.NodesAmountZ;
        sw.WriteLine((nr - 1) * (nz - 1));
        for (int k = 0; k < nz - 1; k++)
            for (int i = 0; i < nr - 1; i++)
            {
                sw.WriteLine($"{k * nr + i} {k * nr + i + 1} " +
                             $"{(k + 1) * nr + i} {(k + 1) * nr + i + 1} " +
                             $"{SelectMuAndSigma(mesh, k, k + 1)}");
            }
    }

    private static string SelectMuAndSigma(Mesh2Dim mesh, int a, int d)
    {
        int index = -1;
        for (int i = 0; i < mesh.NodesZRefs.Count - 1; i++)
            if (mesh.NodesZRefs[i] <= a && d <= mesh.NodesZRefs[i + 1])
                index = i;
        if (index == -1) throw new IndexOutOfRangeException("Failure during selecting mu and sigma");
        return $"{mesh.Elems[index].mu} {mesh.Elems[index].sigma}";
    }

    public static void OutputListOfBorders(Mesh2Dim mesh)
    {
        using var sw = new StreamWriter(_borders2DPath);
        sw.WriteLine(2 * (mesh.NodesAmountR - 1 + mesh.NodesAmountZ - 1));
        foreach (var border in mesh.borders)
        {
            if (border.R0 == border.R1) // border || oR
            {
                int iter = border.R0 == 0 ? 0 : mesh.NodesAmountR - 1;
                for (int i = 0; i < mesh.NodesAmountZ - 1; i++)
                {
                    sw.WriteLine($"{border.BorderType} {border.BorderFormula} {iter} {iter + mesh.NodesAmountR}");
                    iter += mesh.NodesAmountR;
                }
            }
            else if (border.Z0 == border.Z1) // border || oZ
            {
                int iter = border.Z0 == 0 ? 0 : mesh.NodesAmountR * (mesh.NodesAmountZ - 1);
                for (int i = 0; i < mesh.NodesAmountR - 1; i++)
                {
                    sw.WriteLine($"{border.BorderType} {border.BorderFormula} {iter} {iter + 1}");
                    iter++;
                }
            }
        }
    }

    public static void ConstructMesh(ref Mesh3Dim mesh)
    {
        GenerateMesh(ref mesh);
        RemakeBorders(ref mesh);
        var arr = OutputPoints(ref mesh);
        GenerateListOfRibs(ref mesh, arr);
        OutputListOfElems(ref mesh); 
        OutputListOfBorders(mesh);
    }

    public static void ConstructMesh(ref Mesh3Dim mesh, Layer layer, int layerNum)
    {
        GenerateMesh(ref mesh, layer);
        RemakeBorders(ref mesh);
        var arr = OutputPoints(ref mesh, layerNum);
        GenerateListOfRibs(ref mesh, arr, layerNum);
        OutputListOfElems(ref mesh, layerNum);
        OutputListOfBorders(mesh, layerNum);
    }

    public static void ConstructMeshAnomaly(ref Mesh3Dim mesh, string path)
    {
        GenerateNodesAnomaly(ref mesh);
        RemakeBorders(ref mesh);
        var arr = OutputPoints(ref mesh, path);
        var arr1 = GenerateListOfRibs(ref mesh, arr, path);
        OutputListOfElems(ref mesh, arr1, path);
        OutputListOfBorders(mesh, path);
    }

    public static void GenerateNodesAnomaly(ref Mesh3Dim mesh)
    {
        int currentPosition = 0;
        string[] kek = mesh.infoAboutX.Split();
        for (int i = 0; i < mesh.NodesXWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesX[1 + currentPosition] - mesh.nodesX[currentPosition];
            double denominator = 0.0;
            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);
            double x0 = h / denominator;
            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesX.Insert(currentPosition + 1, mesh.nodesX[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesXRefs.Add(currentPosition);
        }
        mesh.AnomalyBorders[0] = mesh.NodesXRefs[mesh.AnomalyBorders[0]];
        mesh.AnomalyBorders[1] = mesh.NodesXRefs[mesh.AnomalyBorders[1]];
        currentPosition = 0;
        kek = mesh.infoAboutY.Split();
        for (int i = 0; i < mesh.NodesYWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesY[1 + currentPosition] - mesh.nodesY[currentPosition];
            double denominator = 0.0;
            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);
            double x0 = h / denominator;
            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesY.Insert(currentPosition + 1, mesh.nodesY[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesYRefs.Add(currentPosition);
        }
        mesh.AnomalyBorders[2] = mesh.NodesYRefs[mesh.AnomalyBorders[2]];
        mesh.AnomalyBorders[3] = mesh.NodesYRefs[mesh.AnomalyBorders[3]];
        
        mesh.nodesZ = [.. mesh.NodesZWithoutFragmentation];
        mesh.NodesZWithoutFragmentation = [.. mesh.nodesZ];
        kek = mesh.infoAboutZ.Split();
        currentPosition = 0;
        for (int i = 0; i < mesh.NodesZWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesZ[1 + currentPosition] - mesh.nodesZ[currentPosition];
            double denominator = 0.0;
            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);
            double x0 = h / denominator;
            for(int j = 0; j < negr - 1 && mesh.nodesZ[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j) < mesh.NodesZWithoutFragmentation[^1]; j++)
            {
                mesh.nodesZ.Insert(currentPosition + 1, mesh.nodesZ[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesZRefs.Add(currentPosition);
        }
        mesh.AnomalyBorders[4] = mesh.NodesZRefs[mesh.AnomalyBorders[4]];
        mesh.AnomalyBorders[5] = mesh.NodesZRefs[mesh.AnomalyBorders[5]];
    }

    public static void GenerateMesh(ref Mesh3Dim mesh)
    {
        int currentPosition = 0;
        string[] kek = mesh.infoAboutX.Split();
        for (int i = 0; i < mesh.NodesXWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesX[1 + currentPosition] - mesh.nodesX[currentPosition];
            double denominator = 0.0;
            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);
            double x0 = h / denominator;
            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesX.Insert(currentPosition + 1, mesh.nodesX[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesXRefs.Add(currentPosition);
        }
        currentPosition = 0;
        kek = mesh.infoAboutY.Split();
        for (int i = 0; i < mesh.NodesYWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesY[1 + currentPosition] - mesh.nodesY[currentPosition];
            double denominator = 0.0;
            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);
            double x0 = h / denominator;
            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesY.Insert(currentPosition + 1, mesh.nodesY[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesYRefs.Add(currentPosition);
        }
        mesh.nodesZ = mesh.NodesZWithoutFragmentation.ToList();
        mesh.NodesZWithoutFragmentation = [.. mesh.nodesZ];
        kek = mesh.infoAboutZ.Split();
        currentPosition = 0;
        for (int i = 0; i < mesh.NodesZWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesZ[1 + currentPosition] - mesh.nodesZ[currentPosition];
            double denominator = 0.0;
            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);
            double x0 = h / denominator;
            for(int j = 0; j < negr - 1 && mesh.nodesZ[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j) < mesh.NodesZWithoutFragmentation[^1]; j++)
            {
                mesh.nodesZ.Insert(currentPosition + 1, mesh.nodesZ[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesZRefs.Add(currentPosition);
        }
    }

    public static void GenerateMesh(ref Mesh3Dim mesh, Layer layer)
    {
        int currentPosition = 0;
        string[] kek = mesh.infoAboutX.Split();
        for (int i = 0; i < mesh.NodesXWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesX[1 + currentPosition] - mesh.nodesX[currentPosition];
            double denominator = 0.0;
            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);
            double x0 = h / denominator;
            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesX.Insert(currentPosition + 1, mesh.nodesX[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesXRefs.Add(currentPosition);
        }
        currentPosition = 0;
        kek = mesh.infoAboutY.Split();
        for (int i = 0; i < mesh.NodesYWithoutFragmentation.Length - 1; i++)
        {
            double h = mesh.nodesY[1 + currentPosition] - mesh.nodesY[currentPosition];
            double denominator = 0.0;
            int negr = int.Parse(kek[2 * i]);
            for (int j = 0; j < negr; j++)
                denominator += Math.Pow(double.Parse(kek[2 * i + 1]), j);
            double x0 = h / denominator;
            for(int j = 0; j < negr - 1; j++)
            {
                mesh.nodesY.Insert(currentPosition + 1, mesh.nodesY[currentPosition] + x0 * Math.Pow(double.Parse(kek[2 * i + 1]), j));
                currentPosition++;
            }
            currentPosition++;
            mesh.NodesYRefs.Add(currentPosition);
        }
        List<double> ndsz = [mesh.nodesZ[0], layer.Z0, layer.Z1, mesh.nodesZ[^1]];
        ndsz.Sort();
        mesh.nodesZ = ndsz.Distinct().ToList();
        mesh.NodesZWithoutFragmentation = [.. mesh.nodesZ];
        currentPosition = 0;
        mesh.NodesZRefs = [0];
        for (int i = 0; i < mesh.NodesZWithoutFragmentation.Length - 1; i++)
        {
            // Генерируем внутри слоя.
            if (layer.Z0 == mesh.NodesZWithoutFragmentation[i] && layer.Z1 == mesh.NodesZWithoutFragmentation[i + 1])
            {
                double h = mesh.nodesZ[1 + currentPosition] - mesh.nodesZ[currentPosition];
                double denominator = 0.0;
                int negr = layer.Z1 - layer.Z0 > 30.0D ? (int)(layer.Z1 - layer.Z0) / 30 : 1;
                double k = 1.0D;
                for (int j = 0; j < negr; j++)
                    denominator += Math.Pow(k, j);
                double x0 = h / denominator;
                for(int j = 0; j < negr - 1; j++)
                {
                    mesh.nodesZ.Insert(currentPosition + 1, mesh.nodesZ[currentPosition] + x0 * Math.Pow(k, j));
                    currentPosition++;
                }
            }
            // Генерируем вне слоя.
            else
            {
                double h = mesh.nodesZ[1 + currentPosition] - mesh.nodesZ[currentPosition];
                double denominator = 0.0;
                int negr = (int)(h / (mesh.nodesZ[^1] - mesh.nodesZ[0]) * 100);
                double k = mesh.NodesZWithoutFragmentation[i] < layer.Z0 ? 1.0D / 1.05D : 1.2D;
                for (int j = 0; j < negr; j++)
                    denominator += Math.Pow(k, j);
                double x0 = mesh.NodesZWithoutFragmentation[i] == layer.Z1 ? 30.0D : h / denominator;
                for(int j = 0; j < negr - 1 && mesh.nodesZ[currentPosition] + x0 * Math.Pow(k, j) < mesh.NodesZWithoutFragmentation[^1]; j++)
                {
                    mesh.nodesZ.Insert(currentPosition + 1, mesh.nodesZ[currentPosition] + x0 * Math.Pow(k, j));
                    currentPosition++;
                }
            }
            currentPosition++;
            mesh.NodesZRefs.Add(currentPosition);
        }
    }

    public static void RemakeBorders(ref Mesh3Dim mesh)
    {
        for (int i = 0; i < mesh.borders.Count; i++)
        {
            var newBorder = new Border3D(mesh.borders[i].BorderType, mesh.borders[i].BorderFormula,
                                         mesh.NodesXRefs[mesh.borders[i].X0], mesh.NodesXRefs[mesh.borders[i].X1],
                                         mesh.NodesYRefs[mesh.borders[i].Y0], mesh.NodesYRefs[mesh.borders[i].Y1],
                                         mesh.NodesZRefs[mesh.borders[i].Z0], mesh.NodesZRefs[mesh.borders[i].Z1]);
            mesh.borders[i] = newBorder;
        }
    }

    public static ArrayOfPoints3D OutputPoints(ref Mesh3Dim mesh, string path)
    {
        using var sw = new StreamWriter(path + "Points.poly");
        var arrPnt = new ArrayOfPoints3D(mesh.NodesAmountTotal);
        sw.WriteLine(mesh.NodesAmountTotal);
        int i = 0;
        foreach (var Z in mesh.nodesZ)
            foreach (var Y in mesh.nodesY)
                foreach (var X in mesh.nodesX)
                {
                    var pnt = SetPointType3D(new Point3D(X, Y, Z), mesh);
                    sw.WriteLine($"{pnt}");
                    arrPnt.Append(pnt);
                    i++;
                }
        sw.Close();
        return arrPnt;
    }

    public static ArrayOfPoints3D OutputPoints(ref Mesh3Dim mesh, int layerNum = -1)
    {
        string path = _3dValuesPath + (layerNum == -1 ? $"AfterConvertation\\Points.poly" : $"Field{layerNum}/Points.poly");
        using var sw = new StreamWriter(path);
        var arrPnt = new ArrayOfPoints3D(mesh.NodesAmountTotal);
        sw.WriteLine(mesh.NodesAmountTotal);
        int i = 0;
        foreach (var Z in mesh.nodesZ)
            foreach (var Y in mesh.nodesY)
                foreach (var X in mesh.nodesX)
                {
                    var pnt = SetPointType3D(new Point3D(X, Y, Z), mesh);
                    sw.WriteLine($"{pnt}");
                    arrPnt.Append(pnt);
                    i++;
                }
        sw.Close();
        return arrPnt;
    }
    
    private static Point3D SetPointType3D(Point3D pnt, Mesh3Dim mesh)
    {
        double xMin = mesh.nodesX[0];
        double xMax = mesh.nodesX[^1];
        double yMin = mesh.nodesY[0];
        double yMax = mesh.nodesY[^1];
        double zMin = mesh.nodesZ[0];
        double zMax = mesh.nodesZ[^1];
        if (pnt.Z == zMax)
            pnt.Type = Location.BoundaryI;
        else if (pnt.Z == zMin || pnt.X == xMin ||pnt.X == xMax || pnt.Y == yMin || pnt.Y == yMax)
            pnt.Type = Location.BoundaryI;
        else
            pnt.Type = Location.Inside;
        return pnt;
    }
    
    public static void OutputListOfElems(ref Mesh3Dim mesh, ArrayOfRibs arrRibs, string path)
    {
        if (mesh.arrayOfRibs is null) throw new ArgumentNullException("array of ribs didn't generated");
        using var sw = new StreamWriter(path + "Elems.poly");
        int rx = mesh.NodesAmountX - 1;
        int ry = mesh.NodesAmountY - 1;
        int nx = mesh.NodesAmountX;
        int ny = mesh.NodesAmountY;
        int nz = mesh.NodesAmountZ;
        int rxy = rx * ny + ry * nx;
        int nxy = nx * ny;
        sw.WriteLine((nx - 1) * (ny - 1) * (nz - 1));
        for (int k = 0; k < nz - 1; k++)
            for (int j = 0; j < ny - 1; j++)
                for (int i = 0; i < nx - 1; i++)
                {
                    int curr = i + j * (nx + rx) + k * (rxy + nxy);
                    double mui = 0.0D;
                    double sigmai = 0.0D;
                    (mui, sigmai) = SelectMuAndSigma(mesh,
                                                     mesh.arrayOfRibs[curr].a.X, mesh.arrayOfRibs[curr].b.X,
                                                     mesh.arrayOfRibs[curr + rx].a.Y, mesh.arrayOfRibs[curr + rx].b.Y,
                                                     mesh.arrayOfRibs[curr + rxy - j * rx].a.Z, mesh.arrayOfRibs[curr + rxy - j * rx].b.Z);
                    List<int> arr_i = [                curr,               curr + rx,             curr + rx + 1,               curr + rx + nx,
                                        curr + rxy - j * rx, curr + rxy + 1 - j * rx,  curr + rxy + nx - j * rx, curr + rxy + nx + 1 - j * rx,
                                           curr + rxy + nxy,   curr + rxy + nxy + rx, curr + rxy + nxy + rx + 1,   curr + rxy + nxy + rx + nx];
                    sw.WriteLine($"{0} {arr_i[0]} {arr_i[1]} {arr_i[2]} {arr_i[3]}" + 
                                    $" {arr_i[4]} {arr_i[5]} {arr_i[6]} {arr_i[7]}" + 
                                    $" {arr_i[8]} {arr_i[9]} {arr_i[10]} {arr_i[11]} {mui} {sigmai}");
                }
        sw.Close();
    }

    public static void OutputListOfElems(ref Mesh3Dim mesh, int layerNum = -1)
    {
        if (mesh.arrayOfRibs is null) throw new ArgumentNullException("array of ribs didn't generated");
        string path = _3dValuesPath + (layerNum == -1 ? "AfterConvertation\\Elems.poly" : $"Field{layerNum}/Elems.poly");
        using var sw = new StreamWriter(path);
        int rx = mesh.NodesAmountX - 1;
        int ry = mesh.NodesAmountY - 1;      
        int nx = mesh.NodesAmountX;
        int ny = mesh.NodesAmountY;
        int nz = mesh.NodesAmountZ;
        int rxy = rx * ny + ry * nx;
        int nxy = nx * ny;
        sw.WriteLine((nx - 1) * (ny - 1) * (nz - 1));
        for (int k = 0; k < nz - 1; k++)
            for (int j = 0; j < ny - 1; j++)
                for (int i = 0; i < nx - 1; i++)
                {
                    int curr = i + j * (nx + rx) + k * (rxy + nxy);
                    double mui = 0.0D;
                    double sigmai = 0.0D;
                    (mui, sigmai) = SelectMuAndSigma(mesh, mesh.arrayOfRibs[curr + rxy - j * rx].a.Z, mesh.arrayOfRibs[curr + rxy - j * rx].b.Z);
                    sw.WriteLine($"{0} {curr} {curr + rx} {curr + rx + 1} {curr + rx + nx}" + 
                                 $" {curr + rxy - j * rx} {curr + rxy + 1 - j * rx} {curr + rxy + nx - j * rx} {curr + rxy + nx + 1 - j * rx}" + 
                                 $" {curr + rxy + nxy} {curr + rxy + nxy + rx} {curr + rxy + nxy + rx + 1} {curr + rxy + nxy + rx + nx} {mui} {sigmai}");
                }
        sw.Close();
    }

    private static (double, double) SelectMuAndSigma(Mesh3Dim mesh, double x0, double x1,
                                                                    double y0, double y1,
                                                                    double z0, double z1)
    {
        double mu = 0.0D;
        double sigma = 0.0D;
        for (int i = 0; i < mesh.Elems.Count; i++)
            if (mesh.NodesXWithoutFragmentation[mesh.Elems[i].Arr[0]] <= x0 && x1 <= mesh.NodesXWithoutFragmentation[mesh.Elems[i].Arr[1]] && 
                mesh.NodesYWithoutFragmentation[mesh.Elems[i].Arr[2]] <= y0 && y1 <= mesh.NodesYWithoutFragmentation[mesh.Elems[i].Arr[3]] && 
                mesh.NodesZWithoutFragmentation[mesh.Elems[i].Arr[4]] <= z0 && z1 <= mesh.NodesZWithoutFragmentation[mesh.Elems[i].Arr[5]])
            {
                mu = mesh.Elems[i].mu;
                sigma = mesh.Elems[i].sigma;
                return (mu, sigma);
            }
        throw new Exception("Out of ranges");
    }

    private static (double, double) SelectMuAndSigma(Mesh3Dim mesh, double z0, double z1)
    {
        int index = -1;
        if (mesh.NodesZRefs.Count > 2)
            if (z1 <= 0.0D)
                return (4.0D * Math.PI * Math.Pow(10, -7), 0.01D);
            else
                return (4.0D * Math.PI * Math.Pow(10, -7), 0.0D);
        for (int i = 0; i < mesh.NodesZRefs.Count - 1; i++)
            if (mesh.nodesZ[mesh.NodesZRefs[i]] <= z0 && z1 <= mesh.nodesZ[mesh.NodesZRefs[i + 1]])
                index = i;
        if (index == -1) throw new IndexOutOfRangeException("Failure during selecting mu and sigma");
        return (mesh.Elems[index].mu, mesh.Elems[index].sigma);
    }

    public static void SelectRibs(ref ArrayOfRibs arrRibs, ref ArrayOfElems arrEl)
    {
        int ii = 0;
        while (ii < arrRibs.Count)
        {
            Console.WriteLine($"{ii * 100.0D / arrRibs.Count:E5}\r");
            if (arrRibs[ii].typeOfRib == TypeOfRib.BoundaryI)
            {
                foreach (Elem elem in arrEl)
                    foreach (int item in elem)
                    {
                        if (item > ii) break;
                        if (item == ii)
                        {
                            elem.Remove(item);
                            break;
                        }
                    }
                arrRibs.Remove(ii);
                for (int i = 0; i < arrEl.Length; i++)
                    for (int j = 0; j < arrEl[i].Count; j++) 
                        if (arrEl[i][j] > ii)
                            arrEl[i].Arr[j] = arrEl[i].Arr[j] - 1;
            }
            else
                ii++;
        }
    }
    
    public static void OutputListOfBorders(Mesh3Dim mesh, string path)
    {
        using var sw = new StreamWriter(path + "Borders.poly");
        int nx = mesh.NodesAmountX;
        int ny = mesh.NodesAmountY;
        int nz = mesh.NodesAmountZ;
        int nxny = nx * ny;
        int rxy = (nx - 1) * ny + (ny - 1) * nx;
        sw.WriteLine(2 * ((nx - 1) * (ny - 1) + (nx - 1) * (nz - 1) + (ny - 1) * (nz - 1)));
        // XY0
        for (int i = 0; i < ny - 1; i++)
            for (int j = 0; j < nx - 1; j++)
                sw.WriteLine($"{1} {1} {i * (2 * nx - 1) + j} {i * (2 * nx - 1) + j + nx - 1} {i * (2 * nx - 1) + j + nx} {i * (2 * nx - 1) + j + nx + nx - 1}");
        // X0Z
        for (int i = 0; i < nz - 1; i++)
            for (int j = 0; j < nx - 1; j++)
                sw.WriteLine($"{1} {2} {i * (rxy + nxny) + j} {i * (rxy + nxny) + j + rxy} {i * (rxy + nxny) + j + rxy + 1} {i * (rxy + nxny) + j + rxy + nxny}");
        // 0YZ
        for (int i = 0; i < nz - 1; i++)
            for (int j = 0; j < ny - 1; j++)
                sw.WriteLine($"{1} {3} {nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1)} {rxy + i * (rxy + nxny) + j * nx} {rxy + nx + i * (rxy + nxny) + j * nx} {rxy + nxny + nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1)}");
        // XY1
        for (int i = 0; i < ny - 1; i++)
            for (int j = 0; j < nx - 1; j++)
                sw.WriteLine($"{1} {4} {(nz - 1) * (rxy + nxny) + j + i * (2 * nx - 1)} {(nz - 1) * (rxy + nxny) + nx - 1 + j + i * (2 * nx - 1)} {(nz - 1) * (rxy + nxny) + nx + j + i * (2 * nx - 1)} {(nz - 1) * (rxy + nxny) + nx + nx - 1 + j + i * (2 * nx - 1)}");        
        // X1Z
        for (int i = 0; i < nz - 1; i++)
            for (int j = 0; j < nx - 1; j++)
                sw.WriteLine($"{1} {5} {(ny - 1) * nx + (ny - 1) * (nx - 1) + j + i * (rxy + nxny)} {(ny - 1) * nx + (ny - 1) * (nx - 1) + nxny - 1 + j + i * (rxy + nxny)} {(ny - 1) * nx + (ny - 1) * (nx - 1) + nxny + j + i * (rxy + nxny)} {(ny - 1) * nx + (ny - 1) * (nx - 1) + nxny + j + rxy + i * (rxy + nxny)}");
        // 1YZ
        for (int i = 0; i < nz - 1; i++)
            for (int j = 0; j < ny - 1; j++)
                sw.WriteLine($"{1} {6} {nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1) + nx - 1} {rxy + i * (rxy + nxny) + j * nx + nx - 1} {rxy + nx + i * (rxy + nxny) + j * nx + nx - 1} {rxy + nxny + nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1) + nx - 1}");
    }

    public static void OutputListOfBorders(Mesh3Dim mesh, int layerNum = -1)
    {
        string path = _3dValuesPath + (layerNum == -1 ? "AfterConvertation\\Borders.poly" : $"Field{layerNum}/Borders.poly");
        using var sw = new StreamWriter(path);
        int nx = mesh.NodesAmountX;
        int ny = mesh.NodesAmountY;
        int nz = mesh.NodesAmountZ;
        int nxny = nx * ny;
        int rxy = (nx - 1) * ny + (ny - 1) * nx;
        sw.WriteLine(2 * ((nx - 1) * (ny - 1) + (nx - 1) * (nz - 1) + (ny - 1) * (nz - 1)));
        // XY0
        for (int i = 0; i < ny - 1; i++)
            for (int j = 0; j < nx - 1; j++)
                sw.WriteLine($"{1} {1} {i * (2 * nx - 1) + j} {i * (2 * nx - 1) + j + nx - 1} {i * (2 * nx - 1) + j + nx} {i * (2 * nx - 1) + j + nx + nx - 1}");
        // X0Z
        for (int i = 0; i < nz - 1; i++)
            for (int j = 0; j < nx - 1; j++)
                sw.WriteLine($"{1} {2} {i * (rxy + nxny) + j} {i * (rxy + nxny) + j + rxy} {i * (rxy + nxny) + j + rxy + 1} {i * (rxy + nxny) + j + rxy + nxny}");
        // 0YZ
        for (int i = 0; i < nz - 1; i++)
            for (int j = 0; j < ny - 1; j++)
                sw.WriteLine($"{1} {3} {nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1)} {rxy + i * (rxy + nxny) + j * nx} {rxy + nx + i * (rxy + nxny) + j * nx} {rxy + nxny + nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1)}");
        // XY1
        for (int i = 0; i < ny - 1; i++)
            for (int j = 0; j < nx - 1; j++)
                sw.WriteLine($"{1} {4} {(nz - 1) * (rxy + nxny) + j + i * (2 * nx - 1)} {(nz - 1) * (rxy + nxny) + nx - 1 + j + i * (2 * nx - 1)} {(nz - 1) * (rxy + nxny) + nx + j + i * (2 * nx - 1)} {(nz - 1) * (rxy + nxny) + nx + nx - 1 + j + i * (2 * nx - 1)}");
        // X1Z
        for (int i = 0; i < nz - 1; i++)
            for (int j = 0; j < nx - 1; j++)
                sw.WriteLine($"{1} {5} {(ny - 1) * nx + (ny - 1) * (nx - 1) + j + i * (rxy + nxny)} {(ny - 1) * nx + (ny - 1) * (nx - 1) + nxny - 1 + j + i * (rxy + nxny)} {(ny - 1) * nx + (ny - 1) * (nx - 1) + nxny + j + i * (rxy + nxny)} {(ny - 1) * nx + (ny - 1) * (nx - 1) + nxny + j + rxy + i * (rxy + nxny)}");
        // 1YZ
        for (int i = 0; i < nz - 1; i++)
            for (int j = 0; j < ny - 1; j++)
                sw.WriteLine($"{1} {6} {nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1) + nx - 1} {rxy + i * (rxy + nxny) + j * nx + nx - 1} {rxy + nx + i * (rxy + nxny) + j * nx + nx - 1} {rxy + nxny + nx - 1 + i * (rxy + nxny) + j * (2 * nx - 1) + nx - 1}");
    }
    
    public static ArrayOfRibs GenerateListOfRibs(ref Mesh3Dim mesh, ArrayOfPoints3D arrPt, string path)
    {
        ArgumentNullException.ThrowIfNull(arrPt);
        int nx = mesh.NodesAmountX;
        int ny = mesh.NodesAmountY;
        int nz = mesh.NodesAmountZ;
        int nxny = nx * ny;
        mesh.arrayOfRibs = new(3 * nx * ny * nz - nx * ny - nx * nz - ny * nz);
        // Генерируем список всех ребер.
        for (int k = 0; k < nz; k++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx - 1; i++) 
                    mesh.arrayOfRibs.Add(new Rib(arrPt[k * nxny + nx * j + i], arrPt[k * nxny + nx * j + i + 1]));
                if (j != ny - 1)
                    for (int i = 0; i < nx; i++)
                        mesh.arrayOfRibs.Add(new Rib(arrPt[k * nxny + nx * j + i], arrPt[k * nxny + nx * (j + 1) + i]));
            }
            if (k != nz - 1)
                for (int j = 0; j < ny; j++)
                    for (int i = 0; i < nx; i++)
                        mesh.arrayOfRibs.Add(new Rib(arrPt[k * nxny + nx * j + i], arrPt[(k + 1) * nxny + nx * j + i]));
        }
        using var sw = new StreamWriter(path + "Ribs.poly");
        sw.WriteLine($"{mesh.arrayOfRibs.Count}");
        for (int i = 0; i < mesh.arrayOfRibs.Count; i++)
            sw.WriteLine(mesh.arrayOfRibs[i]);
        sw.Close();
        return mesh.arrayOfRibs;
    }

    public static void GenerateListOfRibs(ref Mesh3Dim mesh, ArrayOfPoints3D arrPt, int layerNum = -1)
    {
        ArgumentNullException.ThrowIfNull(arrPt);
        int nx = mesh.NodesAmountX;
        int ny = mesh.NodesAmountY;
        int nz = mesh.NodesAmountZ;
        int nxny = nx * ny;
        mesh.arrayOfRibs = new(3 * nx * ny * nz - nx * ny - nx * nz - ny * nz);
        // Генерируем список всех ребер.
        for (int k = 0; k < nz; k++)
        {
            for (int j = 0; j < ny; j++)
            {
                for (int i = 0; i < nx - 1; i++) 
                    mesh.arrayOfRibs.Add(new Rib(arrPt[k * nxny + nx * j + i], arrPt[k * nxny + nx * j + i + 1]));
                if (j != ny - 1)
                    for (int i = 0; i < nx; i++)
                        mesh.arrayOfRibs.Add(new Rib(arrPt[k * nxny + nx * j + i], arrPt[k * nxny + nx * (j + 1) + i]));
            }
            if (k != nz - 1)
                for (int j = 0; j < ny; j++)
                    for (int i = 0; i < nx; i++)
                        mesh.arrayOfRibs.Add(new Rib(arrPt[k * nxny + nx * j + i], arrPt[(k + 1) * nxny + nx * j + i]));
        }
        string path = _3dValuesPath + (layerNum == -1 ? "AfterConvertation\\Ribs.poly" : $"Field{layerNum}\\Ribs.poly");
        using var sw = new StreamWriter(path);
        for (int i = 0; i < mesh.arrayOfRibs.Count; i++)
            sw.WriteLine(mesh.arrayOfRibs[i]);
        sw.Close();
    }
}
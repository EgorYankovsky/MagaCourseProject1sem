using System.Collections.Immutable;
using DataStructs;

namespace Grid;

public class Mesh2Dim(List<double> nodesR, string infoAboutR,
                      List<double> nodesZ, string infoAboutZ,
                      List<Elem> elems, double lastR) : Mesh(elems)
{
    public override int NodesAmountTotal => NodesAmountR * NodesAmountZ;
    public List<double> nodesR = [.. nodesR, lastR];
    public List<double> nodesZ = nodesZ;
    public override int ElemsAmount
    {
        get => (NodesAmountR - 1) * (NodesAmountZ - 1);
        set => ElemsAmount = value;
    }
    public int NodesAmountR => nodesR.Count;
    public int NodesAmountZ => nodesZ.Count;
    public List<int> NodesR_Refs = [0];
    public List<int> NodesZRefs = [0];
    internal ImmutableArray<double> NodesRWithoutFragmentation { get; set; } = [.. nodesR, lastR];
    public ImmutableArray<double> NodesZWithoutFragmentation { get; set; } = [.. nodesZ];
    internal string infoAboutR = infoAboutR;
    internal string infoAboutZ = infoAboutZ;
    internal List<Border2D> borders = [];
    public void SetLastR(double val)
    {
        int i = 1;
        if (val <= nodesR[^i])
        {
            while (val < nodesR[^i]) i++; 
            nodesR[^i] = val;
            nodesR = nodesR.Take(i + 1).ToList();
        }
        else
            nodesR.Add(val);
    }

    public void SetBorders(List<Border3D> borders3D)
    {
        // Set lower border.
        borders.Add(new Border2D(borders3D[0].BorderType, borders3D[0].BorderFormula,
                                 0, NodesRWithoutFragmentation.Length - 1, 0, 0));
        // Set left border.
        borders.Add(new Border2D(borders3D[1].BorderType, borders3D[1].BorderFormula,
                                 0, 0, 0, NodesZWithoutFragmentation.Length - 1));
        // Set right border.
        borders.Add(new Border2D(borders3D[3].BorderType, borders3D[3].BorderFormula,
                                 NodesRWithoutFragmentation.Length - 1, 
                                 NodesRWithoutFragmentation.Length - 1,
                                 0, NodesZWithoutFragmentation.Length - 1));
        // Set upper border.
        borders.Add(new Border2D(borders3D[^1].BorderType, borders3D[^1].BorderFormula,
                                 0, NodesRWithoutFragmentation.Length - 1,
                                 NodesZWithoutFragmentation.Length - 1, 
                                 NodesZWithoutFragmentation.Length - 1));
    }

    public override Mesh2Dim Clone() => (Mesh2Dim)MemberwiseClone();
}

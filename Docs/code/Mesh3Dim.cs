using System.Collections.Immutable;
using DataStructs;

namespace Grid;

public class Mesh3Dim(List<double> nodesX, string infoAboutX,
                      List<double> nodesY, string infoAboutY,
                      List<double> nodesZ, string infoAboutZ,
                      List<Elem> elems, List<Border3D> borders) : Mesh(elems)
{
    public override int NodesAmountTotal 
    { 
        get => NodesAmountX * NodesAmountY * NodesAmountZ;
    }

    public override int ElemsAmount
    {
        get => (NodesAmountX - 1) * (NodesAmountY - 1) * (NodesAmountZ - 1);
        set => ElemsAmount = value;
    }

    public int NodesAmountX => nodesX.Count;
    public int NodesAmountY => nodesY.Count;
    public int NodesAmountZ => nodesZ.Count;
    public List<double> nodesX = nodesX;
    public List<double> nodesY = nodesY;
    public List<double> nodesZ = nodesZ;
    public void CommitFieldBorders(List<int> ints) => FieldBorders = ints;
    public void CommitAnomalyBorders(List<int> ints) => AnomalyBorders = ints;
    public void CommitSecondAnomalyBorders(List<int> ints) => SecondAnomaly = ints;
    public List<int> FieldBorders = [];
    public List<int> AnomalyBorders = [];
    public List<int> SecondAnomaly = [];
    internal ImmutableArray<double> NodesXWithoutFragmentation { get; set; } = [.. nodesX];
    internal ImmutableArray<double> NodesYWithoutFragmentation { get; set; } = [.. nodesY];
    public ImmutableArray<double> NodesZWithoutFragmentation { get; set; } = [.. nodesZ];
    internal string infoAboutX = infoAboutX;
    internal string infoAboutY = infoAboutY;
    internal string infoAboutZ = infoAboutZ;
    public List<int> NodesXRefs = [0];
    public List<int> NodesYRefs = [0];
    public List<int> NodesZRefs = [0];
    public List<Border3D> borders = borders;
    public ArrayOfRibs? arrayOfRibs;
    public override Mesh3Dim Clone() => (Mesh3Dim)MemberwiseClone();
}
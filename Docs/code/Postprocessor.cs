using System.Diagnostics;

namespace Processor;

public static class Postprocessor
{
    private static readonly string _sourceToDrawPlot = "..\\..\\..\\..\\Drawer\\PictureDrawer2D.py";
    private static readonly string _readToDraw2_dim_A = Path.GetFullPath("..\\..\\..\\..\\Data\\Output\\ToDraw\\2_dim\\Aphi\\");   
    private static readonly string _readToDraw2_dim_E = Path.GetFullPath("..\\..\\..\\..\\Data\\Output\\ToDraw\\2_dim\\Ephi\\");
    private static readonly string _outputDrawn2_dim_A = Path.GetFullPath("..\\..\\..\\..\\Drawer\\Pictures\\A_phi\\");
    private static readonly string _outputDrawn2_dim_E = Path.GetFullPath("..\\..\\..\\..\\Drawer\\Pictures\\E_phi\\");
    private static readonly string _sourceGrapher2D = "..\\..\\..\\..\\Drawer\\2D_graphicsDrawer.py";
    private static readonly string _readDataForA = "..\\..\\..\\..\\Data\\Output\\ToDraw\\2_dim\\Receivers\\A.txt";
    private static readonly string _readDataForE = "..\\..\\..\\..\\Data\\Output\\ToDraw\\2_dim\\Receivers\\E.txt";
    private static readonly string _writePictures = "..\\..\\..\\..\\Drawer\\Graphics\\";
    private static readonly string _arg4 = Path.GetFullPath("../../../../Data/Output/E_phi/Answer/");
    private static readonly string _arg5 = Path.GetFullPath("../../../../Drawer/Pictures/E_phi/");
    private static readonly string _source3D = "..\\..\\..\\..\\Drawer\\3DplotVectors.py";
    private static readonly string _source3DRead = "..\\..\\..\\..\\Data\\Output\\E_phi\\ToDraw\\ConvertedTo3D\\";
    private static readonly string _source3DWrite = "..\\..\\..\\..\\Drawer\\Pictures\\E_phi3D\\";
    private static readonly string _fileName = "python";
    private static Process? process;

    public static int DrawA_phi()
    {
        process = new();
        process.StartInfo.Arguments = $"{_sourceToDrawPlot} {_readToDraw2_dim_A} {_outputDrawn2_dim_A}";
        process.StartInfo.FileName = _fileName;
        process.Start();
        process.WaitForExit();
        return process.ExitCode;
    }

    public static int DrawE_phi()
    {
        process = new();
        process.StartInfo.Arguments = $"{_sourceToDrawPlot} {_readToDraw2_dim_E} {_outputDrawn2_dim_E}";
        process.StartInfo.FileName = _fileName;
        process.Start();
        process.WaitForExit();
        return process.ExitCode;
    }

    public static int DrawE3D(string _readToDraw3_dim_E, string _outputDrawn3_dim_E)
    {
        process = new();
        process.StartInfo.Arguments = $"{_sourceToDrawPlot} {_readToDraw3_dim_E} {_outputDrawn3_dim_E}";
        process.StartInfo.FileName = _fileName;
        process.Start();
        process.WaitForExit();
        return process.ExitCode;
    }

    public static int DrawA3D(string _readToDraw3_dim_A, string _outputDrawn3_dim_A)
    {
        process = new();
        process.StartInfo.Arguments = $"{_sourceToDrawPlot} {_readToDraw3_dim_A} {_outputDrawn3_dim_A}";
        process.StartInfo.FileName = _fileName;
        process.Start();
        process.WaitForExit();
        return process.ExitCode;
    }

    public static int DrawGraphics2D()
    {
        process = new();
        process.StartInfo.Arguments = $"{_sourceGrapher2D} {_readDataForA} {_readDataForE} {_writePictures}";
        process.StartInfo.FileName = _fileName;
        process.Start();
        process.WaitForExit();
        return process.ExitCode;
    }
}
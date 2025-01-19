using System.Collections;

namespace Manager;

public static class FolderManager
{
    public static void ClearFolder(string path)
    {
        DirectoryInfo di = new(path);
        foreach (var dir in di.GetDirectories())
            ClearFolder(path + dir.Name + "\\");
        foreach (FileInfo file in di.GetFiles())
            file.Delete();
    }

    public static void ClearFolders(IEnumerable paths)
    {
        foreach (string path in paths)
            ClearFolder(path);
    }

    public static int CountFilesAmount(string path)
    {
        DirectoryInfo di = new (path);
        return di.GetFiles().Length;
    }
}
